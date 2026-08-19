"""
ERA5 Climate Matching to Parthenium Observations    
==========================================

Workflow:
1. Load ERA5 hourly data 
2. For each observation (date + location):
   - Extract ERA5 window: [observation_date - 30 days] to [observation_date]
   - Truncated linear GDD model: Tbase=7.2°C, Tmax=42.8°C
   - Compute 7-day rolling averages from this window
   - Compute 30-day rolling sums from this window
   - Extract computed values at observation_date
3. Save matched dataset with climate variables

Variable set:
  temp_mean_7d      — recent thermal conditions
  gdd_7d            — short-term heat accumulation (truncated)
  gdd_calendar      — seasonal heat budget: Jan 1 → obs date (truncated)
  precip_sum_30d    — monthly moisture total
  precip_freq_30d   — rain burst frequency

USAGE:
   python climate_match_v3.py           
"""

import sys
import numpy as np
import pandas as pd
import xarray as xr
import glob
import os
import warnings
# This only silences the "this feature will change in the future" messages
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# ==============================================================================
# CONFIG
# ==============================================================================

ERA5_PATH        = "/home/anaga-ambady/Documents/DISSERT/DATA/climate/individual_files/*.nc"
PARTHENIUM_FILE  = "/home/anaga-ambady/Documents/DISSERT/DATA/pakistan_clean_with_datetime.csv"
OUTPUT_FILE      = "/results/parthenium_with_climate_matched.csv"

# Climate parameters for Parthenium
PARTHENIUM_TBASE  = 7.2   # Base temperature for GDD (°C)
PARTHENIUM_TCEIL  = 42.8  # Ceiling temperature (cap heat accumulation)
RAIN_THRESHOLD_MM = 2.0   # Minimum daily rain to count as "rainy day"
LOOKBACK_DAYS     = 60    # Window for rolling climate calculations

CLIMATE_COLS = [
    'temp_mean_7d', 'precip_sum_7d', 'precip_freq_7d', 'gdd_7d',
    'temp_mean_30d', 'precip_sum_30d', 'precip_freq_30d', 'gdd_30d',
]


# ==============================================================================
# SECTION 1: LOAD ERA5 DATA
# ==============================================================================

def load_era5(file_pattern):
    """
    Load all ERA5 files and prepare dataset.
    Handles both 'valid_time' and 'time' coordinate names.
    Sorts latitude to prevent nearest-neighbor bugs.
    """
    files = sorted(glob.glob(file_pattern))
    # Try each supported engine to open the dataset (prefer netcdf4, fallback to h5netcdf)
    print(f"\n✓ Found {len(files)} ERA5 files")

    if len(files) == 0:
        raise FileNotFoundError(f"No files matching: {file_pattern}")

    ds = None
    for engine in ('netcdf4', 'h5netcdf'):
        try:
            ds = xr.open_mfdataset(files, engine=engine, combine='by_coords')
            print(f"✓ Loaded with engine: {engine}")
            break
        except Exception as e:
            print(f"  {engine} failed: {e}")

    if ds is None:
        raise RuntimeError("Could not load ERA5 files with any engine")

    if 'valid_time' in ds.dims:
        ds = ds.rename({'valid_time': 'time'})
        print("✓ Renamed 'valid_time' → 'time'")

    if 'expver' in ds.coords:
        ds = ds.drop_vars('expver')
        print("✓ Dropped 'expver' coordinate")

    if not ds.indexes["latitude"].is_monotonic_increasing:
        ds = ds.sortby("latitude")
        print("✓ Sorted latitude (was not monotonic)")

    if 't2m' in ds.data_vars:
        ds['t2m'] = ds['t2m'] - 273.15
        print("✓ Converted temperature K → °C")

    print(f"✓ ERA5 time range: {ds.time.values[0]} → {ds.time.values[-1]}")
    print(f"✓ ERA5 shape: {ds.dims}")

    return ds


# ==============================================================================
# SECTION 2: LOAD OBSERVATIONS
# ==============================================================================

def load_observations(filepath):
    """Load field observations and prepare time columns."""
    df = pd.read_csv(filepath)

    # Ensure a usable datetime column exists and normalize it for date-based joins

    if 'datetime' not in df.columns:
        if 'date' in df.columns:
            df['datetime'] = pd.to_datetime(df['date'])
        else:
            raise ValueError("CSV must contain 'datetime' or 'date' column")

    df['datetime'] = pd.to_datetime(df['datetime'])
    df = df.sort_values('datetime').reset_index(drop=True)
    df['obs_date'] = df['datetime'].dt.normalize()

    print(f"\n✓ Loaded {len(df)} observations")
    print(f"✓ Unique observation dates: {df['obs_date'].nunique()}")
    print(f"✓ Parthenium presence: {df['presence'].sum()}/{len(df)} ({100*df['presence'].mean():.1f}%)")

    return df


# ==============================================================================
# SECTION 3: COVERAGE DIAGNOSTICS
# ==============================================================================

def coverage_check(ds, obs_df):
    """Diagnose temporal coverage and warn about at-risk observations."""
    print(f"\n{'='*70}")
    print("COVERAGE DIAGNOSTICS")
    print(f"{'='*70}")

    # Compare ERA5 time coverage with observation dates
    era5_start = pd.Timestamp(ds.time.values[0])
    era5_end   = pd.Timestamp(ds.time.values[-1])
    obs_start  = obs_df["datetime"].min()
    obs_end    = obs_df["datetime"].max()

    print(f"ERA5 coverage: {era5_start.date()} → {era5_end.date()}")
    print(f"Obs coverage:  {obs_start.date()} → {obs_end.date()}")

    safe_start = era5_start + pd.Timedelta(days=LOOKBACK_DAYS)
    at_risk = (obs_df["datetime"] < safe_start).sum()

    if at_risk > 0:
        print(f"\n⚠️  WARNING: {at_risk} observations before {safe_start.date()}")
        print(f"   These will have INCOMPLETE rolling window values (<{LOOKBACK_DAYS}d)")
        print(f"   Recommendation: Download earlier ERA5 data or accept NaN for early obs")
    else:
        print(f"\n✓ OK: All {len(obs_df)} observations have full {LOOKBACK_DAYS}-day lookback coverage")

    print(f"{'='*70}\n")


# ==============================================================================
# SECTION 4: CLIMATE EXTRACTION (OPTIMIZED)
# ==============================================================================

def compute_climate(point_ds):
    """
    Compute all climate variables for one location's time window.

    OPTIMIZATIONS:
      - .clip() replaces xr.where() for temperature capping
      - .isel(time=-1) extracts the final rolling value in one shot
        (replaces the per-timestamp for-loop + .sel() that was the main bottleneck)
      - Caller is expected to pass an already-.load()'ed DataArray
        so no repeated lazy IO happens here

    Returns:
      dict of {climate_column: scalar_value}
    """
    # --- Daily aggregates ---
    daily_t  = point_ds['t2m'].resample(time='D').mean()
    daily_t_capped = daily_t.clip(max=PARTHENIUM_TCEIL)           

    # ERA5-Land 'tp' accumulates from 01:00 UTC to 00:00 UTC of the next day,
    # i.e. each hourly value is a running total since the bucket reset, and
    # the value AT 00:00 (start of the next calendar day) is the true total
    # for the day that just ended.
    #
    # xarray's resample(time='D') floors timestamps to midnight, so the
    # 00:00 reading gets grouped into the NEXT day's bin rather than the
    # day whose bucket it actually closes out. Shifting timestamps back
    # 1 hour before resampling corrects this label, so resample('D').max()
    # then attributes each day's true total to the correct calendar day.
    tp_shifted = point_ds['tp'].copy()
    tp_shifted['time'] = tp_shifted['time'] - pd.Timedelta(hours=1)
    daily_tp = tp_shifted.resample(time='D').max() * 1000
    

    daily_gdd   = (daily_t_capped - PARTHENIUM_TBASE).clip(min=0) 
    
    #--- Rainy day frequency: 1 if daily precipitation >= threshold, else 0
    rainy_day = (daily_tp >= RAIN_THRESHOLD_MM).astype(int)

    # --- Rolling series: grab last value only (.isel(time=-1)) ---
    def roll_mean(da, n): return da.rolling(time=n, min_periods=n).mean().isel(time=-1)
    def roll_sum(da, n):  return da.rolling(time=n, min_periods=n).sum().isel(time=-1)

    return {
        'temp_mean_7d':    float(roll_mean(daily_t,   7)),
        'precip_sum_7d':   float(roll_sum(daily_tp,   7)),
        'precip_freq_7d':  float(roll_sum(rainy_day,  7)),
        'gdd_7d':          float(roll_sum(daily_gdd,  7)),
        'temp_mean_30d':   float(roll_mean(daily_t,  30)),
        'precip_sum_30d':  float(roll_sum(daily_tp,  30)),
        'precip_freq_30d': float(roll_sum(rainy_day, 30)),
        'gdd_30d':         float(roll_sum(daily_gdd, 30)),
    }


# ==============================================================================
# SECTION 5: MAIN MATCHING LOOP (OPTIMIZED)
# ==============================================================================

def match_observations(ds, obs_df):
    """
    Match observations to climate data.

    OPTIMIZATIONS:
      - Outer loop groups by unique DATE (not unique point-date)
        → time window sliced once per date, not once per (date × location)
      - Inner loop iterates unique (lat, lon) within that date
      - .load() is called after spatial selection so only the needed
        location is pulled into memory (avoids repeated lazy IO)
    """
    # Initialize climate columns with NaN; they'll be filled per matched point
    for col in CLIMATE_COLS:
        obs_df[col] = np.nan

    unique_dates = sorted(obs_df['obs_date'].unique())
    n_dates = len(unique_dates)

    print(f"\n{'='*70}")
    print(f"CLIMATE MATCHING")
    print(f"  {len(obs_df)} observations | {n_dates} unique dates")
    print(f"{'='*70}\n")

    successful_matches = 0
    failed_matches     = 0
    report_every       = max(1, n_dates // 20)

    # Iterate over each unique observation DATE (slice time once per date)
    for i, obs_date in enumerate(unique_dates):
        obs_date = pd.Timestamp(obs_date)

        if (i + 1) % report_every == 0:
            pct = 100 * (i + 1) / n_dates
            print(f"  {i+1}/{n_dates} ({pct:.0f}%) | Matched: {successful_matches}")

        # --- OPTIMIZATION: slice time window ONCE per date ---
        window_start = obs_date - pd.Timedelta(days=LOOKBACK_DAYS)
        window_end = obs_date + pd.Timedelta(hours=23)
        try:
            time_window = ds.sel(time=slice(window_start, window_end))
        except Exception as e:
            print(f"  Time slice failed for {obs_date.date()}: {e}")
            date_mask = obs_df['obs_date'] == obs_date
            failed_matches += date_mask.sum()
            continue

        if len(time_window.time) == 0:
            date_mask = obs_df['obs_date'] == obs_date
            failed_matches += date_mask.sum()
            continue

        # --- Inner loop: unique locations within this date ---
        date_rows = obs_df[obs_df['obs_date'] == obs_date]

        for (lat, lon), _ in date_rows.groupby(['latitude', 'longitude']):
            try:
                # Spatial selection then .load() → pulls only this point into RAM
                point = time_window.sel(
                    latitude=lat,
                    longitude=lon,
                    method='nearest'
                ).load()

                climate_vals = compute_climate(point)

                mask = (
                    (obs_df['obs_date'] == obs_date) &
                    (obs_df['latitude'] == lat) &
                    (obs_df['longitude'] == lon)
                )
                n_matched = mask.sum()

                for col in CLIMATE_COLS:
                    obs_df.loc[mask, col] = climate_vals[col]

                successful_matches += n_matched

            except Exception as e:
                print(f"    Error at {obs_date.date()}, ({lat:.3f}, {lon:.3f}): {e}")
                mask = (
                    (obs_df['obs_date'] == obs_date) &
                    (obs_df['latitude'] == lat) &
                    (obs_df['longitude'] == lon)
                )
                failed_matches += mask.sum()

    print(f"\n✓ Matching complete:")
    print(f"  Successful: {successful_matches}/{len(obs_df)}")
    print(f"  Failed:     {failed_matches}/{len(obs_df)}")
    print(f"  Success rate: {100*successful_matches/len(obs_df):.1f}%\n")

    return obs_df


# ==============================================================================
# SECTION 6: FINALIZE & SAVE
# ==============================================================================

def finalize(obs_df, output_file):
    """Clean up and save results."""
    # Remove helper columns used for matching
    obs_df = obs_df.drop(columns=['obs_date'], errors='ignore')

    complete = obs_df[CLIMATE_COLS].notna().all(axis=1).sum()
    pct = 100 * complete / len(obs_df)

    print(f"{'='*70}")
    print("FINAL STATISTICS")
    print(f"{'='*70}")
    print(f"Complete cases (all variables filled): {complete}/{len(obs_df)} ({pct:.1f}%)")

    if pct < 70:
        print(f"⚠️  WARNING: <70% complete — check ERA5 date range coverage")
    else:
        print(f"✓ Data quality good for modeling\n")

    obs_df.to_csv(output_file, index=False)
    print(f"✓ Saved: {output_file}")
    print(f"  Rows: {len(obs_df)}")
    print(f"  Columns: {len(obs_df.columns)}")
    print(f"\nNew climate columns added:")
    for col in CLIMATE_COLS:
        print(f"  • {col}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ds = load_era5(ERA5_PATH)
    obs_df = load_observations(PARTHENIUM_FILE)
    coverage_check(ds, obs_df)
    obs_df = match_observations(ds, obs_df)
    finalize(obs_df, OUTPUT_FILE)


if __name__ == "__main__":
    main()
