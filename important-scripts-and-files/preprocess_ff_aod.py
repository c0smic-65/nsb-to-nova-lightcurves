#!/usr/bin/env python3
"""
Preprocess H.E.S.S. NSB files:
- per telescope (TEL_ID from CLI)
- apply flatfield + AOD correction
- save corrected per-pixel CSVs: Time;Pixel;NSB_corr

Usage:
    python3 preprocess_ff_aod.py --tel 1
    python3 preprocess_ff_aod.py -t 2
"""
import sys
sys.path.append('/lfs/l7/hess/users/sghosh/projects/Satellite-trails')
import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path
import datetime as _dt

# ===== NOVA-SPECIFIC CONFIG =====
NOVA_NAME = "MGAB-V207"   # purely documentary, but useful
BASE_DIR = Path(
    "/lfs/l7/hess/users/sghosh/projects/Satellite-trails/MGAB-V207/nsb"
)
OUTPUT_DIR = Path(
    "/lfs/l7/hess/users/sghosh/projects/Satellite-trails/MGAB-V207/nsb_corrected"
)
AERONET_FILE = Path(
    "/lfs/l7/hess/users/sghosh/projects/Satellite-trails/MGAB-V207/important-scripts-and-files/20200717_20200729_HESS.lev20"
)
AOD_WAVELENGTH     = 380               # target nm
USE_CLOSEST_AVAIL  = True              # if exact nm empty, use nearest AOD_*nm with >0 data
MATCH_STRAT        = "bracket_mean"    # "nearest" | "daily" | "bracket_mean"
NEAR_TOL           = "6H"              # for "nearest" / "daily"
BIN_INTERVAL       = None              # keep native resolution

# ===== FLAT-FIELD (FF) CONFIG =====
# We will use Gerrit's ff_corr.npy (3, 960, 2) instead of the HESS DB flat-field.
#
# Two modes:
#   1) "gerrit_scale_like_hess": treat a[pix]=corr[pix,0] like Coefficient_hg,
#      build factors = a / median(a), and MULTIPLY raw NSB by factors (same *formula style* as before).
#   2) "gerrit_affine": apply Gerrit's actual affine correction per pixel:
#         NSB_corr = (NSB_raw - b) / a, where a=corr[:,0], b=corr[:,1].
#
# Pick one below.
FF_MODE = "gerrit_affine"  # "gerrit_affine" | "gerrit_scale_like_hess"

# Path to Gerrit's correction file
GERRIT_FF_CORR_FILE = Path(
    "/lfs/l7/hess/users/sghosh/projects/Satellite-trails/ff_corr.npy"
)

# Run-number boundaries used to select which of the 3 correction blocks to use
GERRIT_RUN_BOUNDARIES = (163380, 166200)  # (<b0) -> idx0, [b0,b1) -> idx1, (>=b1) -> idx2

# ── Zenith angles per run (same for all telescopes) ──
NOVA_CONFIG_ZENITH = {
    161096: 53.0,
    161097: 48.3,
    161098: 44.9,
    161131: 53.9,
    161132: 49.2,
    161133: 45.1,
    161163: 54.7,
    161164: 49.9,
    161165: 45.3,
    161195: 58.0,
    161196: 54.2,
    161197: 49.5,
    161198: 45.3,
    161224: 57.2,
    161225: 53.5,
    161226: 48.8,
    161227: 44.2,
    161268: 57.3,
    161269: 52.6,
    161270: 48.1,
    161271: 44.5,
    161306: 52.4,
    161307: 47.7,
    161308: 43.2,
    161334: 51.8,
    161335: 48.1,
    161336: 43.1,
    161366: 51.2,
    161367: 46.5,
    161368: 43.0,
    161426: 51.9,
    161427: 48.2,
    161428: 43.7,
    161476: 51.5,
    161477: 48.1,
    161478: 42.9,
}


# ──────────────── COLORS / PRETTY OUTPUT ────────────────
RESET = "\033[0m"
COLORS = {
    "blue": "\033[94m",
    "green": "\033[92m",
    "yellow": "\033[93m",
    "red": "\033[91m",
    "cyan": "\033[96m",
}

def c(text: str, color: str) -> str:
    """Color helper (no external deps)."""
    return COLORS.get(color, "") + str(text) + RESET


# ───────────────── HELPERS ──────────────────────────────
def airmass_kasten_young(z_deg: float) -> float:
    """Kasten–Young airmass formula."""
    z = np.clip(float(z_deg), 0.0, 89.999)
    return 1.0 / (np.cos(np.deg2rad(z)) + 0.50572 * (96.07995 - z)**-1.6364)


def detect_aeronet_header(filepath: str) -> int:
    """Return the line index where the AERONET data header (Date()...) starts."""
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            if line.strip().lower().startswith("date("):
                return i
    raise RuntimeError("AERONET header not found")


def available_aod_columns(df: pd.DataFrame):
    """Return list of (nm:int, col:str) for AOD_*nm columns present."""
    out = []
    for ccol in df.columns:
        m = re.fullmatch(r"AOD_(\d+)nm", str(ccol))
        if m:
            out.append((int(m.group(1)), ccol))
    return sorted(out)


def choose_best_aod_column(df: pd.DataFrame, target_nm: int):
    """
    Choose the closest AOD_*nm column that has (>0) data.
    Returns (chosen_nm:int, chosen_col:str).
    """
    candidates = []
    for nm, col in available_aod_columns(df):
        s = pd.to_numeric(df[col], errors="coerce")
        pos = (s > 0).sum()
        if pos > 0:
            candidates.append((abs(nm - int(target_nm)), nm, col, pos))
    if not candidates:
        raise RuntimeError("No AOD_*nm columns with positive data found in AERONET file.")
    candidates.sort()
    _, nm, col, pos = candidates[0]
    print("  " + c(f"[AOD] Using column {col} ({nm} nm), rows>0={pos}", "cyan"))
    return nm, col


def _daily_bracket_mean(aero_df: pd.DataFrame, aod_col: str) -> pd.Series:
    """
    For each calendar day D, compute mean(last AOD on D, first AOD on D+1).
    Returns a Series indexed by date (datetime.date).
    """
    if not isinstance(aero_df.index, pd.DatetimeIndex):
        raise ValueError("aero_df must be indexed by DatetimeIndex")

    aero = aero_df.sort_index()
    by_day = dict(tuple(aero.groupby(aero.index.date)))
    all_days = sorted(by_day.keys())

    out = {}
    for d in all_days:
        v_last = None
        v_next_first = None

        g = by_day.get(d, None)
        if g is not None and not g.empty:
            vals = pd.to_numeric(g[aod_col], errors="coerce").dropna()
            v_last = None if vals.empty else float(vals.iloc[-1])

        d_next = d + _dt.timedelta(days=1)
        g2 = by_day.get(d_next, None)
        if g2 is not None and not g2.empty:
            vals2 = pd.to_numeric(g2[aod_col], errors="coerce").dropna()
            v_next_first = None if vals2.empty else float(vals2.iloc[0])

        if v_last is None and v_next_first is None:
            continue
        if v_last is None:
            out[d] = v_next_first
        elif v_next_first is None:
            out[d] = v_last
        else:
            out[d] = 0.5 * (v_last + v_next_first)

    return pd.Series(out)


def load_and_prepare_aeronet(aeronet_file: str, target_nm: int, use_closest: bool = True):
    """Load AERONET file once, pick appropriate AOD column, return (aero_df, aod_col)."""
    print(c("[AERONET] Loading AERONET file...", "blue"))
    skip = detect_aeronet_header(aeronet_file)
    aero = pd.read_csv(aeronet_file, skiprows=skip, na_values=[-999, -999.0, -9999])

    date_col = next(cc for cc in aero.columns if str(cc).startswith("Date("))
    time_col = next(cc for cc in aero.columns if str(cc).startswith("Time("))
    date_str = aero[date_col].astype(str).str.replace(":", "-", regex=False)
    time_str = aero[time_col].astype(str)
    dt = (pd.to_datetime(date_str + " " + time_str, dayfirst=True, utc=True, errors="coerce")
            .dt.tz_convert(None))

    if use_closest:
        chosen_nm, aod_col = choose_best_aod_column(aero, target_nm)
        print("  " + c(f"[AOD] Target {target_nm} nm → using closest at {chosen_nm} nm", "cyan"))
    else:
        needle = f"AOD_{int(target_nm)}"
        match = [cc for cc in aero.columns if needle in str(cc)]
        if not match:
            raise KeyError(f"No exact AOD column found for {target_nm} nm")
        aod_col = match[0]
        print("  " + c(f"[AOD] Using exact column {aod_col}", "cyan"))

    aero[aod_col] = pd.to_numeric(aero[aod_col], errors="coerce")
    aero = pd.DataFrame({"Datetime": dt, aod_col: aero[aod_col]})
    aero = aero.dropna(subset=["Datetime", aod_col]).set_index("Datetime").sort_index()
    print("  " + c(f"[AERONET] Time range: {aero.index.min()} → {aero.index.max()} | points: {len(aero)}", "blue"))
    return aero, aod_col


def aod_for_times(time_index: pd.DatetimeIndex,
                  aero: pd.DataFrame,
                  aod_col: str,
                  match_strategy: str,
                  nearest_tolerance: str):
    """
    Return an AOD Series indexed like time_index using the chosen strategy.
    """
    if time_index.empty:
        return pd.Series([], dtype=float, index=time_index)

    tmp_hess = pd.DataFrame({"Time": time_index}).sort_values("Time")
    tmp_aod  = aero[[aod_col]].reset_index().rename(columns={"Datetime": "AOD_Time"}).sort_values("AOD_Time")

    if match_strategy == "nearest":
        merged = pd.merge_asof(
            tmp_hess, tmp_aod,
            left_on="Time", right_on="AOD_Time",
            direction="nearest",
            tolerance=pd.Timedelta(nearest_tolerance)
        )
        used = merged[aod_col]

    elif match_strategy == "daily":
        aod_daily = aero[[aod_col]].resample("1D").median().dropna()
        tmp_aod = aod_daily.reset_index().rename(columns={"Datetime": "AOD_Time"})
        merged = pd.merge_asof(
            tmp_hess, tmp_aod,
            left_on="Time", right_on="AOD_Time",
            direction="nearest",
            tolerance=pd.Timedelta(nearest_tolerance)
        )
        used = merged[aod_col]

    elif match_strategy == "bracket_mean":
        # compute daily bracket means: mean(last of D, first of D+1)
        aod_daily = _daily_bracket_mean(aero[[aod_col]], aod_col)  # index: date
        aod_daily_df = aod_daily.rename("AOD_bmean").to_frame()
        aod_daily_df["HESS_Date"] = aod_daily_df.index

        tmp_hess["HESS_Date"] = tmp_hess["Time"].dt.date
        merged = tmp_hess.merge(aod_daily_df[["HESS_Date","AOD_bmean"]],
                                on="HESS_Date", how="left")
        used = merged["AOD_bmean"]
        merged = merged.drop(columns=["HESS_Date"])

    else:
        raise ValueError(f"Unknown MATCH_STRAT='{match_strategy}'. Use 'nearest' | 'daily' | 'bracket_mean'.")

    used.index = time_index        # align back to the original index
    return used.astype(float)


"""Flat-field correction source: Gerrit's ff_corr.npy"""

def _load_gerrit_ff_corr(ff_path: Path) -> np.ndarray:
    arr = np.load(str(ff_path))
    if not (isinstance(arr, np.ndarray) and arr.ndim == 3 and arr.shape[1] == 960 and arr.shape[2] == 2):
        raise ValueError(
            f"Unexpected shape for {ff_path}. Expected (3, 960, 2), got {getattr(arr, 'shape', None)}"
        )
    return arr


_GERRIT_ARR = _load_gerrit_ff_corr(GERRIT_FF_CORR_FILE)


def gerrit_corr_for_run(run_id: int) -> tuple:
    """Return (a_scale, b_offset) as pandas Series indexed by Pixel for the given run."""
    b0, b1 = GERRIT_RUN_BOUNDARIES
    if run_id < b0:
        idx = 0
    elif run_id >= b1:
        idx = 2
    else:
        idx = 1

    corr = _GERRIT_ARR[idx]  # shape (960, 2)
    pix_idx = pd.Index(np.arange(corr.shape[0]), name="Pixel")
    a = pd.Series(corr[:, 0], index=pix_idx, dtype=float)  # scale
    b = pd.Series(corr[:, 1], index=pix_idx, dtype=float)  # offset
    return a, b


# ───────────────── MAIN PROCESSOR ──────────────────────
def process_telescope(tel_id: int):
    base = BASE_DIR / f"tel_{tel_id}"
    out_dir = OUTPUT_DIR / f"tel_{tel_id}"

    if not base.exists():
        print(c(f"[ERROR] Base directory does not exist: {base}", "red"))
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    pattern = f"nsb_file_{tel_id}_*.csv"
    files = sorted(base.glob(pattern)) or sorted(base.rglob(pattern))
    if not files:
        print(c(f"[WARN] No files found under {base} matching {pattern}", "yellow"))
        return

    run_ids = sorted({
        int(re.search(rf"nsb_file_{tel_id}_(\d+)\.csv$", p.name).group(1))
        for p in files if re.search(rf"nsb_file_{tel_id}_(\d+)\.csv$", p.name)
    })

    print()
    print(c(f"=== Preprocessing NSB files for CT{tel_id} ===", "green"))
    print(c(f"Base dir: {base}", "blue"))
    print(c(f"Out dir : {out_dir}", "blue"))
    print(c(f"Runs    : {run_ids}", "blue"))

    # Load AERONET once
    aero_df, aod_col = load_and_prepare_aeronet(AERONET_FILE, AOD_WAVELENGTH, USE_CLOSEST_AVAIL)

    n_files = len(files)
    n_ok = 0

    for i, f in enumerate(files, start=1):
        # progress header
        pct = 100.0 * i / n_files
        print()
        print(c(f"[CT{tel_id}] [{i}/{n_files}] ({pct:5.1f}%)", "green"), c(f"{f.name}", "cyan"))

        m = re.search(rf"nsb_file_{tel_id}_(\d+)\.csv$", f.name)
        if not m:
            print("  " + c("[skip] Could not parse run_id from filename", "yellow"))
            continue
        run_id = int(m.group(1))

        # load raw NSB
        df = pd.read_csv(f, sep=";", header=None, names=["Pixel","NSB","Time"])
        df["Time"] = (pd.to_datetime(df["Time"].str.replace("UTC: ", ""), utc=True, errors="coerce")
                        .dt.tz_convert(None))
        df = df.dropna(subset=["Time"]).sort_values("Time")

        if df.empty:
            print("  " + c("[warn] No valid times in file, skipping.", "yellow"))
            continue

        # guard duplicates and pivot to wide: Time x Pixel
        df = df.groupby(["Time","Pixel"], as_index=False)["NSB"].mean()
        df_wide = df.pivot(index="Time", columns="Pixel", values="NSB")

        # No binning here (BIN_INTERVAL=None), but keep code flexible
        if BIN_INTERVAL:
            df_wide = df_wide.resample(BIN_INTERVAL).mean()

        # ── Flat-field (Gerrit ff_corr.npy) ──
        # corr[:,0] = a (scale) ; corr[:,1] = b (offset)
        a_scale_all, b_off_all = gerrit_corr_for_run(run_id)
        denom = float(np.nanmedian(a_scale_all.values))

        # align to whatever pixels exist in this file
        cols = df_wide.columns
        a_scale = a_scale_all.reindex(cols).astype(float)
        b_off   = b_off_all.reindex(cols).astype(float)
        a_scale = a_scale.fillna(1.0)
        b_off   = b_off.fillna(0.0)
        a_scale = a_scale.replace(0.0, np.nan)

        if FF_MODE == "gerrit_scale_like_hess":
            # Same "style" as your old HESS FF block: factors = a/median(a), then multiply.
            factors = (a_scale / denom).fillna(1.0)
            df_flat = df_wide.multiply(factors, axis=1)
            print("  " + c(f"[FF] Gerrit scale-only (HESS-like): multiplied by a/median(a) (denom={denom:.6g}).", "blue"))

        elif FF_MODE == "gerrit_affine":
            # Gerrit's actual correction: (img - b) / a
            df_flat = df_wide.subtract(b_off, axis=1).divide(a_scale, axis=1)
            print("  " + c("[FF] Gerrit affine: (NSB_raw - b) / a", "blue"))

        else:
            raise ValueError(f"Unknown FF_MODE='{FF_MODE}'. Use 'gerrit_affine' | 'gerrit_scale_like_hess'.")

        print("  " + c(f"[FF] Applied to {df_flat.shape[1]} pixels.", "blue"))

        # AOD for each time stamp
        if run_id not in NOVA_CONFIG_ZENITH:
            print("  " + c(f"[ERROR] No zenith available for run {run_id}. Add it to NOVA_CONFIG_ZENITH.", "red"))
            continue

        zenith_deg = float(NOVA_CONFIG_ZENITH[run_id])
        m_airmass  = airmass_kasten_young(zenith_deg)
        print("  " + c(f"[zenith] {zenith_deg:.2f} deg → airmass = {m_airmass:.3f}", "blue"))

        aod_series = aod_for_times(df_flat.index, aero_df, aod_col,
                                   MATCH_STRAT, NEAR_TOL)
        matched_frac = aod_series.notna().mean()
        print("  " + c(f"[AOD] matched {100*matched_frac:.1f}% of time bins.", "cyan"))

        # build correction factor; drop times without AOD
        corr_factor = np.exp(aod_series * m_airmass)
        valid_mask = corr_factor.notna()

        df_flat = df_flat.loc[valid_mask]
        corr_factor = corr_factor.loc[valid_mask]

        if df_flat.empty:
            print("  " + c("[warn] No rows after AOD matching, skipping.", "yellow"))
            continue

        df_corr = df_flat.multiply(corr_factor, axis=0)

        # Save as long-format CSV: Time;Pixel;NSB_corr
        df_long = df_corr.reset_index().melt(id_vars="Time",
                                             var_name="Pixel",
                                             value_name="NSB_corr")
        df_long = df_long.dropna(subset=["NSB_corr"])

        out_file = out_dir / f.name.replace(".csv", "_ff_aod.csv")
        df_long.to_csv(out_file, sep=";", index=False)
        print("  " + c(f"[ok] Saved corrected file to {out_file}", "green"))
        n_ok += 1

    print()
    print(c(f"Done for CT{tel_id}. Wrote {n_ok} corrected files to {out_dir}", "green"))


# ───────────────── CLI ENTRYPOINT ──────────────────────
def parse_args():
    parser = argparse.ArgumentParser(
        description="Flatfield + AOD preprocess H.E.S.S. NSB files per telescope."
    )
    parser.add_argument(
        "-t", "--tel",
        type=int,
        required=True,
        help="Telescope ID (e.g. 1, 2, 3, 4)."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    process_telescope(args.tel)
