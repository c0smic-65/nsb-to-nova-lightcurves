#!/usr/bin/env python3
import csv
import subprocess
import sys
from pathlib import Path

CSV_PATH = Path("hess_run_list_V572-Vel.csv")   # <-- updated
NOVA_NAME = "V572-Vel"                          # <-- updated
LLIM = 0
TEL_IDS = [1, 2, 3, 4]

def main():
    if not CSV_PATH.exists():
        print(f"ERROR: CSV not found: {CSV_PATH}", file=sys.stderr)
        sys.exit(1)

    runs = []
    with CSV_PATH.open(newline="") as f:
        reader = csv.DictReader(f)
        if "Run Number" not in reader.fieldnames:
            print(
                f"ERROR: Expected header 'Run Number'. Found: {reader.fieldnames}",
                file=sys.stderr,
            )
            sys.exit(1)

        for row in reader:
            val = (row.get("Run Number") or "").strip()
            if not val:
                continue
            runs.append(int(val))

    print(f"Loaded {len(runs)} runs from {CSV_PATH}")

    failed = []
    for run in runs:
        for tel in TEL_IDS:
            cmd = [
                "python3",
                "convert_nsb.py",
                str(tel),
                str(run),
                str(LLIM),
                NOVA_NAME,
            ]
            print("Running:", " ".join(cmd))
            res = subprocess.run(cmd)
            if res.returncode != 0:
                failed.append((tel, run, res.returncode))

    if failed:
        print("\nFAILED jobs:", file=sys.stderr)
        for tel, run, code in failed:
            print(f"  tel={tel} run={run} exit={code}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()
