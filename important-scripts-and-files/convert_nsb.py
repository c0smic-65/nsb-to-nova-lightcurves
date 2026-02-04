#!/usr/bin/env python3

import sys
import numpy as np
import os
from subprocess import Popen, PIPE

try:
    from StringIO import StringIO  # noqa: F401
except ImportError:
    from io import StringIO

if len(sys.argv) != 5:
    print("Usage: python convert_nsb.py <telescope> <run> <llim> <nova_name>")
    sys.exit(1)

itel = int(sys.argv[1])
nrun = int(sys.argv[2])
llim = int(sys.argv[3])
nova_name = sys.argv[4]

root_exec = "root"
macro_cmd = f"get_run_nsb.C+({itel},{nrun},{llim})"
command = [root_exec, "-l", "-b", "-q", macro_cmd]

process = Popen(command, stdout=PIPE, stderr=PIPE)
output, error = process.communicate()

# Print ROOT stderr (warnings, compile errors, etc.)
if error:
    sys.stderr.write(error.decode(errors="replace"))

# Bail out early if ROOT failed
if process.returncode != 0:
    sys.stderr.write(f"\nERROR: ROOT exited with code {process.returncode}\n")
    sys.exit(process.returncode)

text = output.decode(errors="replace")

# Robust parsing
if "[run]" not in text:
    sys.stderr.write(
        "\nERROR: Expected marker '[run]' not found in ROOT stdout.\n"
        "Macro likely did not run or output format changed.\n"
        "---- ROOT STDOUT (first 1200 chars) ----\n"
        f"{text[:1200]}\n"
        "--------------------------------------\n"
    )
    sys.exit(2)

after_run = text.split("[run]", 1)[1]
cleaned = after_run.split("Inside", 1)[0] if "Inside" in after_run else after_run

rec_arr = np.genfromtxt(
    StringIO(cleaned),
    delimiter=";",
    skip_header=1,
    skip_footer=1,
    dtype=[("Pixel", "i8"), ("NSB", "f8"), ("Time", "U26")]
)

base_dir = os.path.join(nova_name, "nsb")
tel_dir = os.path.join(base_dir, f"tel_{itel}")
os.makedirs(tel_dir, exist_ok=True)

csv_path = os.path.join(tel_dir, f"nsb_file_{itel}_{nrun}.csv")

np.savetxt(
    csv_path,
    rec_arr,
    fmt=["%d", "%.6f", "%s"],
    delimiter=";",
    comments=""
)

print(f"Wrote {csv_path}")
