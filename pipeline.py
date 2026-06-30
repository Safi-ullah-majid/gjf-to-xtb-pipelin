# =====================================================================
# GJF -> GFN2-xTB Geometry Optimization Pipeline (Google Colab)
# =====================================================================
# Usage in Colab:
#   1. Run the SETUP cell once (installs xtb).
#   2. Upload your .gjf file when prompted.
#   3. The script parses atoms/coordinates from the .gjf, runs a
#      GFN2-xTB optimization, and writes:
#        - <name>_optimized.xyz   (optimized geometry, XYZ format)
#        - <name>_optimized.gjf   (new Gaussian input with optimized coords)
#        - xtbopt.log / xtb output kept for inspection
# =====================================================================

# --------------------- CELL 1: SETUP (run once, no restart needed) ---------------------
# Uses micromamba instead of condacolab — no kernel restart required.
#
# Run these as separate %%bash cells, IN THIS ORDER:
#
#   %%bash
#   curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
#   mv bin/micromamba /usr/local/bin/
#
#   %%bash
#   micromamba create -y -n xtb -c conda-forge xtb
#
#   %%bash
#   micromamba run -n xtb xtb --version
#
# Confirm a version number prints before proceeding to Cell 2 below.
#
# IMPORTANT: %%bash shell activation does NOT carry over into Python's
# subprocess calls. This script therefore calls xtb via
# `micromamba run -n xtb xtb ...` from Python directly, so it works
# regardless of activation state in any given cell.
# ----------------------------------------------------------------------

import os
import re
import subprocess
from google.colab import files

PERIODIC_ELEMENTS = {
    "h","he","li","be","b","c","n","o","f","ne","na","mg","al","si","p","s",
    "cl","ar","k","ca","sc","ti","v","cr","mn","fe","co","ni","cu","zn","ga",
    "ge","as","se","br","kr","rb","sr","y","zr","nb","mo","tc","ru","rh","pd",
    "ag","cd","in","sn","sb","te","i","xe","cs","ba","la","ce","pr","nd","pm",
    "sm","eu","gd","tb","dy","ho","er","tm","yb","lu","hf","ta","w","re","os",
    "ir","pt","au","hg","tl","pb","bi","po","at","rn"
}

def parse_gjf(path):
    """
    Parse a Gaussian .gjf/.com file and extract:
      - charge, multiplicity
      - list of (element, x, y, z)
      - route section (for re-use in regenerated gjf header)
    """
    with open(path, "r") as f:
        lines = f.readlines()

    # Strip Gaussian link0/route comment lines, find blank-line separated blocks
    blocks = []
    current = []
    for line in lines:
        if line.strip() == "":
            blocks.append(current)
            current = []
        else:
            current.append(line.rstrip("\n"))
    if current:
        blocks.append(current)
    blocks = [b for b in blocks if b]  # drop empty blocks

    # Typical gjf structure (after blank-line splitting):
    # block0: %chk, # route lines
    # block1: title
    # block2: charge/mult + atom coordinates
    route_lines = []
    charge, mult = 0, 1
    atoms = []

    coord_block = None
    for b in blocks:
        # charge/multiplicity + coordinate block looks like:
        # "0 1" followed by "Element  x  y  z" lines
        if re.match(r"^\s*-?\d+\s+\d+\s*$", b[0]):
            coord_block = b
            break

    if coord_block is None:
        raise ValueError("Could not locate charge/multiplicity + coordinate block in .gjf file")

    charge, mult = coord_block[0].split()
    charge, mult = int(charge), int(mult)

    for line in coord_block[1:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        elem = parts[0]
        # Strip possible isotope/freeze flags like "C(Iso=13)" -> "C"
        elem_clean = re.sub(r"[^A-Za-z]", "", elem)
        if elem_clean.lower() not in PERIODIC_ELEMENTS:
            continue
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            continue
        atoms.append((elem_clean, x, y, z))

    # capture route section (first block, lines starting with #)
    for line in blocks[0]:
        if line.strip().startswith("#"):
            route_lines.append(line.strip())

    if not atoms:
        raise ValueError("No atoms parsed from .gjf file — check formatting.")

    return charge, mult, atoms, route_lines


def write_xyz(atoms, path, comment="Generated from GJF"):
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for elem, x, y, z in atoms:
            f.write(f"{elem:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")


def read_xyz(path):
    with open(path, "r") as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        elem = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append((elem, x, y, z))
    return atoms


def write_gjf(atoms, charge, mult, path, route="#p opt", title="Optimized structure"):
    with open(path, "w") as f:
        f.write(f"%chk={os.path.splitext(os.path.basename(path))[0]}.chk\n")
        f.write(f"{route}\n\n")
        f.write(f"{title}\n\n")
        f.write(f"{charge} {mult}\n")
        for elem, x, y, z in atoms:
            f.write(f"{elem:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")
        f.write("\n")


MICROMAMBA = "/usr/local/bin/micromamba"
XTB_ENV = "xtb"

def xtb_cmd(*args):
    """Build a command that runs xtb inside the 'xtb' micromamba env,
    regardless of whether the current shell has it activated."""
    return [MICROMAMBA, "run", "-n", XTB_ENV, "xtb", *args]


def run_xtb_optimization(xyz_path, charge=0, uhf=0, workdir="."):
    """
    Runs: micromamba run -n xtb xtb molecule.xyz --opt --chrg <charge> --uhf <uhf>
    uhf = multiplicity - 1 (number of unpaired electrons)
    Produces xtbopt.xyz in workdir.
    """
    cmd = xtb_cmd(
        os.path.abspath(xyz_path),
        "--opt",
        "--chrg", str(charge),
        "--uhf", str(uhf),
    )
    result = subprocess.run(
        cmd, cwd=workdir, capture_output=True, text=True
    )
    log_path = os.path.join(workdir, "xtb_run.log")
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write("\n----- STDERR -----\n")
        f.write(result.stderr)

    optxyz = os.path.join(workdir, "xtbopt.xyz")
    if result.returncode != 0 or not os.path.exists(optxyz):
        print("xTB STDOUT (tail):\n", result.stdout[-3000:])
        print("xTB STDERR (tail):\n", result.stderr[-2000:])
        raise RuntimeError("xTB optimization failed — see xtb_run.log for details.")

    return optxyz, log_path


def check_xtb_installed():
    from shutil import which

    if which("micromamba") is None and not os.path.exists(MICROMAMBA):
        raise EnvironmentError(
            "micromamba was not found.\n\n"
            "Run the SETUP cell first (separate %%bash cells, in order):\n\n"
            "  %%bash\n"
            "  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba\n"
            "  mv bin/micromamba /usr/local/bin/\n\n"
            "  %%bash\n"
            "  micromamba create -y -n xtb -c conda-forge xtb\n\n"
            "  %%bash\n"
            "  micromamba run -n xtb xtb --version\n\n"
            "Only run this pipeline once `xtb --version` prints a version number."
        )

    result = subprocess.run(xtb_cmd("--version"), capture_output=True, text=True)
    if result.returncode != 0 or "xtb version" not in (result.stdout + result.stderr).lower():
        raise EnvironmentError(
            "xtb env exists but failed to run via micromamba.\n"
            f"STDOUT: {result.stdout[-1000:]}\n"
            f"STDERR: {result.stderr[-1000:]}\n\n"
            "Make sure you ran:\n"
            "  %%bash\n"
            "  micromamba create -y -n xtb -c conda-forge xtb\n"
            "and that it completed without errors."
        )
    print("xtb found:", result.stdout.strip().splitlines()[0] if result.stdout.strip() else result.stderr.strip().splitlines()[0])


def pipeline(gjf_path, workdir="xtb_workdir"):
    check_xtb_installed()
    os.makedirs(workdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(gjf_path))[0]

    print(f"[1/4] Parsing {gjf_path} ...")
    charge, mult, atoms, route_lines = parse_gjf(gjf_path)
    uhf = mult - 1
    print(f"   Charge={charge}, Multiplicity={mult}, Atoms={len(atoms)}")

    print("[2/4] Writing starting .xyz ...")
    start_xyz = os.path.join(workdir, f"{base}_start.xyz")
    write_xyz(atoms, start_xyz, comment=f"From {gjf_path}")

    print("[3/4] Running GFN2-xTB optimization ...")
    optxyz, log_path = run_xtb_optimization(start_xyz, charge=charge, uhf=uhf, workdir=workdir)
    print(f"   Done. Log: {log_path}")

    print("[4/4] Writing optimized outputs ...")
    opt_atoms = read_xyz(optxyz)

    final_xyz = f"{base}_optimized.xyz"
    final_gjf = f"{base}_optimized.gjf"
    write_xyz(opt_atoms, final_xyz, comment=f"GFN2-xTB optimized geometry of {base}")
    write_gjf(
        opt_atoms, charge, mult, final_gjf,
        route=(route_lines[0] if route_lines else "#p opt"),
        title=f"GFN2-xTB optimized: {base}",
    )

    print(f"\nFinished. Outputs:\n  {final_xyz}\n  {final_gjf}")
    return final_xyz, final_gjf


# --------------------- CELL 2: RUN (upload + execute) ---------------------
if __name__ == "__main__":
    print("Upload your .gjf file:")
    uploaded = files.upload()
    gjf_file = list(uploaded.keys())[0]

    final_xyz, final_gjf = pipeline(gjf_file)

    files.download(final_xyz)
    files.download(final_gjf)
