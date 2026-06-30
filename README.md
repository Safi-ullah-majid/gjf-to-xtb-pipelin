# GJF → GFN2-xTB Optimization Pipeline

A lightweight Google Colab pipeline that takes a Gaussian input file (`.gjf`) and produces a geometry-optimized structure using GFN2-xTB — fast, free, and good enough for most starting-geometry needs before higher-level DFT.

## What it does

1. Parses atoms, charge, and multiplicity from a `.gjf` file
2. Converts the geometry to `.xyz` format
3. Runs a GFN2-xTB optimization via [xtb](https://github.com/grimme-lab/xtb)
4. Writes back:
   - `<name>_optimized.xyz` — optimized geometry, viewable in Avogadro/VMD/PyMOL
   - `<name>_optimized.gjf` — ready to drop into Gaussian for follow-up DFT work

## Why GFN2-xTB

GFN2-xTB is a semi-empirical tight-binding method — much faster than DFT, and accurate enough for reliable bond lengths/angles on typical organic and organometallic structures. It's not a substitute for DFT when you need precise energies, reaction barriers, or non-covalent interaction accuracy — but it's an excellent first pass for getting a clean, reasonable starting geometry in seconds rather than minutes or hours.

## Setup (Google Colab)

xtb is installed via [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) (no kernel restart required, unlike condacolab).

Run these as separate cells, in order:

```bash
%%bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest -o micromamba.tar.bz2
mkdir -p micromamba_extract
tar -xvjf micromamba.tar.bz2 -C micromamba_extract
find micromamba_extract -name micromamba -exec cp {} /usr/local/bin/micromamba \;
chmod +x /usr/local/bin/micromamba
/usr/local/bin/micromamba --version
```

```bash
%%bash
micromamba create -y -n xtb -c conda-forge xtb
micromamba run -n xtb xtb --version
```

Once `xtb --version` prints a version number, you're ready to run the pipeline.

## Usage

```python
# paste the contents of pipeline.py into a Colab cell, then run it.
# It will prompt you to upload a .gjf file, run the optimization,
# and auto-download the optimized .xyz and .gjf files.
```

Or, if running locally instead of Colab (with `xtb` and `micromamba` available on your system):

```bash
python pipeline.py
```

## Notes

- Calls to `xtb` are routed through `micromamba run -n xtb xtb ...` directly from Python, so the script doesn't depend on shell activation state — it works the same whether run interactively or as a script.
- Multiplicity is converted to xtb's `--uhf` flag automatically (`uhf = multiplicity - 1`).
- Full xtb output is saved to `xtb_run.log` in the working directory for troubleshooting non-converging optimizations.

## License

MIT (or replace with your preferred license)
