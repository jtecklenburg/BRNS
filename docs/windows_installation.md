# Windows installation guide

This page gives a beginner-friendly setup guide for BRNS on Windows. The same overall workflow as on Linux applies, but the installation details are adjusted to the Windows environment.

The setup is:

- use Python in a virtual environment,
- install a Fortran compiler,
- clone the repository,
- install BRNS and `macrofor`,
- run the minimal example,
- inspect the results in the notebook.

---

## 1. Install Python

### Option A: Install Python with the official installer

Download Python from:

https://www.python.org/downloads/windows/

When installing, make sure you check the box:

- "Add Python to PATH"
- "Install launcher for all users"

Then open a new PowerShell window and verify:

```powershell
python --version
python -m pip --version
```

If both commands work, Python is ready.

### Option B: install Python with winget

```powershell
winget install --id Python.Python.3.12 -e
```

Then reopen the terminal and check:

```powershell
python --version
```

---

## 2. Prepare a Python environment

It is strongly recommended to install BRNS in a dedicated virtual environment.

Open PowerShell in the repository folder after cloning the project, or create the folder first:

```powershell
cd C:\path\to\BRNS
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

If PowerShell blocks the activation script, run:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

You should now see `(.venv)` in the terminal prompt.

---

## 3. Install a Fortran compiler

BRNS generates and compiles Fortran source files. You need a Fortran compiler such as `gfortran`.

### Recommended beginner-friendly option: install MSYS2

1. Download and install MSYS2 from:

   https://www.msys2.org/

2. Start the MSYS2 MINGW64 terminal and install the GNU Fortran toolchain,
    LAPACK, and Git:

```bash
pacman -S --needed \
   git \
   mingw-w64-x86_64-gcc-fortran \
   mingw-w64-x86_64-lapack
```

3. Verify the installation:

```bash
gfortran --version
```

This gives you a working Fortran compiler on Windows and is a straightforward setup for new users.

---

## 4. Install Git

If you do not already have Git, install it:

```powershell
winget install --id Git.Git -e
```

Verify:

```powershell
git --version
```

---

## 5. Clone the repository

```powershell
git clone https://github.com/jtecklenburg/BRNS.git
cd BRNS
```

If you already have the repository, just go to the folder and continue.

---

## 6. Install BRNS and `macrofor`

Make sure the virtual environment is active.

```powershell
.\.venv\Scripts\Activate.ps1
```

Install the Python dependency from the README:

```powershell
pip install "macrofor @ git+https://github.com/jtecklenburg/macrofor.git"
```

Then install BRNS itself:

```powershell
pip install -e .
```

This installs the project in editable mode so you can work with the local source from the repository.

---

## 7. Run the minimal example

The minimal example in the README uses the script `./build_python.sh`.

### Use Git Bash

Open Git Bash in the repository root and run:

```bash
chmod +x ./build_python.sh
./build_python.sh -c ./models/single_species/single_species_example.yaml -i ./models/single_species
```

### What this does

The script does the following:

1. reads the YAML model,
2. runs the Python-based ACG pipeline,
3. generates Fortran files,
4. prepares the build directory,
5. compiles the generated code,
6. runs the simulation and writes the result files.

If the example succeeds, the script prints the result directory. Typical output locations are under a build folder such as:

```bash
./build_output/single_species_example/results
```

The exact directory name depends on the model and output options.

---

## 8. Inspect the results in Jupyter

After the model run finishes, install Jupyter in the environment if needed:

```powershell
python -m pip install notebook jupyter nbconvert matplotlib
```

Then start Jupyter:

```powershell
jupyter notebook notebooks\plot_results.ipynb
```
Open the notebook and update the result-path cell if the output folder is different from the default path. To use another result folder without editing the notebook, set `BRNS_RESULT_DIR` before executing it:

```powershell
$env:BRNS_RESULT_DIR = Join-Path $PWD "build_output\single_species_example\results"
jupyter nbconvert `
   --to notebook `
   --execute `
   --ExecutePreprocessor.timeout=120 `
   --output executed_plot_results.ipynb `
   notebooks/plot_results.ipynb
```

This creates `notebooks/executed_plot_results.ipynb`, which contains the figures and cell output. The selected directory must contain the simulation `.dat` files.

---

## Troubleshooting

### `python` is not recognized

Close and reopen the terminal after installation, or check whether Python was added to PATH.

### PowerShell blocks activation

Run:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

### `gfortran` is not found

Install MSYS2 and then check that `gfortran --version` works from a new terminal session.

### `Permission denied` on `build_python.sh`

Run in Git Bash or a Linux shell:

```bash
chmod +x ./build_python.sh
```

### Notebook cannot find result files

Set `$env:BRNS_RESULT_DIR` to the exact generated result folder before opening or executing the notebook.
