# Linux and WSL installation guide and first steps

This page is a beginner-friendly installation guide for BRNS on Linux and on Windows Subsystem for Linux (WSL). The goal is to get from a fresh system to a working BRNS example as quickly and clearly as possible.

The workflow is:

1. install Python if needed,
2. prepare a Python environment,
3. install the Fortran compiler,
4. install BRNS and its Python dependency `macrofor`,
5. run the minimal demonstration model,
6. inspect the generated results in the Jupyter notebook.

> BRNS is designed for a Linux-like environment. If you run Windows, WSL is the easiest way to use it reliably.

---

## 1. Install Python (if needed)

### On Ubuntu / Debian / WSL

Open a terminal and run:

```bash
sudo apt update
sudo apt install -y python3 python3-pip python3-venv git build-essential
```

This installs the basic software needed for Python and compilation.

### Check that Python is available

```bash
python3 --version
python3 -m pip --version
```

If these commands work, Python is ready.

### Optional: create a dedicated virtual environment

It is strongly recommended to install BRNS in its own virtual environment, so your Python packages stay isolated.

From the repository root, create a project environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

You will use this environment for all BRNS installation and execution steps.

---

## 2. Install the Fortran compiler

BRNS uses a generated Fortran code workflow, so a Fortran compiler is required.

### Install gfortran

```bash
sudo apt install -y gfortran
```

### Verify the installation

```bash
gfortran --version
```

If the version information is printed, the compiler is installed correctly.

---

## 3. Clone the BRNS repository

```bash
git clone https://github.com/jtecklenburg/BRNS.git
cd BRNS
```

If you are already in the project folder, just continue there.

---

## 4. Install BRNS in the Python environment

BRNS requires the Python package `macrofor` as well as the project itself.

### Activate the environment

```bash
source .venv/bin/activate
```

### Install `macrofor`

This is the dependency from the project README:

```bash
pip install "macrofor @ git+https://github.com/jtecklenburg/macrofor.git"
```

### Install BRNS itself

From the repository root:

```bash
pip install -e .
```

This installs the project in editable mode, which is useful during development and testing.

---

## 5. Run the minimal example

The README includes a small example workflow that uses the generated Python pipeline and the build script.

### Check the example files

The repository includes the example model under the project folder:

```bash
ls ./models/single_species
```

You should see files such as:

- `single_species_example.yaml`
- `diss_a.inp`

### Run the minimal model

Use the build script from the repository root:

```bash
./build_python.sh -c ./models/single_species/single_species_example.yaml -i ./models/single_species
```

### What this does

The script performs the following:

1. reads the YAML model,
2. calls the Python-based code generator,
3. generates Fortran source files,
4. prepares the build directory,
5. compiles the generated code,
6. runs the model and writes result files.

If the command succeeds, it will print the output directory and the result path. The generated model results are typically placed under a build folder such as:

```bash
./build_output/single_species_example/results
```

The exact location depends on your output folder and model name, but the script prints the path in the terminal.

---

## 6. Inspect the result in Jupyter

After the model run finishes, open the plotting notebook from the repository:

```bash
jupyter notebook notebooks/plot_results.ipynb
```

If `jupyter` is not installed yet, install it with:

```bash
python -m pip install notebook jupyter
```

Then run again.

### Adjust the notebook to your result folder

Open [notebooks/plot_results.ipynb](notebooks/plot_results.ipynb) and change the example path in the first code cell if needed. The notebook expects a result directory with `.dat` files and then plots them by depth or time snapshot.