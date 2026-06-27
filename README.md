# Sensor Modelling

Simulation and analysis code for evaluating chemosensor array designs, with a focus on indicator displacement assays (IDAs), mixed-host sensor arrays, and PLS-DA classification efficiency.

The project generates candidate formation-constant configurations, simulates spectral responses under chemical-equilibrium models, ranks sensor designs by classification efficiency, and produces correlation/heatmap visualizations for the resulting parameter spaces.

## Project Contents

### MATLAB simulation code

- `Data_D.m`, `data_S.m`, `values_mixed.mat`, `piers_values.mat` - input data and parameter values used by the simulations.
- `IDAa_f.m`, `IDAs_f.m`, `Ms_f.m`, `Ma_f.m` - sensor response models for IDA and mixed-host architectures.
- `NewtonRaphson.m`, `gauss.m` - numerical solvers used for equilibrium calculations.
- `spectra_IDA.m`, `spectra_M.m`, `piers_spectra_IDA.m`, `piers_spectra_m.m` - simulated spectral basis functions.
- `PLSDA_train.m`, `PLSDA_pred.m` - PLS-DA training and prediction utilities.
- `*_poolgen*_ffd.m` - parameter-pool generators for full-factorial or sampled formation-constant designs.
- `ffd_a.m`, `ffd_s.m`, `IDAa_ffdU.m`, `Ma_ffdU.m`, `Ms_ffdU.m` - long-running evaluation scripts that compute classification efficiency and save sorted results.
- `plot*.m`, `GUI_piers.m` - plotting and exploratory utilities.

### Python post-processing code

- `Ms_tables.py`, `IDAs_table.py` - generate correlation tables and heatmaps from CSV results.
- `maximafinder.py` - interpolates gridded response surfaces and finds local maxima.
- `convexchecker.py` - checks convexity or concavity of regular-grid response data.
- `filter.py` - small CSV filtering helper.

### Reports and outputs

- `IDA_report.tex`, `IDA_report.pdf` - report source and compiled report.
- `*_res*.csv`, `*_res*.mat`, `*_checkpoint.mat` - generated result and checkpoint files.
- `*_outputs*`, `sorted_plots`, `Et_comp` - generated plots and comparison outputs.

Many result files are large and are produced by long-running simulations. Treat the MATLAB and Python scripts as the main reproducible source files; treat the `.mat`, `.csv`, `.png`, `.rar`, and log files as generated artifacts unless you intentionally need to preserve a specific run.

## Requirements

### MATLAB

The simulation workflow expects MATLAB with:

- Parallel Computing Toolbox, for `parpool`, `parfor`, and chunked long-running evaluations.
- Statistics and Machine Learning Toolbox, for functions such as `lhsdesign`.

Large runs can require substantial memory. Several scripts use checkpointing and `matfile` to reduce memory pressure, but the full datasets may still be heavy.

### Python

Python post-processing uses:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`

Install them with:

```bash
pip install numpy pandas matplotlib scipy
```

## Typical Workflow

### 1. Generate parameter pools

Run the appropriate pool generator from MATLAB. For example:

```matlab
A = Ma_poolgenU_ffd(meanHD, deltaD1, deltaK1, deltaC1, deltaG1, ...
                    deltaHD, deltaD2, deltaK2, deltaC2, deltaG2, 5000);
save('Ma_cons9U.mat', 'A', '-v7.3');
```

The generated matrix usually contains formation constants followed by the input parameters used to construct them.

### 2. Evaluate sensor designs

Run the matching evaluation script:

```matlab
ffd_a
```

or:

```matlab
ffd_s
```

These scripts process parameter rows in chunks, train/test PLS-DA classifiers, save checkpoints, and produce sorted result matrices. If a run stops early, rerunning the script will resume from the checkpoint when compatible.

### 3. Export or analyze results

Use the generated `.mat` or `.csv` files for analysis. The common result columns are:

- formation constants, such as `K1`, `K2`, ...
- derived spacing parameters, such as `deltaD`, `deltaK`, `deltaC`, `deltaG`
- total classification efficiency, usually named `Et`
- per-class efficiencies, usually named `Class1` and `Class2`

### 4. Generate figures

Run the Python table/plot scripts:

```bash
python Ms_tables.py
python IDAs_table.py
```

These scripts read CSV result files and write correlation tables, mean-efficiency heatmaps, and stratified-correlation plots into the corresponding output folders.

### 5. Find local maxima

Use `maximafinder.py` for response-surface interpolation and local maximum detection:

```bash
python maximafinder.py
```

Edit `csv_path`, `knob_cols`, and `value_col` in the script to target a different result file or response column.

## Naming Notes

The project uses short prefixes for sensor families and analysis variants:

- `IDAa` - IDA array workflow.
- `IDAs` - IDA sensor workflow.
- `Ma` - mixed-host array workflow.
- `Ms` - mixed-host sensor workflow.
- `P`, `N`, `U` - positive, negative, and union subsets used in efficiency/correlation analysis.

## Checkpoints and Large Files

Long-running MATLAB scripts write checkpoint files such as:

- `Ma_res9U_checkpoint.mat`
- `Ms_res9U_checkpoint.mat`
- `IDAa_res4_checkpoint.mat`

Keep these files if you want to resume interrupted runs. Delete or archive them only when you are sure the corresponding final result file is complete.

## Reproducibility Notes

- Some scripts seed the random number generator for repeatable sampling or model evaluation.
- Evaluation scripts shuffle train/test rows during each run, so changes to random seeding or chunking can change exact results.
- The MATLAB scripts are designed for Windows paths and local generated files in the project root.
- Before rerunning a large simulation, check available disk space and memory.

## Suggested Cleanup Before Sharing

If publishing or sharing the repository, consider excluding generated artifacts with `.gitignore` rules for:

```gitignore
*.mat
*.csv
*.xlsx
*.png
*.rar
*.log
*_outputs*/
sorted_plots/
Et_comp/
```

Keep any final report files or curated figures that are needed for the paper/report.
