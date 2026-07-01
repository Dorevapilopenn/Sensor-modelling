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
- `IDAs_sobol_experiment.m`, `IDAa_sobol_experiment.m`, `Ms_sobol_experiment.m`, `Ma_sobol_experiment.m` - Sobol-sampled optimization scripts that evaluate classification efficiency across chemical/design variables and classifier settings.
- `plot*.m`, `GUI_piers.m` - plotting and exploratory utilities.

### Python post-processing code

- `Ms_tables.py`, `IDAs_tables.py`, `IDAa_tables.py`, `Ma_tables.py` - generate correlation tables and heatmaps from CSV results.
- `maximafinder.py` - interpolates sampled response surfaces and finds local maxima.
- `convexchecker.py` - legacy helper for checking convexity or concavity of regular-grid response data.
- `filter.py` - small CSV filtering helper.

### Reports and outputs

- `IDA_report.tex`, `IDA_report.pdf` - report source and compiled report.
- `*_sobol_optimal*.json`, `*_sobol_checkpoint.mat` - generated Sobol result and checkpoint files.
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

### 1. Run Sobol optimization

Run the appropriate Sobol experiment from MATLAB. For example:

```matlab
Ma_sobol_experiment
```

The Sobol scripts process sampled chemical/design rows, evaluate PLS-DA classifier settings, save checkpoints, and write top-fraction JSON summaries. If a run stops early, rerunning the script will resume from the checkpoint when compatible.

### 2. Export or analyze results

Use the generated `.mat` or `.csv` files for analysis. The common result columns are:

- formation constants, such as `K1`, `K2`, ...
- derived spacing parameters, such as `deltaD`, `deltaK`, `deltaC`, `deltaG`
- total classification efficiency, usually named `Et`
- per-class efficiencies, usually named `Class1` and `Class2`

### 3. Generate figures

Run the Python table/plot scripts:

```bash
python Ms_tables.py
python IDAs_tables.py
```

These scripts read CSV result files and write correlation tables, mean-efficiency heatmaps, and stratified-correlation plots into the corresponding output folders.

### 4. Find local maxima

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

- `Ma_sobol_checkpoint.mat`
- `Ms_sobol_checkpoint.mat`
- `IDAa_sobol_checkpoint.mat`
- `IDAs_sobol_checkpoint.mat`

Keep these files if you want to resume interrupted runs. Delete or archive them only when you are sure the corresponding final result file is complete.

## Reproducibility Notes

- Some scripts seed the random number generator for repeatable sampling or model evaluation.
- Evaluation scripts shuffle train/test rows during each run, so changes to random seeding or chunking can change exact results.
- Sobol scripts write generated checkpoints and JSON summaries in the project root.
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
