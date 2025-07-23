# Industrial Wireless Analysis (rev2)

## Overview

This repository provides MATLAB tools for analyzing complex impulse responses from RF measurements in industrial wireless systems, developed by the National Institute of Standards and Technology (NIST). The `rev2` directory contains scripts and classes for processing RF measurement data, computing metrics such as path gain, K-factor, delay spread, and tap estimation, and generating reports. The primary script, `doallMobile.m`, orchestrates analysis in production or test modes, using a manifest file (`manifest.xlsx`) to specify measurement parameters. The code is organized into subdirectories for reporting, measurement run management, and metric computation.

**Author**: Rick Candell, NIST  
**Contact**: rick.candell@nist.gov

## Prerequisites

- **MATLAB**: Version compatible with the code (e.g., R2019b or later recommended).
- **Input Data**: `.mat` files containing RF measurement data (e.g., `pdp` structure with `delay` in seconds and `pow_dB` in dB), stored in `rev2/data/` or `rev2/dataCloud/`.
- **Manifest File**: An Excel file (`manifest.xlsx`) in the repository root, with columns for file paths, run IDs, and metadata (string, string, double, double, string, double, double).
- **Dependencies**:
  - MATLAB Signal Processing Toolbox (for CIR analysis, if used).
  - `MeasurementRun` class (in `rev2/@MeasurementRun/`).
  - `MeasurementRunMetric` class (in `rev2/@MeasurementRunMetric/`, if used).

## File Structure

- `rev2/`: Core directory for RF analysis tools.
  - `doallMobile.m`: Main script to process RF measurement data in production or test mode.
  - `reporting/`: Scripts for generating reports or visualizations (e.g., plots, text summaries).
  - `@MeasurementRun/`: MATLAB class directory containing `MeasurementRun.m`, defining an object for managing measurement runs and channel estimation.
  - `@MeasurementRunMetric/`: MATLAB class directory containing `MeasurementRunMetric.m`, for computing RF metrics (e.g., mean delay, RMS delay spread).
  - `data/`: Directory for production `.mat` files (e.g., `AAPlantD1_2GHz_TX1_vpol_run3_pp.mat`).
  - `dataCloud/`: Directory for test or cloud-based measurement data.
- `manifest.xlsx`: Excel file specifying measurement metadata.
- `stats.dat`: Output file for computed statistics.

## Usage

1. **Prepare Data**:
   - Place `.mat` files in `rev2/data/` (production) or `rev2/dataCloud/` (test).
   - Ensure `manifest.xlsx` is in the repository root with correct column types (string, string, double, double, string, double, double).

2. **Run the Script**:
   - Open MATLAB and navigate to the `rev2/` directory.
   - Execute the main script:
     ```matlab
     doallMobile()
     ```
   - **Production Mode** (`TEST=false`): Processes all `.mat` files in `rev2/data/` using `manifest.xlsx`.
   - **Test Mode** (`TEST=true`): Uses a test dataset (`rev2/data/AAPlantD1_2GHz_TX1_vpol_run3_pp.mat`) in `rev2/dataCloud/`.

3. **Analysis Options**:
   - The `analyze_cwd` function in `doallMobile.m` accepts an `OPTS` parameter to control tasks:
     - `1`: Compute path gain.
     - `2`: Compute K-factor.
     - `3`: Compute delay spread.
     - `4`: Compute average CIR.
     - `5`: Estimate number of taps (NTAP).
     - `6`: Write statistics to `stats.dat`.
     - `7`: Generate plots (disabled by default).
   - Use `OPTS='all'` to enable all options except NTAP approximation and plotting, or specify a custom vector:
     ```matlab
     OPTS = [1; 1; 1; 1; 1; 1; 1]; % Enable all options, including plots
     analyze_cwd(OPTS, 'manifest.xlsx');
     ```

4. **Output**:
   - Statistics are written to `stats.dat` in the current directory.
   - Plots (if enabled via `OPTS(7)=1`) are displayed as docked MATLAB figures.
   - Console messages indicate progress or errors (e.g., missing manifest file).

## Example

```matlab
% Run in production mode
cd rev2
doallMobile()

% Run in test mode (edit doallMobile.m to set TEST=true)
doallMobile()

% Custom analysis with plots enabled
OPTS = [1; 1; 1; 1; 1; 1; 1];
analyze_cwd(OPTS, 'manifest.xlsx');
```

## Notes

- **File Paths**: Update paths in `doallMobile.m` and `manifest.xlsx` to match your local setup. The default path (`C:\Users\rick\...`) is hardcoded for NIST’s environment.
- **Data Structure**: `.mat` files must contain a `pdp` structure with `delay` (seconds) and `pow_dB` (dB) fields.
- **Class Dependencies**: The script requires the `MeasurementRun` class in `rev2/@MeasurementRun/`, which defines methods like `estimate_channel_cwd`. Ensure this class is in MATLAB’s path.
- **Test Mode**: Uses a global `TEST_DATA` variable with a predefined dataset for debugging.
- **Subdirectories**:
  - `reporting/`: Likely contains scripts for generating plots or text reports (specific files not provided).
  - `@MeasurementRun/`: Contains `MeasurementRun.m`, implementing channel estimation and data management.
  - `@MeasurementRunMetric/`: Likely contains `MeasurementRunMetric.m` for computing metrics like mean delay or RMS delay spread.
- **Limitations**: The `reporting` directory’s contents are not specified; additional scripts may be required for full functionality.

## Citation

If you use this code, please cite:
- Candell, R. (2025). Industrial Wireless Analysis. National Institute of Standards and Technology. https://doi.org/10.18434/T4359D

## Contact

For questions or contributions, contact Rick Candell at rick.candell@nist.gov.
