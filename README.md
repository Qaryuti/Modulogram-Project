# Modulogram Project

This repository contains MATLAB scripts used for generating modulograms and for computing statistics on phase-amplitude coupling (PAC).

## Repository Structure

- `run_modulogram_pipeline_wrapper.m` – main entry point to create modulograms for a collection of subjects and sessions.
- `config_stats.m` and `run_pac_stats.m` – utilities to compute PAC statistics from modulogram output.
- `render_single_from_struct_file.m` – helper for rendering a modulogram from a saved data structure.
- `results_stats/` – example results including summary plots and statistics.

Most scripts expect data files that are not part of this repository (e.g., `loss-1`, `win-1`, or subject session files). See the comments inside each script for details about required paths and variables.

## Usage

1. Edit `run_modulogram_pipeline_wrapper.m` to configure subject IDs, session numbers, and sampling rates.
2. Run the wrapper script from MATLAB:
   ```matlab
   run_modulogram_pipeline_wrapper(config);
   ```
   where `config` contains any parameters you wish to override.
3. After modulograms are generated, use `config_stats.m` followed by `run_pac_stats.m` to compute summary statistics and plots.

## Requirements

These scripts require MATLAB with the Statistics and Signal Processing Toolboxes.

## License

This project is provided as-is for research purposes. See individual file headers for author information and usage restrictions.

