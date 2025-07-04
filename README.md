# Modulogram Project

This repository contains MATLAB scripts used for generating modulograms and for computing statistics on phase-amplitude coupling (PAC).

## Repository Structure

- `modulogram_pipeline.m` – orchestrates the modulogram computation for one session using helpers in `src/`.
- `execute2.m` – example script showing how to configure and run the pipeline.
- `config_stats.m` and `run_pac_stats.m` – utilities to compute PAC statistics from modulogram output.
- `render_single_from_struct_file.m` – helper for rendering a modulogram from a saved data structure.
- `results_stats/` – example results including summary plots and statistics.

Most scripts expect data files that are not part of this repository (e.g., `loss-1`, `win-1`, or subject session files). See the comments inside each script for details about required paths and variables.

## Usage

1. Edit `execute2.m` to set up the configuration structure. Key fields include:
   - `subjects`       – cell array of subject IDs
   - `sessionNum`     – session number to process
   - `alignments`     – alignment conditions (e.g., `{'win','loss'}`)
   - `dataDir` and `resultsDir` – paths to your formatted data and where results should be written
   - frequency and window parameters
2. From MATLAB, run the script:
   ```matlab
   execute2
   ```
   which will call `modulogram_pipeline(config)`.
  Modulogram images for each channel are written to
  `resultsDir/<Subject>/Session_<N>/<alignment>/CH###_<alignment>_Modulogram.png`.
  A log file named `pipeline_<timestamp>.log` is also written to `resultsDir` and
  a consolidated MATLAB structure containing all computed data and metadata is
 saved as `allSubjects_session<N>_<timestamp>.mat` for downstream analysis.
3. After modulograms are generated, use `config_stats.m` followed by `run_pac_stats.m` to compute summary statistics and plots.
4. To analyze modulograms for subject EMU024 across three sessions, use
   `analyze_EMU024_sessions(baseDir, outDir)`. The optional `outDir`
   argument allows writing results to a custom directory (defaults to

## Requirements

These scripts require MATLAB with the Statistics and Signal Processing Toolboxes.

## License

This project is provided as-is for research purposes. See individual file headers for author information and usage restrictions.

