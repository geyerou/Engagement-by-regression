# Engagement by Regression

A MATLAB pipeline for analyzing brain engagement using ridge regression encoding models, functional connectivity, and gradient-based parcellation.

## Overview

This repository contains the analysis pipeline for the "Engagement by Regression" project, which investigates how brain networks contribute to cognitive engagement through:

- **Ridge regression encoding models** (`s004`, `s012`)
- **Brain space gradient analysis** (`s009a`, `s009b`)
- **Functional connectivity baselines** (`s011`, `s011b`)
- **Network specificity and contribution** (`s008`, `s010`)
- **Stability and null model analyses** (`s014`, `s015`)
- **Local extensions and validations** (`s021`–`s030`)
- **Behavioral prediction** (`s017`, `s026`)

## Pipeline Structure

Scripts are numbered in execution order (`sXXX_*.m`). The main entry points are:

| Script | Description |
|--------|-------------|
| `s000_config_project.m` | Project configuration |
| `s001_check_inputs_and_masks.m` | Input validation |
| `s002_extract_timeseries_with_wm_smoothing.m` | Time series extraction |
| `s003_build_model_matrices.m` | Design matrix construction |
| `s004_run_ridge_main_model.m` | Main ridge regression model |
| `s005_write_r2_maps.m` | R² map output |
| `s006_extract_beta_profiles.m` | Beta weight profiles |
| `s007_write_roi_beta_maps.m` | ROI beta maps |
| `RunAllFigures.m` | Generate all main figures |
| `RunSupplementFigures.m` | Generate supplementary figures |

## Helper Modules

- `functions/` — Core analysis functions
- `figure_helpers/` — Visualization utilities

## Requirements

- MATLAB (R2019b or later recommended)
- Required toolboxes: Statistics and Machine Learning Toolbox, Parallel Computing Toolbox (optional)

## Usage

1. Configure project paths in `s000_config_project.m`
2. Run scripts sequentially from `s001` onward
3. Generate figures with `RunAllFigures.m` and `RunSupplementFigures.m`
