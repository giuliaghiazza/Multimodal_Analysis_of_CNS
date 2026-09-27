# AGENTS.md — EEG Asymmetry Analysis Pipeline

## Overview

MATLAB pipeline for group-level EEG asymmetry (AI), GFP, and HFP analysis
comparing rest (eyes-closed) vs. task (memory) using EEGLAB 2026.0.0 on BIDS
BrainVision data (64 channels, 500 Hz).

Subproject of the `Multimodal` repo — see `../README.md` for the other
projects (`fMRI/`, `NIRS/`).

## Pipeline execution order (strict)

1. `processing.m` — preprocesses raw .vhdr files → PROCESSED_DATA/*.set
2. `GFP_HFP.m` and `AI.m` — feature extraction (can run in parallel) → RESULTS/
3. `statisticalTest.m` — Wilcoxon + FDR-BH + max-t permutation → RESULTS/
   (`statisticalTestHedges_maxTfixed.m` is a variant also adding
   temporal-stability metrics)
4. `visualize.m` — bar plots, heatmap, topoplots from Group_Statistics.mat

Step 2 depends on step 1 output. Step 3 depends on step 2 output.
Step 4 depends on step 3 output.

## Environment

- EEGLAB is vendored in `../eeglab2026.0.0/` (repo root, gitignored).
  The scripts auto-add it to the MATLAB path via a guarded `addpath` —
  no manual setup needed if the folder exists. Otherwise install EEGLAB
  separately and `addpath` it once.
- `GFP_HFP.m` uses `bandpass()` → requires **MATLAB Signal Processing Toolbox**.

## Data

- Raw data: `DATA/` (BIDS 1.2, BrainVision .vhdr). Sub-01 through sub-09 + sub-43.
  Only **session 1** (`ses-session1`) and the target tasks are processed
  (see the `is_target_task` filter in `processing.m`).
- `DATA.zip` is a compressed archive of the same; extract if `DATA/` is empty.
- `PROCESSED_DATA_mathematic.zip` — archive of the mathematic-task processed run.
- `DATA/DATA` may contain the same dataset structure nested; `processing.m`
  uses recursive `dir(**/*.vhdr)` so both locations work.
- All data dirs are gitignored — nothing heavy is pushed.

## Key conventions

- Scripts use `fileparts(mfilename('fullpath'))` for base_dir — run from the
  script's own directory or expect paths relative to the script location.
- Output directories (`PROCESSED_DATA/`, `RESULTS/`) are auto-created if missing.
- Group D = 10 subjects (sub-01 through sub-09, sub-43).
- Frequency bands: Delta [1–4], Theta [4–8], Alpha [8–13], Beta [13–30], Gamma [30–45] Hz.

## Folders

- `experiments/` — past analysis runs (10 s vs 20 s windows, memory vs
  matematica tasks, STD variants). Local only, not on git.
- `archive/` — phase-1 work: single-subject asymmetry scripts (`asymmetry*.m`,
  `GFP_HFP_BandPowerAnalysis.m`), `assignment/` (sub-43 BIDS), `lab_exercises/`,
  `datasets/`, `immagini/`, `scratch/` (old `prova.m` + test datasets).
  Local only, not on git.
- `references/` — reference papers (Smith et al. 2017 frontal asymmetry primer).

## Gotchas

- The `bh_fdr` local function in `statisticalTest.m` sorts p-values; ensure
  MATLAB's `sort` with two outputs behaves as expected for your version.
- `run_wilcoxon` at `statisticalTest.m:259` omits a trailing semicolon on
  the validity-index assignment — harmless but noisy in output.
- `eeglab2026.0.0/` has its own docs and plugin READMEs — those are
  EEGLAB-internal, not project-specific.
