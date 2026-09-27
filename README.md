# Multimodal

MATLAB pipelines for multimodal neuroimaging — three lab projects in one repo.

| Folder | Project | Description |
|--------|---------|-------------|
| `EEG_Asymmetry/` | EEG asymmetry | Group-level frontal asymmetry index (AI), GFP/HFP features, Wilcoxon + FDR + max-t permutation stats. BIDS BrainVision data, 64 ch, 500 Hz |
| `fMRI/` | fMRI | SPM batch preprocessing + ROI-to-ROI sliding-window Pearson connectivity |
| `NIRS/` | fNIRS | Homer3 GLM pipeline (`Codes/`) + TD-NIRS motor / DPF processing (`Motor/`) |

## Data policy

Raw and processed data stay **local only** (see `.gitignore`): `DATA.zip`,
`PROCESSED_DATA/`, `experiments/`, `archive/`, `*.nii`, `*.fdt`, `*.mat`, …

On git: code, docs, and small figures (`.png`, `.csv`, `.txt`).

To reproduce the EEG analysis, extract `EEG_Asymmetry/DATA.zip` to
`EEG_Asymmetry/DATA/` (BIDS 1.2, BrainVision `.vhdr`).

## Requirements

- **MATLAB** + Signal Processing Toolbox (`GFP_HFP.m` uses `bandpass`)
- **EEGLAB 2026.0.0** — vendored in `eeglab2026.0.0/` (gitignored, not pushed).
  The EEG scripts auto-add it to the MATLAB path; if missing, reinstall EEGLAB
  and run `addpath('eeglab2026.0.0')` once.
- **SPM12** for `fMRI/`
- **Homer3** for `NIRS/Codes/` (edit `HOMER3_DIR` in `PipiLineHomer.m`)

## EEG pipeline order

Run from `EEG_Asymmetry/`:

1. `processing.m` — preprocess `.vhdr` → `PROCESSED_DATA/*.set`
2. `GFP_HFP.m` and `AI.m` — feature extraction → `RESULTS/` (can run in parallel)
3. `statisticalTest.m` (or `statisticalTestHedges_maxTfixed.m`) — group stats
4. `visualize.m` — bar plots, heatmaps, topoplots → `RESULTS/FIGURES/`

See `EEG_Asymmetry/AGENTS.md` for the full reference.

## Repo layout

```
Multimodal/
├── EEG_Asymmetry/        # EEG asymmetry project (phases 1+2)
│   ├── *.m               # pipeline scripts
│   ├── DATA.zip          # raw BIDS data (local)
│   ├── PROCESSED_DATA/   # preprocessed .set/.fdt (local)
│   ├── RESULTS/          # .mat results (local) + FIGURES/*.png (git)
│   ├── references/       # reference papers
│   ├── experiments/      # past runs: 10s/20s windows, memory/matematica (local)
│   └── archive/          # phase-1 scripts & single-subject data (local)
├── fMRI/
│   ├── corr_analysis.m   # sliding-window ROI connectivity
│   ├── spm_batch/        # SPM batch .mat files
│   ├── roi/              # ROI definitions & locations
│   └── data/, output/, outputanalysis/, cnn_sliding_window_output/ (local)
├── NIRS/
│   ├── Codes/            # Homer3 pipeline + GLM + ROI configs
│   ├── Motor/            # TD-NIRS / DPF scripts & data
│   └── TDNIRS_output/    # results (local)
└── eeglab2026.0.0/       # vendored EEGLAB (local only)
```
