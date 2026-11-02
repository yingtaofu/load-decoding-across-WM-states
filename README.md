# Load decoding across WM states/data and scripts

This repository contains behavioral data, MATLAB analysis scripts, and documentation for the load-decoding-across-WM-states study. The preprocessed EEG files are hosted on Zenodo and are not included in this GitHub repository.

## Contents

- `behav/DataAll.xlsx`: behavioral data from 33 participants.
- `behav/statistics_plot.m`: summarizes probe accuracy by load and attended/unattended state, performs repeated-measures statistics and paired comparisons, and plots accuracy.
- `EEG/cleaned_data/README.md`: instructions and link for obtaining the preprocessed EEG data from Zenodo.
- `EEG/*.m`: EEG analysis scripts and helper functions.
- `EEG/biosemi64_layout.mat`: channel layout file.

## Data availability

The preprocessed EEG dataset for 33 participants is available on Zenodo: [Download the EEG dataset](https://doi.org/10.5281/zenodo.23764504).

After downloading the dataset, place each participant's matching `.set` and `.fdt` files in `EEG/cleaned_data/` before running the EEG analysis scripts. Keep each pair together.

## Behavioral data fields

| Columns | Meaning |
| --- | --- |
| `Name`, `Id`, `Sex`, `Age` | Participant label and demographics. `Name` is not a real name. |
| `Block`, `Trial` | Block and trial numbers. |
| `AtLoad`, `UnatLoad` | Attended and unattended memory load. |
| `AtPos` | Cued hemifield: 1 = left, 2 = right. |
| `CueColor`, `nonCueColor` | Color labels for attended and unattended items. |
| `IfChange1`, `Probe1Color` | Change status (1 = change, 0 = no change) and color for probe 1, relative to cued items. |
| `IfChange2`, `Probe2Color` | Corresponding information for probe 2, relative to uncued items. |
| `Probe1Key`, `Probe1ACC`, `Probe1RT` | Probe 1 response key (J = change, F = no change), accuracy, and response time. |
| `Probe2Key`, `Probe2ACC`, `Probe2RT` | Corresponding measures for probe 2. |

## EEG bin coding

The preprocessed EEG data are hosted on Zenodo (link above). Epochs are time-locked to memory-array onset (time 0). A and U refer to attended and unattended load; the number is the number of items.

| Bin | Condition | Cued hemifield |
| --- | --- | --- |
| 1 | A1U1 | Left |
| 2 | A1U1 | Right |
| 3 | A1U3 | Left |
| 4 | A3U1 | Right |
| 5 | A3U1 | Left |
| 6 | A1U3 | Right |

The `.set` files can be opened in MATLAB with EEGLAB. Sampling rate, channel information, epoch timing, and trial information are stored in each dataset.

## Analysis scripts

| Script | Purpose |
| --- | --- |
| `TF_HilbertTransform_6to16hz.m` | 6-16 Hz Hilbert-transform power analysis. |
| `decoding_timeresolved.m` | Time-resolved MVPA decoding of attended and unattended load. |
| `searchlight_chantime_remap.m` | Channel-by-time searchlight analysis|
| `timegenralization_parietooccipital.m` | Time-generalization decoding using parieto-occipital electrodes. |
| `crossdecoding_parietooccipital_AtoU.m` | Cross-decoding from attended to unattended representations. |
| `crossdecoding_parietooccipital_UtoA.m` | Cross-decoding from unattended to attended representations. |
| `mirror_trials_lr_WML.m` | Helper function that swaps left/right hemisphere channel data for right-cued trials. |
| `cluster_based_permutation_test.m` | Cluster-based permutation test helper. |

## Requirements and use

The scripts require MATLAB. EEG analyses call EEGLAB, FieldTrip, and MVPA-light functions. The behavioral script uses MATLAB table, repeated-measures, plotting, and density-estimation functions. Analyses were run using MATLAB R2023a, EEGLAB [2025.0.0.], FieldTrip [20230926], and MVPA-light-master. Add required packages to the MATLAB path and set the current folder to the appropriate analysis folder before running scripts, because paths are relative to the current folder.

## Data notes

- Block 6 was excluded for participant 5 in both behavioral and EEG data.
- Block numbers may not be consecutive because some participants exited and restarted the experiment. Participants completed eight blocks, except when extra blocks were completed to replace unusable trials.
- There were two participants with `Id = 1`. Distinguish participants using the combination of `Id` and `Name` when joining or summarizing data.

## Contact

For questions, please contact Yingtao Fu at yingtaofu@zju.edu.cn.
