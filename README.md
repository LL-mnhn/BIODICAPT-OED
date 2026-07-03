# BIODICAPT-OED

## Overview
This repository focuses on **Optimal Experimental Designs (OEDs)** for the BIODICAPT project. 

> [!WARNING]
> **Early Version (V0).** This repository is still under active development and might not run as intented on external machines.

## Context
### BIODICAPT
BIODICAPT is a French initiative that aims at monitoring biodiversity of agricultural lands on a large (national) scale through the use of various recording devices.

The project consists of two phases:
1. Data collection in a "research" network of agricultural plots 
2. Data collection in the 500 ENI network (voluntary participation of farmers)

The transition from the first phase to the second involves a significant shift in scale for data collection, which is where OED will be needed (more details [here](resources/README.md)).

## Description
### Structure
```
├── scripts                 # Core scripts to run
├── R                       # Functions
├── data/
│   ├── config              # Configuration files
│   ├── raw_data            # Original datasets
│   └── preprocessed_data   # Preprocessed datasets
├── outputs/
│   ├── figures             # Visualizations and plots
│   └── results             # Statistical outputs and reports
├── renv                    # Information about R environment (packages)
├── DESCRIPTION             # Standard DESCRIPTION file for R packages
└── README.md               # This file
```

### Progress tracker
- [X] Anonymize agricultural plots locations
- [X] Load STOC dataset
- [X] Pre-processing
- [X] Training of HMSC models
- [ ] Implementation of custom cost function
- [ ] OED using HMSC on STOC
- [ ] Add other models (GJAM, RF, DNN,...)
- [ ] Results and comparison


## Getting Started
### How to use
1\. Clone this repository on your machine
```bash
cd /your/local/folder
git clone https://github.com/LL-mnhn/BIODICAPT-OED.git
```

2\. Install dependencies

Open `BIODICAPT-OED` as a new session in R (either with [Rstudio](https://docs.posit.co/ide/user/) or [Positron](https://positron.posit.co/welcome.html)). 

Install `renv` if not already installed on your machine. Then run:
```R
install.packages("renv")
renv::restore()
```

3\. Environment is ready, local scripts can be run.

### Usage Notes
Important results and figures are already saved in the `outputs` folder.

When running on an external machine: 
- `0-verify_datasets.R` can not be run (raw dataset are only accessible by authors)
- `1-pre_processing.R` can only be run partially (it processes raw datasets, which are unavailable; but it can show and save basic figures).
- Other files in `scripts` can be run without restrictions.

## Contact
For inquiries, please contact: <loic.lehnhoff@mnhn.fr>
