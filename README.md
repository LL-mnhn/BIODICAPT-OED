# BIODICAPT-OED

## Overview
This repository focuses on **Optimal Experimental Designs (OEDs)** for the BIODICAPT project. 

> [!WARNING]
> **Early Version (V0).** This repository is still under active development and might not run as intented on external machines.

## Context
### BIODICAPT
BIODICAPT is a French initiative that aims at monitoring biodiversity of agricultural lands on a large (national) scale through the use of various recording devices.

The project will consist of two phases:
1. Data collection in a "research" network of agricultural plots 
2. Data collection in the 500 ENI network (voluntary participation of farmers)

The transition from the first phase to the second involves a significant shift in scale for data collection, which is where OED will be needed.

### OED
Optimal Experimental Design (OED) is a method for designing experiments to maximise information gain while minimising uncertainty, cost, time,...

The key concept is that **not all observations are equally informative**. Applied to BIODICAPT, this concept can be particularly interesting when working with 'participartory science' and other fields. Here, we want to find the best argicultural plots to conduct our observation, with the constraints of the project (participatory & reduced budget).

## Description
### Repository Structure
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
├── README.md               # This file
```

### Progress Tracker
- [X] Anonymize agricultural plots locations
- [X] Auto-load datasets
- [X] Pre-processing of datasets
- [X] Implementation of various models (HMSC, GJAM, RF, DNN,...)
- [ ] Implementation of custom cost function (functions?)
- [ ] Implementation of exchange algorithm
- [ ] Results and comparison


## Getting Started
### Clone the Repo
Clone this repository to your machine
```bash
cd /your/local/folder
git clone https://github.com/LL-mnhn/BIODICAPT-OED.git
```

### Install dependencies
> [!NOTE]
> Add renv usage

### Usage Notes
> [!NOTE]
> Add tips


## Contact
For questions and/or inquiries, please contact: <loic.lehnhoff@mnhn.fr>
