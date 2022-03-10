# MoonNMR 
A Shiny application for high-throughput Nuclear Magnetic Resonance workflows

MoonNMR is a Shiny application for the optimization of high-throughput NMR and 1H Metabolomics workflows. 

---
## Installation
Install directly from github using devtools (package)
```
devtools::install_github("funkam/MoonNMR")
```
---
## Start-Up
There is only single function in the package -> run_moonnmr()
To start:
```
library(MoonNMR)
run_moonnmr()
```
---
## Overview of modules and manual
---
## Preparation Module
### IVDr Template Creator
Create ICON template for use with IVDr protocols and Bruker SampleJet, either by upload or by sample-by-sample creation.
Includes automatic filling of sample spaces in racks (1-96), automatic change of rack spaces (1-5) and interactive experiment selection.
Allows quick creation of large sample tables ready for submission.

### ICON Template Creator
An additinal Tool for ICON is available that allows non-IVDr template creation with an upload of a list of experiments to be performed.

## Proccessing Module
### Import .xml into csv
Upload a zipped folder containing multiple .xml from B.I.QUANT-Methods or B.I.LISA
Create an organised and cleaned-up datatable.

### Clean-Up machine-CSV
Extract data from Bruker created files. Removes unwanted columns and uses HMBD names, when possible.
Due to the uneven text output by Bruker of Metabolic Profiles (;-csv) and Lipidomics Analysis (Tab-csv), it is possible to upload either one.

### Combine tables
Combine two tables when a shared column (UniqueID) is present. The columns can be selected from list.

## Analysis Module
### Normalization
Normalize data in preparation for parametrric statistical analysis. Histogramm and Distrubtion plot interactively change upon option selection.

### Automated Bar Plots
Upload a data table and interactively select columns to plot against a group column. Creation of large amount of plots possible, but are stored in a single PDF.
Suggested grouped plots (e.g. Lipid-Subclasses involving Triglycerides) are also available.

---
## Archive
The main difference to the web version of MoonNMR is that allows archiving of submitted samples, using a file called archive.csv found in the libary folder. A backup is also created for every submission. The 

## Examples
A few example are included to see the formatting for a few of the modules, but generally the format is very simple. 
For example, the user can upload an Experiment list to be used in ICON template creation:
| Experiments|                 
| ------------- | 
| NOESY  | 
| JRES  | 

Or the table of sample names uploaded for automatic ICON template creation
| Name|                 
| ------------- | 
| Sample 1  | 
| Sample 2 | 

## Contact
