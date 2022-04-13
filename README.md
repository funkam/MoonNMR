# MoonNMR 
<b>M</b>etabol<b>O</b>mics <b>O</b>rga<b>N</b>iser <b>NMR</b> -A Shiny application for high-throughput Nuclear Magnetic Resonance workflows
<img src="https://user-images.githubusercontent.com/88379260/157672281-8f3902d3-998e-48cc-a445-25dc17a42fa5.png" width="130" height="130" align="right">



MoonNMR is a Shiny application for the optimization of high-throughput NMR and <sup>1</sup>H Metabolomics workflows. 
The application is designed primarily to be used with the Bruker BioSpin Avance IVDr (and methods), and a SampleJet sample robot. MoonNMR provides for tools for template creation, data extraction and manipulation and a series of analysis tools. However, the modules are designed to work independently and can easily be used by other researchers.


## Installation


#### Option 1: Via Github
```
library(devtools)
devtools::install_github("funkam/MoonNMR",build = TRUE, build_vignettes = TRUE)
```

#### Option 2: Download .tar.gz [here](https://drive.google.com/file/d/1mm_5oNZlgQsQgISwC4KD7OKvrnR0gjLB/view?usp=sharing)
```
Use R Studio or
R CMD INSTALL MoonNMR_1.1.2.tar.gz
```

---

## Start-Up
There are only two single functions in the package. One for creating the archive files to be used and one to start the shiny application.
The archive functions creates an archive in the current directory, it should be run before MoonNMR is started for the first time. This step can also be skipped if the Preparation Module is not going to be used.

Before first time start:
```
moonnmr_archivecreator()
```

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
Create ICON templates for use with IVDr protocols and Bruker SampleJet. There are two sub-functions: Either upload a file with sample names or create a template sample by sample.
Both functions includes automatic filling of sample spaces in racks (1-96), automatic change of rack spaces (1-5) and interactive experiment selection.
The templates can be directly uploaded into ICON-NMR. In ICON-NMR it might be necessary to set up the header names the first time they are imported.
#### Instructions
There are 2 sets of tools available for template creation for import into ICON. An uploader function and a manual creation function. The first set are the IVDr tools, here the template is geared towards the use of the IVDr methods by Bruker BioSpin (e.g. B.I.Methods, B.I.Lisa). The ICON version allows manual input of solvent and the upload of a file for experiments. This tool can be used for the creation of any ICON template.
##### Uploader
In the uploader function a table can be uploaded with sample names and the input of SampleJet parameters (Rack, Position) as well as experimental parameters can be selected dynamically.
There is a example table in ~/app/examples, but it follows a very simple format:
| Name|                 
| ------------- | 
| Sample 1  | 
| Sample 2 | 

The samples are automatically archived. The final table can be downloaded as .xlsx.
#### Creator
The creator function allows creation of a table sample by sample. Here a sample can be scanned or entered and the table will be built sample by sample.
The final table can be downloaded as .xlsx.

#### ICON Version
The ICON version simply allows a direct input of solvents and experiments instead of pre-selected options. A list of experiments can be uploaded via the button, it requires a simple single column .csv file in the following style. Technically the entire experiment library can be downloaded from TOPSPIN and uploaded here.
| Experiments|                 
| ------------- | 
| NOESY  | 
| zg30 | 
#### Archiving
All submitted samples are saved in the Archive.csv (+archive_backup.csv) that was created in the beginning.

## Proccessing Module
### Import .xml into csv
Upload a zipped folder containing multiple .xml from B.I.QUANT-Methods or B.I.LISA
Create an organised and cleaned-up datatable, which can be downloaded as .csv.

### Clean-Up machine-CSV
Extract data from Bruker created files. Removes unwanted columns and uses HMBD names, when possible.
Due to the uneven text output by Bruker of Metabolic Profiles (B.I.Methods uses ;-csv) and Lipidomics Subclass-Analysis (B.I.Lisa uses Tab-csv), it is possible to upload either one.
### Delete
An interface to allow deletion of specific columns (parameters), rows (samples) or whole groups interactively.
### Combine 
Combine two tables sharing a common column. Designed for unification of sample data information and sample grouping information. Or combining outputs from different panels.

## Analysis Module
### Normalization
Normalize data in preparation for parametric statistical analysis (e.g. log10, MeanCenter). Histogramm and Distrubtion plot interactively change upon option selection. Plots can be downloaded directly.

### Box Plots
Interactively select the grouping column and the parameters to be plotted. A single plot is visible on the right for a quick reference. A small selection of plots can be grouped together in a single TIFF-file. A selection of paramters can also be plotted in a PDF-File containing all plots individually.
### PCA Analysis
Principal component analysis can be performed and visualized by scores and loadings/biplots. The use plotly and can be interactively changed by the user.
