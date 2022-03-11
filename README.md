# MoonNMR 
<img src="https://user-images.githubusercontent.com/88379260/157672281-8f3902d3-998e-48cc-a445-25dc17a42fa5.png" width="130" height="130" align="right">
A Shiny application for high-throughput Nuclear Magnetic Resonance workflows

MoonNMR is a Shiny application for the optimization of high-throughput NMR and <sup>1</sup>H Metabolomics workflows. 
The application is designed primarily to be used with the Bruker BioSpin Avance IVDr (and methods), and a SampleJet sample robot. MoonNMR provides for tools for easier template creation, as well as data extraction tools and a series of analysis tools. However, the tools can also be used by other users

---
## Installation


#### Option 1: Via Github
```
library(devtools)
devtools::install_github("funkam/MoonNMR")
```

#### Option 2: Download .tar.gz [here](https://drive.google.com/file/d/1iAUxgL9HdxZ7EBpT8dv8WB1tFJFtOCJe/view?usp=sharing)
```
Use R Studio or
R CMD INSTALL MoonNMR_1.0.0.tar.gz
```

---

## Start-Up
There are two single functions in the package. One for creating the archive files to be used and one to start the shiny application.
The archive functions creates an archive in the current directory, it should be run before MoonNMR is started for the first time. It has a simple format. It can then be pre-populated with previous samples.

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
Create ICON template for use with IVDr protocols and Bruker SampleJet, either by upload or by sample-by-sample creation.
Includes automatic filling of sample spaces in racks (1-96), automatic change of rack spaces (1-5) and interactive experiment selection.
Allows quick creation of large sample tables ready for submission.

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


## Proccessing Module
### Import .xml into csv
Upload a zipped folder containing multiple .xml from B.I.QUANT-Methods or B.I.LISA
Create an organised and cleaned-up datatable, which can be downloaded as .csv.

### Clean-Up machine-CSV
Extract data from Bruker created files. Removes unwanted columns and uses HMBD names, when possible.
Due to the uneven text output by Bruker of Metabolic Profiles (B.I.Methods uses ;-csv) and Lipidomics Subclass-Analysis (B.I.Lisa uses Tab-csv), it is possible to upload either one.

### Combine tables
A simple tool to combine two tables sharing a common column. Designed for unification of sample data information and sample grouping information. Or combining outputs from different panels.

## Analysis Module
### Normalization
Normalize data in preparation for parametric statistical analysis (e.g. log10, MeanCenter). Histogramm and Distrubtion plot interactively change upon option selection.

### Automated Bar Plots
Upload a data table and interactively select columns to plot against a group column. Creation of large amount of plots possible, but are stored in a single PDF.
Suggested grouped plots (e.g. Lipid-Subclasses involving Triglycerides) are also available.

---
## Differences to Web-Version
### Archive
The main difference to the web version of MoonNMR is that allows archiving of submitted samples, using a file called archive.csv found in the libary folder. A backup is also created for every submission. There is an overview plot embedded on the main page detailing amount of samples measured and contribution of each project to the overall count.




