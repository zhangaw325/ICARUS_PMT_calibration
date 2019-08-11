# ICARUS_PMT_calibration
PMT calibration scripts

Data is processed in a series of steps. First, preprocessing is performed using one of the **Python** scripts. Charge distributions are fit using one of the **FitChargeDistribution** scripts. Dark rate analysis can be performed separately using scripts in the **dark_rate** directory, while afterpulse analysis can be performed with scripts in the **afterpulse** directory.

Bash scripts used to read analysis outputs and write them to text in a format convenient for spreadsheet entry are located in the **spreadsheet_compilation_scripts** directory. Bash scripts for creating and running Condor analysis jobs are located in the **condor_analysis_scripts** directory.

Finally, old code is located in the **old_code** directory for archival reasons.

## Preprocessing

The ProcessData_v1.py is the script processes the raw waveform data taken with a tektronix MSO64 scope. 
  -  4 channels,
  -  20 us time window for each waveform
  -  each waveform is saved in a single file
  -  a filename list is processed to loop through all waveform data files
  -  a root file will be output

## Charge distribution fits

The FitChargeDistributions.C reads in the root file produced from previous step, and fit charge distributions with Poisson function convoluted with a Guassian function.
  -  All 4 PMTs' distributions and all 3 HV data points are fit
  -  a text file will be output to give the fit parameters

## Dark rate analysis

## Afterpulse analysis
