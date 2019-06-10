# ICARUS_PMT_calibration
PMT calibration scripts

The ProcessData_v1.py is the script processes the raw waveform data taken with a tektronix MSO64 scope. 
  -  4 channels,
  -  20 us time window for each waveform
  -  each waveform is saved in a single file
  -  a filename list is processed to loop through all waveform data files
  -  a root file will be output
  
The FitChargeDistributions.C reads in the root file produced from previous step, and fit charge
 distributions with Poisson function convoluted with a Guassian function.
  -  All 4 PMTs' distributions and all 3 HV data points are fit
  -  a text file will be output to give the fit parameters
