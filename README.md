# ICARUS_PMT_calibration
PMT calibration scripts

1) The ProcessData_v1.py is the script processes the raw waveform data taken with a tektronix MSO64 scope.
  -  4 channels,
  -  20 us time window for each waveform
  -  each waveform is saved in a single file
  -  a filename list is processed to loop through all waveform data files
  -  a root file will be output

2) The FitChargeDistributions.C reads in the root file produced from ProcessData_v1.py (1), and fits charge
 distributions with Poisson function convoluted with a Gaussian function. The fit is for high charges, so a
 a combined fit is performed using a  single chi2 value and constant mu across three voltages.
  - Processes four PMTs at a time
  - Combined fitting code adapted from ROOT tutorials: combinedFit.C
  - Output:
     - Root file
     - .pdf of the canvas contents
     - .txt of comma separated values including:
          - final parameters: voltage, channel number
                              mu, mu error, q, q error, sigma, sigma error,
                              amplitude, amplitude error, chi2, ndf, fit probability
          - initial parameters: channel id,
                                start fit, end fit, rebin factor,
                                mu, q, sigma, amplitude
  - Uses parameter limits, which results in very poor errors
     - To improve the errors, save the initial parameters and rerun using
       FitChargeDistributions_InitParam.C
  - The *LedOn*result.root files must be included in the same folder that this macro is stored in
  - Do not have any *LedOff*result.root files in the folder
  - Processes three root files in a batch (same four PMTs at three different voltages)
  - Example function call in root:
       .x FitChargeDistributions.C("A10", 5, 6, 7, 8, 1400, 1430, 1460)
       - voltages as parameters must match the nominal values matching the *result.root files

3) The FitChargeDistributions_LC_batch.C reads in the root file produced from ProcessData_v1.py (1), and fits
 charge distributions with a Poisson function convoluted with a Gaussian function. The fit is for low charges,
 so the model equation of a Poisson function convoluted with a Gaussian function is fit to each distribution
 individually (each PMT, each voltage); contrast to FitChargeDistributions.C (2)
  - Processes four PMTs at a time (contrast to FitChargeDistributions_LC_ind.C (4))
  - Output:
    - Root file
    - .pdf of the canvas contents
    - .txt of comma separated values including:
       - final parameters: voltage, channel number
                           mu, mu error, q, q error, sigma, sigma error,
                           amplitude, amplitude error, chi2, ndf, fit probability
       - initial parameters: channel id,
                             start fit, end fit, rebin factor,
                             mu, q, sigma, amplitude
  - Uses parameter limits, which results in very poor errors
    - To improve the errors, save the initial parameters and rerun using
      FitChargeDistributions_InitParam.C
  - The *LedOn*result.root files must be included in the same folder that this macro is stored in
  - Do not have any *LedOff*result.root files in the folder
  - Processes three root files in a batch (same four PMTs at three different voltages)
  - Example function call in root:
      .x FitChargeDistributions_LC_batch.C("A10", 5, 6, 7, 8, 1400, 1430, 1460)
      - Voltages as parameters must match the nominal values matching the *result.root files

4) The FitChargeDistributions_LC_ind.C reads in the root file produced from ProcessData_v1.py (1), and fits
 charge distributions with a Poisson function convoluted with a Gaussian function. The fit is for low charges,
 so the model equation of a Poisson function convoluted with a Gaussian function is fit to each distribution
 individually (each PMT, each voltage); contrast to FitChargeDistributions.C (2)
  - Processes a single PMT at a time (contrast to FitChargeDistributions_LC_batch.C (4))
  - Output:
    - Root file
    - .pdf of the canvas contents
    - .txt of comma separated values including:
       - final parameters: voltage, channel number
                           mu, mu error, q, q error, sigma, sigma error,
                           amplitude, amplitude error, chi2, ndf, fit probability
       - initial parameters: channel id,
                             start fit, end fit, rebin factor,
                             mu, q, sigma, amplitude
  - Uses parameter limits, which results in very poor errors
    - To improve the errors, save the initial parameters and rerun using
      FitChargeDistributions_InitParam.C
  - The *LedOn*result.root files must be included in the same folder that this macro is stored in
  - Do not have any *LedOff*result.root files in the folder
  - Processes three root files in a batch (same four PMTs at three different voltages)
  - Example function call in root:
      .x FitChargeDistributions_LC_ind.C("A10", 5, 6, 7, 8, 1400, 1430, 1460, 0)
      - Last parameter is the index of the desired PMT out of the four, ie. 0, 1, 2, or 3
      - Voltages as parameters must match the nominal values matching the *result.root files

*** NOTE: It is necessary when fitting using (2 - 4) to go into the code to adjust initial parameters if desired

5) The FitChargeDistributions_InitParam.C is used to take in a .csv value of initial parameters and performs fits
  identical to those seen in (2 - 4), with the exception of the fact that parameter limits are no longer in place.
 - Example table (to export as .csv) can be found here:
        https://docs.google.com/spreadsheets/d/19rzmrZNi8R2X_QgjkQUhmoAZPIGufCa9ip4rDvazvKE/edit?usp=sharing
 - The macro utilizes a .csv file containing initial parameters for
   each PMT at three voltages in order to recreate the original fits.
    - Each PMT should have three rows dedicated to it
 - The .csv file must contain two headers of the following form:
   ,,,,Initial parameters,,,,,,,,Flags
   ROOT File Name,Chimney,PMT channel,Voltage,Fit Begin,Fit End, Rebin Factor,Y Max,Mu (NPE),q (SPE),Sigma (SPE),Amplitude,Charge
    - There should not be any empty rows or columns in the .csv file
 - The *LedOn*result.root files must be included in the same folder that this macro is stored in
 - Do not have any *LedOff*result.root files in the folder
 - Output (for each PMT):
     - Root file
     - .pdf of the canvas contents
     - .txt of comma separated values including:
          - final parameters: voltage, channel number
                              mu, mu error, q, q error, sigma, sigma error,
                              amplitude, amplitude error, chi2, ndf, fit probability
          - initial parameters: channel id,
                                start fit, end fit, rebin factor,
                                mu, q, sigma, amplitude
 - Example function call in root:
   FitChargeDistributions_InitParam("data.csv")