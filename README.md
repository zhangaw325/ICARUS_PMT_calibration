# ICARUS_PMT_calibration
PMT calibration scripts

Data is processed in a series of steps. First, preprocessing is performed using one of the **Python** scripts. Charge distributions are fit using one of the **FitChargeDistribution** scripts. Dark rate analysis can be performed separately using scripts in the **dark_rate** directory, while afterpulse analysis can be performed with scripts in the **afterpulse** directory.

Bash scripts used to read analysis outputs and write them to text in a format convenient for spreadsheet entry are located in the **spreadsheet_compilation_scripts** directory. Bash scripts for creating and running Condor analysis jobs are located in the **condor_analysis_scripts** directory.

Finally, old code is located in the **old_code** directory for archival purposes.

## Preprocessing

To preprocess the raw data, **PreprocessData.py** is run. It assumes that raw waveform data was taken with a Tektronix MSO64 scope with the following settings: 

  -  4 channels,
  -  20 us time window for each waveform
  -  each waveform is saved in a single file

The script does the following:

  -  a ROOT file is output with various analysis histograms and selected raw waveforms
      +  the histograms include charge distributions
  -  the following information is printed:
      +  the directory name with the waveforms
      +  the number of waveforms for each PMT in this run
      +  the number of pulses counted
      +  the number of afterpulse candidates counted

To run the preprocessing analysis over an entire data set, a batch job can be submitted as follows: 
    
 1. a file named **WaveformDirectories.txt** must be created with the names of the waveform directories containing data. Only the directory names should be included, not the path to the directories. E.g. the directory at `/media/disk_a/ICARUS/PMT/Data201906/B10_PMT_5_6_7_8_1440V_LedOff` should be labeled as `B10_PMT_5_6_7_8_1440V_LedOff`.
    - this can be done, for example, using `ls > WaveformDirectories.txt` from within the directory of choice containing all of the data to be analyzed (e.g. `/media/disk_a/ICARUS/PMT/Data201906`)

 2. change the `data_path` variable on line 3 of **analyze_cl.sh** to the proper directory (e.g. `/media/disk_a/ICARUS/PMT/Data201906`)

 3. run the **create_waveform_analysis_files.sh** script
    - ensure that all the necessary files are present, including **fix_led_names.sh**
    - this will create individual folders for each group of data, to aid in parallelizing the analysis job

 4. ensure that the `input` line in **run_full_analysis_batch.job** is correct. It should point to the directory you are working in. Also change the number after `Queue` at the end of the file to equal the number of directories created by **create_waveform_analysis_files.sh**

 5. submit **run_full_analysis_batch.job** as a Condor job using `condor_submit run_full_analysis_batch.job`

For each directory containing data, the printed results of the preprocessing script will be written to a text file labeled as **\*output.txt**, and the ROOT file will be output as usual.

## Charge distribution fits

The FitChargeDistributions.C reads in the root file produced from previous step, and fit charge distributions with Poisson function convoluted with a Guassian function.
  -  All 4 PMTs' distributions and all 3 HV data points are fit
  -  a text file will be output to give the fit parameters
  
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

### Gain vs voltage analysis

The script **GainVoltage.C** performs analysis to determine the PMT gain as a function of voltage. It does the following:

 - **input:** a text file containing the following 4 space-separated columns:
    1. PMT number
    2. PMT voltage
    3. gain
    4. gain error

    Each row is information for a single PMT at a single voltage. The order of the rows should not matter, nor should it matter if a PMT number is missing. Problems will arise if more than 6 rows contain information for a single PMT.

    The function input is the name of the text file, without `.txt` at the end. E.g. if the input file is `A10.txt`, then `GainVoltage("A10")` should be called.

 - **output:**
     + PDF files with a gain-voltage fit for each PMT on both log-log and linear plots
     + ROOT file containing all gain-voltage fits
     + comma-separated text file containing fit parameters for each PMT, in increasing numerical PMT order

## Dark rate analysis

Dark pulses are counted using the number of pulses counted in LED-off data. These pulses are counted in preprocessing. Currently, dark rates themselves are calculated externally, in a spreadsheet. The script **dark_rate/HistogramDarkRate.C** exists solely to produce histograms of the dark rates.

## Afterpulse analysis

Most of the afterpulse probability determination is performed during preprocessing. The number of afterpulse candidate events is recorded then. In order to convert to a probability with dark rate subtraction, the Bash script **afterpulse/calculate_afterpulse_probabilities.sh** is used. The script assumes that three data sets are collected for each chimney --- the **\*.output** files for the chimney must all be present in the same directory as the script. As an input, the script takes a chimney name, e.g. `A10`. It then searches for the output files with `A10` at the front, and proceeds from there.

**afterpulse/AnalyzeAfterpulse.C** contains some preliminary code to do some basic afterpulse time-structure analysis.
