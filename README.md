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

  -  a filename list is processed to loop through all waveform data files
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