# ICARUS_PMT_calibration
PMT calibration scripts

Data is processed in a series of steps. First, preprocessing is performed using one of the **Python** scripts. Charge distributions are fit using one of the **FitChargeDistribution** scripts. Dark rate analysis can be performed separately using scripts in the **dark_rate** directory, while afterpulse analysis can be performed with scripts in the **afterpulse** directory.

Bash scripts used to read analysis outputs and write them to text in a format convenient for spreadsheet entry are located in the **spreadsheet_compilation_scripts** directory. Bash scripts for creating and running Condor analysis jobs are located in the **condor_analysis_scripts** directory.

Finally, old code is located in the **old_code** directory for archival reasons.

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

## Dark rate analysis

## Afterpulse analysis
