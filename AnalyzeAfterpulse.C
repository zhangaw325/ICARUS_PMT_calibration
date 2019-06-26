#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TStyle.h"
#include "TSpectrum.h"

// function declarations
Double_t* findPeaksInRange(TH1 *hist, Double_t start, Double_t end, Double_t sigma = 2, Double_t threshold = 0.05, Int_t maxPeaks = 2);
Double_t countAfterpulse(TH1 *hist, Double_t startTime);
Double_t calculateAfterpulseProb(TH1 *hist, Double_t startTime);
Double_t getXMeanInRange(TH1 *hist, Double_t start, Double_t end);
Double_t getXMeanErrorInRange(TH1 *hist, Double_t start, Double_t end);
void scaleXAxis(TH1 *hist, Double_t scaleFactor);

void AnalyzeAfterpulse(string pmtRow,
			    char pmt1, char pmt2, char pmt3, char pmt4,
			    int volt1, int volt2, int volt3, bool led){

  const int NCH = 4; // 4 PMTs
  // histogram and fit options
  int rbf_0 = 2;
  double fbc_0 = 0.0;
  double fec_0 = 90.0;
  int rebinfactor[NCH]={rbf_0, rbf_0,rbf_0,rbf_0}; // rebin histograms
  double fitbeginch[NCH]={fbc_0,fbc_0,fbc_0,fbc_0}; // fit start locations
  double fitendch[NCH] = {fec_0,fec_0,fec_0,fec_0}; // fit end locations

  // the function to be used to do fit
  gStyle->SetOptFit(1111);

  double par[6];
  double parerr[6];

  const int NVOLTAGES = 1;

  string rtfilenames[NVOLTAGES];
  string resultnames[4];
  double channelnames[4] = {0,1,2,3};
  string strchimney = pmtRow + "_PMT_";
  string strpmt = to_string(pmt1) + "_" + to_string(pmt2) + "_" + to_string(pmt3) + "_" + to_string(pmt4) + "_";
  string voltagestr[3] = {to_string(volt1), to_string(volt2), to_string(volt3)};
  string ledstr;
  if(led)
    ledstr = "On";
  else
    ledstr = "Off";
  
  // Setup input files; ONLY considering highest voltage
  for(int i=0; i < NVOLTAGES; i++){
    rtfilenames[i]  = strchimney + strpmt + voltagestr[i+2] + "V_Led" + ledstr + "_result.root";
    cout << rtfilenames[i] << endl;
  }
  
  // Setup output files
  for(int i=0;i<4; i++){
    resultnames[i]  = strchimney + strpmt + "CH" + channelnames[i] + "_" + "Led" + ledstr + "_afterpulse.pdf";
    cout << resultnames[i] <<endl;
  }

  string outnameroot = strchimney + strpmt + "afterpulse.root";
  string outnametxt = strchimney + strpmt + "afterpulse.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TH1F* hPulseDist[NVOLTAGES]; // histograms for each canvas
  TCanvas* c[4];
  char tempname[100];

  TFile* files[NVOLTAGES];
  // Read each of the three data ROOT files and save them to files[]
  for(int i = 0; i < NVOLTAGES; i++){
    files[i] = new TFile(rtfilenames[i].c_str(),"read");
  }


  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);
  
  // Generate 4 canvases and plot the PMT histograms on them
  for(int i = 0; i < 4; i++){
    sprintf(tempname, "c_%d",i);
    string canvasTitle = strchimney + to_string(i+1);
    //sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
    c[i] = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
    c[i]->Divide(NVOLTAGES); // divide canvas into 3 pads along the width
    c[i]->SetLogy(); // set log scale on y axis

    Double_t afterpulseProb[NVOLTAGES];
    Double_t firstPulseMean[NVOLTAGES];
    Double_t firstPulseMeanError[NVOLTAGES];
    Double_t secondPulseMean[NVOLTAGES];
    Double_t secondPulseMeanError[NVOLTAGES];

    // generate histograms and calculate values
    sprintf(tempname, "Results/PulseTimeDist_%d",i);
    for(int j = 0; j < NVOLTAGES; j++){
      hPulseDist[j] = (TH1F*)files[j]->Get(tempname);
      hPulseDist[j]->SetTitle((voltagestr[j] + "V").c_str());

      scaleXAxis(hPulseDist[j], 1.6/1000); // scale X axis to be a time, rather than a bin number

      hPulseDist[j]->Rebin(rebinfactor[i]);
      hPulseDist[j]->SetXTitle("Pulse time (us; CHECK LATER)");

      // calculate afterpulse probability
      afterpulseProb[j] = calculateAfterpulseProb(hPulseDist[j], 0.48); // 480 ns to cut out LED pulse

      // cout << afterpulseProb[j] << endl;

      // switch to proper pad to begin fitting
      c[i]->cd(j+1);

      // find pulses
      TF1 *pulse1 = new TF1("pulse1", "gaus", 0.5, 3.5);
      TF1 *pusle2 = new TF1("pulse2", "gaus", 4.5, 8.5);

      hPulseDist[j]->Fit(pulse1, "R");
      hPulseDist[j]->Fit(pulse2, "R+"); // fit and add to list of fitted functions

      // Double_t* peakLocs = findPeaksInRange(hPulseDist[j], 0.48, 4, 1.25, 0.05, 1); // find first two peaks between 0.48 us and 4 us
      // firstPulseMean[j] = peakLocs[0];

      // Double_t binWidth = hPulseDist[j]->GetBinWidth(0);
      // firstPulseMeanError[j] = 2 * binWidth; // two times bin width for the error, for now

      //TF1 *fitSecond = (TF1*)hPulseDist[j]->GetListOfFunctions()->FindObject("gaus");
      firstPulseMean[j] = pulse1->GetParameter(1);
      firstPulseMeanError[j] = pulse1->GetParError(1);
      secondPulseMean[j] = pulse2->GetParameter(1);
      secondPulseMeanError[j] = pusle2->GetParError(1);

      // cout << secondPulseMean[j] << "\t" << secondPulseMeanError[j] << endl;

      // draw histogram
      hPulseDist[j]->Draw();
      gPad->SetLogy();

    }

    // write to output txt file
    // foutFit << "********** afterpulse probabilities **********" << endl;

    // write to output txt file
    foutFit << "Afterpulse Prob:\tchID\t" << i;
    for(int j = 0; j < NVOLTAGES; j++){
    	foutFit << "\t" << afterpulseProb[j];
    }
    foutFit << endl;

    // foutFit << "********** pulse locations and errors (first then second) **********" << endl;

    for(int j = 0; j < NVOLTAGES; j++){
    	foutFit << "Pulse Locs:\tchID\t" << i 
    			<< "\t" << firstPulseMean[j]
    			<< "\t" << firstPulseMeanError[j]
    			<< "\t" << secondPulseMean[j]
    			<< "\t" << secondPulseMeanError[j]
    			<< endl;
    }

    // write results to output ROOT file
    outROOTfile->cd();
    c[i]->Write();
    c[i]->Print(resultnames[i].c_str(),"pdf");
  }

  // close output ROOT file
  outROOTfile->Close();
}

Double_t* findPeaksInRange(TH1 *hist, Double_t start, Double_t end, Double_t sigma = 2, Double_t threshold = 0.05, Int_t maxPeaks = 2){
	if (!hist) return NULL; // precaution

	TAxis *xAxis = hist->GetXaxis();

	Double_t originalStart = xAxis->GetXmin();
	Double_t originalEnd = xAxis->GetXmax();

	// set range within which to calculate mean
	xAxis->SetRangeUser(start, end);

	TSpectrum *spec = new TSpectrum(maxPeaks);
	spec->Search(hist, sigma, "", threshold);

	// set original range
	xAxis->SetRangeUser(originalStart, originalEnd);

	return spec->GetPositionX();
}

// Count number of pulses due to afterpulsing
Double_t countAfterpulse(TH1 *hist, Double_t startTime){
	if (!hist) return -1.0; // precaution

	TAxis *xAxis = hist->GetXaxis();

	// store axis range for future reference
	Double_t originalMin = xAxis->GetXmin();
	Double_t originalMax = xAxis->GetXmax();

	// change range to cut out LED pulse and count pulses in the range
	xAxis->SetRangeUser(startTime, originalMax);
	Double_t afterpulseCounts = hist->Integral();

	// set histogram to original range
	xAxis->SetRangeUser(originalMin, originalMax);

	// return
	return afterpulseCounts;
}

Double_t calculateAfterpulseProb(TH1 *hist, Double_t startTime){
	if (!hist) return -1.0; // precaution

	Double_t totalEntries = hist->GetEntries();
	Double_t afterpulseCounts = countAfterpulse(hist, startTime);

	Double_t prob = afterpulseCounts/totalEntries;

	return prob;
}

Double_t getXMeanInRange(TH1 *hist, Double_t start, Double_t end){
	if (!hist) return -0.0; // precaution

	TAxis *xAxis = hist->GetXaxis();

	Double_t originalStart = xAxis->GetXmin();
	Double_t originalEnd = xAxis->GetXmax();

	// set range within which to calculate mean
	xAxis->SetRangeUser(start, end);

	// calculate mean
	Double_t mean = hist->GetMean(1);

	// set original range
	xAxis->SetRangeUser(originalStart, originalEnd);

	// return
	return mean;
}

Double_t getXMeanErrorInRange(TH1 *hist, Double_t start, Double_t end){
	if (!hist) return -0.0; // precaution

	TAxis *xAxis = hist->GetXaxis();

	Double_t originalStart = xAxis->GetXmin();
	Double_t originalEnd = xAxis->GetXmax();

	// set range within which to calculate mean
	xAxis->SetRangeUser(start, end);

	// calculate mean
	Double_t meanError = hist->GetMeanError(1);

	// set original range
	xAxis->SetRangeUser(originalStart, originalEnd);

	// return
	return meanError;
}

// Code to scale X axis
void scaleXAxis(TH1 *hist, Double_t scaleFactor){
	if (!hist) return; // precautionary

	TAxis *axis = hist->GetXaxis();

	Double_t newXmin = axis->GetXmin()*scaleFactor;
	Double_t newXmax = axis->GetXmax()*scaleFactor;

	axis->Set(axis->GetNbins(), newXmin, newXmax); // new Xmax

  	return;
}
