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

// function declarations
double IdealResponse(double *x,double *par);
Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3);

void AnalyzeAfterpulse(string pmtRow,
			    char pmt1, char pmt2, char pmt3, char pmt4,
			    int volt1, int volt2, int volt3, bool led){

  const int NCH = 4; // 4 PMTs
  // histogram and fit options
  int rbf_0 = 5;
  double fbc_0 = 0.0;
  double fec_0 = 90.0;
  int rebinfactor[NCH]={rbf_0, rbf_0,rbf_0,rbf_0}; // rebin histograms
  double fitbeginch[NCH]={fbc_0,fbc_0,fbc_0,fbc_0}; // fit start locations
  double fitendch[NCH] = {fec_0,fec_0,fec_0,fec_0}; // fit end locations

  // the function to be used to do fit
  gStyle->SetOptFit(1111);

  double par[6];
  double parerr[6];

  // process 3 files in a batch
  string rtfilenames[3];
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
  
  // Setup input files
  for(int i=0; i<3; i++){
    rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_Led" + ledstr + "_result.root";
    cout << rtfilenames[i] << endl;
  }
  
  // Setup output files
  for(int i=0;i<4; i++){
    resultnames[i]  = trchimney + strpmt + "CH" + channelnames[i] + "_" + "Led" + ledstr + "_afterpulse.pdf";
    cout << resultnames[i] <<endl;
  }

  string outnameroot = strchimney + strpmt + "afterpulse.root";
  string outnametxt = strchimney + strpmt + "afterpulse.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TH1F* hCharge[3]; // histograms for each canvas
  TCanvas* c[4];
  char tempname[100];

  TFile* files[3];
  // Read each of the three data ROOT files and save them to files[]
  for(int i = 0; i < 3; i++){
    files[i] = new TFile(rtfilenames[i].c_str(),"read");
  }


  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);
  
  // Generate 4 canvases and plot the PMT histograms on them
  for(int i = 0; i < 4; i++){
    sprintf(tempname, "c_%d",i);
    string canvasTitle = strchimney + to_string(i+1);
    //sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
    c[i] = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
    c[i]->Divide(3); // divide canvas into 3 pads along the width
    c[i]->SetLogy(); // set log scale on y axis

    // generate histograms
    sprintf(tempname, "Results/PulseTimeDist_%d",i);
    for(int j = 0; j < 3; j++){
      hCharge[j] = (TH1F*)files[j]->Get(tempname);
      hCharge[j]->SetTitle((voltagestr[j] + "V").c_str());
      hCharge[j]->Rebin(rebinfactor[i]);
      hCharge[j]->SetXTitle("Pulse time bin (16 ns/bin)");
    }

    // write results to output ROOT file
    outROOTfile->cd();
    c[i]->Write();
    c[i]->Print(resultnames[i].c_str(),"pdf");
  }

  // close output ROOT file
  outROOTfile->Close();
}
