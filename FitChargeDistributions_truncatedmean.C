#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"

//#include "Fit_Spe.C" // the fit function is defined in this file

double PI = TMath::Pi();

double IdealResponse(double *x,double *par);
Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3);

void FitChargeDistributions_truncatedmean(string pmtRow,
			    char pmt1, char pmt2, char pmt3, char pmt4,
			    int volt1, int volt2, int volt3, bool led){

  const int NCH = 4; // 4 PMTs
  // histogram and fit options
  int rebinfactor[NCH]={5, 5, 5, 5}; // rebin histograms
  double fitbeginch[NCH]={1.0,1.0,1.0,1.0}; // fit start locations
  double fitendch[NCH] = {60,60,60,60}; // fit end locations

  // the function to be used to do fit
  gStyle->SetOptFit(1111);
  TF1* Fideal = new TF1("Fideal",IdealResponse, 0, 500, 4);
  Fideal->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
  Fideal->SetLineColor(2); Fideal->SetLineStyle(1);
  double par[4];
  double parerr[4];

  // process 3 files in a batch
  string rtfilenames[3];
  string strchimney = pmtRow + "_PMT_";
  string strpmt = to_string(pmt1) + "_" + to_string(pmt2) + "_" + to_string(pmt3) + "_" + to_string(pmt4) + "_";
  //string voltagestr[3]={"1440","1470","1500"};
  string voltagestr[3] = {to_string(volt1), to_string(volt2), to_string(volt3)};
  string ledstr;
  if(led)
    ledstr = "On";
  else
    ledstr = "Off";
  for(int i=0; i<3; i++){
    rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_Led" + ledstr + "_result.root";
    cout << rtfilenames[i] << endl;
  }

  string outnameroot = strchimney + strpmt + "gain.root";
  string outnametxt = strchimney + strpmt + "gain_fit.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TH1F* hCharge[3]; // histograms for each canvas
  TCanvas* c[3];
  char tempname[100];

  TFile* files[3];
  // Read each of the three data ROOT files and save them to files[]
  for(int i = 0; i < 3; i++){
    files[i] = new TFile(rtfilenames[i].c_str(),"read");
  }

  // Generate 4 canvases and plot the PMT histograms on them
  for(int i = 0; i < 4; i++){
    sprintf(tempname, "c_%d",i);
    string canvasTitle = strchimney + to_string(i+1);
    //sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
    c[i] = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
    c[i]->Divide(3); // divide canvas into 3 pads along the width

    // Plot the histograms
    for(int j = 0; j<3; j++){
      c[i]->cd(j+1); // switch pads
      sprintf(tempname,"Results/FinalCharge_%d",i); // store string "Results/FinalCharge_%d" in tempname
      hCharge[j] = (TH1F*)files[j]->Get(tempname); // read histogram data from the ROOT file
      hCharge[j]->SetTitle((voltagestr[j] + "V").c_str());
      hCharge[j]->Rebin(rebinfactor[i]);
      hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
      hCharge[j]->Draw();

      // Set initial fit parameters
      //Fideal->SetParameter(1,1.6);
      //Fideal->SetParameter(2,1.6*0.4);
      //Fideal->SetParameter(3,hCharge[j]->Integral());
      Fideal->SetParLimits(0,0.1,100);
      Fideal->SetParLimits(1,0.5,10);
      Fideal->SetParLimits(2,0.1,10);
      Fideal->SetParLimits(3,0.1,20000);

      Double_t hist_mean = truncatedMean(hCharge[j],10);
      Fideal->SetParameter(1,hist_mean);

      // Iteratively fit more than once
      for(int k=0; k<2;k++){
        hCharge[j]->Fit("Fideal","RQ","",fitbeginch[i],fitendch[i]); // Fit the histogram
        Fideal->GetParameters(par);
        Fideal->SetParameters(par); // Set fit parameters for next iteration
      }
      
      //hCharge[j]->Fit("expo","","",fitbeginch[i],fitendch[i]);
      hCharge[j]->Fit("Fideal","","",fitbeginch[i],fitendch[i]);
      
      hCharge[j]->GetXaxis()->SetRangeUser(0, fitendch[i]); // Set axes

      // write parameters to output txt file
      Fideal->GetParameters(par);
      //parerr = Fideal->GetParErrors();
      foutFit<<"voltage\t"<<voltagestr[j]<<"\tchID\t"<<i
             <<"\t"<<par[0]<<"\t"<<Fideal->GetParError(0)
             <<"\t"<<par[1]<<"\t"<<Fideal->GetParError(1)
             <<"\t"<<par[2]<<"\t"<<Fideal->GetParError(2)
             <<"\t"<<Fideal->GetChisquare()
             <<"\t"<<Fideal->GetNDF()
             <<"\t"<<Fideal->GetProb()
             <<endl;
    }

    // write results to output ROOT file
    outROOTfile->cd();
    c[i]->Write();
  }

  // close output ROOT file
  outROOTfile->Close();
}


double IdealResponse(double *x,double *par){
    double mu = par[0];
    double q = par[1];
    double sigma = par[2];
    double amplitude = par[3];
    double sum=0;
    for(Int_t n=1; n<50; n++){
        sum += TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)*TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*PI*n));
    }
    return sum*amplitude;
}

Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3){
  Double_t mean = 0.;
  Double_t stddev = 0.;

  // Zoom out
  hist->GetXaxis()->UnZoom();

  // Calculate truncated mean
  for(int i = 0; i < n_iterations; i++){
    mean = hist->GetMean(1);
    stddev = hist->GetMeanError();

    cout << mean;

    // Truncate
    Double_t new_start = mean - n_rejection_stddevs * stddev;
    Double_t new_end = mean + n_rejection_stddevs * stddev;

    hist->GetXaxis()->SetRangeUser(new_start, new_end);
  }

  // Zoom histogram back out
  hist->GetXaxis()->UnZoom();

  return mean;
}

