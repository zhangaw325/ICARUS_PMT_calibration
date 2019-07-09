#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"

#include <string>
#include <iostream>

//#include "Fit_Spe.C" // the fit function is defined in this file

double PI = TMath::Pi();

double IdealResponse(double *x,double *par);
double RealResponse(double *x, double *par);
Double_t iGnE(Double_t x, Double_t q0, Double_t sigma0, 
              Double_t q1, Double_t sigma1, Double_t a, Double_t n);
Double_t Gn(Double_t x, Double_t q0, Double_t sigma0, Double_t q1, 
            Double_t sigma1, Double_t n);
Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3);

void FitChargeDistributions_background(std::string pmtRow,
			    char pmt1, char pmt2, char pmt3, char pmt4,
			    int volt1, int volt2, int volt3, bool led){

  const int NCH = 4; // 4 PMTs
  const int NPAR = 8;
  // histogram and fit options
  int rebinfactor[NCH]={1, 1, 1, 1}; // rebin histograms
  double fitbeginch[NCH]={-0.75,-0.75,-0.75,-0.75}; // fit start locations
  double fitendch[NCH] = {30,30,30,30}; // fit end locations

  // the function to be used to do fit
  gStyle->SetOptFit(1111);
  TF1* Freal = new TF1("Freal",RealResponse, 0, 500, NPAR);
  Freal->SetParNames("meanNpe","pedestalPeak","pedestalSigma","spePeak","speSigma","w","alpha","amplitude");
  Freal->SetLineColor(2); Freal->SetLineStyle(1);
  double par[NPAR];
  double parerr[NPAR];

  // process 3 files in a batch
  std::string rtfilenames[3];
  std::string strchimney = pmtRow + "_PMT_";
  std::string strpmt = std::to_string(pmt1) + "_" + std::to_string(pmt2) + "_" + std::to_string(pmt3) + "_" + std::to_string(pmt4) + "_";
  //std::string voltagestr[3]={"1440","1470","1500"};
  std::string voltagestr[3] = {std::to_string(volt1), std::to_string(volt2), std::to_string(volt3)};
  std::string ledstr;
  if(led)
    ledstr = "On";
  else
    ledstr = "Off";
  for(int i=0; i<3; i++){
    rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_Led" + ledstr + "_result.root";
    std::cout << rtfilenames[i] << std::endl;
  }

  std::string outnameroot = strchimney + strpmt + "gain.root";
  std::string outnametxt = strchimney + strpmt + "gain_fit.txt";
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

  // set call limits
  //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);

  // Generate 4 canvases and plot the PMT histograms on them
  for(int i = 0; i < 4; i++){
    sprintf(tempname, "c_%d",i);
    std::string canvasTitle = strchimney + std::to_string(i+1);
    //sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
    c[i] = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
    c[i]->Divide(3); // divide canvas into 3 pads along the width

    // Plot the histograms
    for(int j = 0; j<3; j++){
      c[i]->cd(j+1); // switch pads
      sprintf(tempname,"Results/FinalCharge_%d",i); // store std::string "Results/FinalCharge_%d" in tempname
      hCharge[j] = (TH1F*)files[j]->Get(tempname); // read histogram data from the ROOT file
      hCharge[j]->SetTitle((voltagestr[j] + "V").c_str());
      hCharge[j]->Rebin(rebinfactor[i]);
      hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
      hCharge[j]->Draw();

      // Set initial fit parameters
      //Freal->SetParameter(1,1.6);
      //Freal->SetParameter(2,1.6*0.4);
      //Freal->SetParameter(3,hCharge[j]->Integral());
      Freal->SetParLimits(0,1,30);
      Freal->SetParLimits(1,-1,1);
      Freal->SetParLimits(2,0.0001,1);
      Freal->SetParLimits(3,0.1, 4);
      Freal->SetParLimits(4, 0.1, 4);
      Freal->SetParLimits(5, 0.00001, 1);
      Freal->SetParLimits(6, 0.00000001, 1);
      Freal->SetParLimits(7, 1, 100000);

      Double_t hist_mean = truncatedMean(hCharge[j],10);
      Freal->SetParameter(0,hist_mean);
      Freal->SetParameters(hist_mean, 0.005, 0.08, 1.2, 1, 0.5, 0.003, 2000);

      // Iteratively fit more than once
      // for(int k=0; k<2;k++){
      //   hCharge[j]->Fit("Freal","RQ","",fitbeginch[i],fitendch[i]); // Fit the histogram
      //   Freal->GetParameters(par);
      //   Freal->SetParameters(par); // Set fit parameters for next iteration
      // }
      
      //hCharge[j]->Fit("expo","","",fitbeginch[i],fitendch[i]);
      hCharge[j]->Fit("Freal","LM","",fitbeginch[i],fitendch[i]);
      hCharge[j]->GetXaxis()->SetRangeUser(fitbeginch[i], fitendch[i]); // Set axes

      // write parameters to output txt file
      Freal->GetParameters(par);
      //parerr = Freal->GetParErrors();
      //foutFit<<"voltage\t"<<voltagestr[j]<<"\tchID\t"<<i
            // <<"\t"<<par[0]<<"\t"<<Freal->GetParError(0)
            // <<"\t"<<par[1]<<"\t"<<Freal->GetParError(1)
            // <<"\t"<<par[2]<<"\t"<<Freal->GetParError(2)
            // <<"\t"<<Freal->GetChisquare()
            // <<"\t"<<Freal->GetNDF()
            // <<"\t"<<Freal->GetProb()
            // <<std::endl;
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

double RealResponse(double *x, double *par){
  double mu = par[0];
  double q0 = par[1];
  double sigma0 = par[2];
  double q1 = par[3];
  double sigma1 = par[4];
  double w = par[5];
  double a = par[6];
  double amplitude = par[7];

  Double_t sum = 0;

  for(Int_t n = 0; n < 50; n++){
    sum += (TMath::Power(mu, n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n))
          * ((1 - w)*Gn(x[0], q0, sigma0, q1, sigma1, n) 
          + w*iGnE(x[0], q0, sigma0, q1, sigma1, a, n));
  }

  return amplitude*sum;
}

Double_t iGnE(Double_t x, Double_t q0, Double_t sigma0, 
              Double_t q1, Double_t sigma1, Double_t a, Double_t n){

  Double_t qn = q0 + n*q1;
  Double_t sigman = TMath::Sqrt(sigma0*sigma0 + n*sigma1*sigma1);
  Double_t sign = 1;
  if(x - qn - sigman*sigman*a > 0){
    sign = 1;
  } else if(x - qn - sigman*sigman*a == 0){
    sign = 0;
  } else {
    sign = -1;
  }

  Double_t iGnE = (a/2)*TMath::Exp(-1.0*a*(x - qn - a*sigman*sigman))
            *(TMath::Erf(TMath::Abs(q0 - qn - sigman*sigman*a)/(sigman*TMath::Sqrt2()))
            + sign*TMath::Erf(TMath::Abs(x - qn - sigman*sigman*a)/(sigman*TMath::Sqrt2())));

  return iGnE;
}

Double_t Gn(Double_t x, Double_t q0, Double_t sigma0, Double_t q1, Double_t sigma1, Double_t n){
  Double_t Gn = 0;
  if(n == 0){
    Gn = TMath::Exp(-1.0*(x - q0)*(x - q0)/(2*sigma0*sigma0))
           /(sigma0*TMath::Sqrt(2*TMath::Pi()));
  } else {
    Gn = TMath::Exp(-1.0*(x - n*q1)*(x - n*q1)/(2*n*sigma1*sigma1))
            /(sigma1*TMath::Sqrt(2*TMath::Pi()*n));
  }

  return Gn;
}

Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs){
  Double_t mean = 0.;
  Double_t stddev = 0.;

  // Zoom out
  hist->GetXaxis()->UnZoom();

  // Calculate truncated mean
  for(int i = 0; i < n_iterations; i++){
    mean = hist->GetMean(1);
    stddev = hist->GetMeanError();

    // std::cout << mean;

    // Truncate
    Double_t new_start = mean - n_rejection_stddevs * stddev;
    Double_t new_end = mean + n_rejection_stddevs * stddev;

    hist->GetXaxis()->SetRangeUser(new_start, new_end);
  }

  // Zoom histogram back out
  hist->GetXaxis()->UnZoom();

  return mean;
}

