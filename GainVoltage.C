//*******************************************************************
// Plots gain versus voltage and performs a power log fit on data
//*******************************************************************

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

double power(double* x, double* par);

void GainVoltage(){
  //Variables
  double fbegin = 1000;   // begin fit
  double fend = 2000;     // end fit
  int NPAR = 2;        // number of parameters
  
  // Define power law fit function
  TF1 * fit = new TF1("fit",power,fbegin,fend,NPAR);
  double par[NPAR];
  double parerr[NPAR];
  fit->SetParNames("Amplitude Constant", "Exponent");
  fit->SetLineColor(2);
  fit->SetLineStyle(1);
  
  // Create canvas
  auto c1 = new TCanvas("c1","c1",200,10,600,400);
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);

  // Gathers data from file: a .txt file with format:
  // <voltage>      <gain>      <gain error>
  TGraph* data = new TGraphErrors("./input.txt","%lg %lg %lg");
  data->SetTitle("PMT Gain versus Voltage;"
		 "Voltage [V];"
		 "Gain");
  data->SetMarkerStyle(kCircle);
  // Fit
  fit->SetParameters(0.003,0.33);
  data->Fit("fit");
  gStyle->SetOptFit();
  // fit->GetParameters(par);
  // cout<<"\t"<<par[0]<<"\t"<<fit->GetParError(0)
  //     <<"\t"<<par[1]<<"\t"<<fit->GetParError(1)
  //     <<"\t"<<par[2]<<"\t"<<fit->GetParError(2)
  //     <<"\t"<<fit->GetChisquare()
  //     <<"\t"<<fit->GetNDF()
  //     <<"\t"<<fit->GetProb()
  //     <<endl;

  // Plot gain and voltage
  data->Draw("ap");
}

// Power law function definition: a * (x ^ k)
double power(double* x, double* par){
  double k = par[0];
  double a = par[1];
  return k*TMath::Power(x[0],a);
}
