//*******************************************************************
// Plots gain versus voltage and performs a power log fit on data
//
// Input must be CHIMNEY.txt
// Input format must be PMT# Voltage Gain GainError
//*******************************************************************

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

#include <fstream>

double power(double* x, double* par);

void GainVoltage(string chimney, Int_t pmt_num){
  //Variables
  double fbegin = 1000;   // begin fit
  double fend = 2000;     // end fit
  int NPAR = 2;        // number of parameters
  
  // read gain and voltage values for the PMT in question
  string input_file_name = chimney + ".txt";
  std::ifstream input_file(input_file_name);

  Double_t voltage_raw[6];
  Double_t gain_raw[6];
  Double_t gain_error_raw[6];

  // read line by line
  double p, v, g, ge;
  int num_data_points = 0;

  while (input_file >> p >> v >> g >> ge){
    if (p == pmt_num){
      voltage_raw[num_data_points] = v;
      gain_raw[num_data_points] = g;
      gain_error_raw[num_data_points] = ge;
      num_data_points++;
    }
  }

  // create properly sized arrays
  Double_t voltage[num_data_points];
  Double_t gain[num_data_points];
  Double_t gain_error[num_data_points];

  for(int i = 0; i < num_data_points; i++){
    voltage[i] = voltage_raw[i];
    gain[i] = gain_raw[i];
    gain_error[i] = gain_error_raw[i];
  }

  if(num_data_points != 3 && num_data_points != 6){
    cout << "Improper number of data points";
    return;
  }

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

  // Create graph of data
  TGraph* data = new TGraphErrors(num_data_points, voltage, gain, 0, gain_error);
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
