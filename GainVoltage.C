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

const int NPMT = 10;

void GainVoltage(string chimney){
  // Variables
  double fbegin = 1000;   // begin fit
  double fend = 2000;     // end fit
  int NPAR = 2;        // number of parameters
  
  // output files
  string outnametxt[10];
  string outnamepdf[10];

  for(int i = 0; i < NPMT; i++){
    outnametxt[i] = chimney + "_" + to_string(i + 1) + "_gainvsvoltage.txt";
    outnamepdf[i] = chimney + "_" + to_string(i + 1) + "_gainvsvoltage.pdf";
  }

  string outnameroot = chimney + "_gainvsvoltage.root";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream fout(outnametxt.c_str(),ios::out);

  for(int pmt_num = 0; pmt_num < NPMT; pmt_num++){

    // read gain and voltage values for the PMT in question
    string input_file_name = chimney + ".txt";
    std::ifstream input_file(input_file_name);

    Double_t voltage_raw[6];
    Double_t voltage_error_raw[6] = {2,2,2,2,2,2};
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
    Double_t voltage_error[num_data_points];
    Double_t gain[num_data_points];
    Double_t gain_error[num_data_points];

    for(int i = 0; i < num_data_points; i++){
      voltage[i] = TMath::Log(voltage_raw[i]);
      gain[i] = TMath::Log(gain_raw[i]);
      gain_error[i] = gain_error_raw[i]/gain_raw[i];
      voltage_error[i] = voltage_error_raw[i]/voltage_raw[i];
    }

    if(num_data_points != 3 && num_data_points != 6){
      cout << "Improper number of data points";
      return;
    }

    // Define power law fit function

    // Perform linear fit on log-log plot
    TF1 *fit = new TF1("fit", "pol1", fbegin, fend);
    double par[NPAR];
    double parerr[NPAR];
    fit->SetParNames("Constant", "Exponent");
    fit->SetLineColor(2);
    fit->SetLineStyle(1);
    
    // Create canvas
    auto c1 = new TCanvas("c1","c1",200,10,600,400);
    c1->SetGrid();
    c1->GetFrame()->SetBorderSize(12);

    // Create graph of data
    TGraph* data = new TGraphErrors(num_data_points, voltage, gain, voltage_error, gain_error);
    string title =  "PMT " + chimney + "_" + to_string(pmt_num) + " gain vs voltage;"
                    + "log(voltage [V]);"
                    + "log(gain)";
    data->SetTitle(title.c_str());
    data->SetMarkerStyle(kCircle);
    // Fit
    fit->SetParameters(0.003,0.33);
    data->Fit("fit");
    gStyle->SetOptFit();
    fit->GetParameters(par);
    fout<<"\t"<<par[0]<<"\t"<<fit->GetParError(0)
         <<"\t"<<par[1]<<"\t"<<fit->GetParError(1)
         <<"\t"<<par[2]<<"\t"<<fit->GetParError(2)
         <<"\t"<<fit->GetChisquare()
         <<"\t"<<fit->GetNDF()
         <<"\t"<<fit->GetProb()
         <<endl;

    // Plot gain and voltage
    data->Draw("ap");

    outROOTfile->cd();
    c1->Write();
    c1->Print(outnamepdf.c_str(),"pdf");
  }
}

// Power law function definition: a * (x ^ k)
double power(double* x, double* par){
  double k = par[0];
  double a = par[1];
  return k*TMath::Power(x[0],a);
}
