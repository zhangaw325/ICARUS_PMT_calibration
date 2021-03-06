//*******************************************************************
// Combined fitting code adapted from ROOT tutorials: combinedFit.C
// Result root files saved as PDF
//*******************************************************************

#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"

//#include "Fit_Spe.C" // the fit function is defined in this file

// constants
double PI = TMath::Pi();
// define which parameters correspond to which indices
const int NPAR_i = 4;
const int NPAR = NPAR_i*3-2;
int ipar1[NPAR_i] = {  0, // common mu
                  1, // q for hist 1
                  2, // sigma for hist 1
                  3, // amplitude for hist 1
};

int ipar2[NPAR_i] = {  0, // common mu
                  1, // q for hist 2
                  2, // sigma for hist 2
                  3  // amplitude for hist 2
};

int ipar3[NPAR_i] = {  0,  // common mu
                  1, // q for hist 3
                  2, // sigma for hist 3
                  3 // amplitude for hist 3
};

// function declarations
double IdealResponse(double *x,double *par);
Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3);

/*
// create combined chi2 structure
struct GlobalChi2 {
  
  GlobalChi2( ROOT::Math::IMultiGenFunction & f1,
              ROOT::Math::IMultiGenFunction & f2,
              ROOT::Math::IMultiGenFunction & f3) : 
    fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

  double operator() (const double *par) const {
    double p1[NPAR_i];
    double p2[NPAR_i];
    double p3[NPAR_i];

    for (int i = 0; i < NPAR_i; ++i){
      p1[i] = par[ipar1[i]];
      p2[i] = par[ipar2[i]];
      p3[i] = par[ipar3[i]];
    }

    return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
  }

  const  ROOT::Math::IMultiGenFunction * fChi2_1;
  const  ROOT::Math::IMultiGenFunction * fChi2_2;
  const  ROOT::Math::IMultiGenFunction * fChi2_3;
};
*/


void FitChargeDistributions_lowcharge(string pmtRow,
			    char pmt1, char pmt2, char pmt3, char pmt4,
			    int volt1, int volt2, int volt3, bool led){

  const int NCH = 4; // 4 PMTs
  // histogram and fit options
  int rbf_0 = 4;
  double fbc_0 = 2.0;
  double fec_0 = 90.0;
  int rebinfactor[NCH]={rbf_0, rbf_0, rbf_0, rbf_0}; // rebin histograms
  double fitbeginch[NCH]={fbc_0, fbc_0, fbc_0, fbc_0}; // fit start locations
  double fitendch[NCH] = {fec_0, fec_0, fec_0, fec_0}; // fit end locations
  string initparam[4];

  // the function to be used to do fit
  gStyle->SetOptFit(1111);
  
  /*
  TF1* Fideal = new TF1("Fideal",IdealResponse, 0, 500, 6);
  Fideal->SetParNames("meanNpe","spePeak","speWidth","Amplitude","expAmp","expCoeff");
  Fideal->SetLineColor(2); Fideal->SetLineStyle(1);
  */
  double par[NPAR_i];
  double parerr[NPAR_i];

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
  
  for(int i=0; i<3; i++){
    rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_Led" + ledstr + "_result.root";
    cout << rtfilenames[i] << endl;
  }
  
  for(int i=0;i<4; i++){
    resultnames[i]  = strchimney + strpmt + "CH" + channelnames[i] + "_" + "Led" + ledstr + "_lowcharge.pdf";
    cout << resultnames[i] <<endl;
  }

  string outnameroot = strchimney + strpmt + "_lowcharge_gain.root";
  string outnametxt = strchimney + strpmt + "_lowcharge_gain_fit.txt";
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

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);

  // Print header
  foutFit << "**************************** PARAMETER VALUES ****************************" << endl;
  
  // Generate 4 canvases and plot the PMT histograms on them
  for(int i = 0; i < 4; i++){
    sprintf(tempname, "c_%d",i);
    string canvasTitle = strchimney + to_string(i+1);
    //sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
    c[i] = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
    c[i]->Divide(3); // divide canvas into 3 pads along the width

    //*************************
    // Begin fit
    //*************************

    // generate histograms
    sprintf(tempname, "Results/FinalCharge_%d",i);
    for(int j = 0; j < 3; j++){
      hCharge[j] = (TH1F*)files[j]->Get(tempname);
      hCharge[j]->SetTitle((voltagestr[j] + "V").c_str());
      hCharge[j]->Rebin(rebinfactor[i]);
      hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
    }

    // generate fit functions
    TF1* fit_ideal_1 = new TF1("fit_ideal_1",IdealResponse, fitbeginch[i], fitendch[i], NPAR_i);
    fit_ideal_1->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
    fit_ideal_1->SetLineColor(2);
    fit_ideal_1->SetLineStyle(1);
    TF1* fit_ideal_2 = new TF1("fit_ideal_2",IdealResponse, fitbeginch[i], fitendch[i], NPAR_i);
    fit_ideal_2->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
    fit_ideal_2->SetLineColor(2);
    fit_ideal_2->SetLineStyle(1);
    TF1* fit_ideal_3 = new TF1("fit_ideal_3",IdealResponse, fitbeginch[i], fitendch[i], NPAR_i);
    fit_ideal_3->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
    fit_ideal_3->SetLineColor(2);
    fit_ideal_3->SetLineStyle(1);

    // initial parameters
    Double_t q_0 = 1.6;
    Double_t sigma_0 = 1.6*0.4;
    Double_t expAmp_0 = 2.0;
    Double_t expCoeff_0 = 2.0;

    Double_t hist_mean_1 = truncatedMean(hCharge[1], 10);

    // set ranges on fit parameters
    fit_ideal_1->SetParameters(hist_mean_1,q_0,sigma_0,hCharge[0]->Integral()/2);
    fit_ideal_2->SetParameters(hist_mean_1,q_0,sigma_0,hCharge[0]->Integral()/2);
    fit_ideal_3->SetParameters(hist_mean_1,q_0,sigma_0,hCharge[0]->Integral()/2); 

    // mu
    fit_ideal_1->SetParLimits(0,5,20);
    fit_ideal_2->SetParLimits(0,5,20);
    fit_ideal_3->SetParLimits(0,5,20);

    // q
    fit_ideal_1->SetParLimits(1,0.01, 10);
    fit_ideal_2->SetParLimits(1,0.01, 10);
    fit_ideal_3->SetParLimits(1,0.01, 10);
    
    // sigma
    fit_ideal_1->SetParLimits(2,0.1, 3.1);
    fit_ideal_2->SetParLimits(2,0.1, 3.1);
    fit_ideal_3->SetParLimits(2,0.1, 3.1);
    
    // amplitude
    fit_ideal_1->SetParLimits(3,0.01, 20000);
    fit_ideal_2->SetParLimits(3,0.01, 20000);
    fit_ideal_3->SetParLimits(3,0.01, 20000);
    /*
    for(int j = 0; j < 3; j++){
      if (j == 0)
	hCharge[j]->Fit("fit_ideal_1");
      else if (j == 1)
	   hCharge[j]->Fit("fit_ideal_2");
      else
	hCharge[j]->Fit("fit_ideal_3");
    }
    */
    //hCharge[0]->Fit("fit_ideal_1");
    //hCharge[1]->Fit("fit_ideal_2");
    //hCharge[2]->Fit("fit_ideal_3");
    
    for(int j = 0; j < 3; j++){
      c[i]->cd(j+1);

       if (j == 0)
	hCharge[j]->Fit("fit_ideal_1");
      else if (j == 1)
	   hCharge[j]->Fit("fit_ideal_2");
      else
	hCharge[j]->Fit("fit_ideal_3");

      hCharge[j]->GetXaxis()->SetRangeUser(0, fitendch[i]);
      hCharge[j]->GetYaxis()->SetRangeUser(0, 400);
      gStyle->SetOptFit();

      
      hCharge[j]->Draw();
    }

    // write parameters to output txt file
    initparam[i]="chID\t"+to_string(i)+"\t"       //channel id
      +fitbeginch[i]+"\t"         //start fit
      +fitendch[i]+"\t"           //end fit
      +rebinfactor[i]+"\t"        //rebin factor
      +hist_mean_1+"\t"           //mu
      +q_0+"\t"                   //q
      +sigma_0+"\t"               //sigma
      +hCharge[0]->Integral()/2;  //amplitude
    
    // Voltage 1
    fit_ideal_1->GetParameters(par);
    foutFit<<"voltage\t"<<voltagestr[0]<<"\tchID\t"<<i
	   <<"\t"<<par[0]<<"\t"<<fit_ideal_1->GetParError(0)
	   <<"\t"<<par[1]<<"\t"<<fit_ideal_1->GetParError(1)
	   <<"\t"<<par[2]<<"\t"<<fit_ideal_1->GetParError(2)
	   <<"\t"<<par[3]<<"\t"<<fit_ideal_1->GetParError(3)
	   <<"\t"<<fit_ideal_1->GetChisquare()
	   <<"\t"<<fit_ideal_1->GetNDF()
	   <<"\t"<<fit_ideal_1->GetProb()
	   <<endl;
    
    // Voltage 2
    fit_ideal_2->GetParameters(par);
    foutFit<<"voltage\t"<<voltagestr[1]<<"\tchID\t"<<i
	   <<"\t"<<par[0]<<"\t"<<fit_ideal_2->GetParError(0)
	   <<"\t"<<par[1]<<"\t"<<fit_ideal_2->GetParError(1)
	   <<"\t"<<par[2]<<"\t"<<fit_ideal_2->GetParError(2)
	   <<"\t"<<par[3]<<"\t"<<fit_ideal_2->GetParError(3) 
	   <<"\t"<<fit_ideal_2->GetChisquare()
	   <<"\t"<<fit_ideal_2->GetNDF()
	   <<"\t"<<fit_ideal_2->GetProb()
	   <<endl;

    // Voltage 3
    fit_ideal_3->GetParameters(par);
    foutFit<<"voltage\t"<<voltagestr[2]<<"\tchID\t"<<i
	   <<"\t"<<par[0]<<"\t"<<fit_ideal_3->GetParError(0)
	   <<"\t"<<par[1]<<"\t"<<fit_ideal_3->GetParError(1)
	   <<"\t"<<par[2]<<"\t"<<fit_ideal_3->GetParError(2)
	   <<"\t"<<par[3]<<"\t"<<fit_ideal_3->GetParError(3)
	   <<"\t"<<fit_ideal_3->GetChisquare()
	   <<"\t"<<fit_ideal_3->GetNDF()
	   <<"\t"<<fit_ideal_3->GetProb()
	   <<endl;

    //*************************
    // End fit
    //*************************

    // write results to output ROOT file
    outROOTfile->cd();
    c[i]->Write();
    c[i]->Print(resultnames[i].c_str(),"pdf");
  }

  // print initial parameters at end of root file
  foutFit << "**************************** INITIAL PARAMETERS ****************************" << endl;
  for(int i = 0; i < 4; i++){
    foutFit << initparam[i] << endl;
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
    sum += (TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)*TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*PI*n)));
  }
  return amplitude * sum;
}

Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3){
  Double_t mean = 0;
  Double_t stddev = 0.;

  // Zoom out
  hist->GetXaxis()->UnZoom();

  // Calculate truncated mean
  for(int i = 0; i < n_iterations; i++){
    mean = hist->GetMean(1);
    stddev = hist->GetStdDev(1);

    //cout << mean << endl;

    // Truncate
    Double_t new_start = mean - n_rejection_stddevs * stddev;
    Double_t new_end = mean + n_rejection_stddevs * stddev;

    hist->GetXaxis()->SetRangeUser(new_start, new_end);
  }

  // Zoom histogram back out
  hist->GetXaxis()->UnZoom();

  return mean;
}
