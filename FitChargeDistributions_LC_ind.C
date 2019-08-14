//*********************************************************************************************/
// - For LOW-CHARGE distributions
//     - Fits equation to each distribution (each PMT, each voltage)
// - Processes a single PMTs at a time (contrast to FitChargeDistributions_LC_batch.C)
// - Output:
//     - Root file
//     - .pdf of the canvas contents
//     - .txt of comma separated values including:
//          - final parameters: voltage, channel number
//                              mu, mu error, q, q error, sigma, sigma error,
//                              amplitude, amplitude error, chi2, ndf, fit probability
//          - initial parameters: channel id,
//                                start fit, end fit, rebin factor,
//                                mu, q, sigma, amplitude
// - Uses parameter limits, which results in very poor errors
//     - To improve the errors, save the initial parameters and rerun using
//       FitChargeDistributions_InitParam.C
// - The *LedOn*result.root files must be included in the same folder that this macro is stored in
// - Do not have any *LedOff*result.root files in the folder
// - Processes three root files in a batch (same four PMTs at three different voltages)
// - Example function call in root:
//       .x FitChargeDistributions_LC_ind.C("A10", 5, 6, 7, 8, 1400, 1430, 1460, 0)
//       - Last parameter is the index of the desired PMT out of the four, ie. 0, 1, 2, or 3
//       - Voltages as parameters must match the nominal values matching the *result.root files
//**********************************************************************************************/

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

void FitChargeDistributions_LC_ind(string pmtRow,
					     char pmt1, char pmt2, char pmt3, char pmt4,
					     int volt1, int volt2, int volt3, int pmt){

  const int NCH = 4; // 4 PMTs
  // histogram and fit options
  int rbf_0 = 2;
  double fbc_0 = 0.5;
  double fec_0 = 40.0;
  gStyle->SetOptFit(1111); // formatting
  //******* CHANGE VALUES HERE TO BE WRITTEN TO OUTPUT FILE *******
  // rebin histogram
  int rebinfactor[NCH][3]={{rbf_0, rbf_0, rbf_0}, // pmt 0
			   {rbf_0, rbf_0, rbf_0}, // pmt 1
			   {rbf_0, rbf_0, rbf_0}, // pmt 2
			   {rbf_0, rbf_0, rbf_0}};// pmt 3
  // fit start locations
  double fitbeginch[NCH][3]={{fbc_0, fbc_0, fbc_0}, // pmt 0
			     {fbc_0, fbc_0, fbc_0}, // pmt 1
			     {fbc_0, fbc_0, fbc_0}, // pmt 2
			     {fbc_0, fbc_0, fbc_0}};// pmt 3
  // fit end locations
  double fitendch[NCH][3] = {{fec_0, fec_0, fec_0}, // pmt 0
			     {fec_0, fec_0, fec_0}, // pmt 1
			     {fec_0, fec_0, fec_0}, // pmt 2
			     {fec_0, fec_0, fec_0}};// pmt 3
  // ***************************************************************
    
  // arrays storing parameter information for a particular function
  string initparam[NCH][3]; // for outputting initial parameters
  double par[NPAR_i];
  double parerr[NPAR_i];

  // file handling
  string rtfilenames[3];
  string pdfname[4];
  double channelnames[4] = {0,1,2,3};
  string strchimney = pmtRow + "_PMT_";
  string strpmt = to_string(pmt1) + "_" + to_string(pmt2) + "_" + to_string(pmt3) + "_" + to_string(pmt4) + "_";
  string voltagestr[3] = {to_string(volt1), to_string(volt2), to_string(volt3)};
  for(int i=0; i<3; i++){
    rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_LedOn_result.root";
    cout << rtfilenames[i] << endl;
  }
  pdfname[pmt]  = strchimney + strpmt + "CH" + channelnames[pmt] + "_"
    + "LedOn_lowcharge.pdf";
  cout << pdfname[pmt] <<endl;
  string outnameroot = strchimney + strpmt + "lowcharge_CH"+to_string(pmt)+"_gain.root";
  string outnametxt = strchimney + strpmt + "lowcharge_CH"+to_string(pmt)+"_gain_fit.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TFile* files[3];
  // Read each of the three data ROOT files and save them to files[]
  for(int i = 0; i < 3; i++){
    files[i] = new TFile(rtfilenames[i].c_str(),"read");
  }

  // canvas handling
  TH1F* hCharge[3]; // histograms for each canvas
  TCanvas* c[4];
  char tempname[100];
  
  // minimize a max of 5000 times
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);

  // Print header
  foutFit << "**************************** PARAMETER VALUES ****************************" << endl;

  //*************************
  // Begin fit
  //*************************

  // generate histograms
  sprintf(tempname, "Results/FinalCharge_%d",pmt);
  for(int j = 0; j < 3; j++){
    hCharge[j] = (TH1F*)files[j]->Get(tempname);
    hCharge[j]->SetTitle((voltagestr[j] + "V").c_str());
    hCharge[j]->Rebin(rebinfactor[pmt][j]);
    hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
  }
  
  // generate fit functions
  TF1* fit_ideal_1 = new TF1("fit_ideal_1",IdealResponse, fitbeginch[pmt][0], fitendch[pmt][0], NPAR_i);
  fit_ideal_1->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
  fit_ideal_1->SetLineColor(2);
  fit_ideal_1->SetLineStyle(1);
  TF1* fit_ideal_2 = new TF1("fit_ideal_2",IdealResponse, fitbeginch[pmt][1], fitendch[pmt][1], NPAR_i);
  fit_ideal_2->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
  fit_ideal_2->SetLineColor(2);
  fit_ideal_2->SetLineStyle(1);
  TF1* fit_ideal_3 = new TF1("fit_ideal_3",IdealResponse, fitbeginch[pmt][2], fitendch[pmt][2], NPAR_i);
  fit_ideal_3->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
  fit_ideal_3->SetLineColor(2);
  fit_ideal_3->SetLineStyle(1);
  
  // default parameters
  //Double_t mu_def = hist_mean_1;
  Double_t q_def = 1.6;
  Double_t sigma_def = 1.6*0.4;

  // initial parameters for each voltage
  // ****** CHANGE THESE ARRAY VALUES TO RECORD TO OUTPUT *******
  Double_t mu_0[3], q_0[3], sigma_0[3], amp_0[3];
  for(int k = 0; k < 3; k++){
    mu_0[k] = truncatedMean(hCharge[k],10);
    q_0[k] = q_def;
    sigma_0[k] = sigma_def;
    amp_0[k] = hCharge[k]->Integral()/2;
  
  // write initial parameters to output txt file
  initparam[pmt][k]=
    "voltage\t"+voltagestr[k]+"\t"
    +"chID\t"+to_string(pmt)+"\t"       //channel id
    +fitbeginch[pmt][k]+"\t"          //start fit
    +fitendch[pmt][k]+"\t"            //end fit
    +rebinfactor[pmt][k]+"\t"         //rebin factor
    +mu_0[k]+"\t"                  //mu
    +q_0[k]+"\t"                   //q
    +sigma_0[k]+"\t"               //sigma
    +amp_0[k];                     //amplitude
  }
  
  // set ranges on fit parameters
  fit_ideal_1->SetParameters(mu_0[0],q_0[0],sigma_0[0],amp_0[0]); 
  fit_ideal_2->SetParameters(mu_0[1],q_0[1],sigma_0[1],amp_0[1]); 
  fit_ideal_3->SetParameters(mu_0[2],q_0[2],sigma_0[2],amp_0[2]); 

  // mu
  fit_ideal_1->SetParLimits(0,1,30);
  fit_ideal_2->SetParLimits(0,1,30);
  fit_ideal_3->SetParLimits(0,1,30);

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
 
  // Generate canvas and plot the PMT histograms on them
  sprintf(tempname, "c_%d",pmt);
  string canvasTitle = strchimney + to_string(pmt);
  //sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
  c[pmt] = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
  c[pmt]->Divide(3); // divide canvas into 3 pads along the width
  for(int j = 0; j < 3; j++){
    c[pmt]->cd(j+1);
    hCharge[j]->GetXaxis()->SetRangeUser(0, fitendch[pmt][j]);
    hCharge[j]->GetYaxis()->SetRangeUser(0, 400);
    gStyle->SetOptFit();
    if(j==0)
      hCharge[0]->Fit("fit_ideal_1","","",fitbeginch[pmt][j],fitendch[pmt][j]);
    else if(j==1)
      hCharge[1]->Fit("fit_ideal_2","","",fitbeginch[pmt][j],fitendch[pmt][j]);
    else
      hCharge[2]->Fit("fit_ideal_3","","",fitbeginch[pmt][j],fitendch[pmt][j]);
  }
  
 
  // Voltage 1
  fit_ideal_1->GetParameters(par);
  foutFit<<"voltage\t"<<voltagestr[0]<<"\tchID\t"<<pmt
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
  foutFit<<"voltage\t"<<voltagestr[1]<<"\tchID\t"<<pmt
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
  foutFit<<"voltage\t"<<voltagestr[2]<<"\tchID\t"<<pmt
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
  c[pmt]->Write();
  c[pmt]->Print(pdfname[pmt].c_str(),"pdf");

  // print initial parameters at end of root file
  foutFit << "**************************** INITIAL PARAMETERS ****************************" << endl;
  for(int k = 0; k < 3; k++){
    foutFit << initparam[pmt][k] << endl;
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
