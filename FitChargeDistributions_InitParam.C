//***********************************************************************************************************************************************
// - Code adapts both FitChargeDistributions.C (high charge, simultaneous fit) and FitChargeDistributions_LC_ind.C (low charge fit).
//
// - The macro utilizes a .csv file containing initial parameters for
//   each PMT at three voltages in order to recreate the original fits.
//    - Each PMT should have three rows dedicated to it
//
// - The .csv file must contain two headers of the following form:
//   ,,,,Initial parameters,,,,,,,,Flags
//   ROOT File Name,Chimney,PMT channel,Voltage,Fit Begin,Fit End, Rebin Factor,Y Max,Mu (NPE),q (SPE),Sigma (SPE),Amplitude,Charge
//    - There should not be any empty rows or columns in the .csv file
//
// - The *LedOn*result.root files must be included in the same folder that this macro is stored in
//
// - Do not have any *LedOff*result.root files in the folder
//
// - Output (for each PMT):
//     - Root file
//     - .pdf of the canvas contents
//     - .txt of comma separated values including:
//          - final parameters: voltage, channel number
//                              mu, mu error, q, q error, sigma, sigma error,
//                              amplitude, amplitude error, chi2, ndf, fit probability
//          - initial parameters: channel id,
//                                start fit, end fit, rebin factor,
//                                mu, q, sigma, amplitude
//
// - Example function call in root:
//   FitChargeDistributions_InitParam("data.csv");
//***********************************************************************************************************************************************

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
                  4, // q for hist 2
                  5, // sigma for hist 2
                  6  // amplitude for hist 2
};

int ipar3[NPAR_i] = {  0,  // common mu
                  7, // q for hist 3
                  8, // sigma for hist 3
                  9 // amplitude for hist 3
};

// function declarations
double IdealResponse(double *x,double *par); // equation to be fit
Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3); // calculates truncated mean of the distribution
bool fileExist (const std::string& name); // returns whether a given file exists in the folder
string fit_low(const string rootFile1, const string rootFile2, const string rootFile3, // fits a single PMT using the low-charge method
	     const string pmtRow, const int pmt,
	     const int volt1, const int volt2, const int volt3,
	     const double * fitbeginch, const double * fitendch,
	     const int * rebinfactor, const double ymax,
	     const double * mu_0, const double * q_0,
	     const double * sigma_0, const double * amp_0
	     );
string fit_high(const string rootFile1, const string rootFile2, const string rootFile3, // fits a single PMT using the high-charge method
	      const string pmtRow, const int pmt,
	      const int volt1, const int volt2, const int volt3,
	      const double fbc_0, const double fec_0,
	      const int rbf_0, const double ymax,
	      const double mu_0, const double * q_0,
	      const double * sigma_0, const double * amp_0
	      );
void FitChargeDistributions_InitParam(string csvFile);

// Initial parameter storage structure
struct InitParam {
  string root; // root file name
  string chim; // chimney
  int pmt;   // channel #
  double volt; // voltage
  double fb;   // fit begin
  double fe;   // fit end
  int rb;      // rebin factor
  double y;    // maximum y value
  double m;    // mu
  double q;    // q
  double s;    // sigma
  double a;    // amplitude
};

// Combined chi2 structure
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
    // returns sum of the three chi2 of the individual functions
    return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
  }
  const  ROOT::Math::IMultiGenFunction * fChi2_1;
  const  ROOT::Math::IMultiGenFunction * fChi2_2;
  const  ROOT::Math::IMultiGenFunction * fChi2_3;
};

string fit_low(const string rootFile1, const string rootFile2, const string rootFile3,
	     const string pmtRow, const int pmt,
	     const int volt1, const int volt2, const int volt3,
	     const double * fitbeginch, const double * fitendch,
	     const int * rebinfactor, const double ymax,
	     const double * mu_0, const double * q_0,
	     const double * sigma_0, const double * amp_0
	     ){
  
  // arrays storing parameter information for a particular function
  string initparam[3];
  double par[NPAR_i];
  double parerr[NPAR_i];

  // formatting
  gStyle->SetOptFit(1111);

  // generate strings to hold information about the pmt and files
  string rtfilenames[3] = {rootFile1, rootFile2, rootFile3};
  string pdfname;
  string strchimney = pmtRow + "_PMT_";
  string strpmt = to_string(pmt) + "_";
  string voltagestr[3] = {to_string(volt1), to_string(volt2), to_string(volt3)};

  // get the index of the PMT given the group of four PMT's
  regex pmtVoltRgx(".*PMT_(\\w+)V.*");
  smatch pmtVolt;
  regex_search(rootFile1.begin(), rootFile1.end(), pmtVolt, pmtVoltRgx);
  istringstream s(pmtVolt.str(1));
  string field;
  int pmtIndex; // contains the index of the PMT (0, 1, 2, 3)
  for(int i = 0; i < 4; i++){
    getline(s, field, '_');
    if(stoi(field) == pmt)
      pmtIndex = i;
  }


  // if .pdf output already exists for the pmt, add a version number to the PMT
  struct stat buffer;
  int version = 1;
  string pmtVersion = "v" + to_string(version) + "_";
  pdfname = strchimney + strpmt + pmtVersion + "lowcharge.pdf";
  while(fileExist(pdfname) || fileExist(strchimney + strpmt + pmtVersion + "highcharge.pdf")){
    version ++;
    pmtVersion = "v" + to_string(version) + "_";
    pdfname = strchimney + strpmt + pmtVersion + "lowcharge.pdf";
  }
  cout << pdfname <<endl;
  string outnameroot = strchimney + strpmt + pmtVersion + "lowcharge_gain.root";
  string outnametxt = strchimney + strpmt + pmtVersion + "lowcharge_gain_fit.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TH1F* hCharge[3]; // histograms for each canvas
  TCanvas* c;
  char tempname[100];
  TFile* files[3];
  
  // Read each of the three data ROOT files and save them to files[]
  for(int i = 0; i < 3; i++){
    files[i] = new TFile(rtfilenames[i].c_str(),"read");
  }

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);

  // Print header
  foutFit << "**************************** PARAMETER VALUES ****************************" << endl;

  //*************************
  // Begin fit
  //*************************

  // generate histograms
  sprintf(tempname, "Results/FinalCharge_%d",pmtIndex);
  for(int j = 0; j < 3; j++){
    hCharge[j] = (TH1F*)files[j]->Get(tempname);
    hCharge[j]->SetTitle((strchimney + strpmt + pmtVersion + voltagestr[j] + "V").c_str());
    hCharge[j]->Rebin(rebinfactor[j]);
    hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
  }
  
  // generate fit functions
  TF1* fit_ideal_1 = new TF1("fit_ideal_1",IdealResponse, fitbeginch[0], fitendch[0], NPAR_i);
  fit_ideal_1->SetParNames("meanNpe","spePeak","speSigma","Amplitude");
  fit_ideal_1->SetLineColor(2);
  fit_ideal_1->SetLineStyle(1);
  TF1* fit_ideal_2 = new TF1("fit_ideal_2",IdealResponse, fitbeginch[1], fitendch[1], NPAR_i);
  fit_ideal_2->SetParNames("meanNpe","spePeak","speSigma","Amplitude");
  fit_ideal_2->SetLineColor(2);
  fit_ideal_2->SetLineStyle(1);
  TF1* fit_ideal_3 = new TF1("fit_ideal_3",IdealResponse, fitbeginch[2], fitendch[2], NPAR_i);
  fit_ideal_3->SetParNames("meanNpe","spePeak","speSigma","Amplitude");
  fit_ideal_3->SetLineColor(2);
  fit_ideal_3->SetLineStyle(1);

  // initial parameters for each voltage
  for(int k = 0; k < 3; k++){  
    // write initial parameters to output txt file
    initparam[k]=
      rtfilenames[k]+","    //root file
      +pmtRow+","           //pmt row
      +to_string(pmt)+","   //channel id
      +voltagestr[k]+","    //voltage
      +fitbeginch[k]+","    //start fit
      +fitendch[k]+","      //end fit
      +rebinfactor[k]+","   //rebin factor
      +mu_0[k]+","          //mu
      +q_0[k]+","           //q
      +sigma_0[k]+","       //sigma
      +amp_0[k]+","         //amplitude
      +"0";                 //lowcharge
  }
  
  fit_ideal_1->SetParameters(mu_0[0],q_0[0],sigma_0[0],amp_0[0]); 
  fit_ideal_2->SetParameters(mu_0[1],q_0[1],sigma_0[1],amp_0[1]); 
  fit_ideal_3->SetParameters(mu_0[2],q_0[2],sigma_0[2],amp_0[2]);
  
  // Generate canvas and plot the PMT histograms on them
  sprintf(tempname, "c_%d",pmt);
  string canvasTitle = strchimney + strpmt + "v" + to_string(version);
  c = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
  c->Divide(3); // divide canvas into 3 pads along the width
  for(int j = 0; j < 3; j++){
    c->cd(j+1);
    hCharge[j]->GetXaxis()->SetRangeUser(0, fitendch[j]);
    hCharge[j]->GetYaxis()->SetRangeUser(0, ymax);
    gStyle->SetOptFit();
    if(j==0)
      hCharge[0]->Fit("fit_ideal_1","","",fitbeginch[j],fitendch[j]);
    else if(j==1)
      hCharge[1]->Fit("fit_ideal_2","","",fitbeginch[j],fitendch[j]);
    else
      hCharge[2]->Fit("fit_ideal_3","","",fitbeginch[j],fitendch[j]);
  }
  
 
  // Voltage 1
  fit_ideal_1->GetParameters(par);
  foutFit<<"PMT "<<pmt<<","
	 <<"Voltage "<<voltagestr[0]
	 <<","<<par[0]<<","<<fit_ideal_1->GetParError(0)
	 <<","<<par[1]<<","<<fit_ideal_1->GetParError(1)
	 <<","<<par[2]<<","<<fit_ideal_1->GetParError(2)
	 <<","<<par[3]<<","<<fit_ideal_1->GetParError(3)
	 <<","<<fit_ideal_1->GetChisquare()
	 <<","<<fit_ideal_1->GetNDF()
	 <<","<<fit_ideal_1->GetProb()
	 <<endl;
    
  // Voltage 2
  fit_ideal_2->GetParameters(par);
  foutFit<<"PMT "<<pmt<<","
	 <<"Voltage "<<voltagestr[1]
	 <<","<<par[0]<<","<<fit_ideal_1->GetParError(0)
	 <<","<<par[1]<<","<<fit_ideal_1->GetParError(1)
	 <<","<<par[2]<<","<<fit_ideal_1->GetParError(2)
	 <<","<<par[3]<<","<<fit_ideal_1->GetParError(3)
	 <<","<<fit_ideal_1->GetChisquare()
	 <<","<<fit_ideal_1->GetNDF()
	 <<","<<fit_ideal_1->GetProb()
	 <<endl;

  // Voltage 3
  fit_ideal_3->GetParameters(par);
  foutFit<<"PMT "<<pmt<<","
	 <<"Voltage "<<voltagestr[2]
	 <<","<<par[0]<<","<<fit_ideal_1->GetParError(0)
	 <<","<<par[1]<<","<<fit_ideal_1->GetParError(1)
	 <<","<<par[2]<<","<<fit_ideal_1->GetParError(2)
	 <<","<<par[3]<<","<<fit_ideal_1->GetParError(3)
	 <<","<<fit_ideal_1->GetChisquare()
	 <<","<<fit_ideal_1->GetNDF()
	 <<","<<fit_ideal_1->GetProb()
	 <<endl;

  // Output string - lists parameters in a .csv style 
  string output =to_string(pmt)
    +","+pdfname
    +","+to_string(fit_ideal_1->GetParameter(0))+","+to_string(fit_ideal_1->GetParError(0))
    +","+to_string(fit_ideal_1->GetParameter(1))+","+to_string(fit_ideal_1->GetParError(1))
    +","+to_string(fit_ideal_1->GetParameter(2))+","+to_string(fit_ideal_1->GetParError(2))
    +","+to_string(fit_ideal_1->GetParameter(3))+","+to_string(fit_ideal_1->GetParError(3))
    +","+to_string(fit_ideal_1->GetChisquare())
    +","+to_string(fit_ideal_1->GetNDF())
    +","+to_string(fit_ideal_1->GetProb())
    +"\n"
    +to_string(pmt)
    +","+pdfname
    +","+to_string(fit_ideal_2->GetParameter(0))+","+to_string(fit_ideal_2->GetParError(0))
    +","+to_string(fit_ideal_2->GetParameter(1))+","+to_string(fit_ideal_2->GetParError(1))
    +","+to_string(fit_ideal_2->GetParameter(2))+","+to_string(fit_ideal_2->GetParError(2))
    +","+to_string(fit_ideal_2->GetParameter(3))+","+to_string(fit_ideal_2->GetParError(3))
    +","+to_string(fit_ideal_2->GetChisquare())
    +","+to_string(fit_ideal_2->GetNDF())
    +","+to_string(fit_ideal_2->GetProb())
    +"\n"
    +to_string(pmt)
    +","+pdfname
    +","+to_string(fit_ideal_3->GetParameter(0))+","+to_string(fit_ideal_3->GetParError(0))
    +","+to_string(fit_ideal_3->GetParameter(1))+","+to_string(fit_ideal_3->GetParError(1))
    +","+to_string(fit_ideal_3->GetParameter(2))+","+to_string(fit_ideal_3->GetParError(2))
    +","+to_string(fit_ideal_3->GetParameter(3))+","+to_string(fit_ideal_3->GetParError(3))
    +","+to_string(fit_ideal_3->GetChisquare())
    +","+to_string(fit_ideal_3->GetNDF())
    +","+to_string(fit_ideal_3->GetProb())
    +"\n";
  
  //*************************
  // End fit
  //*************************

  // write results to output ROOT file
  outROOTfile->cd();
  c->Write();
  c->Print(pdfname.c_str(),"pdf");

  // print initial parameters at end of root file
  // delete input ROOT files
  foutFit << "**************************** INITIAL PARAMETERS ****************************" << endl;
  for(int k = 0; k < 3; k++){
    foutFit << initparam[k] << endl;
    delete files[k];
  }
  
  // close output files
  outROOTfile->Close();
  foutFit.close();

  // return string with three lines of final parameters
  return output;
}

string fit_high(const string rootFile1, const string rootFile2, const string rootFile3,
	      const string pmtRow, const int pmt,
	      const int volt1, const int volt2, const int volt3,
	      const double fbc_0, const double fec_0,
	      const int rbf_0, const double ymax,
	      const double mu_0, const double * q_0,
	      const double * sigma_0, const double * amp_0
	      ){
  
  // histogram and fit options
  int rebinfactor = rbf_0; // rebin histograms
  double fitbeginch = fbc_0; // fit start locations
  double fitendch = fec_0; // fit end locations
  gStyle->SetOptFit(1111); // formatting

  // arrays storing parameter information for a particular function
  string initparam[4]; // for outputting initial parameters
  double par[NPAR_i];
  double parerr[NPAR_i];

  // process 3 files in a batch
  string rtfilenames[3] = {rootFile1, rootFile2, rootFile3};
  string pdfname;
  string strchimney = pmtRow + "_PMT_";
  string strpmt = to_string(pmt) + "_" ;
  string voltagestr[3] = {to_string(volt1), to_string(volt2), to_string(volt3)};

  // get the index of the PMT given the group of four PMT's
  regex pmtVoltRgx(".*PMT_(\\w+)V.*");
  smatch pmtVolt;
  regex_search(rootFile1.begin(), rootFile1.end(), pmtVolt, pmtVoltRgx);
  istringstream s(pmtVolt.str(1));
  string field;
  int pmtIndex; // contains the index of the PMT (0, 1, 2, 3)
  for(int i = 0; i < 4; i++){
    getline(s, field, '_');
    if(stoi(field) == pmt)
      pmtIndex = i;
  }

  //OUTPUT
  // if .pdf output already exists for the pmt, add a version number to the PMT
  struct stat buffer;
  int version = 1;
  string pmtVersion = "v" + to_string(version) + "_";
  pdfname = strchimney + strpmt + pmtVersion + "highcharge.pdf";
  while(fileExist(pdfname) || fileExist(strchimney + strpmt + pmtVersion + "lowcharge.pdf")){
    version ++;
    pmtVersion = "v" + to_string(version) + "_";
    pdfname = strchimney + strpmt + pmtVersion + "highcharge.pdf";
  }
  cout << pdfname <<endl;
  string outnameroot = strchimney + strpmt + pmtVersion + "highcharge_gain.root";
  string outnametxt = strchimney + strpmt + pmtVersion + "highcharge_gain_fit.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TFile* files[3];
  // Read each of the three data ROOT files and save them to files[]
  for(int i = 0; i < 3; i++){
    files[i] = new TFile(rtfilenames[i].c_str(),"read");
  }

  // canvas handling
  TH1F* hCharge[3]; // histograms for each canvas
  TCanvas* c;
  char tempname[100];

  // minimize a max of 5000 times
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(5000);

  // Print header
  foutFit << "**************************** PARAMETER VALUES ****************************" << endl;
  
  // Generate canvas and plot the PMT histogram
  sprintf(tempname, "c_%d",pmtIndex);
  string canvasTitle = strchimney + strpmt + "v" + to_string(version);
  c = new TCanvas(tempname,canvasTitle.c_str(),1400,600); // generate canvas
  c->Divide(3); // divide canvas into 3 pads along the width
  
  //*************************
  // Begin fit
  //*************************
  
  // generate histograms
  sprintf(tempname, "Results/FinalCharge_%d",pmtIndex);
  for(int j = 0; j < 3; j++){
      hCharge[j] = (TH1F*)files[j]->Get(tempname);
      hCharge[j]->SetTitle((strchimney + strpmt + pmtVersion + voltagestr[j] + "V").c_str());
      hCharge[j]->Rebin(rebinfactor);
      hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
  }
  
  // generate fit functions
  TF1* fit_ideal_1 = new TF1("fit_ideal_1",IdealResponse, 0, 500, NPAR_i);
  fit_ideal_1->SetParNames("meanNpe","spePeak","speSigma","Amplitude");
  fit_ideal_1->SetLineColor(2);
  fit_ideal_1->SetLineStyle(1);
  TF1* fit_ideal_2 = new TF1("fit_ideal_2",IdealResponse, 0, 500, NPAR_i);
  fit_ideal_2->SetParNames("meanNpe","spePeak","speSigma","Amplitude");
  fit_ideal_2->SetLineColor(2);
  fit_ideal_2->SetLineStyle(1);
  TF1* fit_ideal_3 = new TF1("fit_ideal_3",IdealResponse, 0, 500, NPAR_i);
  fit_ideal_3->SetParNames("meanNpe","spePeak","speSigma","Amplitude");
  fit_ideal_3->SetLineColor(2);
  fit_ideal_3->SetLineStyle(1);
  
  ROOT::Math::WrappedMultiTF1 wf1(*fit_ideal_1,1);
  ROOT::Math::WrappedMultiTF1 wf2(*fit_ideal_2,1);
  ROOT::Math::WrappedMultiTF1 wf3(*fit_ideal_3,1);
  
  // set data range
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range;
  
  range.SetRange(fitbeginch, fitendch);
  
  // get bin data
  ROOT::Fit::BinData chargeData[3]; 
  
  for(int j = 0; j < 3; j++){
    ROOT::Fit::FillData(chargeData[j], hCharge[j]);
  }
  
    // create chi2 function
  ROOT::Fit::Chi2Function chi2_1(chargeData[0], wf1);
  ROOT::Fit::Chi2Function chi2_2(chargeData[1], wf2);
  ROOT::Fit::Chi2Function chi2_3(chargeData[2], wf3);
  GlobalChi2 globalChi2(chi2_1, chi2_2, chi2_3);
  
  // create fitter
  ROOT::Fit::Fitter fitter;
  Double_t hist_mean_1 = truncatedMean(hCharge[1], 10);
  Double_t par0[NPAR] = { mu_0,
			  q_0[0],
			  sigma_0[0],
			  amp_0[0],
			  q_0[1],
			  sigma_0[1],
			  amp_0[1],
			  q_0[2],
			  sigma_0[2],
			  amp_0[2]}; // starting values
  
  // set ranges on fit parameters
  fitter.Config().SetParamsSettings(NPAR, par0)

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  
  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(NPAR,globalChi2,0, chargeData[0].Size() + chargeData[1].Size() + chargeData[2].Size(), true);
  ROOT::Fit::FitResult result;
  for(int j = 0; j < 3; j++){
    result = fitter.Result();
    // fitter updates fit parameters after fitting
  }
  result.Print(std::cout);
  
  // display results
  fit_ideal_1->SetFitResult(result,ipar1);
  fit_ideal_1->SetRange(range().first, range().second);
  hCharge[0]->GetListOfFunctions()->Add(fit_ideal_1);
  
  fit_ideal_2->SetFitResult(result,ipar2);
  fit_ideal_2->SetRange(range().first, range().second);
  hCharge[1]->GetListOfFunctions()->Add(fit_ideal_2);
  
  fit_ideal_3->SetFitResult(result, ipar3);
  fit_ideal_3->SetRange(range().first, range().second);
  hCharge[2]->GetListOfFunctions()->Add(fit_ideal_3);
  
  for(int j = 0; j < 3; j++){
    c->cd(j+1);
    hCharge[j]->GetXaxis()->SetRangeUser(0, fitendch);
    hCharge[j]->GetYaxis()->SetRangeUser(0, ymax);
    hCharge[j]->Draw();
  }
  
  // write initial parameters to output txt file
  for(int k = 0; k < 3; k++){  
    // write initial parameters to output txt file
    initparam[k]=
      rtfilenames[k]+","    //root file
      +pmtRow+","           //pmt row
      +to_string(pmt)+","   //channel id
      +voltagestr[k]+","    //voltage
      +fitbeginch+","       //start fit
      +fitendch+","         //end fit
      +rebinfactor+","      //rebin factor
      +mu_0+","             //mu
      +q_0[k]+","           //q
      +sigma_0[k]+","       //sigma
      +amp_0[k]+","         //amplitude
      +"1";                 //highcharge
  }
  
  // Voltage 1
  fit_ideal_1->GetParameters(par);
  foutFit<<"PMT "<<pmt<<","
	 <<"Voltage "<<voltagestr[0]
	 <<","<<par[0]<<","<<fit_ideal_1->GetParError(0)
	 <<","<<par[1]<<","<<fit_ideal_1->GetParError(1)
	 <<","<<par[2]<<","<<fit_ideal_1->GetParError(2)
	 <<","<<par[3]<<","<<fit_ideal_1->GetParError(3)
	 <<","<<fit_ideal_1->GetChisquare()
	 <<","<<fit_ideal_1->GetNDF()
	 <<","<<fit_ideal_1->GetProb()
	 <<endl;
    
  // Voltage 2
  fit_ideal_2->GetParameters(par);
  foutFit<<"PMT "<<pmt<<","
	 <<"Voltage "<<voltagestr[1]
	 <<","<<par[0]<<","<<fit_ideal_1->GetParError(0)
	 <<","<<par[1]<<","<<fit_ideal_1->GetParError(1)
	 <<","<<par[2]<<","<<fit_ideal_1->GetParError(2)
	 <<","<<par[3]<<","<<fit_ideal_1->GetParError(3)
	 <<","<<fit_ideal_1->GetChisquare()
	 <<","<<fit_ideal_1->GetNDF()
	 <<","<<fit_ideal_1->GetProb()
	 <<endl;

  // Voltage 3
  fit_ideal_3->GetParameters(par);
  foutFit<<"PMT "<<pmt<<","
	 <<"Voltage "<<voltagestr[2]
	 <<","<<par[0]<<","<<fit_ideal_1->GetParError(0)
	 <<","<<par[1]<<","<<fit_ideal_1->GetParError(1)
	 <<","<<par[2]<<","<<fit_ideal_1->GetParError(2)
	 <<","<<par[3]<<","<<fit_ideal_1->GetParError(3)
	 <<","<<fit_ideal_1->GetChisquare()
	 <<","<<fit_ideal_1->GetNDF()
	 <<","<<fit_ideal_1->GetProb()
	 <<endl;

  // Output string
  string output =to_string(pmt)
    +","+pdfname
    +","+to_string(fit_ideal_1->GetParameter(0))+","+to_string(fit_ideal_1->GetParError(0))
    +","+to_string(fit_ideal_1->GetParameter(1))+","+to_string(fit_ideal_1->GetParError(1))
    +","+to_string(fit_ideal_1->GetParameter(2))+","+to_string(fit_ideal_1->GetParError(2))
    +","+to_string(fit_ideal_1->GetParameter(3))+","+to_string(fit_ideal_1->GetParError(3))
    +","+to_string(fit_ideal_1->GetChisquare())
    +","+to_string(fit_ideal_1->GetNDF())
    +","+to_string(fit_ideal_1->GetProb())
    +"\n"
    +to_string(pmt)
    +","+pdfname
    +","+to_string(fit_ideal_2->GetParameter(0))+","+to_string(fit_ideal_2->GetParError(0))
    +","+to_string(fit_ideal_2->GetParameter(1))+","+to_string(fit_ideal_2->GetParError(1))
    +","+to_string(fit_ideal_2->GetParameter(2))+","+to_string(fit_ideal_2->GetParError(2))
    +","+to_string(fit_ideal_2->GetParameter(3))+","+to_string(fit_ideal_2->GetParError(3))
    +","+to_string(fit_ideal_2->GetChisquare())
    +","+to_string(fit_ideal_2->GetNDF())
    +","+to_string(fit_ideal_2->GetProb())
    +"\n"
    +to_string(pmt)
    +","+pdfname
    +","+to_string(fit_ideal_3->GetParameter(0))+","+to_string(fit_ideal_3->GetParError(0))
    +","+to_string(fit_ideal_3->GetParameter(1))+","+to_string(fit_ideal_3->GetParError(1))
    +","+to_string(fit_ideal_3->GetParameter(2))+","+to_string(fit_ideal_3->GetParError(2))
    +","+to_string(fit_ideal_3->GetParameter(3))+","+to_string(fit_ideal_3->GetParError(3))
    +","+to_string(fit_ideal_3->GetChisquare())
    +","+to_string(fit_ideal_3->GetNDF())
    +","+to_string(fit_ideal_3->GetProb())
    +"\n";
  
  //*************************
  // End fit
  //*************************
  
  // write results to output ROOT file
  outROOTfile->cd();
  c->Write();
  c->Print(pdfname.c_str(),"pdf");

  // print initial parameters at end of root file
  // delete input ROOT files
  foutFit << "**************************** INITIAL PARAMETERS ****************************" << endl;
  for(int k = 0; k < 3; k++){
    foutFit << initparam[k] << endl;
    delete files[k];
  }

  // close output files
  outROOTfile->Close();
  foutFit.close();

  // return output
  return output;
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
  Double_t mean = 0.;
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

bool fileExist (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

void FitChargeDistributions_InitParam(string csvFile)
{  
  // Names of branches
  char branchName[13][20] = {"ROOTFileName",
			     "Channel",
			     "PMT",
			     "Voltage",
			     "FitBegin",
			     "FitEnd",
			     "Rebin",
			     "Ymax",
			     "Mu",
			     "q",
			     "Sigma",
			     "Amplitude",
			     "ChargeStatus"};
  
  // Variables used to fill branches
  Char_t chargeDataFile[500];
  Char_t chimney[500];
  Int_t pmtChannel, voltage, rebin;
  Double_t begin, end, ymax, mu, q, sigma, amplitude;
  Bool_t chargeStatus;
  Int_t nlines = 0;
  
  // READ .CSV FILE
  // File handling
  ifstream inputFile;
  inputFile.open(csvFile.c_str());
  string treeName = csvFile.substr(0,csvFile.find(".csv"));
  string rootOutputName = treeName + ".root";
  TFile *rootOutputFile = new TFile(rootOutputName.c_str(), "RECREATE");
  TTree *tree =  new TTree(treeName.c_str(),treeName.c_str());
  
  // Instantiate branches
  tree->Branch(branchName[0],&chargeDataFile,"chargeDataFile/C");
  tree->Branch(branchName[1],&chimney,"chimney/C");
  tree->Branch(branchName[2],&pmtChannel,"pmtChannel/I");
  tree->Branch(branchName[3],&voltage,"voltage/I");
  tree->Branch(branchName[4],&begin,"begin/D");
  tree->Branch(branchName[5],&end,"end/D");
  tree->Branch(branchName[6],&rebin,"rebin/I");
  tree->Branch(branchName[7],&ymax,"ymax/D");
  tree->Branch(branchName[8],&mu,"mu/D");
  tree->Branch(branchName[9],&q,"q/D");
  tree->Branch(branchName[10],&sigma,"sigma/D");
  tree->Branch(branchName[11],&amplitude,"amplitude/D");
  tree->Branch(branchName[12],&chargeStatus,"chargeStatus/O");
  
  string line, field;
  // Grabs first two lines of headers
  getline(inputFile, line);
  getline(inputFile, line);
  while(getline(inputFile, line)){
    istringstream s(line);
    while (getline(s, field, ',')){
      strlcpy(chargeDataFile, field.c_str(), sizeof(chargeDataFile));
      getline(s, field, ',');
      strlcpy(chimney, field.c_str(), sizeof(chimney));
	getline(s, field, ',');
	pmtChannel = stoi(field);
	getline(s, field, ',');
	voltage = stoi(field);
	getline(s, field, ',');
	begin = stod(field);
	getline(s, field, ',');
	end = stod(field);
	getline(s, field, ',');
	rebin = stoi(field);
	getline(s, field, ',');
	ymax = stod(field);
	getline(s, field, ',');
	mu = stod(field);
	getline(s, field, ',');
	q = stod(field);
	getline(s, field, ',');
	sigma = stod(field);
	getline(s, field, ',');
	amplitude = stod(field);
	getline(s, field, ',');
	chargeStatus = field.compare("0");
	tree->Fill();
    }
  }
  
  rootOutputFile->Write();

  // Output file with all final parameters written
  string txtOutputFile = treeName + ".txt";
  fstream fout(txtOutputFile,ios::out);
  fout << "pmt,pdf_name,mean_npe,mean_npe_err,spe_mean,spe_mean_err,spe_sigma,spe_sigma_err,amplitude,amplitude_err,chi2,ndf,fit_probability"<<endl;
  
  // Initial parameter storage structures used to pass data into fit function
  InitParam ip[3];

  // Tracks current chimney; sets initial value
  Char_t currentChim[500];
  tree->GetEvent(0);
  memset(currentChim, 0, 500);
  strncpy(currentChim, chimney, 500);
  fout << currentChim << ",,,,,,,,,,,," << endl;
  
  // Loop over all entries of the TTree accessing data in groups of 3.
  int count = 0;
  while (count < tree->GetEntries()) {
    // Store information on initial parameters
    for(int i = 0; i < 3; i++){
      memset(chargeDataFile, 0, 500);
      memset(chimney, 0, 500);
      tree->GetEvent(count);
      ip[i].root = chargeDataFile;
      ip[i].chim = chimney;
      if(strcmp(currentChim,chimney) != 0){
	strncpy(currentChim, chimney, 500);
	fout << currentChim << ",,,,,,,,,,,," << endl;
      }
      ip[i].pmt = pmtChannel;
      ip[i].volt = voltage;
      ip[i].fb = begin;
      ip[i].fe = end;
      ip[i].rb = rebin;
      ip[i].y = ymax;
      ip[i].m = mu;
      ip[i].q = q;
      ip[i].s = sigma;
      ip[i].a = amplitude;
      count++;
      if(count >= tree->GetEntries())
	break;
    }
    // If charge is 1, perform high charge fit
    if(chargeStatus){
      double q_list[3] = {ip[0].q, ip[1].q,ip[2].q};
      double s_list[3] = {ip[0].s, ip[1].s,ip[2].s};
      double a_list[3] = {ip[0].a, ip[1].a,ip[2].a};
      cout << ip[0].root << " " << ip[1].root << " " << ip[2].root << endl;
      fout << fit_high(ip[0].root, ip[1].root, ip[2].root,
		       ip[2].chim, ip[2].pmt,
		       ip[0].volt, ip[1].volt, ip[2].volt,
		       ip[2].fb, ip[2].fe,
		       ip[2].rb, ip[2].y,
		       ip[2].m,
		       q_list, s_list, a_list
		       );
    }
    // If charge is 0, perform low charge fit
    else{
      double fb_list[3] = {ip[0].fb, ip[1].fb, ip[2].fb};
      double fe_list[3] = {ip[0].fe, ip[1].fe, ip[2].fe};
      int rb_list[3] = {ip[0].rb, ip[1].rb, ip[2].rb};
      double m_list[3] = {ip[0].m, ip[1].m, ip[2].m};
      double q_list[3] = {ip[0].q, ip[1].q,ip[2].q};
      double s_list[3] = {ip[0].s, ip[1].s,ip[2].s};
      double a_list[3] = {ip[0].a, ip[1].a,ip[2].a}; 
      fout << fit_low(ip[0].root, ip[1].root, ip[2].root,
		      ip[2].chim, ip[2].pmt,
		      ip[0].volt, ip[1].volt, ip[2].volt,
		      fb_list,
		      fe_list,
		      rb_list,
		      ip[2].y,
		      m_list,
		      q_list,
		      s_list,
		      a_list
		      );
    }
  }
  
  fout.close();

  // Print statements
  cout << "Analysis on " << csvFile << " complete" << endl;
  cout << "Output written to:"<< endl;
  cout << "\t" << rootOutputName << endl;
  cout << "\t" << txtOutputFile << endl;
}
