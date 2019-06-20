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
const int NPAR = 10;
int ipar1[4] = {  0, // common mu
                  1, // q for hist 1
                  2, // sigma for hist 1
                  3  // amplitude for hist 1
};

int ipar2[4] = {  0, // common mu
                  4, // q for hist 2
                  5, // sigma for hist 2
                  6  // amplitude for hist 2
};

int ipar3[4] = {  0, // common mu
                  7, // q for hist 3
                  8, // sigma for hist 3
                  9  // amplitude for hist 3
};

// function declarations
double IdealResponse(double *x,double *par);
Double_t truncatedMean(TH1 *hist, int n_iterations, int n_rejection_stddevs = 3);

// create combined chi2 structure
struct GlobalChi2 {

  GlobalChi2( ROOT::Math::IMultiGenFunction & f1,
              ROOT::Math::IMultiGenFunction & f2,
              ROOT::Math::IMultiGenFunction & f3) : 
    fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

  double operator() (const double *par) const {
    double p1[4];
    double p2[4];
    double p3[4];

    for (int i = 0; i < 4; ++i){
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

void FitChargeDistributions(string pmtRow,
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
  string resultnames[3];
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
    resultnames[i]  = strchimney + strpmt + voltagestr[i] + "V_Led" + ledstr + ".pdf";
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
    TF1* fit_ideal_1 = new TF1("fit_ideal_1",IdealResponse, 0, 500, 4);
    fit_ideal_1->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
    fit_ideal_1->SetLineColor(2);
    fit_ideal_1->SetLineStyle(1);
    TF1* fit_ideal_2 = new TF1("fit_ideal_2",IdealResponse, 0, 500, 4);
    fit_ideal_2->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
    fit_ideal_2->SetLineColor(2);
    fit_ideal_2->SetLineStyle(1);
    TF1* fit_ideal_3 = new TF1("fit_ideal_1",IdealResponse, 0, 500, 4);
    fit_ideal_3->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
    fit_ideal_3->SetLineColor(2);
    fit_ideal_3->SetLineStyle(1);

    ROOT::Math::WrappedMultiTF1 wf1(*fit_ideal_1,1);
    ROOT::Math::WrappedMultiTF1 wf2(*fit_ideal_2,1);
    ROOT::Math::WrappedMultiTF1 wf3(*fit_ideal_3,1);

    // set data range
    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange range;

    range.SetRange(fitbeginch[i], fitendch[i]);

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
    Double_t par0[NPAR] = { hist_mean_1,
                            1.6,
                            1.6*0.4,
                            hCharge[0]->Integral(),
                            1.6,
                            1.6*0.4,
                            hCharge[1]->Integral(),
                            1.6,
                            1.6*0.4,
                            hCharge[2]->Integral()}; // starting values

    // set ranges on fit parameters
    fitter.Config().SetParamsSettings(NPAR, par0);

    fitter.Config().ParSettings(0).SetLimits(0.1, 100);
    // q
    for(int j = 1; j < 8; j += 3){
      fitter.Config().ParSettings(j).SetLimits(0.1, 10);
    }
    // sigma
    for(int j = 2; j < 9; j += 3){
      fitter.Config().ParSettings(j).SetLimits(0.1, 10);
    }
    // amplitude
    for(int j = 3; j < 10; j += 3){
      fitter.Config().ParSettings(j).SetLimits(0.1, 20000);
    }

    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");

    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(NPAR,globalChi2,0, chargeData[0].Size() + chargeData[1].Size() + chargeData[2].Size() ,true);
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
      c[i]->cd(j+1);
      hCharge[j]->GetXaxis()->SetRangeUser(0, fitendch[i]);
      hCharge[j]->Draw();
    }

    // write parameters to output txt file

    // Voltage 1
    fit_ideal_1->GetParameters(par);
    foutFit<<"voltage\t"<<voltagestr[0]<<"\tchID\t"<<i
            <<"\t"<<par[0]<<"\t"<<fit_ideal_1->GetParError(0)
            <<"\t"<<par[1]<<"\t"<<fit_ideal_1->GetParError(1)
            <<"\t"<<par[2]<<"\t"<<fit_ideal_1->GetParError(2)
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
