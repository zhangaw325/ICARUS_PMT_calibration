#include "TH1.h"
#include "TCanvas.h"

//#include "Fit_Spe.C" // the fit function is defined in this file

double PI = TMath::Pi();

double IdealResponse(double *x,double *par);

void FitChargeDistributions(string pmtRow,
			    char pmt1, char pmt2, char pmt3, char pmt4,
			    int volt1, int volt2, int volt3, bool led){
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
  const int NCH = 4; // 4 PMTs
 
  int rebinfactor[NCH]={1, 1, 1, 1}; // rebin histograms
  double fitbeginch[NCH]={1.0,1.0,1.0,1.0};
  double fitend = 60;

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
    sprintf(canvasTitle,strchimney + "%d",i); // generate canvas title
    c[i] = new TCanvas(tempname,canvasTitlev.c_str(),1400,600); // generate canvas
    c[i]->Divide(3); // divide canvas into 3 pads along the width

    // Plot the histograms
    for(int j = 0; j<3; j++){
      c[i]->cd(j+1); // switch pads
      sprintf(tempname,"Results/FinalCharge_%d",i); // store string "Results/FinalCharge_%d" in tempname
      hCharge[j] = (TH1F*)f[j]->Get(tempname); // read histogram data from the ROOT file
      hCharge[j]->Rebin(rebinfactor[i]);
      hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
      hCharge[j]->Draw();

      // Set initial fit parameters
      Fideal->SetParameter(1,1.6);
      Fideal->SetParameter(2,1.6*0.4);
      Fideal->SetParameter(3,hCharge[j]->Integral());
      Fideal->SetParLimits(0,0.1,100);
      Fideal->SetParLimits(1,0.5,10);
      Fideal->SetParLimits(2,0.1,10);
      Fideal->SetParLimits(3,0.1,20000);

      // Iteratively fit more than once
      for(int k=0; k<2;k++){
        hCharge[j]->Fit("Fideal","RQ","",fitbeginch[j],fitend); // Fit the histogram
        Fideal->GetParameters(par);
        Fideal->SetParameters(par); // Set fit parameters for next iteration
      }

      hCharge[j]->GetXaxis()->SetRangeUser(0, fitend); // Set axes

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

  // Read each of the three output ROOT files and generate 3 canvases
  // for(int i=0; i<3; i++){
  //   TFile* f = new TFile(rtfilenames[i].c_str(),"read");
  //   sprintf(tempname,"c_%d",i);
  //   c[i] = new TCanvas(tempname,rtfilenames[i].c_str(),1200,900); // generate canvas
  //   c[i]->Divide(2,2); // divide canvas into four panels

  //   // Generate histograms and fits for each of the four PMTs in each ROOT file
  //   for(int j=0; j<NCH; j++){
  //     c[i]->cd(j+1); // switch to the proper panel
  //     sprintf(tempname,"Results/FinalCharge_%d",j); // store string "Results/FinalCharge_%d" in tempname
  //     hCharge[j] = (TH1F*)f->Get(tempname);
  //     hCharge[j]->Rebin(rebinfactor[j]);
  //     hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
  //     hCharge[j]->Draw(); 

  //     // Set initial fit parameters
  //     Fideal->SetParameter(1,1.6);
  //     Fideal->SetParameter(2,1.6*0.4);
  //     Fideal->SetParameter(3,hCharge[j]->Integral());
  //     Fideal->SetParLimits(0,0.1,100);
  //     Fideal->SetParLimits(1,0.5,10);
  //     Fideal->SetParLimits(2,0.1,10);
  //     Fideal->SetParLimits(3,0.1,20000);

  //     // Iteratively fit more than once
  //     for(int k=0; k<5;k++){
  //       hCharge[j]->Fit("Fideal","RQ","",fitbeginch[j],fitend); // Fit the histogram
  //       Fideal->GetParameters(par);
  //       Fideal->SetParameters(par); // Set fit parameters for next iteration
  //     }

  //     hCharge[j]->GetXaxis()->SetRangeUser(0, fitend); // Set axes

  //     Fideal->GetParameters(par);
  //     //parerr = Fideal->GetParErrors();
  //     foutFit<<"voltage\t"<<voltagestr[i]<<"\tchID\t"<<j
  //            <<"\t"<<par[0]<<"\t"<<Fideal->GetParError(0)
  //            <<"\t"<<par[1]<<"\t"<<Fideal->GetParError(1)
  //            <<"\t"<<par[2]<<"\t"<<Fideal->GetParError(2)
  //            <<"\t"<<Fideal->GetChisquare()
  //            <<"\t"<<Fideal->GetNDF()
  //            <<"\t"<<Fideal->GetProb()
  //            <<endl;
  //   }

  //   //c[i]->Update();
  //   outROOTfile->cd();
  //   c[i]->Write();
  //   //f->Close();
  //   //if(i>0) break;
  // }

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

