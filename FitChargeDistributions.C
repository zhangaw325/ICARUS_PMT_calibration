//#include "Fit_Spe.C" // the fit function is defined in this file
double PI = TMath::Pi();

double IdealResponse(double *x,double *par);

void FitChargeDistributions(){
  // the function to be used to do fit
  gStyle->SetOptFit(1111);
  TF1* Fideal = new TF1("Fideal",IdealResponse, 0, 500, 4);
  Fideal->SetParNames("meanNpe","spePeak","speWidth","Amplitude");
  Fideal->SetLineColor(2); Fideal->SetLineStyle(1);
  double par[4];
  double parerr[4];

  // process 3 files in a batch
  string rtfilenames[3];
  string strchimney = "A18_PMT_";
  string strpmt = "5_6_7_8_";
  string voltagestr[3]={"1460","1490","1520"};
  for(int i=0; i<3; i++){
     rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_LedOn_result.root";
  }
  const int NCH = 4; // 4 PMTs
  /*
  string PMTIDstr[NCH]={
                    "A12_PMT1",
                    "A12_PMT2",
                    "A12_PMT3",
                    "A12_PMT4"
                    };
                    */
  int rebinfactor[NCH]={5, 1, 5, 5}; // histograms need some rebin
  double fitbeginch[NCH]={2,1.5,1.5,0.5};

  string outnameroot = strchimney + strpmt + "gain.root";
  string outnametxt = strchimney + strpmt + "gain_fit.txt";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream foutFit(outnametxt.c_str(),ios::out);
  TH1F* hCharge[NCH];
  TCanvas* c[3];
  char tempname[100];
  for(int i=0; i<3; i++){
    TFile* f = new TFile(rtfilenames[i].c_str(),"read");
    sprintf(tempname,"c_%d",i);
    c[i] = new TCanvas(tempname,rtfilenames[i].c_str(),1200,900);
//    c[i]->SetCanvasSize(1200,900);
    c[i]->Divide(2,2);
    for(int j=0; j<NCH; j++){
      c[i]->cd(j+1);
      sprintf(tempname,"Results/FinalCharge_%d",j);
      hCharge[j] = (TH1F*)f->Get(tempname);
      hCharge[j]->Rebin(rebinfactor[j]);
      hCharge[j]->SetXTitle("Charge in pC, (10^{7} electrons = 1.6 pC)");
      hCharge[j]->Draw(); 

      //Fideal->SetParameter(0,0.1);
      Fideal->SetParameter(1,1.6);
      Fideal->SetParameter(2,1.6*0.4);
      Fideal->SetParameter(3,hCharge[j]->Integral());
      Fideal->SetParLimits(0,0.1,100);
      Fideal->SetParLimits(1,0.5,10);
      Fideal->SetParLimits(2,0.1,10);
      Fideal->SetParLimits(3,0.1,20000);

      if(j==1) Fideal->SetParameter(3,1000);

      if(j==1 && i==2) fitbeginch[j] = 1.0;
      
      //hCharge[j]->GetXaxis()->SetRangeUser(0, hCharge[j]->GetMean()*5.0);
      //double fitbegin = 0;
      //if(j==2) fitbegin = 0.8;
      //if(j==3) fitbegin = 4;
      double fitend = 10;
      for(int k=0; k<2;k++){
        if(j!=1)
           hCharge[j]->Fit("Fideal","RQ","",fitbeginch[j],hCharge[j]->GetMean()*5.0);
        else hCharge[j]->Fit("Fideal","RQ","",fitbeginch[j],fitend);
        Fideal->GetParameters(par);
        Fideal->SetParameters(par);
      }

      if(j!=1) 
          hCharge[j]->GetXaxis()->SetRangeUser(0, hCharge[j]->GetMean()*5.0);
      else{
          hCharge[j]->GetXaxis()->SetRangeUser(0, fitend);
          hCharge[j]->GetYaxis()->SetRangeUser(0, 200);
      }

      Fideal->GetParameters(par);
      //parerr = Fideal->GetParErrors();
      foutFit<<"voltage\t"<<voltagestr[i]<<"\tchID\t"<<j
             <<"\t"<<par[0]<<"\t"<<Fideal->GetParError(0)
             <<"\t"<<par[1]<<"\t"<<Fideal->GetParError(1)
             <<"\t"<<par[2]<<"\t"<<Fideal->GetParError(2)
             <<"\t"<<Fideal->GetChisquare()
             <<"\t"<<Fideal->GetNDF()
             <<"\t"<<Fideal->GetProb()
             <<endl;
    }

    //c[i]->Update();
    outROOTfile->cd();
    c[i]->Write();
    //f->Close();
    //if(i>0) break;
  }
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

