//#include "Fit_Spe.C" // the fit function is defined in this file
double PI = TMath::Pi();

double IdealResponse(double *x,double *par);

double RealResponse(double *x, double *par);

void FitChargeDistributions_v1(){
  // the function to be used to do fit
  gStyle->SetOptFit(1111);
  TF1* FReal = new TF1("FReal",RealResponse, -10, 500, 9);
  FReal->SetParNames("meanNpe","spePeak","speWidth", "amp","w","Q0","#sigma_{0}","#alpha","amp0");
  FReal->SetLineColor(2); FReal->SetLineStyle(1);
  double par[9];
  double parerr[9];

  // process 3 files in a batch
  string rtfilenames[3];
  string strchimney = "A18_PMT_";
  string strpmt = "5_6_7_8_";
  string voltagestr[3]={"1460","1490","1520"};
  for(int i=0; i<3; i++){
     rtfilenames[i]  = strchimney + strpmt + voltagestr[i] + "V_LedOn_result.root";
  }
  const int NCH = 4; // 4 PMTs
  int rebinfactor[NCH]={1, 1, 1, 1}; // histograms need some rebin
  double fitbeginch[NCH]={-1,-1,-1,-1};

  string outnameroot = strchimney + strpmt + "gain_realResponse.root";
  string outnametxt = strchimney + strpmt + "gain_fit_realResponse.txt";
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

      FReal->SetParameter(0,0.1);  FReal->SetParLimits(0,0,50);  // mean npe
      FReal->SetParameter(1,1.6);  FReal->SetParLimits(1,0.5,10);  // spe peak
      FReal->SetParameter(2,1.6*0.4); FReal->SetParLimits(2,0.1,10); // spe width
      FReal->SetParameter(3,hCharge[j]->Integral(56,2550)); FReal->SetParLimits(3,0.0,10.0*hCharge[j]->Integral(56,2550)); // scale factor
      cout<<hCharge[j]->Integral(56,2550)<<"\t"<<hCharge[j]->Integral(1,55)<<endl;

      FReal->SetParameter(4,0.0); FReal->SetParLimits(4,0,0.5); // w
      FReal->SetParameter(5,0.0); FReal->SetParLimits(5,-1.0,1.0); // Q0
      FReal->SetParameter(6,0.1); FReal->SetParLimits(6,0,1); // sigma0
      FReal->SetParameter(7,0.1); FReal->SetParLimits(7,0,1); // alpha
      FReal->SetParameter(8,hCharge[i]->Integral(1,55)); FReal->SetParLimits(8,0,10.0*hCharge[i]->Integral(1,55)); // scale factor, the numbers here depend on bin size, x-axis range and number of bins, 
      FReal->SetParameter(8,100); FReal->SetParLimits(8,0,10000);
      
      //if(j==1) FReal->SetParameter(3,1000);

      //if(j==1 && i==2) fitbeginch[j] = 1.0;
      
      //hCharge[j]->GetXaxis()->SetRangeUser(0, hCharge[j]->GetMean()*5.0);
      //double fitbegin = 0;
      //if(j==2) fitbegin = 0.8;
      //if(j==3) fitbegin = 4;
      double fitend = 10;
      for(int k=0; k<2;k++){
        if(j!=1)
           hCharge[j]->Fit("FReal","RQ","",fitbeginch[j],hCharge[j]->GetMean()*5.0);
        else hCharge[j]->Fit("FReal","RQ","",fitbeginch[j],fitend);
        FReal->GetParameters(par);
        FReal->SetParameters(par);
      }

      if(j!=1) 
          hCharge[j]->GetXaxis()->SetRangeUser(-10, hCharge[j]->GetMean()*5.0);
      else{
          hCharge[j]->GetXaxis()->SetRangeUser(-10, fitend);
          //hCharge[j]->GetYaxis()->SetRangeUser(0, 200);
      }

      FReal->GetParameters(par);
      //parerr = Fideal->GetParErrors();
      foutFit<<"voltage\t"<<voltagestr[i]<<"\tchID\t"<<j
             <<"\t"<<par[0]<<"\t"<<FReal->GetParError(0)
             <<"\t"<<par[1]<<"\t"<<FReal->GetParError(1)
             <<"\t"<<par[2]<<"\t"<<FReal->GetParError(2)
             <<"\t"<<par[3]<<"\t"<<FReal->GetParError(3)
             <<"\t"<<par[4]<<"\t"<<FReal->GetParError(4)
             <<"\t"<<par[5]<<"\t"<<FReal->GetParError(5)
             <<"\t"<<par[6]<<"\t"<<FReal->GetParError(6)
             <<"\t"<<par[7]<<"\t"<<FReal->GetParError(7)
             <<"\t"<<par[8]<<"\t"<<FReal->GetParError(8)
             <<"\t"<<FReal->GetChisquare()
             <<"\t"<<FReal->GetNDF()
             <<"\t"<<FReal->GetProb()
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


double RealResponse(double *x,double *par){
    double mu = par[0];
    double q = par[1];
    double sigma = par[2];
    double amplitude = par[3];

    double w = par[4];
    double Q0 = par[5];
    double sigma0 = par[6];
    double alpha = par[7];
    double Qsh = w/alpha;
    double amplitude0 = par[8];
    
    double sum=0;
    for(Int_t n=1; n<50; n++){
        sum += TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)*TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*PI*n));
    }

    double bkgterm=0;
    double stepvalue = 0;
    
    if(x[0]>=Q0){
        stepvalue = 1.0;
    }
    else stepvalue = 0.0;

    bkgterm = ( (1.0-w)/sigma0/TMath::Sqrt(2.0*PI)*TMath::Exp(-1.0*(x[0]-Q0)*(x[0]-Q0)/(2.0*sigma0*sigma0)) + w*stepvalue*alpha*TMath::Exp(-1.0*alpha*(x[0]-Q0)) )*TMath::Exp(-1.0*mu);

    return sum*amplitude + bkgterm*amplitude0;
    //return sum*amplitude;
}

void testRealResponse(){
  TF1* FReal = new TF1("FReal",RealResponse, -10, 100, 9);
  FReal->SetParNames("meanNpe","spePeak","speWidth", "amp","w","Q0","#sigma_{0}","#alpha","amp0");
  FReal->SetParameter(0,11.89);
  FReal->SetParameter(1,1.83);
  FReal->SetParameter(2,1.13);
  FReal->SetParameter(3,965);
  FReal->SetParameter(4,0.2451);
  FReal->SetParameter(5,-0.233);
  FReal->SetParameter(6,0.216);
  FReal->SetParameter(7,1);
  FReal->SetParameter(8,7.716e6);    

  FReal->Draw();
}
