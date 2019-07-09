#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"

#include <string>
#include <iostream>

//#include "Fit_Spe.C" // the fit function is defined in this file

double PI = TMath::Pi();

double IdealResponse(double *x,double *par);
double RealResponse(double *x, double *par);
Double_t iGnE(Double_t x, Double_t q0, Double_t sigma0, 
              Double_t q1, Double_t sigma1, Double_t a, Double_t n);
Double_t Gn(Double_t x, Double_t q0, Double_t sigma0, Double_t q1, Double_t sigma1, Double_t n);

void test_background(){
  TF1 *f1 = new TF1("RealResponse",RealResponse,0,30,8);
  f1->SetParameters(1.68, 23.26, 0.192, 35.04, 11.73, 0.383, 0.034, 1);
  std::cout << f1->Eval(23) << std::endl;
  f1->Draw();
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

double RealResponse(double *x, double *par){
	double mu = par[0];
	double q0 = par[1];
	double sigma0 = par[2];
	double q1 = par[3];
	double sigma1 = par[4];
	double w = par[5];
	double a = par[6];
  double amplitude = par[7];

	Double_t sum = 0;

	for(Int_t n = 0; n < 50; n++){
		sum += (TMath::Power(mu, n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n))
				  * ((1 - w)*Gn(x[0], q0, sigma0, q1, sigma1, n) 
          + w*iGnE(x[0], q0, sigma0, q1, sigma1, a, n));
	}

	return amplitude*sum;
}

Double_t iGnE(Double_t x, Double_t q0, Double_t sigma0, 
              Double_t q1, Double_t sigma1, Double_t a, Double_t n){

  Double_t qn = q0 + n*q1;
  Double_t sigman = TMath::Sqrt(sigma0*sigma0 + n*sigma1*sigma1);
  Double_t sign = 1;
  if(x - qn - sigman*sigman*a > 0){
    sign = 1;
  } else if(x - qn - sigman*sigman*a == 0){
    sign = 0;
  } else {
    sign = -1;
  }

  Double_t iGnE = (a/2)*TMath::Exp(-1.0*a*(x - qn - a*sigman*sigman))
            *(TMath::Erf(TMath::Abs(q0 - qn - sigman*sigman*a)/(sigman*TMath::Sqrt2()))
            + sign*TMath::Erf(TMath::Abs(x - qn - sigman*sigman*a)/(sigman*TMath::Sqrt2())));

  return iGnE;
}

Double_t Gn(Double_t x, Double_t q0, Double_t sigma0, Double_t q1, Double_t sigma1, Double_t n){
  Double_t Gn = 0;
  if(n == 0){
    Gn = TMath::Exp(-1.0*(x - q0)*(x - q0)/(2*sigma0*sigma0))
           /(sigma0*TMath::Sqrt(2*TMath::Pi()));
  } else {
    Gn = TMath::Exp(-1.0*(x - n*q1)*(x - n*q1)/(2*n*sigma1*sigma1))
            /(sigma1*TMath::Sqrt(2*TMath::Pi()*n));
  }

  return Gn;
}

