#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TF1.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

void AfterpulseProbPlot(std::string input_file_name, std::string title="", std::string x_title="", 
						std::string y_title="", Int_t nBins=100/*, Double_t start=0, Double_t end=100*/){
	// prep output files
	std::string out_name = input_file_name.substr(0, input_file_name.size() - 4);
	std::string out_name_pdf = out_name + ".pdf";
	std::string out_name_root = out_name + ".root";

	TFile *output_file_root = new TFile(out_name_root.c_str(), "recreate");

	// load input file
	std::ifstream input_file(input_file_name);

	// prep plot to display
	TH2 *plot = new TH2F("afterpulse_plot",title.c_str(),nBins,0,65,nBins,-50,400);
	plot->SetXTitle(x_title.c_str());
	plot->SetYTitle(y_title.c_str());
	//plot->SetLineColor(1); // make lines black
	plot->SetMarkerStyle(20);
	plot->SetMarkerSize(0.4);

	// prep normalized histogram
	TH1 *normHist = new TH1F("normalized_hist","Normalized afterpulse probabilities",150,-10,50);
	normHist->SetXTitle(y_title.c_str());
	normHist->SetYTitle("counts");
	normHist->SetLineColor(1);

	// read input file line by line and fill histogram
	std::string line;
	Double_t PMT, nPE, gain, prob, current;

	Double_t current_cut_lower = 0;
	Double_t current_cut_upper = 65;

	int i = 0;
	while (std::getline(input_file, line)){
		std::istringstream iss(line);
		do {
			if(!(iss >> PMT >> nPE >> gain >> prob >> current)){continue;};

			current = nPE * gain;

			if(current > current_cut_lower && current < current_cut_upper){
				plot->Fill(current, prob*100);
			}

			normHist->Fill(prob*100/current);
		}
		while (iss);
	}

	// create canvas
	TCanvas *canvas = new TCanvas("c1", title.c_str(), 1200, 600);
	canvas->Divide(2,1);

	canvas->cd(1);
	canvas->GetPad(1)->SetGrid();
	gStyle->SetOptStat(0);
	// TF1 *linearFit = new TF1("linear","[0]*x",0,5);
	// gStyle->SetOptStat(0);
	// plot->Fit("pol1");
	// gStyle->SetOptFit(111);
	// plot->Fit("linear","R+");
	// TF1 *f = (TF1*) plot->GetListOfFunctions()->FindObject("linear");
	// f->SetLineColor(3);
	// f->SetLineWidth(2);
	// gStyle->SetOptFit(111);

	// draw histogram
	//gStyle->SetGridStyle(0);
	plot->Draw("scat");

	canvas->cd(2);
	normHist->Draw();

	// print mean probability
	std::cout << "Mean normalized afterpulsing probability is: " << normHist->GetMean() << std::endl;

	// save scatterplot PDF
	std::string out_name_scat_pdf = out_name + "_scatterplot.pdf";
	canvas->GetPad(1)->Print(out_name_scat_pdf.c_str(), "pdf");

	// save ROOT file and PDF
	output_file_root->cd();
	plot->Write();
	normHist->Write();
	canvas->Write();

	canvas->Print(out_name_pdf.c_str(), "pdf");

	// save histogram PDF
	TCanvas *histCanvas = new TCanvas("cHist",title.c_str(),600,600);
	histCanvas->cd();

	normHist->Draw();

	std::string out_name_hist_pdf = out_name + "_histogram.pdf";
	histCanvas->Print(out_name_hist_pdf.c_str(), "pdf");

	// close this secondary canvas
	if (histCanvas) { histCanvas->Close(); gSystem->ProcessEvents(); }
}