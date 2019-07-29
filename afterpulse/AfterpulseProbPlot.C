#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TProfile.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

void AfterpulseProbPlot(std::string input_file_name, std::string title="", std::string x_title="", 
						std::string y_title="", Int_t nBins=100/*, Double_t start=0, Double_t end=100*/){
	// prep output files
	std::string out_name = input_file_name.substr(0, input_file_name.size() - 4) + "_histogram";
	std::string out_name_pdf = out_name + ".pdf";
	std::string out_name_root = out_name + ".root";

	TFile *output_file_root = new TFile(out_name_root.c_str(), "recreate");

	// load input file
	std::ifstream input_file(input_file_name);

	// prep plot to display
	TH2 *plot = new TH2F("afterpulse_plot",title.c_str(),nBins,0,35,nBins,0,60);
	plot->SetXTitle(x_title.c_str());
	plot->SetYTitle(y_title.c_str());
	//plot->SetLineColor(1); // make lines black
	plot->SetMarkerStyle(20);
	plot->SetMarkerSize(0.4);

	// read input file line by line and fill histogram
	std::string line;
	Double_t PMT, nPE, gain, prob, current;
	int i = 0;
	while (std::getline(input_file, line)){
		std::istringstream iss(line);
		do {
			if(!(iss >> PMT >> nPE >> gain >> prob >> current)){continue;};

			current = nPE * gain;

			if(current > 2){
					plot->Fill(current, prob*100);
			}	
			// if (i < 20){
			// 	cout << PMT << " " << nPE << " " << gain << " " << prob << endl;
			// 	i++;
			// }
		}
		while (iss);
	}

	// create canvas
	TCanvas *canvas = new TCanvas("c1", title.c_str(), 600, 600);

	//plot->Fit("pol1");
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);

	// draw histogram
	plot->Draw("scat");

	// save ROOT file and PDF
	output_file_root->cd();
	plot->Write();

	canvas->Print(out_name_pdf.c_str(), "pdf");
}