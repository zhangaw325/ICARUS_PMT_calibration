#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

void HistogramAfterpulse(std::string input_file_name){
	// prep output files
	std::string out_name = "afterpulse_probabilities";
	std::string out_name_pdf = out_name + ".pdf";
	std::string out_name_root = out_name + ".root";

	TFile *output_file_root = new TFile(out_name_root.c_str(), "recreate");

	// prep histogram
	std::string title = "Afterpulse probabilities";

	// load input file
	std::ifstream input_file(input_file_name);

	// prep histogram to display
	TH1F *hist = new TH1F("probabilities", title.c_str(), 100, 0, 100);
	hist->SetXTitle("afterpulse probability (%)");
	hist->SetYTitle("number of PMTs");
	hist->SetLineColor(1); // make lines black

	// read input file line by line and fill histogram
	std::string line;
	while (std::getline(input_file, line)){
		std::istringstream iss(line);
		Double_t p;
		if (!(iss >> p)) { continue; } // pass lines that aren't a number
	    hist->Fill(p*100);
	}

	// create canvas
	TCanvas *canvas = new TCanvas("c1", title.c_str(), 600, 600);

	// draw histogram
	hist->Draw();

	// save ROOT file and PDF
	output_file_root->cd();
	hist->Write();

	canvas->Print(out_name_pdf.c_str(), "pdf");
}