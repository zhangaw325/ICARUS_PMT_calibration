#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "THStack.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


// output file names
const std::string out_name = "dark_rates";
const std::string out_name_pdf = out_name + ".pdf";
const std::string out_name_root = out_name + ".root";

// histogram name
const std::string title = "Dark rates";

// row by row
void HistogramDarkRate(	std::string input_file_name_1, std::string input_file_name_2,
						std::string input_file_name_3, std::string input_file_name_4,
						bool stackHist=true){
	// prep output file
	TFile *output_file_root = new TFile(out_name_root.c_str(), "recreate");

	// create canvas
	TCanvas *canvas = new TCanvas("c1", title.c_str(), 800, 600);

	// load input files
	std::ifstream input_file[4];
	input_file[0] = std::ifstream(input_file_name_1);
	input_file[1] = std::ifstream(input_file_name_2);
	input_file[2] = std::ifstream(input_file_name_3);
	input_file[3] = std::ifstream(input_file_name_4);

	// prep stack
	THStack *stack = new THStack("stack", "");
	stack->SetTitle(title.c_str());

	// create histograms
	TH1F *hist[4];
	for(int i = 0; i < 4; i++){
		std::string plot_title;
		if(i == 0) plot_title = "Row A";
		if(i == 1) plot_title = "Row B";
		if(i == 2) plot_title = "Row C";
		if(i == 3) plot_title = "Row D";

		hist[i] = new TH1F("dark_rates", plot_title.c_str(), 100, 0, 15);
		hist[i]->SetXTitle("dark rate (kHz)");
		hist[i]->SetYTitle("number of PMTs");
		hist[i]->SetLineColor(i+1);

		// read input file line by line
		std::string line;
		while (std::getline(input_file[i], line)){
			std::istringstream iss(line);

			Double_t dr;
			if(!(iss >> dr) || dr == 0) {continue;} // pass lines that aren't a non-zero number

			hist[i]->Fill(dr);
		}

		// add histogram to stack, or draw accordingly
		hist[i]->Scale(1/(hist[i]->GetEntries()));
		if(stackHist){
			stack->Add(hist[i]);
		} else {
			if(i == 0){
				hist[i]->Draw("hist");
			} else {
				hist[i]->Draw("samehist");
			}
		}
	}

	if(stackHist){
		stack->Draw("hist");
	}

	stack->GetXaxis()->SetTitle("dark rate (kHz)");
	stack->GetYaxis()->SetTitle("proportion of PMTs");

	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	canvas->Update();

	output_file_root->cd();
	stack->Write();

	canvas->Print(out_name_pdf.c_str(), "pdf");
}

// two files
void HistogramDarkRate(std::string input_file_name_1, std::string input_file_name_2){
	// prep output file
	TFile *output_file_root = new TFile(out_name_root.c_str(), "recreate");

	// create canvas
	TCanvas *canvas = new TCanvas("c1", title.c_str(), 600, 600);

	// load input files
	std::ifstream input_file_1(input_file_name_1);
	std::ifstream input_file_2(input_file_name_2);

	// prep histograms
	TH1F *hist1 = new TH1F("dark_rates", title.c_str(), 100, 0, 25);
	hist1->SetXTitle("dark rate (kHz)");
	hist1->SetYTitle("number of PMTs");
	hist1->SetLineColor(4); // make lines blue

	TH1F *hist2 = new TH1F("dark_rates2", title.c_str(), 100, 0, 25);
	hist2->SetLineColor(2); // make lines orange

	// read input file 1 line by line and fill histogram
	std::string line1;
	while (std::getline(input_file_1, line1)){
		std::istringstream iss(line1);

		Double_t dr;
		if (!(iss >> dr)) { continue; } // pass lines that aren't a number

	    hist1->Fill(dr);
	}

	// read input file 2 line by line and fill histogram
	std::string line2;
	while (std::getline(input_file_2, line2)){
		std::istringstream iss(line2);

		Double_t dr;
		if (!(iss >> dr)) { continue; } // pass lines that aren't a number
		//cout << dr << endl;
	    hist2->Fill(dr);
	}

	// draw first histogram
	hist1->Draw();
	canvas->Update();

	// draw second histogram
	hist2->Draw("same");
	canvas->Update();

	// create legend
	TLegend *legend = new TLegend(0.6, 0.5, 0.8, 0.7);
	legend->AddEntry(hist1, "May");
	legend->AddEntry(hist2, "June");
	legend->SetTextSize(0.05);
	legend->Draw();
	canvas->Update();
}

// single file
void HistogramDarkRate(std::string input_file_name){
	// prep output file
	TFile *output_file_root = new TFile(out_name_root.c_str(), "recreate");

	// load input file
	std::ifstream input_file(input_file_name);

	// prep histogram to display
	TH1F *hist = new TH1F("dark_rates", title.c_str(), 100, 0, 25);
	hist->SetXTitle("dark rate (kHz)");
	hist->SetYTitle("number of PMTs");
	hist->SetLineColor(1); // make lines black

	// read input file line by line and fill histogram
	std::string line;
	while (std::getline(input_file, line)){
		std::istringstream iss(line);

		Double_t dr;
		if (!(iss >> dr) || dr == 0) { continue; } // pass lines that aren't a number

	    hist->Fill(dr);
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