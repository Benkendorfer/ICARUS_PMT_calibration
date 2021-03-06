//***************************************************************************
// Plots gain versus voltage and performs a power log fit on data
//
// - Input must be CHIMNEY.txt
// - Input format must be:
//     PMT# Voltage Gain GainError
//
// - Output:
//    - PDF file containing log-log fit (and linear display) for each PMT
//    - ROOT file containing all gain-voltage fits
//    - comma-separated text file containing fit parameters for each PMT, 
//      in increasing numerical PMT order
//***************************************************************************

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TFrame.h"

#include <fstream>

const int NPMT = 10;

void GainVoltage(std::string chimney){

  // some constants
  const double fbegin = 1000;   // begin fit
  const double fend = 2000;     // end fit
  const int NPAR = 2;           // number of parameters
  
  // generate output PDF file names
  std::string outnamepdf[10];

  for(int i = 0; i < NPMT; i++){
    outnamepdf[i] = chimney + "_" + std::to_string(i + 1) + "_gainvsvoltage.pdf";
  }

  // generate other output file names
  std::string outnametxt = chimney + "_gainvsvoltage.txt";
  std::string outnameroot = chimney + "_gainvsvoltage.root";
  TFile* outROOTfile = new TFile(outnameroot.c_str(),"recreate");  
  fstream fout(outnametxt.c_str(),ios::out);

  // prepare canvases for each fit
  TCanvas *c[NPMT];

  // analyze each PMT
  for(int pmt_num = 0; pmt_num < NPMT; pmt_num++){

    // read gain and voltage values for the PMT in question
    std::string input_file_name = chimney + ".txt";
    std::ifstream input_file(input_file_name);

    Double_t voltage_raw[6];
    Double_t voltage_error_raw[6] = {2,2,2,2,2,2};
    Double_t gain_raw[6];
    Double_t gain_error_raw[6];

    // read line by line
    double p, v, g, ge;
    int num_data_points = 0;

    while (input_file >> p >> v >> g >> ge){
      if (p == pmt_num + 1){
        voltage_raw[num_data_points] = v;
        gain_raw[num_data_points] = g*TMath::Power(10,7);
        gain_error_raw[num_data_points] = ge*TMath::Power(10,7);
        num_data_points++;
      }
    }

    // create properly sized arrays
    Double_t voltage[num_data_points];
    Double_t voltage_error[num_data_points];
    Double_t gain[num_data_points];
    Double_t gain_error[num_data_points];

    Double_t voltage_nolog[num_data_points];
    Double_t voltage_error_nolog[num_data_points];
    Double_t gain_nolog[num_data_points];
    Double_t gain_error_nolog[num_data_points];

    // populate properly sized arrays
    for(int i = 0; i < num_data_points; i++){
      voltage[i] = TMath::Log(voltage_raw[i]);
      gain[i] = TMath::Log(gain_raw[i]);
      gain_error[i] = gain_error_raw[i]/gain_raw[i];
      voltage_error[i] = voltage_error_raw[i]/voltage_raw[i];

      voltage_nolog[i] = voltage_raw[i];
      voltage_error_nolog[i] = voltage_error_raw[i];
      gain_nolog[i] = gain_raw[i];
      gain_error_nolog[i] = gain_error_raw[i];
    }

    if(num_data_points != 3 && num_data_points != 6){
      std::cout << "Improper number of data points for PMT " 
                << pmt_num + 1
                << ". SKIPPING" 
                << std::endl;
      continue;
    }

    // Perform linear fit on log-log plot
    TF1 *fit = new TF1("fit", "pol1", fbegin, fend);
    double par[NPAR];
    double parerr[NPAR];
    fit->SetParNames("Constant", "Exponent");
    fit->SetLineColor(2);
    fit->SetLineStyle(1);
    
    // Create canvas
    char tempname[100];
    sprintf(tempname, "c_%d",pmt_num);
    std::string canvasTitle = chimney + "_" + std::to_string(pmt_num + 1);
    c[pmt_num] = new TCanvas(tempname,canvasTitle.c_str(),200,10,600,700);
    c[pmt_num]->Divide(1,2);

    // Create fit graph of data
    c[pmt_num]->cd(1);
    c[pmt_num]->SetGrid();
    c[pmt_num]->GetFrame()->SetBorderSize(12);

    // Create graph of data
    TGraph* data = new TGraphErrors(num_data_points, voltage, gain, voltage_error, gain_error);
    std::string title =  "PMT " + chimney + "_" + std::to_string(pmt_num + 1) + " gain vs voltage (log);"
                    + "log(voltage [V]);"
                    + "log(gain)";
    data->SetTitle(title.c_str());
    data->SetMarkerStyle(kOpenSquare);

    // Fit
    fit->SetParameters(-30,7);
    std::cout << "Fitting " << pmt_num + 1 << std::endl;
    data->Fit("fit","ME");
    for(int j = 0; j < 9; j++){
      fit->GetParameters(par);
      fit->SetParameters(par[0],par[1]);
      data->Fit("fit","ME");
    }

    // extract fit parameters
    gStyle->SetOptFit();
    fit->GetParameters(par);
    Double_t constant = par[0];
    Double_t exponent = par[1];
    Double_t amplitude = TMath::Exp(constant);

    // write fit parameters to output text file
    for(int j = 0; j < 2; j++){
      // fill first two lines for each PMT with dashes
      fout  << "--" << "," << "--"
            << "," << "--" << "," << "--"
            << "," << "--"
            << "," << "--"
            << "," << "--"
            << std::endl;
    }
    // populate third line for each PMT with actual values
    fout << constant <<","<<fit->GetParError(0)
         << "," << exponent <<","<<fit->GetParError(1)
         << "," << fit->GetChisquare()
         << "," << fit->GetNDF()
         << "," << fit->GetProb()
         << std::endl;

    // Plot gain and voltage
    data->Draw("ap");

    // Create linear graph of data
    c[pmt_num]->cd(2);
    c[pmt_num]->SetGrid();
    c[pmt_num]->GetFrame()->SetBorderSize(12);

    // Plot points
    TGraph *dataLinear = new TGraphErrors(  num_data_points, voltage_nolog, 
                                        gain_nolog, voltage_error_nolog, 
                                        gain_error_nolog);
    title = "PMT " + chimney + "_" + std::to_string(pmt_num + 1) + " gain vs voltage (linear);"
                   + "voltage [V];"
                   + "gain";
    dataLinear->SetTitle(title.c_str());
    dataLinear->SetMarkerStyle(kOpenSquare);
    dataLinear->Draw("ap");

    // Draw fit result
    TF1 *gainFunc = new TF1("gainfunc",power,1000,2000,2);
    gainFunc->SetParameter(0,amplitude);
    gainFunc->SetParameter(1,exponent);
    gainFunc->Draw("SAME");

    // Write output files
    outROOTfile->cd();
    c[pmt_num]->Write();
    c[pmt_num]->Print(outnamepdf[pmt_num].c_str(),"pdf");
  }
}
