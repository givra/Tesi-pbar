#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTree.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

void Beamfile_correction(){

    static TH1D *GeneratedBeam[4];
    double genX;
    double genY;
    double slopeX;
    double slopeY;

    double normX;
    double meanX; 
    double sigmaX;
    double normY;
    double meanY; 
    double sigmaY;
    double normdXdZ;
    double meandXdZ;
    double sigmadXdZ;
    double normdYdZ;
    double meandYdZ;
    double sigmadYdZ;

    TFile hfileC("treeC.root","RECREATE");	
    TTree *ctree = new TTree("CTree","CTree");
    ctree->Branch("genX", &genX,"genX/D");
    ctree->Branch("genY", &genY,"genY/D");
    ctree->Branch("slopeX", &slopeX,"slopeX/D");
    ctree->Branch("slopeY", &slopeY, "slopeY/D");

    TCanvas *C = new TCanvas("C","Generated beamfile data",200,10,600,400);
    C->Divide(2,2);

    GeneratedBeam[0] = new TH1D("X", "Generated Gaussian Events; X [cm];", 100, -2., 2.);      // meanX - 5 * sigmaX, meanX + 5 * sigmaX
    GeneratedBeam[1] = new TH1D("Y", "Generated Gaussian Events; Y [cm];", 100, -2., 2.);      // meanY - 5 * sigmaY, meanY + 5 * sigmaY
    GeneratedBeam[2] = new TH1D("slopeX", "Generated Gaussian Events; dXdZ [mrad];", 100, -1., 1. );      // meandXdZ - 5 * sigmadXdZ, meandXdZ + 5 * sigmadXdZ
    GeneratedBeam[3] = new TH1D("slopeY", "Generated Gaussian Events; dYdZ [mrad];", 100, -1., 1.);       // meandYdZ - 5 * sigmadYdZ, meandYdZ + 5 * sigmadYdZ


    ifstream infile("fit_parameters.txt");  // Replace with your file path
    if(!infile){
        cerr << "Error: Could not open file 'fit_parameters.txt'!" << endl;
        return;
    }

    if (infile.is_open()) {
        infile >> normX >> meanX >> sigmaX >> normY >> meanY >> sigmaY >> normdXdZ >> meandXdZ >> sigmadXdZ >> normdYdZ >> meandYdZ >> sigmadYdZ; 
        infile.close();
    }

    //debug
    cout << "Read meanX = " << meanX << ", sigmaX = " << sigmaX << endl;
    cout << "Read meanY = " << meanY << ", sigmaY = " << sigmaY << endl;
    cout << "Read meandXdZ = " << meandXdZ << ", sigmadXdZ = " << sigmadXdZ << endl;
    cout << "Read meandYdZ = " << meandYdZ << ", sigmadYdZ = " << sigmadYdZ << endl;

    // Gaussian function
    TF1 *gausX = new TF1("gausX","gaus",-1.5,1.5);
    gausX->SetParameter(0, normX);       // normalization
    gausX->SetParameter(1, meanX);       // mean
    gausX->SetParameter(2, sigmaX);      // sigma


    TF1 *gausY = new TF1("gausY","gaus",-1.,1.);
    gausY->SetParameter(0, normY);       // normalization
    gausY->SetParameter(1, meanY);       // mean
    gausY->SetParameter(2, sigmaY);      // sigma

    TF1 *gausdXdZ = new TF1("gausdXdZ","gaus",-0.6,0.6);
    gausdXdZ->SetParameter(0, normdXdZ);       // normalization
    gausdXdZ->SetParameter(1, meandXdZ);       // mean
    gausdXdZ->SetParameter(2, sigmadXdZ);      // sigma

    TF1 *gausdYdZ = new TF1("gausdYdZ","gaus",-0.6,0.6);
    gausdYdZ->SetParameter(0, normdYdZ);       // normalization
    gausdYdZ->SetParameter(1, meandYdZ);       // mean
    gausdYdZ->SetParameter(2, sigmadYdZ);      // sigma


    // random events generation using the Gaussian function
    const int nEvents = 1000000;
    
    for (int i = 0; i < nEvents; ++i) {

        genX = gausX->GetRandom();  // Generate a random number from the Gaussian
        GeneratedBeam[0]->Fill(genX);

        genY = gausY->GetRandom();
        GeneratedBeam[1]->Fill(genY);

        slopeX = gausdXdZ->GetRandom();
        GeneratedBeam[2]->Fill(slopeX);

        slopeY = gausdYdZ->GetRandom();
        GeneratedBeam[3]->Fill(slopeY);

        ctree->Fill();
    }

    // Drawing the histogram and the Gaussian curve
    C->cd(1);
    GeneratedBeam[0]->Draw();
    // gausX->SetLineColor(kRed);
    // gausX->Draw("SAME");

    C->cd(2);
    GeneratedBeam[1]->Draw();
    // gausY->SetLineColor(kRed);
    // gausY->Draw("SAME");

    C->cd(3);
    GeneratedBeam[2]->Draw();
    // gausdXdZ->SetLineColor(kRed);
    // gausdXdZ->Draw("SAME");

    C->cd(4);
    GeneratedBeam[3]->Draw();
    // gausdYdZ->SetLineColor(kRed);
    // gausdYdZ->Draw("SAME");

    C->Update();
    C->SaveAs("/eos/home-g/gmeinard/BeamfileFit/GeneratedBeam.png");
    ctree->Write(); 

}