#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"
#include "PaMetaDB.h"
#include "G3part.h"
#include <fstream>
#include <iomanip>
#include <string>

void UserEvent109(PaEvent& e){

       //apertura file di input
        TFile hfileB("treeB.root");
        ofstream outFile("fit_parameters.txt");

        static TH1D* beam[10];
        static double meanX;
        static double sigmaX;
        static double meanY;
        static double sigmaY;
        static double meandXdZ;
        static double sigmadXdZ;
        static double meandYdZ;
        static double sigmadYdZ;
        int Nevent;
        TCanvas *C1 = new TCanvas("C1","X profile",200,10,600,400);
        TCanvas *C2 = new TCanvas("C2","Y profile",200,10,600,400);
        TCanvas *C3 = new TCanvas("C3","dXdZ profile",200,10,600,400);
        TCanvas *C4 = new TCanvas("C4","dYdZ profile",200,10,600,400);
        
       //lettura TTree  e branch


        TTree *tree2 = (TTree*)hfile2.Get("BTree");
        TBranch *b1=tree2->GetBranch("X");
        TBranch *b2=tree2->GetBranch("Y");
        TBranch *b3=tree2->GetBranch("dXdZ");
        TBranch *b4=tree2->GetBranch("dYdZ");

        b1->SetAddress(&beamX);
        b2->SetAddress(&beamy);
        b3->SetAddress(&dXdZ);
        b4->SetAddress(&dYdZ);

        Nevent = tree2->GetEntries();
        cout << "entries tree " << Nevent << endl;

        // histograms
        beam[0] = new TH1D("beamprofile_X","X profile (cm)",nbin,-5.,5.);
        beam[1] = new TH1D("beamprofile_Y","Y profile (cm)",nbin,-5.,5.);
        beam[2] = new TH1D("beamprofile_dXdZ","dXdZ (mrad)",nbin,-5.,5.);
        beam[3] = new TH1D("beamprofile_dYdZ","dYdZ (mrad)",nbin,-5.,5.);

        for(int ev=0; ev<Neventi; ev++){

                // tree2->GetEvent(ev);
                beam[0] -> Fill(X); 
                beam[1] -> Fill(Y);
                beam[2] -> Fill(dXdZ);
                beam[3] -> Fill(dYdZ);
        }


        // fit with gaussian
        C1->cd();
        beam[0] -> Draw("");
        TF1 *gausX = new TF1("gausX","gaus",-1.5,1.5);              
        meanX = gaussX->GetParameter(1);
        sigmaX = gaussX->GetParameter(2);
        beam[0] -> Fit(gausX);
        C1 -> Update();

        C2->cd();
        beam[1] -> Draw("");
        TF1 *gausY = new TF1("gausY","gaus",-1,1);              
        meanY = gaussY->GetParameter(1);
        sigmaY = gaussY->GetParameter(2);
        beam[1] -> Fit(gausX);
        C2 -> Update();

        C3->cd();
        beam[2] -> Draw("");
        TF1 *gausdXdZ = new TF1("gausdXdZ","gaus",-0.6,0.6);              
        meandXdZ = gaussdXdZ->GetParameter(1);
        sigmadXdZ = gaussdXdZ->GetParameter(2);
        beam[2] -> Fit(gausX);
        C3 -> Update();

        C4->cd();
        beam[3] -> Draw("");
        TF1 *gausdYdZ = new TF1("gausdYdZ","gaus",-0.6,0.6);              
        meandYdZ = gaussdYdZ->GetParameter(1);
        sigmadYdZ = gaussdYdZ->GetParameter(2);
        beam[3] -> Fit(gausX);
        C3 -> Update();

        // debug
        cout << " fit results on x: media " << meanX << " sigma " << sigmaX << endl;
        cout << " fit results on y: media " << meanY << " sigma " << sigmaY << endl;
        cout << " fit results on dxdz: media " << meandXdZ << " sigma " << sigmadXdZ << endl;
        cout << " fit results on dydz: media " << meandYdZ << " sigma " << sigmadYdZ << endl;

        outFile << meanX << " "    // Mean
                << sigmaX << " "    // Sigma
                << endl;
                << meanY << " "    // Mean
                << sigmaY << " "    // Sigma
                << endl;
                << meandXdZ << " "    // Mean
                << sigmadXdZ << " "    // Sigma
                << endl;
                << meandYdZ << " "    // Mean
                << sigmadYdZ << " "    // Sigma
                << endl;

        outFile.close();


}