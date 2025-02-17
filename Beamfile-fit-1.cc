#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTree.h"
#include <TCanvas.h>
#include <TFile.h>
#include "TMath.h"
#include <TGraph.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

double bivariate_gaussian(double *x, double *par){

        // par[0] = Amplitude
        // par[1] = Mean X
        // par[2] = Mean Y
        // par[3] = Sigma X
        // par[4] = Sigma Y
        // par[5] = Correlation coefficient

        // par[0] = 1./(2*TMath::Pi()*par[3]*par[4]*par[5]);
        double dx = (x[0] - par[1])/par[3];
        double dy = (x[1] - par[2])/par[4];
        double z = (dx*dx - 2*par[5]*dx*dy + dy*dy)/(2*(1.-par[5]*par[5]));
        return par[0]*TMath::Exp(-z);
}

void Beamfile_fit(){


        TFile hfileB("treeB.root");

        // random trigger data
        static TH1D *beam[10];
        static TH2D *correlation[10];
        static double normX;
        static double meanX;
        static double sigmaX;
        static double normY;
        static double meanY;
        static double sigmaY;
        static double normdXdZ;
        static double meandXdZ;
        static double sigmadXdZ;
        static double normdYdZ;
        static double meandYdZ;
        static double sigmadYdZ;

        // parameters for bivariate gaussian fit
        static double par0X;
        static double par1X;
        static double par2X;
        static double par3X;
        static double par4X;
        static double par5X;

        static double par0Y;
        static double par1Y;
        static double par2Y;
        static double par3Y;
        static double par4Y;
        static double par5Y;

        double beamX;
        double beamY;
        double dXdZ;
        double dYdZ;

        // new generated data from gaussian fit
        static TH1D *GeneratedBeam[10];
        static TH2D *GeneratedBivGauss[10];                     // 2D histos generated from tf2 bivariate_gaussian
        double genX;
        double genY;
        double slopeX;
        double slopeY;
        // used for GetRandom2
        double bivGaussX;
        double bivGaussdXdZ;
        double bivGaussY;
        double bivGaussdYdZ;

        // for storing generated data
        TFile hfileC("treeC.root","RECREATE");	
        TTree *ctree = new TTree("CTree","CTree");
//        ctree->Branch("genX", &genX,"genX/D");
//        ctree->Branch("genY", &genY,"genY/D");
//        ctree->Branch("slopeX", &slopeX,"slopeX/D");
//        ctree->Branch("slopeY", &slopeY, "slopeY/D");
        ctree->Branch("genBivX", &bivGaussX,"genBivX/D");
        ctree->Branch("genBivY", &bivGaussY,"genBivY/D");
        ctree->Branch("genBivSlopeX", &bivGaussdXdZ,"genBivSlopeX/D");
        ctree->Branch("genBivSlopeY", &bivGaussdYdZ, "genBivSlopeY/D");

        TCanvas* C = nullptr;
        TCanvas* CC = nullptr;
        TCanvas* Cg = nullptr;
        TCanvas* Cg1 = nullptr;
        TCanvas* Cg2 = nullptr;

        int Nevent;

static bool first = false; 

if(first == true){
        // canvases referring to real random-triggered 2023 data

        C = new TCanvas("C","variables",200,10,600,400);
        C->Divide(2,2);

        CC = new TCanvas("CC","correlations",200,10,600,400);
        CC -> Divide(3,2);

        // new generated data
        Cg = new TCanvas("Cg","Generated beamfile data",200,10,600,400);
        Cg->Divide(2,2);
        Cg1 = new TCanvas("Cg1","Generated beamfile data bivGauss",200,10,600,400);
        Cg1->Divide(2,2);
        Cg2 = new TCanvas("Cg2","Generated correlations bivGauss",200,10,600,400);
        Cg2->Divide(2,2);
}

        
       // reading the TTree and its branches

        TTree *tree2 = (TTree*)hfileB.Get("BTree");
        TBranch *b1=tree2->GetBranch("X");
        TBranch *b2=tree2->GetBranch("Y");
        TBranch *b3=tree2->GetBranch("dXdZ");
        TBranch *b4=tree2->GetBranch("dYdZ");

        b1->SetAddress(&beamX);
        b2->SetAddress(&beamY);
        b3->SetAddress(&dXdZ);
        b4->SetAddress(&dYdZ);


        Nevent = tree2->GetEntries();

        int nbin = 100;
        // histograms
        beam[0] = new TH1D("beamprofile_X","X profile; X [cm]",nbin,-2.,2.);
        beam[1] = new TH1D("beamprofile_Y","Y profile; Y [cm]",nbin,-2.,2.);
        beam[2] = new TH1D("beamprofile_dXdZ","slope X; dXdZ [mrad]",nbin,-1.,1.);
        beam[3] = new TH1D("beamprofile_dYdZ","slope Y; dYdZ [mrad]",nbin,-1.,1.);

        GeneratedBeam[0] = new TH1D("X", "Generated X profile; X [cm];", 100, -2., 2.);      
        GeneratedBeam[1] = new TH1D("Y", "Generated Y profile; Y [cm];", 100, -2., 2.);      
        GeneratedBeam[2] = new TH1D("slopeX", "Generated dXdZ profile; dXdZ [mrad];", 100, -1., 1.);      
        GeneratedBeam[3] = new TH1D("slopeY", "Generated dYdZ profile; dYdZ [mrad];", 100, -1., 1.);

        GeneratedBeam[4] = new TH1D("XbivGaus", "Generated X profile from bivGaus; X [cm];", 100, -2., 2.);      
        GeneratedBeam[5] = new TH1D("YbivGaus", "Generated Y profile from bivGaus; Y [cm];", 100, -2., 2.);      
        GeneratedBeam[6] = new TH1D("slopeXbivGaus", "Generated dXdZ profile from bivGaus; dXdZ [mrad];", 100, -1., 1.);      
        GeneratedBeam[7] = new TH1D("slopeYbivGaus", "Generated dYdZ profile from bivGaus; dYdZ [mrad];", 100, -1., 1.);

        correlation[0] = new TH2D("X-Y", "correlation X-Y; X [cm]; Y [cm]", nbin, -2.,2., nbin, -2.,2.);
        correlation[0] -> SetOption("colz");
        correlation[1] = new TH2D("X-slopeX", "correlation X-slopeX; X [cm]; slopeX [mrad]", nbin, -2.,2., nbin,-1.,1.);
        correlation[1] -> SetOption("colz");
        correlation[2] = new TH2D("X-slopeY", "correlation X-slopeY; X [cm]; slopeY [mrad]", nbin, -2.,2., nbin,-1.,1.);
        correlation[2] -> SetOption("colz");
        correlation[3] = new TH2D("Y-slopeX", "correlation Y-slopeX; Y [cm]; slopeX [mrad]", nbin, -2.,2., nbin,-1.,1.);
        correlation[3] -> SetOption("colz");
        correlation[4] = new TH2D("Y-slopeY", "correlation Y-slopeY; Y [cm]; slopeY [mrad]", nbin, -2.,2., nbin,-1.,1.);
        correlation[4] -> SetOption("colz");
        correlation[5] = new TH2D("slopeX-slopeY", "correlation slopeX-slopeY; slopeX [mrad]; slopeY [mrad]", nbin, -1.,1., nbin,-1.,1.);
        correlation[5] -> SetOption("colz");

        GeneratedBivGauss[0] = new TH2D("X-slopeX_bivgauss", "generated X-slopeX; X [cm]; dX/dZ [mrad]", nbin, -2.,2., nbin, -1.,1.);
        GeneratedBivGauss[0] -> SetOption("colz");
        GeneratedBivGauss[1] = new TH2D("Y-slopeY_bivgauss", "generated Y-slopeY; Y [cm]; dY/dZ [mrad]", nbin, -2.,2., nbin, -1.,1.);
        GeneratedBivGauss[1] -> SetOption("colz");

        for(int ev=0; ev<Nevent; ev++){

                tree2->GetEvent(ev);
                beam[0] -> Fill(beamX); 
                beam[1] -> Fill(beamY);
                beam[2] -> Fill(dXdZ);
                beam[3] -> Fill(dYdZ);

                correlation[0] -> Fill(beamX,beamY);
                correlation[1] -> Fill(beamX,dXdZ);
                correlation[2] -> Fill(beamX,dYdZ);
                correlation[3] -> Fill(beamY,dXdZ);
                correlation[4] -> Fill(beamY,dYdZ);
                correlation[5] -> Fill(dXdZ,dYdZ);
        }

         gStyle->SetOptFit(1101);

        // ***************************************************** fit with gaussian *****************************************************
if(first) C->cd(1);
        TF1 *gausX = new TF1("gausX","gaus",-0.9,0.7);
        beam[0] -> Fit(gausX,"LER");        
        normX = gausX->GetParameter(0);      
        meanX = gausX->GetParameter(1);
        sigmaX = gausX->GetParameter(2);

if(first) C->cd(2);
        TF1 *gausY = new TF1("gausY","gaus",-0.8,0.7); 
        beam[1] -> Fit(gausY,"LER");       
        normY = gausY->GetParameter(0);       
        meanY = gausY->GetParameter(1);
        sigmaY = gausY->GetParameter(2);

if(first) C->cd(3);
        TF1 *gausdXdZ = new TF1("gausdXdZ","gaus",-0.15,0.2);
        beam[2] -> Fit(gausdXdZ,"LER");    
        normdXdZ = gausdXdZ->GetParameter(0);           
        meandXdZ = gausdXdZ->GetParameter(1);
        sigmadXdZ = gausdXdZ->GetParameter(2);

if(first) C->cd(4);
        TF1 *gausdYdZ = new TF1("gausdYdZ","gaus",-0.15,0.2);   
        beam[3] -> Fit(gausdYdZ,"LER"); 
        normdYdZ = gausdYdZ->GetParameter(0);            
        meandYdZ = gausdYdZ->GetParameter(1);
        sigmadYdZ = gausdYdZ->GetParameter(2);

 if(first) {
        C -> Update();
        C->SaveAs("/eos/home-g/gmeinard/BeamfileFit/variables.png");
        C->Close();
        delete C;
 }
        // __________________________________________________ 2D histos ___________________________________________________
if(first) CC->cd(1);
        correlation[0] -> Draw("");
        gPad->SetLogz();

if(first) CC->cd(2);
        TF2 *gaus2DX = new TF2("bivgauss_fitX", bivariate_gaussian,-0.9,0.7,-0.15,0.2,6);

        gaus2DX->SetParLimits(0, 0, 1000);     // Amplitude
        gaus2DX->SetParLimits(1, -0.09, -0.078);     // Mean X
        gaus2DX->SetParLimits(2, 0.02, 0.025);      // Mean Y
        gaus2DX->SetParLimits(3, 0.4, 0.45);      // Sigma X
        gaus2DX->SetParLimits(4, 0.08, 0.09);     // Sigma Y
        gaus2DX->SetParLimits(5, -1, 1);        // correlation coefficient

        //gaus2DX->SetParameters(0,500);
        gaus2DX->SetParameters(1,-0.078);
        gaus2DX->SetParameters(2,0.023);
        gaus2DX->SetParameters(3,0.44);
        gaus2DX->SetParameters(4,0.080);
        gaus2DX->SetParameters(5,-0.27);

if(first) correlation[1]->Draw("");
        gPad->SetLogz();
        correlation[1]->Fit(gaus2DX,"SLER+");

        par0X = gaus2DX->GetParameter(0);
        par1X = gaus2DX->GetParameter(1);
        par2X = gaus2DX->GetParameter(2);
        par3X = gaus2DX->GetParameter(3);
        par4X = gaus2DX->GetParameter(4);
        par5X = gaus2DX->GetParameter(5);

if(first){ 
        Cg2->cd(1);
        gaus2DX->Draw("surf1");

        CC->cd(3);
        correlation[2] -> Draw("");
        gPad->SetLogz();

        CC->cd(5);
        correlation[3]-> Draw("");
        gPad->SetLogz();

        CC->cd(6);
}
        TF2 *gaus2DY = new TF2("bivgauss_fitY", bivariate_gaussian,-0.8,0.7,-0.15,0.2,6);

        gaus2DY->SetParLimits(0, 0, 1000);              // Amplitude
        gaus2DY->SetParLimits(1, -0.5, 0);              // Mean X
        gaus2DY->SetParLimits(2, 0, 0.05);               // Mean Y
        gaus2DY->SetParLimits(3, 0., 0.5);           // Sigma X
        gaus2DY->SetParLimits(4, 0.06, 0.09);            // Sigma Y
        gaus2DY->SetParLimits(5, -1, 1);              // correlation coefficient

        //gaus2DY->SetParameters(0,563);                  // amplitude
        gaus2DY->SetParameters(1,-0.12);               // mean x
        gaus2DY->SetParameters(2,0.040);                // mean y
        gaus2DY->SetParameters(3,0.33);                 // sigma x
        gaus2DY->SetParameters(4,0.088);                // sigma y
        //gaus2DY->SetParameters(5,-0.05);                // correlation coefficient

if(first) correlation[4] ->Draw("");

        gPad->SetLogz();
        correlation[4] ->Fit(gaus2DY,"SLER+");

        par0Y = gaus2DY->GetParameter(0);
        par1Y = gaus2DY->GetParameter(1);
        par2Y = gaus2DY->GetParameter(2);
        par3Y = gaus2DY->GetParameter(3);
        par4Y = gaus2DY->GetParameter(4);
        par5Y = gaus2DY->GetParameter(5);

if(first){
        Cg2->cd(3);
        gaus2DY->Draw("surf1");

        CC->cd(4);
        correlation[5] -> Draw("");
        gPad->SetLogz();
        CC -> Update();
        CC -> SaveAs("/eos/home-g/gmeinard/BeamfileFit/correlations.png");
        CC->Close();
        delete CC;
}

        // Gaussian functions
        TF1 *gausGenX = new TF1("gausGenX","gaus",-2.,2.);
        gausGenX->SetParameter(0, normX);       // normalization
        gausGenX->SetParameter(1, meanX);       // mean
        gausGenX->SetParameter(2, sigmaX);      // sigma
                    
        TF1 *gausGenY = new TF1("gausGenY","gaus",-2.,2.);
        gausGenY->SetParameter(0, normY);       // normalization
        gausGenY->SetParameter(1, meanY);       // mean
        gausGenY->SetParameter(2, sigmaY);      // sigma
              
        TF1 *gausGenSlopeX = new TF1("gausGenSlopeX","gaus",-0.6,0.6);
        gausGenSlopeX->SetParameter(0, normdXdZ);       // normalization
        gausGenSlopeX->SetParameter(1, meandXdZ);       // mean
        gausGenSlopeX->SetParameter(2, sigmadXdZ);      // sigma
              
        TF1 *gausGenSlopeY = new TF1("gausGenSlopeY","gaus",-0.6,0.6);
        gausGenSlopeY->SetParameter(0, normdYdZ);       // normalization
        gausGenSlopeY->SetParameter(1, meandYdZ);       // mean
        gausGenSlopeY->SetParameter(2, sigmadYdZ);      // sigma

        // Bivariate Gaussian Functions
        TF2 *gausGen2DX = new TF2("bivgaussGen_fitX", bivariate_gaussian,-2.0,2.0,-0.6,0.6,6);
        gausGen2DX->SetParameter(0, par0X);
        gausGen2DX->SetParameter(1, par1X);
        gausGen2DX->SetParameter(2, par2X);
        gausGen2DX->SetParameter(3, par3X);
        gausGen2DX->SetParameter(4, par4X);
        gausGen2DX->SetParameter(5, par5X);

        TF2 *gausGen2DY = new TF2("bivgaussGen_fitY", bivariate_gaussian,-1.5,1.5,-0.5,0.5,6);
        gausGen2DY->SetParameter(0, par0Y);
        gausGen2DY->SetParameter(1, par1Y);
        gausGen2DY->SetParameter(2, par2Y);
        gausGen2DY->SetParameter(3, par3Y);
        gausGen2DY->SetParameter(4, par4Y);
        gausGen2DY->SetParameter(5, par5Y);
        
        
        // random events generation
        const int nEvents = 10000000;
    
        for (int i = 0; i < nEvents; ++i) {

                genX = gausGenX->GetRandom();  
                GeneratedBeam[0]->Fill(genX);

                genY = gausGenY->GetRandom();
                GeneratedBeam[1]->Fill(genY);

                slopeX = gausGenSlopeX->GetRandom();
                GeneratedBeam[2]->Fill(slopeX);

                slopeY = gausGenSlopeY->GetRandom();
                GeneratedBeam[3]->Fill(slopeY);

                gausGen2DX->GetRandom2(bivGaussX,bivGaussdXdZ);            // extrapolates rndm pair from 2d histo
                GeneratedBivGauss[0]->Fill(bivGaussX,bivGaussdXdZ);

                gausGen2DY->GetRandom2(bivGaussY,bivGaussdYdZ);
                GeneratedBivGauss[1]->Fill(bivGaussY,bivGaussdYdZ);

                GeneratedBeam[4] = GeneratedBivGauss[0]->ProjectionX("projX_gen0");                 // X
                GeneratedBeam[5] = GeneratedBivGauss[1]->ProjectionX("projX_gen1");                 // Y
                GeneratedBeam[6] = GeneratedBivGauss[0]->ProjectionY("projY_gen0");                 // slope X
                GeneratedBeam[7] = GeneratedBivGauss[1]->ProjectionY("projY_gen1");                 // slope Y
               
        
                ctree->Fill();
        }

if(first){
        Cg->cd(1);
        GeneratedBeam[0]->Draw("");
        Cg->cd(2);
        GeneratedBeam[1]->Draw("");
        Cg->cd(3);
        GeneratedBeam[2]->Draw("");
        Cg->cd(4);
        GeneratedBeam[3]->Draw("");
        Cg->Update();
        Cg->SaveAs("/eos/home-g/gmeinard/BeamfileFit/GeneratedBeam.png");
        Cg->Close();
        delete Cg;

        // for generated variables bivGauss
        Cg1->cd(1);
        GeneratedBeam[4]->Draw("");
        Cg1->cd(2);
        GeneratedBeam[5]->Draw("");
        Cg1->cd(3);
        GeneratedBeam[6]->Draw("");
        Cg1->cd(4);
        GeneratedBeam[7]->Draw("");
        Cg1->Update();
        Cg1->SaveAs("/eos/home-g/gmeinard/BeamfileFit/GeneratedBeam1.png");
        Cg1->Close();
        delete Cg1;

        // for generated th2d from bivGauss fit
        Cg2->cd(2);
        GeneratedBivGauss[0]->Draw("");
        gPad->SetLogz();
        Cg2->cd(4);
        GeneratedBivGauss[1]->Draw("");
        gPad->SetLogz();
        Cg2->Update();
        Cg2->SaveAs("/eos/home-g/gmeinard/BeamfileFit/GeneratedBivGauss.png");
        Cg2->Close();
        delete Cg2;
}
        ctree->Write();

}