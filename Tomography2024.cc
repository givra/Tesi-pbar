#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <stdio.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TRandom.h"
#include <TCanvas.h>
#include <TFile.h>
#include "TMath.h"
#include <TGraph.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TLegend.h"

// _____________________________________________________________
//
//          This macro reads the tree USR108 from u108
//          to create the overlapping Z_sim/Z_rec distribution
//          and the tomography of the 2024 target
//______________________________________________________________

using namespace std;

void Tomography2024(){

    // TFile hfileT("treeT.root");
    // TFile *hfileT = new TFile("treeT.root", "READ");
    // hfileT->ls();

    TFile *hfileT = TFile::Open("treeT.root");

    TDirectoryFile *dir = (TDirectoryFile*)hfileT->Get("UserEvent108");
    if (!dir) {
        std::cerr << "ERROR: Directory UserEvent108 not found!" << std::endl;
        return;
    }

    //dir->ls();

    TH2D *h_XvsY_rec = (TH2D*)dir->Get("XvsY_rec");
    TH2D *h_XvsZ_rec = (TH2D*)dir->Get("XvsZ_rec");
    TH2D *h_YvsZ_rec = (TH2D*)dir->Get("YvsZ_rec");
    TH2D *h_XvsY_sim = (TH2D*)dir->Get("XvsY_sim");
    TH2D *h_XvsZ_sim = (TH2D*)dir->Get("XvsZ_sim");
    TH2D *h_YvsZ_sim = (TH2D*)dir->Get("YvsZ_sim");
    TH2D *h_XY_beam = (TH2D*)dir->Get("XY_beam");
    TH2D *h_XZ_beam = (TH2D*)dir->Get("XZ_beam");
    TH2D *h_YZ_beam = (TH2D*)dir->Get("YZ_beam");

    
    TH2D *tom_XY_rec = new TH2D("tom_XY_rec", "XY_rec/XY_beam", 100, -5, 5, 100, -5, 5);
    TH2D *tom_XZ_rec = new TH2D("tom_XZ_rec", "XZ_rec/XZ_beam", 400, -200, 200, 100, -5, 5);
    TH2D *tom_YZ_rec = new TH2D("tom_YZ_rec", "YZ_rec/YZ_beam", 400, -200, 200, 100, -5, 5);
    TH2D *tom_XY_sim = new TH2D("tom_XY_sim", "XY_sim/XY_beam", 100, -5, 5, 100, -5, 5);
    TH2D *tom_XZ_sim = new TH2D("tom_XZ_sim", "XZ_sim/XZ_beam", 400, -200, 200, 100, -5, 5);
    TH2D *tom_YZ_sim = new TH2D("tom_YZ_sim", "YZ_sim/YZ_beam", 400, -200, 200, 100, -5, 5);
  //  TH2D *tom_YZ_sim = new TH2D("tom_YZ_sim", "YZ_sim/YZ_beam", h_YvsZ_sim->GetNbinsX(), h_YvsZ_sim->GetXaxis()->GetXmin(), h_YvsZ_sim->GetXaxis()->GetXmax(), h_YvsZ_sim->GetNbinsY(), h_YvsZ_sim->GetYaxis()->GetXmin(), h_YvsZ_sim->GetYaxis()->GetXmax());
   
   
   tom_XY_rec->Divide(h_XvsY_rec,h_XY_beam,1.0,1.0,"B");
   tom_XZ_rec->Divide(h_XvsZ_rec,h_XZ_beam,1.0,1.0,"B");
   tom_YZ_rec->Divide(h_YvsZ_rec,h_YZ_beam,1.0,1.0,"B");
   tom_XY_sim->Divide(h_XvsY_sim,h_XY_beam,1.0,1.0,"B");
   tom_XZ_sim->Divide(h_XvsZ_sim,h_XZ_beam,1.0,1.0,"B");
   tom_YZ_sim->Divide(h_YvsZ_sim,h_YZ_beam,1.0,1.0,"B");

/*    TH2D *tom_XY_rec = (TH2D*)h_XvsY_rec->Clone("tom_XY_rec");
    TH2D *tom_XZ_rec = (TH2D*)h_XvsZ_rec->Clone("tom_XZ_rec");
    TH2D *tom_YZ_rec = (TH2D*)h_YvsZ_rec->Clone("tom_YZ_rec");
    TH2D *tom_XY_sim = (TH2D*)h_XvsY_sim->Clone("tom_XY_sim");
    TH2D *tom_XZ_sim = (TH2D*)h_XvsZ_sim->Clone("tom_XZ_sim");
    TH2D *tom_YZ_sim = (TH2D*)h_YvsZ_sim->Clone("tom_YZ_sim");

    tom_XY_rec->Divide(h_XY_beam);
    tom_XZ_rec->Divide(h_XZ_beam);
    tom_YZ_rec->Divide(h_YZ_beam);
    tom_XY_sim->Divide(h_XY_beam);
    tom_XZ_sim->Divide(h_XZ_beam);
    tom_YZ_sim->Divide(h_YZ_beam);
*/

TFile *outputTom = new TFile("/eos/home-g/gmeinard/2024Target/outputTom.root", "RECREATE");

    // Write histograms to the output file
    h_XvsY_rec ->Write();
    h_XvsZ_rec ->Write();
    h_YvsZ_rec ->Write();
    h_XvsY_sim ->Write();
    h_XvsZ_sim ->Write();
    h_YvsZ_sim ->Write();
    h_XY_beam ->Write();
    h_XZ_beam ->Write();
    h_YZ_beam ->Write();

    
    tom_XY_rec->Write();
    tom_XZ_rec->Write();
    tom_YZ_rec->Write();
    tom_XY_sim->Write();
    tom_XZ_sim->Write();
    tom_YZ_sim->Write();

    outputTom->Close();
    hfileT->Close();

    std::cout << "Histograms successfully saved to outputTom.root" << std::endl;


}