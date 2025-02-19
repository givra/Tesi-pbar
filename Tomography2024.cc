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

// _____________________________________________________________
//
//          This macro reads the tree USR108 from u108
//          to create the overlapping Z_sim/Z_rec distribution
//          and the tomography of the 2024 target
//______________________________________________________________

using namespace std;

void Tomography(){

    TFile hfileT("USR108.root");                // TO BE CHECKED

    static TH1D *h1[10];
    static TH2D *h2[10];

    // reading the TTree and its branches

    TTree *tree = (TTree*)hfileT.Get("USR108");
    TBranch *b1T=tree->GetBranch("Xvertex_rec");
    TBranch *b2T=tree->GetBranch("Yvertex_rec");
    TBranch *b3T=tree->GetBranch("Zvertex_rec");
    TBranch *b4T=tree->GetBranch("Xvertex_sim");
    TBranch *b5T=tree->GetBranch("Yvertex_sim");
    TBranch *b6T=tree->GetBranch("Zvertex_sim");
    TBranch *b7T=tree->GetBranch("Xbeam");
    TBranch *b8T=tree->GetBranch("Ybeam");

    b1T->SetAddress(&Xvertex_rec);
    b2T->SetAddress(&Yvertex_rec);
    b3T->SetAddress(&Zvertex_rec);
    b4T->SetAddress(&Xvertex_sim);
    b5T->SetAddress(&Yvertex_sim);
    b6T->SetAddress(&Zvertex_sim);
    b7T->SetAddress(&Xbeam);
    b8T->SetAddress(&Ybeam);


    Nevent = tree2->GetEntries();


}