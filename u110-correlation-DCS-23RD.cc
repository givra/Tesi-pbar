#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <stdio.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaDetect.h"
#include "PaRich.h"
#include "PaRichDet.h"
#include "PaEvent.h"
#include "PaParticle.h"
#include "PaPid.h"
#include "TMath.h"
#include "PaMCvertex.h"
#include "PaMCtrack.h"

// RD 2023 
// /eos/experiment/amber/pbar/common/data23/RD/productions/W01t0.8rich1.6/mDST/mDST-304***.root.***

// ______________________________________________________ //
//                                                        //
//     This UE is used to study possible correlations     //
//     between 2023 real data and pressure/temperature    //
//     of the He target stored in several TTrees          //
// _______________________________________________________//

static TH1D* h1[20];
static TH2D* Correlation[20];

void UserEvent110(PaEvent& e){

    static TTree* OutTree(NULL);        // output tree
    TFile hfileFM1("He3_flow_rate_(FM1)_2023-05-19_00_00_00_2023-06-26_00_00_00.root");
    TFile hfileTPO4("He4_evaporator_bottom_(TPO4_temperature)_2023-05-19_00_00_00_2023-06-26_00_00_00.root");
    TFile hfileFM60("He4_flow_rate_(FM60)_2023-05-19_00_00_00_2023-06-26_00_00_00.root");
    TFile hfileTTH4("Mxc_upstream_end_(TTH4_temperature)_2023-05-19_00_00_00_2023-06-26_00_00_00.root");
    TFile hfilePV8("Still_pressure_(PV8)_2023-05-19_00_00_00_2023-06-26_00_00_00.root");
    TFile hfilePV8_piezo("Still_pressure_piezo_(PV8_piezo)_2023-05-19_00_00_00_2023-06-26_00_00_00.root");

    // variables to read from the TTrees

    static double FM1;
    static double TPO4;
    static double FM60;
    static double TTH4;
    static double PV8;
    static double PV8piezo;

    TTimeStamp* time_FM1 = nullptr;
    TTimeStamp* time_TPO4 = nullptr;
    TTimeStamp* time_FM60 = nullptr;
    TTimeStamp* time_TTH4 = nullptr;
    TTimeStamp* time_PV8 = nullptr;
    TTimeStamp* time_PV8piezo = nullptr;

    // time of the event
    static double time_event;

    static int NTracks = 0;
    int NVrtx = e.NVertex();

    int SpillNum = e.SpillNum();
    static double radL;
	static double chi2;
	static int ndf;
    static int run;
	static double TiS; 			// time in spill
    int Nevent = e.UniqueEvNum();
    int NVrtx = e.NVertex();
    int Bvtx = -1;


    static bool first(true);
    if(first){ // histograms and Ntupes booking block
    Phast::Ref().HistFileDir("UserEvent110");

    h1[0] = new TH1D("Nvrtx","number of vertices", 100, 0 , 100);
    h1[1] = new TH1D("Ntrk","number of tracks", 100, 0, 50);
    h1[2] = new TH1D("FM1", "He3 flow rate",  50, 9, 11);
    h1[3] = new TH1D("TPO4", "TPO4 temperature",  10, 0, 0.05);
    h1[4] = new TH1D("FM60", "He4 flow rate", 100, 0, 100);
    h1[5] = new TH1D("TTH4", "TTH4 temperature", 50, 0, 1.6);
    h1[6] = new TH1D("PV8","Still pressure",10 ,0 ,10);
    h1[7] = new TH1D("PV8_piezo"," Still pressure piezo",50, 0, 6);
    // Correlation[0] = new TH2D("FM1-tracks", "correlation He3 flow rate and number of vertices; NVertex; He3 flow", 100, 0, 100, 50, 9, 11);
    // Correlation[0] -> SetOption("colz");

    }

    TTree *treeFM1 = (TTree*)hfileFM1.Get("He3 flow rate (FM1)");
    treeFM1->SetBranchAddress("timestamp", &time_FM1);
    treeFM1->SetBranchAddress("fValue",&FM1);

    TTree *treeTPO4 = (TTree*)hfileTPO4.Get("He4 evaporator bottom (TPO4 temperature)");
    treeTPO4->SetBranchAddress("timestamp", &time_TPO4);
    treeTPO4->SetBranchAddress("fValue",&TPO4);

    TTree *treeFM60 = (TTree*)hfileFM60.Get("He4 flow rate (FM60)");
    treeFM60->SetBranchAddress("timestamp", &time_FM60);
    treeFM60->SetBranchAddress("fValue",&FM60);

    TTree *treeTTH4 = (TTree*)hfileTTH4.Get("Mxc upstream end (TTH4 temperature)");
    treeTTH4->SetBranchAddress("timestamp", &time_TTH4);
    treeTTH4->SetBranchAddress("fValue",&TTH4);

    TTree *treePV8 = (TTree*)hfilePV8.Get("Still pressure (PV8)");
    treePV8->SetBranchAddress("timestamp", &time_PV8);
    treePV8->SetBranchAddress("fValue",&PV8);

    TTree *treePV8_piezo = (TTree*)hfilePV8_piezo.Get("Still pressure piezo (PV8 piezo)");
    treePV8_piezo->SetBranchAddress("timestamp", &time_PV8piezo);
    treePV8_piezo->SetBranchAddress("fValue",&PV8piezo);


}