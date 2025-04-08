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
#include "TLegend.h"
#include <TCanvas.h>

// rec data on LH2 @250GeV
// /eos/experiment/amber/pbar/common/data24/MC/rec/productions/W01t0.1_beamfile_SoftQCDinelastic_batch1/mDST-304002.root.001
// LH2 @250GeV w/ new beamfile
// /eos/experiment/amber/pbar/common/data24/MC/rec/productions/W01t0.1_beamfile_SoftQCDinelastic_batch4/mDST-304117.root.001

static TH1D* h1[20];
static TH2D* h2[20];
// static TLegend* legend;
// static TCanvas* C;

void UserEvent108(PaEvent& e){
			        
	static TTree* tree(NULL);
    // static int Nout = 0;
    // static double mom;
	// static int q;
	static double Zvertex_rec;
	static double Yvertex_rec;
	static double Xvertex_rec;
	static double Zvertex_sim;
	static double Yvertex_sim;
	static double Xvertex_sim;
	static double X2;
	// static int run;
	// static double ztargetB = -70;		// z coord before target
	// static double ztargetA = 70;		// z coord after target
	// static double XtrajB; 
	// static double YtrajB;
	// static double XtrajA; 
	// static double YtrajA;
	// static double Xrid;
	// static double chi2;
	// static int ndf;
	static double Xbeam;
	static double Ybeam;
	static double Zbeam;


	int NVrtx = e.NVertex();
	int NMCvrtx = e.NMCvertex();
	int Nevent = e.UniqueEvNum();
	int N_tracks;			// number of incoming tracks

	// int SpillNum = e.SpillNum();
	// int EventSpill = e.EvInSpill();
	//cout << " Evento " << Nevent << " numero spill " << SpillNum << " evento in spill " << EventSpill << endl;
	

    static bool first(true);
    if(first){ // histograms and Ntupes booking block
    Phast::Ref().HistFileDir("UserEvent108");


    h1[1]  = new TH1D("momentum_rec","reconstructed momentum distribution; GeV/c", 200, 0, 200);
	h1[2]  = new TH1D("z_vertex_rec","reconstructed Z primary vertex; cm", 1000, -200, 200);
	h1[2]->SetLineColor(kBlue);
	h1[3]  = new TH1D("y_vertex_rec","reconstructed Y primary vertex; cm", 500, -5, 5);	
	h1[4]  = new TH1D("x_vertex_rec","reconstructed X primary vertex; cm", 500, -5, 5);
	h1[5]  = new TH1D("N_tracks", "Number of beam tracks", 15, 0, 15);
	h1[6]  = new TH1D("X_beam", "X of beam tracks", 500, -5, 5);
	h1[7]  = new TH1D("Y_beam", "Y of beam tracks", 500, -5, 5);
	h1[8]  = new TH1D("Z_beam", "Z of beam tracks", 500, -5, 5);

	h1[9]  = new TH1D("z_vertex_sim","simulated Z primary vertex; cm", 1000, -200, 200);
	h1[9]->SetLineColor(kRed);
	h1[10]  = new TH1D("y_vertex_sim","simulated Y primary vertex; cm", 500, -5, 5);	
	h1[11]  = new TH1D("x_vertex_sim","simulated X primary vertex; cm", 500, -5, 5);	

	h2[1]  = new TH2D("XvsY_rec","reconstructed XY primary vertex; X [cm]; Y [cm]", 500, -5, 5, 500, -5, 5);	
	h2[1] -> SetOption("colz");
	h2[2]  = new TH2D("XvsZ_rec","reconstructed XZ primary vertex; Z [cm]; X [cm]", 1000, -200, 200, 500, -5, 5);	
	h2[2] -> SetOption("colz");
	h2[3]  = new TH2D("YvsZ_rec","reconstructed YZ primary vertex; Z [cm]; Y [cm]", 1000, -200, 200, 500, -5, 5);	
	h2[3] -> SetOption("colz");
//	h2[4]  = new TH2D("SciFivsSi","Beam track hits; Silicon Hits; SciFi Hits", 12, 0, 12, 8, 0, 8);
//	h2[4] -> SetOption("colz");
	h2[5]  = new TH2D("XvsY_sim","simulated XY primary vertex; X [cm]; Y [cm]", 500, -5, 5, 500, -5, 5);
	h2[5] -> SetOption("colz");
	h2[6]  = new TH2D("XvsZ_sim","simulated XZ primary vertex; Z [cm]; X [cm]", 1000, -200, 200, 500, -5, 5);
	h2[6] -> SetOption("colz");
	h2[7]  = new TH2D("YvsZ_sim","simulated YZ primary vertex; Z [cm]; Y [cm]", 1000, -200, 200, 500, -5, 5);
	h2[7] -> SetOption("colz");

	h2[8] = new TH2D("XY_beam", "XY beam tracks; X [cm]; Y [cm]", 500, -5, 5, 500, -5, 5);
	h2[8] -> SetOption("colz");
	h2[9] = new TH2D("XZ_beam", "XZ beam tracks; Z [cm]; X [cm]", 1000, -200, 200, 500, -5, 5);
	h2[9] -> SetOption("colz");
	h2[10] = new TH2D("YZ_beam", "YZ beam tracks; Z [cm]; Y [cm]", 1000, -200, 200, 500, -5, 5);
	h2[10] -> SetOption("colz");

	h2[11] = new TH2D("XY_rec/beam", "tomography at z = 0; X [cm]; Y [cm]", 500, -5, 5, 500, -5, 5);
	h2[11] -> SetOption("colz");
	h2[12] = new TH2D("XZ_rec/beam", "tomography at y = 0; Z [cm]; X [cm]", 1000, -200, 200, 500, -5, 5);
	h2[12] -> SetOption("colz");
	h2[13] = new TH2D("YZ_rec/beam", "tomography at x = 0; Z [cm]; Y [cm]", 1000, -200, 200, 500, -5, 5);
	h2[13] -> SetOption("colz");

	h2[14] = new TH2D("XY_sim/beam", "tomography at z = 0; X [cm]; Y [cm]", 500, -5, 5, 500, -5, 5);
	h2[14] -> SetOption("colz");
	h2[15] = new TH2D("XZ_sim/beam", "tomography at y = 0; Z [cm]; X [cm]", 1000, -200, 200, 500, -5, 5);
	h2[15] -> SetOption("colz");
	h2[16] = new TH2D("YZ_sim/beam", "tomography at x = 0; Z [cm]; Y [cm]", 1000, -200, 200, 500, -5, 5);
	h2[16] -> SetOption("colz");


    tree = new TTree("USR108","User Ntuple example");
	
	//tree->Branch("momentum_rec", &mom, "mom/D");
	tree->Branch("Zvertex_rec", &Zvertex_rec, "Zvertex_rec/D");
	tree->Branch("Yvertex_rec", &Yvertex_rec, "Yvertex_rec/D");
	tree->Branch("Xvertex_rec", &Xvertex_rec, "Xvertex_rec/D");
	tree->Branch("Zvertex_sim", &Zvertex_sim, "Zvertex_sim/D");
	tree->Branch("Yvertex_sim", &Yvertex_sim, "Yvertex_sim/D");
	tree->Branch("Xvertex_sim", &Xvertex_sim, "Xvertex_sim/D");
	tree->Branch("N_tracks", &N_tracks, "N_tracks/I");
	tree->Branch("Ybeam", &Ybeam, "Ybeam/D");
	tree->Branch("Xbeam", &Xbeam, "Xbeam/D");


    first=false;
    }

	int Bvtx 		    = -1;
	bool extrap 		= false;
	// bool extrapB		= false;		// Before target
	// bool extrapA		= false;		// After target
	bool extrap0     	= false;
	double Xmin 		= 100;
	int ivok 			= e.iBestPrimaryVertex();
	int Vcounter = 0;
	int MCVcounter = 0;

	for(int jv = 0; jv < NVrtx ; jv++){ 

		Bvtx  =  e.iBestPrimaryVertex();
		const PaVertex& v1 = e.vVertex(jv);
		X2 = v1.Chi2();
	
		if(! v1.IsPrimary()) continue;
	
		if(Bvtx == -1){
		
			if(X2 < Xmin){
			 	Xmin = X2;
				ivok = jv;  //si seleziona l'indice corrispondente al X2 minore	
			}
		}
	}
	if(ivok == -1) return;
	const PaVertex& v = e.vVertex(ivok);
	Vcounter++;

	// Nout = v.NOutParticles();
	N_tracks = v.InParticle();

	// h1[5] -> Fill(N_tracks + 1);

	// loop over MC vertices
	for(int kv = 0; kv < NMCvrtx; kv++){

		const PaMCvertex& MCv = e.vMCvertex(kv);
		// const PaMCtrack& MCtrack = e.vMCtrack(kv);
					
		if(! MCv.IsPrimary()) continue;					// to select only primary MC vertex
		MCVcounter ++;

		Zvertex_sim = MCv.Pos(2);
		Yvertex_sim = MCv.Pos(1);
		Xvertex_sim = MCv.Pos(0);
	}


	for( int it = 0; it < e.NTrack(); it++){			// loop over tracks in target
	 	
	 	const PaTrack& beam_track = e.vTrack(it);

      	if( beam_track.NTPar() == 0 ) continue;  				//Skip the track if it has no parameters
		if( beam_track.iParticle() == -1 ) continue;
		if( !beam_track.IsBeam() ) continue;

		PaTPar par0;			// parameters at z = 0
		extrap0 = beam_track.Extrapolate(0., par0);
		
		Xbeam = par0.X();
		Ybeam = par0.Y();
		Zbeam = par0.Z();

		h1[6] -> Fill(Xbeam);
		h1[7] -> Fill(Ybeam);
		h1[8] -> Fill(Zbeam);
		h2[8] -> Fill(Xbeam,Ybeam);

			for(int j=1; j<501; j++){
				double contentX = h1[6]->GetBinContent(j);
				double contentY = h1[7]->GetBinContent(j);

				for(int k=1; k<1001;k++){
				h2[9]->SetBinContent(k,j,contentX);
				h2[10]->SetBinContent(k,j,contentY);		
				}
			}

//		hits_FI01 = beam_tr.NHitsFoundInDetect("FI01"); 			// should return n of hits in SciFi FI01,15,02
//		hits_FI02 = beam_tr.NHitsFoundInDetect("FI02");
//		hits_FI15 = beam_tr.NHitsFoundInDetect("FI15");
//		N_hits_FI = hits_FI01 + hits_FI02 + hits_FI15;
//		N_hits_SI = beam_tr.NHitsFoundInDetect("SI");
		// cout << " Evento " << Nevent << " n hits in SciFi 01 "  << hits_FI01 << " n hits in SciFi 02 " << hits_FI02 << " n hits in SciFi 15 " << hits_FI15 << " N hit tot " << N_hits_FI << endl;

	}

	// run = e.RunNum();	
	
	Zvertex_rec = v.Z();
	Yvertex_rec = v.Y();
	Xvertex_rec = v.X();

	

	/*	for(int ip = 0; ip < Nout; ip++){ 			 // loop over outgoing particles from each vertex
		
		   int index = v.iOutParticle(ip);
		   const PaParticle& p = e.vParticle(index); 
		   if(p.iTrack() == -1) continue;

		   const PaTrack& pt = e.vTrack(p.iTrack());
		   const PaTPar& par = pt.vTPar(0);

		   const PaSetup& setup = PaSetup::Ref();
		   const double z_rich = setup.Rich().DetPos(0).Z();
		   PaTPar Hout;
		   extrap = pt.Extrapolate(z_rich, Hout);  //extrapolates trajectory parameters at z = z_rich and Hout is the result
		   PaTPar parB;			// parameter extrapolated before target
		   PaTPar parA;			// parameter extrapolated after target

		   extrapB = pt.Extrapolate(ztargetB, parB);
		   extrapA = pt.Extrapolate(ztargetA, parA);

		   mom = p.ParInVtx(ivok).Mom();
		   q = p.Q();
		   chi2 = pt.Chi2tot();
		   ndf = pt.Ndf();
		   Xrid = chi2/ndf;
		   XtrajB = parB.X();
		   YtrajB = parB.Y();
		   XtrajA = parA.X();
		   YtrajA = parA.Y();

		   double rB = sqrt((XtrajB*XtrajB) + (YtrajB*YtrajB));
		   double rA = sqrt((XtrajA*XtrajA) + (YtrajA*YtrajA));

		   h1[1] -> Fill(mom);
		 
		//  h2[4] -> Fill(N_hits_SI, N_hits_FI);


        }   //end particles loop
	*/
		tree -> Fill();

		if(Zvertex_rec > 100. || Zvertex_rec < -100.) return;
		if(Zvertex_sim > 100. || Zvertex_sim < -100.) return;

		h1[2] -> Fill(Zvertex_rec);
		h1[3] -> Fill(Yvertex_rec);
		h1[4] -> Fill(Xvertex_rec);		  

		h2[1] -> Fill(Xvertex_rec, Yvertex_rec);
		h2[2] -> Fill(Zvertex_rec, Xvertex_rec);
		h2[3] -> Fill(Zvertex_rec, Yvertex_rec);
	  

	  h1[9] -> Fill(Zvertex_sim);
	  h1[10] -> Fill(Yvertex_sim);
	  h1[11] -> Fill(Xvertex_sim);

	  h2[5] -> Fill(Xvertex_sim, Yvertex_sim);
	  h2[6] -> Fill(Zvertex_sim, Xvertex_sim);
	  h2[7] -> Fill(Zvertex_sim, Yvertex_sim);
		
}

void UserJobEnd108(){

//	C = new TCanvas("C","Z vertex",200,10,600,400);
//	legend = new TLegend(0.15, 0.7, 0.35, 0.9);
//
//	h1[2]->Draw();
//	h1[9]->Draw("SAME");
//	
//	legend->AddEntry(h1[2], "Reconstructed", "l");
//	legend->AddEntry(h1[9], "Simulated", "l");
//	legend->Draw();
//	
//	C->Update();
//	C->SaveAs("/eos/home-g/gmeinard/2024Target/Zvertex.png");
//	C->Close();

	h2[11]->Divide(h2[1],h2[8],1.0,1.0,"B");		// XY
	h2[12]->Divide(h2[2],h2[9],1.0,1.0,"B");		// XZ
	h2[13]->Divide(h2[3],h2[10],1.0,1.0,"B");		// YZ

	h2[14]->Divide(h2[5],h2[8],1.0,1.0,"B");		// XY
	h2[15]->Divide(h2[6],h2[9],1.0,1.0,"B");		// XZ
	h2[16]->Divide(h2[7],h2[10],1.0,1.0,"B");		// YZ
}