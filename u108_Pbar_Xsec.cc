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
//#include "CEDAR_Helper.h"
#include "TMath.h"
#include "PaMCvertex.h"
#include "PaMCtrack.h"

// rec data on LH2 @250GeV
// /eos/experiment/amber/pbar/common/data24/MC/rec/productions/W01t0.1_beamfile_SoftQCDinelastic_batch1/mDST-304002.root.001

void UserEvent108(PaEvent& e){
			        
    static TH1D* h1[20];
	static TH2D* h2[10];
	static TTree* tree(NULL);
    static int Nout = 0;
    static double mom;
	static int q;
	static double Zvertex_rec;
	static double Yvertex_rec;
	static double Xvertex_rec;
	static double Zvertex_sim;
	static double Yvertex_sim;
	static double Xvertex_sim;
	static double X2;
	static int run;
	static double TiS; 			// time in spill
	static double Mtime;
	static double zfirst;
    static double zlast;
	static double ztargetB = -70;		// z coord before target
	static double ztargetA = 70;		// z coord after target
	static double XtrajB; 
	static double YtrajB;
	static double XtrajA; 
	static double YtrajA;
	static double Xrid;
	static double theta;
	static double radL;
	static double chi2;
	static int ndf;

//    static CEDAR* cedar = &CEDAR::Instance();

	int NVrtx = e.NVertex();
	int NMCvrtx = e.NMCvertex();
	int Nevent = e.UniqueEvNum();
	int N_tracks;			// number of incoming tracks
//	int N_hits_FI;
//	int hits_FI01;
//	int hits_FI02;
//	int hits_FI15;
//	int N_hits_SI;
	bool RICH; 				// true if RICH info is != 0

	int SpillNum = e.SpillNum();
	int EventSpill = e.EvInSpill();
	//cout << " Evento " << Nevent << " numero spill " << SpillNum << " evento in spill " << EventSpill << endl;
	

    static bool first(true);
    if(first){ // histograms and Ntupes booking block
    Phast::Ref().HistFileDir("UserEvent108");

    h1[1]  = new TH1D("momentum_rec","reconstructed momentum distribution; GeV/c", 200, 0, 200);
	h1[2]  = new TH1D("z_vertex_rec","reconstructed Z primary vertex; cm", 400, -200, 200);
	h1[3]  = new TH1D("y_vertex_rec","reconstructed Y primary vertex; cm", 50, -4, 4);
	h1[4]  = new TH1D("x_vertex_rec","reconstructed X primary vertex; cm", 50, -4, 4);
//	h1[5]  = new TH1D("N_tracks", "Number of beam tracks", 15, 0, 15);
//	h1[6]  = new TH1D("TiS", "Event Time in Spill; Time [s]; Events", 100, 0, 6);
//	h1[7]  = new TH1D("Mtime", "Mean Time; [ns]; Events", 50, -8, 8);
	h1[9]  = new TH1D("z_vertex_sim","simulated Z primary vertex; cm", 400, -200, 200);
	h1[10]  = new TH1D("y_vertex_sim","simulated Y primary vertex; cm", 50, -4, 4);
	h1[11]  = new TH1D("x_vertex_sim","simulated X primary vertex; cm", 50, -4, 4);

	h2[1]  = new TH2D("XvsY_rec","reconstructed XY primary vertex; X [cm]; Y [cm]", 100, -5, 5, 100, -5, 5);
	h2[1] -> SetOption("colz");
	h2[2]  = new TH2D("XvsZ_rec","reconstructed XZ primary vertex; Z [cm]; X [cm]", 400, -200, 200, 80, -5, 5);
	h2[2] -> SetOption("colz");
	h2[3]  = new TH2D("YvsZ_rec","reconstructed YZ primary vertex; Z [cm]; Y [cm]", 400, -200, 200, 80, -5, 5);
	h2[3] -> SetOption("colz");
//	h2[4]  = new TH2D("SciFivsSi","Beam track hits; Silicon Hits; SciFi Hits", 12, 0, 12, 8, 0, 8);
//	h2[4] -> SetOption("colz");
	h2[5]  = new TH2D("XvsY_sim","simulated XY primary vertex; X [cm]; Y [cm]", 100, -5, 5, 100, -5, 5);
	h2[5] -> SetOption("colz");
	h2[6]  = new TH2D("XvsZ_sim","simulated XZ primary vertex; Z [cm]; X [cm]", 400, -200, 200, 80, -5, 5);
	h2[6] -> SetOption("colz");
	h2[7]  = new TH2D("YvsZ_sim","simulated YZ primary vertex; Z [cm]; Y [cm]", 400, -200, 200, 80, -5, 5);
	h2[7] -> SetOption("colz");

    tree = new TTree("USR108","User Ntuple example");
	
	tree->Branch("momentum_rec", &mom, "mom/D");
	tree->Branch("Zvertex_rec", &Zvertex_rec, "Zvertex_rec/D");
	tree->Branch("Yvertex_rec", &Yvertex_rec, "Yvertex_rec/D");
	tree->Branch("Xvertex_rec", &Xvertex_rec, "Xvertex_rec/D");
	tree->Branch("Zvertex_sim", &Zvertex_sim, "Zvertex_sim/D");
	tree->Branch("Yvertex_sim", &Yvertex_sim, "Yvertex_sim/D");
	tree->Branch("Xvertex_sim", &Xvertex_sim, "Xvertex_sim/D");
//	tree->Branch("N_tracks", &N_tracks, "N_tracks/I");
//	tree->Branch("TiS", &TiS, "TiS/D");
//	tree->Branch("Mtime", &Mtime, "Mtime/D");
	tree->Branch("zfirst", &zfirst, "zfirst/D");
   	tree->Branch("zlast", &zlast, "zlast/D");


    first=false;
    }

    //pbar 				= 0;
	int Bvtx 		    = -1;
	bool extrap 		= false;
	bool extrapB		= false;		// Before target
	bool extrapA		= false;		// After target
	double Xmin 		= 100;
	int ivok 			= e.iBestPrimaryVertex();
	int Vcounter = 0;
	int MCVcounter = 0;

	for(int jv = 0; jv < NVrtx ; jv++){ 

		Bvtx  =  e.iBestPrimaryVertex();
		TiS = e.TimeInSpill();
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
	//cout << "indice vertice scelto " << ivok << endl; 

	Nout = v.NOutParticles();
	N_tracks = v.InParticle();


	// loop over MC vertices
	for(int kv = 0; kv < NMCvrtx; kv++){

		const PaMCvertex& MCv = e.vMCvertex(kv);
		const PaMCtrack& MCtrack = e.vMCtrack(kv);
					
		if(! MCv.IsPrimary()) continue;					// to select only primary MC vertex
		MCVcounter ++;

		Zvertex_sim = MCv.Pos(2);
		Yvertex_sim = MCv.Pos(1);
		Xvertex_sim = MCv.Pos(0);

	}

	// cout << " evento " << Nevent << " vertici simulati " << NMCvrtx << " PV simulati " << MCVcounter << " vertici ricostruiti " << NVrtx << " PV ricostruiti " << Vcounter << endl;

	for(int in_p = 0; in_p < N_tracks; in_p++){				// loop over incoming particles in target
	 	
	 	int inc_index = v.InParticle();
		const PaParticle& beam_p = e.vParticle(inc_index);			// creates a beam particle for every incoming track

		if(beam_p.iTrack() == -1) {
			cout << " Evento " << Nevent<< " incoming particle's index -1" << endl;
			continue;
		}
		const PaTrack& beam_tr = e.vTrack(beam_p.iTrack());

//		hits_FI01 = beam_tr.NHitsFoundInDetect("FI01"); 			// should return n of hits in SciFi FI01,15,02
//		hits_FI02 = beam_tr.NHitsFoundInDetect("FI02");
//		hits_FI15 = beam_tr.NHitsFoundInDetect("FI15");
//		N_hits_FI = hits_FI01 + hits_FI02 + hits_FI15;
//		N_hits_SI = beam_tr.NHitsFoundInDetect("SI");
		// cout << " Evento " << Nevent << " n hits in SciFi 01 "  << hits_FI01 << " n hits in SciFi 02 " << hits_FI02 << " n hits in SciFi 15 " << hits_FI15 << " N hit tot " << N_hits_FI << endl;

	}

	run = e.RunNum();
//  n_rich = 1 + PaMetaDB::Ref().NminusOne(run)/1e6; //refraction index
//  if(PaMetaDB::Ref().NminusOne(run) < 0) n_rich = 1 + PaSetup::Ref().NminusOne()/1e6;
	
	
	Zvertex_rec = v.Z();
	Yvertex_rec = v.Y();
	Xvertex_rec = v.X();

	// incoming particle identification
	// cedar->Search(e);
	
/*	int CE1m =0, CE2m = 0;
        std::vector<PaDigit> rawDigits = e.RawDigits();
      unsigned int hitMask[2] = {0};
      for(int idig = 0; idig < (signed) rawDigits.size(); ++idig) { 

         bool is_ce1 = rawDigits[idig].DecodeMapName() == "CE01P1__";
         bool is_ce2 = rawDigits[idig].DecodeMapName() == "CE02P1__";
         int cedar = is_ce1 ? 1 : (is_ce2 ? 2 : 0);

         if(cedar){

            double hit_time = rawDigits[idig].DigInfo(3);
            int ipm = std::fabs(rawDigits[idig].IWire())-2;

            // calibration by "hand" for W01, looking at time plots
            // if( is_ce1 && hit_time < -3300 && hit_time > -3330 ) hitMask[0] |= 1 << ipm;
            // if( is_ce2 && hit_time < -3350 && hit_time > -3390 ) hitMask[1] |= 1 << ipm;

            if( is_ce1 && hit_time < -3310 && hit_time > -3330 ) hitMask[0] |= 1 << ipm;
                           if( is_ce2 && hit_time < -3355 && hit_time > -3375 ) hitMask[1] |= 1 << ipm;


         }
      }

           for(int ipmt = 0; ipmt < 8; ipmt++){
                        for( int ipad = 0; ipad < 4; ipad++){
                                if( (hitMask[0] >> (ipad + ipmt*4) & 0xF) > 0 ){
                                        CE1m += 1;
                                        break;
                                }
                        }
                        for( int ipad = 0; ipad < 4; ipad++){
                                if( (hitMask[1] >> (ipad + ipmt*4) & 0xF) > 0 ){
                                        CE2m += 1;
                                        break;
                                }
                        }
                }
	*/
		
		for(int ip = 0; ip < Nout; ip++){ 			 // loop over outgoing particles from each vertex
		
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

		   RICH = pt.hasRich();		   


		   mom = p.ParInVtx(ivok).Mom();
		   q = p.Q();
		   Mtime = pt.MeanTime();
		   zfirst = pt.ZFirst();
		   zlast = pt.ZLast();
		   chi2 = pt.Chi2tot();
		   ndf = pt.Ndf();
		   radL = pt.XX0();
		   Xrid = chi2/ndf;
		   theta = par.Theta();
		   XtrajB = parB.X();
		   YtrajB = parB.Y();
		   XtrajA = parA.X();
		   YtrajA = parA.Y();

		   double rB = sqrt((XtrajB*XtrajB) + (YtrajB*YtrajB));
		   double rA = sqrt((XtrajA*XtrajA) + (YtrajA*YtrajA));

		   
		   
	  
		/*  if(Nout < 2 && theta < 0.008 && mom > 188) continue;	// quasi-elastic event

		  if(rB > 1.75 || rA > 1.75) continue;					// beam track extrapolated through the target
																
		  if(RICH != 0) continue;								// RICH info at least 1 hit

		  if(CE1m < 6 || CE2m < 6) continue;					// CEDAR multiplicity

		  if(zfirst < -350 && zlast > 350) continue; 			//	cut for good momentum reconstruction	
		  
		  if(Mtime < -3 || Mtime > 3) continue;					// mean time cut

		  if(Zvertex < -65 || Zvertex > 65) continue;  			// PV outside target region

		  if(Xrid > 10) continue;								// X/ndf too large
			
		  if(abs(Yvertex) > 1.75) continue;

	      if(abs(Xvertex) > 1.75) continue;

		  if(radL > 10) continue; 								// no muons

		  if(mom < 10 || mom > 45) continue;
		  
		  if((theta < 0.01 || theta > 0.18) && (mom < 10 || mom > 60)) continue;			// phase space cut between 10 and 60 GeV

		  if(Mtime > 1e4) continue;								// defined mean time
		  
		// if(tm & 0 << 0 ) continue;							// tm in posizione 1 diverso da 0

		  if(mom > 190) continue;								// anelastic event

//		  if(N_hits_FI < 4 && N_hits_SI < 10) continue;			// >2 beam tracks

//		  if(TiS < 1.2 || TiS > 5.4) continue; 					// time in spill cut
		*/

	// ******************** filling histos **********************
		  h1[1] -> Fill(mom);
		  h1[2] -> Fill(Zvertex_rec);
		  h1[3] -> Fill(Yvertex_rec);
		  h1[4] -> Fill(Xvertex_rec);
		  //h1[5] -> Fill(N_tracks + 1);
		  //h1[6] -> Fill(TiS);
		  //h1[7] -> Fill(Mtime);

		  h2[1] -> Fill(Xvertex_rec, Yvertex_rec);
		  h2[2] -> Fill(Zvertex_rec, Xvertex_rec);
		  h2[3] -> Fill(Zvertex_rec, Yvertex_rec); 
		//  h2[4] -> Fill(N_hits_SI, N_hits_FI); 

		h1[9] -> Fill(Zvertex_sim);
		h1[10] -> Fill(Yvertex_sim);
		h1[11] -> Fill(Xvertex_sim);

		h2[5] -> Fill(Xvertex_sim, Yvertex_sim);
		h2[6] -> Fill(Zvertex_sim, Xvertex_sim);
		h2[7] -> Fill(Zvertex_sim, Yvertex_sim);
		 
		  
/*		  proton = pid.IsProton(mom, mom_p, mom_k, mom_pi, like, count, q);


		  switch(proton)
		  { case 0:    //pion-
				h1[5] -> Fill(proton);
				numPi++;
				h2[1] -> Fill(mom, theta);
				break;
		    case 1:    //kaon-
				h1[5] -> Fill(proton);
				countproton1++;
				numK++;
				h2[2] -> Fill(mom, theta);
			break;
		    case 2:    //pbar
				h1[5] -> Fill(proton);
				pbar++;
				countAproton++;
				numPbar++;
				h2[3] -> Fill(mom, theta);
				break;
		     case 3:     //noID
				h1[5] -> Fill(proton);
				break;
		     case 4:      //pion+
				h1[5] -> Fill(proton);
				break;
		     case 5:	  //kaon+
				h1[5] -> Fill(proton);
				break;
		     case 6:	  //proton
				h1[5] -> Fill(proton);
				break;
		    //default: 
			//cout << " caso default # event " << (Nevent - 1048580)/20 << " proton " << proton << endl;	
			//break;
            }
	*/		
		    tree -> Fill();

	//	    if(proton >= 3) continue; 

        }   //end particles loop
		

}