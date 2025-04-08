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
#include "TTimeStamp.h"
#include "TFile.h"

// RD 2023 
// /eos/experiment/amber/pbar/common/data23/RD/productions/W01t0.8rich1.6/mDST/mDST-300537.root.***

// ______________________________________________________ //
//                                                        //
//     This UE is used to study possible correlations     //
//     between 2023 real data and pressure/temperature    //
//     of the He target stored in several TTrees          //
// _______________________________________________________//


// static int NTracks = 0;
// static double R;            // rate = NVrtx/Nbeam

void UserEvent110(PaEvent& e){

    // static TH1D* h1[20];
    static bool first(true);
    static TTree * OutTree = NULL;        // output tree
    static int Nout = 0;
    int NVrtx = 0;
    int Nbeam = 0;        // # beam particles
//	static double chi2;
//	static int ndf;
	static double TiS; 			// time in spill
    static int Nevent;
    static int TgM;
    static double X2;
    static double UnixSeconds;
    static unsigned int CE1m;
    static unsigned int CE2m;
	static unsigned int hitMask[2];
    static double ztargetB = -70.;		// z coord before target
	static double ztargetA = 70.;		// z coord after target
    static double rtarget = 2.;         // LH2 volume radius
//  static double XtrajB; 
//	static double YtrajB;
//	static double XtrajA; 
//	static double YtrajA;
    static double X;
    static double Y;
    static double Z;

//    static double chi2tot;
//    static double chi2ndf;

    int N_hits_FI;
    int hits_FI01;
    int hits_FI02;
    int hits_FI15;
    int N_hits_SI;

    int Bvtx = -1;
    double Xmin = 100;
    int ivok = e.iBestPrimaryVertex();
    

    if(first){ // histograms and Ntupes booking block
    // Phast::Ref().HistFileDir("UserEvent110");


    OutTree = new TTree("OutTree","OutTree");
	OutTree->Branch("TiS", &TiS, "TiS/D");
	OutTree->Branch("UnixSeconds", &UnixSeconds, "UnixSeconds/D");
	// OutTree->Branch("NTracks", &NTracks, "NTracks/I");
	OutTree->Branch("NVrtx", &NVrtx, "NVrtx/I");
	OutTree->Branch("Nbeam", &Nbeam, "Nbeam/I");
	// OutTree->Branch("R", &R, "R/I");
	OutTree->Branch("Nevent", &Nevent, "Nevent/I");
    
    first = false;
    
    }
    
    Nevent = e.UniqueEvNum();
    UnixSeconds = e.UnixSeconds();         // time of the event
    TiS = e.TimeInSpill();

//  bool extrapB		= false;		// Before target
//	bool extrapA		= false;        // After target

    // if(NVrtx == 0) return;

    if(TiS < 1.2 || TiS > 5.4) return;

    // CEDAR cut

    CE1m = 0;
	CE2m = 0;
        std::vector<PaDigit> rawDigits = e.RawDigits();
      hitMask[0] = 0;
      hitMask[1] = 0;
      for(int idig = 0; idig < (signed) rawDigits.size(); ++idig) { 

         bool is_ce1 = rawDigits[idig].DecodeMapName() == "CE01P1__";
         bool is_ce2 = rawDigits[idig].DecodeMapName() == "CE02P1__";
         int cedar = is_ce1 ? 1 : (is_ce2 ? 2 : 0);

         if(cedar){

            double hit_time = rawDigits[idig].DigInfo(3);
            int ipm = std::fabs(rawDigits[idig].IWire())-1;

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
	
	if( !e.IsMC() && (CE1m < 6 || CE2m < 6) ) return;

    // SciFi and SI cut

    for( int it = 0; it < e.NTrack(); it++){			// loop over tracks in target
	 	
        const PaTrack& beam_track = e.vTrack(it);

       if( beam_track.NTPar() == 0 ) continue;  				//Skip the track if it has no parameters
       if( beam_track.iParticle() == -1 ) continue;
       if( !beam_track.IsBeam() ) continue;

       //Chi2
    /*  chi2tot = beam_track.Chi2tot(); //Chi2
      int ndf = beam_track.Ndf();
      chi2ndf = chi2tot/double(ndf); //Reduced Chi2

      if( chi2ndf > 10. ) continue;  //Skipped tracks with chi2 > 10

      if( track.iParticle() < e.NParticle() ) {

        PaParticle part = e.vParticle(track.iParticle());
        if( part.NVertex() > 0 ) continue;                      //Let's skip particles with associated vertices
      }
     */  
       // cut on SciFis and SI
       hits_FI01 = beam_track.NHitsFoundInDetect("FI01"); 			
       hits_FI02 = beam_track.NHitsFoundInDetect("FI02");
       hits_FI15 = beam_track.NHitsFoundInDetect("FI15");
       N_hits_FI = hits_FI01 + hits_FI02 + hits_FI15;
       N_hits_SI = beam_track.NHitsFoundInDetect("SI");

       if(N_hits_FI < 4 && N_hits_SI < 10) continue;

       // beam track has to be inside target

       PaTPar tar_up, tar_down;
      if( !beam_track.Extrapolate(ztargetA,tar_up) || !beam_track.Extrapolate(ztargetB,tar_down) ) continue;  
      double x_up = tar_up.Pos(0);
      double y_up = tar_up.Pos(1);
      double x_down = tar_down.Pos(0);
      double y_down = tar_down.Pos(1);

      if( sqrt( x_up*x_up + y_up*y_up ) > rtarget ) continue;
      if( sqrt( x_down*x_down + y_down*y_down ) > rtarget ) continue;

      Nbeam ++;      // count beam tracks reconstructed correctly 
    }

    TgM = e.TrigMask();
    if( !e.IsMC() && !(TgM & 1 << 1) )  return;   // physics trigger


   Nout = e.NVertex();

   for(int jv = 0; jv < Nout ; jv++){ 

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

    // for(int ip = 0; ip < Nout; ip++){ 			 // loop over vertices in event 

	/*	const PaParticle& p = e.vParticle(ip); 
		if(p.iTrack() == -1) continue;

        const PaTrack& pt = e.vTrack(p.iTrack());
        PaTPar parB;			// parameter extrapolated before target
		PaTPar parA;			// parameter extrapolated after target

		extrapB = pt.Extrapolate(ztargetB, parB);
		extrapA = pt.Extrapolate(ztargetA, parA);

        XtrajB = parB.X();
		YtrajB = parB.Y();
		XtrajA = parA.X();
		YtrajA = parA.Y();


        if( sqrt((XtrajB*XtrajB) + (YtrajB*YtrajB)) > rtarget ) continue;
        if( sqrt((XtrajA*XtrajA) + (YtrajA*YtrajA)) > rtarget ) continue;

        if(Zvertex_rec < -70 || Zvertex_rec > 70) continue;  			// PV outside target region
        */

        // const PaVertex& v = e.vVertex(ip);
        X = v.X();
        Y = v.Y();
        Z = v.Z();

        // check on vertex position inside target volume

        // cout << " Nevent " << Nevent << " vertex position before cut: Z = " << Z << " X = " << X << " Y = " << Y;

        if(Z < -70 || Z > 70) return;

        if(sqrt((X*X)+(Y*Y)) > rtarget) return;
        
        // cout << " vertex position after cut: Z = " << Z << " X = " << X << " Y = " << Y << endl;

        NVrtx++;
   //  }

   // cout << " Nevent " << Nevent << " nbeam " << Nbeam << " nvrtx " << NVrtx << endl;
        
    // num interazioni in target/num particelle fascio

    OutTree->Fill();
    
    cout << "OutTree successfully save on /eos/home-g/gmeinard/2023DCS-correlation/" << endl;

}

