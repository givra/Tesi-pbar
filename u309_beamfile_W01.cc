// Riccardo Longo's UE for DY beam file creation
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

// /eos/experiment/amber/pbar/common/data23/RD/productions/W01t0.8rich1.6/mDST/mDST-300537.root.001

// double Rndm();
/* modified by Davide G. and Tomas Klasek for pbar use */

double tar_radius = 4; // cm
double tarUp = -70; 
double tarDown = 70; 


void UserEvent309(PaEvent& e){

    // Histograms and Trees pointers
    static TTree * btree = NULL;
    static TH1D* beam_i[10];
    static TH1D* beam_v[10];
    static TH2D* td_beam_i[10];
    static TH2D* td_beam_v[10];

    static TH1D* hced1 = NULL;
    static TH1D* hced2 = NULL;
    static TH1D* stats = NULL;
    static TH1D* htis = NULL;
    static TH1D* hnhits = NULL;

    // static TH2F* h2[10];
    static TH2F* h2xyTar[2];
    static TH1D* htime =NULL;
    static TH1D* hetime =NULL;
    static double dXdZ, dYdZ, beamX, beamY, beamZ, beamMomentum;
    static double beamParticleFlag; 
    static double X1,Y1; // dX1dZ,dY1dZ;

    // static string  sOut = "out.dat";
    // static ofstream out ( sOut.c_str() );

    //Debug flag
    int debug = 0;
    int Nevent = e.UniqueEvNum();

    //Parameters at z = -760 cm
    static double iPos[3];
    static double iSlope[2];


    static int TgM;
    static int Bin_Trig;
    static unsigned int hitMask[2];

    static unsigned int CE1m;
    static unsigned int CE2m;
    static double Zmin, Zmax, lastZ, firstZ;
    static double chi2tot;
    static double chi2ndf;

    //Histogram variables
    string histoname;
    string cutname;
    int nbin = 100;

    //Histograms booking
    static bool first(true);
    if(first){ // histograms and Ntupes booking block


        // h2[0] = new TH2F("hCE1_t_ch","CE1 time vs channel; channel; time",34,0,34,1000,-4000,-3000);
        // h2[1] = new TH2F("hCE2_t_ch","CE2 time vs channel; channel; time",34,0,34,1000,-4000,-3000);

    
    
    /*    cutname = "beam_profile_SCIFI01";
        histoname = cutname + "_0";
        beam_i[0] = new TH1D(histoname.c_str(),"X profile (cm)",nbin,-5.,5.);
        histoname = cutname + "_1";
        beam_i[1] = new TH1D(histoname.c_str(),"Y profile (cm)",nbin,-5.,5.);
        histoname = cutname + "_2";
        beam_i[2] = new TH1D(histoname.c_str(),"X slope profile (mrad)",nbin,-20.,20.);
        histoname = cutname + "_3";
        beam_i[3] = new TH1D(histoname.c_str(),"Y slope profile (mrad)",nbin,-20.,20.);

        histoname = "TD_" + cutname + "_0";
        td_beam_i[0] = new TH2D(histoname.c_str(),"X:Y distribution (z = -760.)",nbin,-5.,5.,nbin,-5.,5.);
        histoname = "TD_" + cutname + "_1";
        td_beam_i[1] = new TH2D(histoname.c_str(),"#theta_{X}:#theta_{Y} distribution (z = -760.)",nbin,-5.,5.,nbin,-5.,5.);

        cutname = "beam_profile_vertex";
        histoname = cutname + "_0";
        beam_v[0] = new TH1D(histoname.c_str(),"X profile (cm)",nbin,-5.,5.);
        histoname = cutname + "_1";
        beam_v[1] = new TH1D(histoname.c_str(),"Y profile (cm)",nbin,-5.,5.);
        histoname = cutname + "_2";
        beam_v[2] = new TH1D(histoname.c_str(),"X slope profile (mrad)",nbin,-20.,20.);
        histoname = cutname + "_3";
        beam_v[3] = new TH1D(histoname.c_str(),"Y slope profile (mrad)",nbin,-20.,20.);

        histoname = "TD_" + cutname + "_0";
        td_beam_v[0] = new TH2D(histoname.c_str(),"X:Y distribution (closerr point to spec.)",nbin,-5.,5.,nbin,-5.,5.);
        histoname = "TD_" + cutname + "_1";
        td_beam_v[1] = new TH2D(histoname.c_str(),"#theta_{X}:#theta_{Y} distribution (closer point to spec.)",nbin,-5.,5.,nbin,-5.,5.);
*/
	      h2xyTar[0] = new TH2F("hxyTar0","X:Y distribution (z = -70. cm)",nbin,-5.,5.,nbin,-5.,5.);
        h2xyTar[1] = new TH2F("hxyTar1","X:Y distribution (z = 70. cm)",nbin,-5.,5.,nbin,-5.,5.);


        stats = new TH1D("stats","Event selection", 14, 0, 14);
        htis = new TH1D("htis","Time in spill", 1000, 0, 6);
        hced1 = new TH1D("hced1","CE01P1__", 33, 0, 33);
        hced2 = new TH1D("hced2","CE02P1__", 33, 0, 33);
        htime = new TH1D("htime","Beam time",400,-100,100);
        hetime = new TH1D("hetime","Beam time error",100,0,5);
        hnhits = new TH1D("hnhitsBT","Hits in beam telescope; n hits; Tracks",21,0,21);

        btree = new TTree("BTree","BTree");
        // btree->Branch("x_FI", &iPos[0],"x_FI/D");
        // btree->Branch("y_FI", &iPos[1],"y_FI/D");
        // btree->Branch("z_FI", &iPos[2],"z_FI/D");
        // btree->Branch("x_slope_FI", &iSlope[0],"x_slope_FI/D");
        // btree->Branch("y_slope_FI", &iSlope[1],"y_slope_FI/D");
        // btree->Branch("x_vtx", &vPos[0],"x_vtx/D");
        // btree->Branch("y_vtx", &vPos[1],"y_vtx/D");
        // btree->Branch("z_vtx", &vPos[2],"z_vtx/D");
        // btree->Branch("x_slope_vtx", &vSlope[0],"x_slope_vtx/D");
        // btree->Branch("y_slope_vtx", &vSlope[1],"y_slope_vtx/D");

        // Needed for beam file creation
        btree->Branch("X", &beamX,"X/D");   // X info using Extrapolate
        btree->Branch("Y", &beamY,"Y/D");
        btree->Branch("Z", &beamZ,"Z/D");
        btree->Branch("X1", &X1, "X/D");    // X info using PaHit
        btree->Branch("Y1", &Y1, "Y/D");    // Y info using PaHit
        btree->Branch("dXdZ", &dXdZ,"dXdZ/D");
        btree->Branch("dYdZ", &dYdZ, "dYdZ/D");
        // btree->Branch("dX1dZ", &dX1dZ,"dX1dZ/D");
        // btree->Branch("dY1dZ", &dY1dZ, "dY1dZ/D");
        btree->Branch("P", &beamMomentum,"P/D");
        btree->Branch("particleFlag", &beamParticleFlag,"particleFlag/D");
        btree->Branch("Chi2", &chi2tot, "Chi2/D");
        btree->Branch("Chi2_Red", &chi2ndf, "chi2ndf/D");
        btree->Branch("TgM", &TgM,"TgM/I");
        btree->Branch("Bin_Trig", &Bin_Trig,"Bin_Trig/I");
        first = false;

    }

    beamParticleFlag = 1;

    const PaSetup& setup = PaSetup::Ref();
    double z_FI = -760.;// not FI01.Z() since in TGEANT we set it by hand
    
    int ind = 0;
    stats->Fill(ind++);

    TgM = e.TrigMask();
    if( !e.IsMC() && !(TgM & 1 << 10) )  return;    // random trigger
    //if( !e.IsMC() && !(TgM & 1 << 1) )  return;   // physics trigger

    double SpillTime = e.TimeInSpill();
    htis->Fill(SpillTime);
    stats->Fill(ind++);

    if( !e.IsMC() && ( SpillTime < 1.2 || SpillTime > 5.4) ) return;
    stats->Fill(ind++);

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
            // if( is_ce1 ) {
            //     h[9]->Fill(hit_time);
            //     h2[3]->Fill(ipm,hit_time);
            // }	
            // if( is_ce2 ) {
            //     h[10]->Fill(hit_time);
            //     h2[4]->Fill(ipm,hit_time);
            // }

            // calibration by "hand" for W01, looking at time plots

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



    hced1->Fill(CE1m);
    hced2->Fill(CE2m);
  
    if( !e.IsMC() && (CE1m < 6 || CE2m < 6) ) return;
    stats->Fill( ind++);

    for( int it = 0; it < e.NTrack(); it++){
  
      ind = 4;
      stats->Fill(ind++); // 4

      const PaTrack& track = e.vTrack(it);

      if( track.NTPar() == 0 ) continue;  //Skip the track if it has no parameters

      stats->Fill(ind++); //5

      const PaTPar& par = track.vTPar(0); //Extract first measured track parameters. Be careful: if it is beam this are parameters of the track approaching vtx


      //Zmin and Zmax
      Zmin = track.Zmin();
      Zmax = track.Zmax();
      firstZ   = track.ZFirst();
      lastZ   = track.ZLast();

      //Chi2
      chi2tot = track.Chi2tot(); //Chi2
      int ndf = track.Ndf();
      chi2ndf = chi2tot/double(ndf); //Reduced Chi2

      if( chi2ndf > 10. ) continue;  //Skipped tracks with chi2 > 10

      stats->Fill(ind++); // 6
     

      // Cut on time
      double time    = track.MeanTime();    // track time
      double etime   = track.SigmaTime();   // errore of track time


      if( track.iParticle() == -1 ) continue;

      stats->Fill(ind++); // 7

      if( track.iParticle() < e.NParticle() ) {

        PaParticle part = e.vParticle(track.iParticle());
        if( part.NVertex() > 0 ) continue; //Let's skip particles with associated vertices
      }

      stats->Fill(ind++); // 8
      
      if( !track.IsBeam() ) continue;

      stats->Fill(ind++); // 9
      

      int nHitsFI = track.NHitsFoundInDetect("FI01") + track.NHitsFoundInDetect("FI02") + track.NHitsFoundInDetect("FI15");
	    int nHitsSI = track.NHitsFoundInDetect("SI");
      //int nHitsBT = track.NHitsFoundInDetect("FI01") + track.NHitsFoundInDetect("FI02") + track.NHitsFoundInDetect("FI15");
      hnhits->Fill(nHitsFI + nHitsSI);
    	if ( nHitsFI < 4 || nHitsSI < 10 ) return; // at least 6 hits in scifi + 9 in 

    //  if( nHitsBT < 6 ) continue;

      stats->Fill(ind++); // 10

      htime->Fill(time);
      hetime->Fill(etime);
      // Time cuts
      if( etime <= 0.) continue;
      if( time >= 10000 ) continue;

      //Check if the time is into 3sigma in a time window defined  as +-6ns
      if( time - 3.* etime > 3. ) continue;
      if( time + 3.* etime < -3. ) continue;

      stats->Fill(ind++); // 11


      // parameters in front of the target
      PaTPar tar_up, tar_down;
      if( !track.Extrapolate(tarUp,tar_up) || !track.Extrapolate(tarDown,tar_down) ) continue;  
      h2xyTar[0]->Fill(tar_up(1),tar_up(2));
      h2xyTar[1]->Fill(tar_down(1),tar_down(2));

      double x_up = tar_up.Pos(0);
      double y_up = tar_up.Pos(1);
      double x_down = tar_down.Pos(0);
      double y_down = tar_down.Pos(1);

      if( sqrt( x_up*x_up + y_up*y_up ) > tar_radius ) continue;
      if( sqrt( x_down*x_down + y_down*y_down ) > tar_radius ) continue;

      // parameters at FI01
      PaTPar first_par;
      if( !track.Extrapolate(z_FI, first_par) ) continue;

     // *********************************************** trying to get Hits from SciFi *******************************************************
      const vector<int>& HitRef = track.vHitRef();      // indices of which hits from "hits" belong to this specific track 
      const vector<PaHit>& AllHits = e.Hits();          // total number of hits registered in a single event     

      for(int k=0; k<HitRef.size();k++){
        const PaHit& hitBeam = AllHits[HitRef[k]];
        // for(int ii=0;ii<hitBeam.size();ii++){
                X1 = hitBeam.X();
                Y1 = hitBeam.Y();
        // }
      }
      
      // X1 = beamPar.X();
      // Y1 = beamPar.Y();
      // dX1dZ = beamPar.dXdZ();
      // dY1dZ = beamPar.dYdZ();

      stats->Fill(ind++); // 12

      iPos[0] = first_par(1);
      iPos[1] = first_par(2);
      iPos[2] = first_par(0);
      iSlope[0] = first_par(3)*1000; // mrad
      iSlope[1] = first_par(4)*1000; // mrad
      float ibmom = first_par.Mom();

      beamX = first_par(1); // cm
      beamY = first_par(2);  // cm
      beamZ = first_par(0); // cm
      dXdZ = first_par(3)*1000; // mrad
      dYdZ = first_par(4)*1000; // mrad
      beamMomentum = first_par.Mom(); // Gev/c

    
      //Filling Histograms at z = -760
    //  beam_i[0]->Fill(iPos[0]);
    //  beam_i[1]->Fill(iPos[1]);
    //  beam_i[2]->Fill(iSlope[0]);
    //  beam_i[3]->Fill(iSlope[1]);
    //  td_beam_i[0]->Fill(iPos[0],iPos[1]);
    //  td_beam_i[1]->Fill(iSlope[0],iSlope[1]);

      //Filling tree
      btree->Fill();

    }//End of the loop on the tracks
    

} 
