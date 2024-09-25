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


void UserEvent108(PaEvent& e){
			        
    static TH1D* h1[10];
	static TH2D* h2[10];
	static TTree* tree(NULL);
    static int Nout = 0;
    static double mom;
	static int q;
	static double Zvertex;
	static double Yvertex;
	static double Xvertex;
	static double X2;
	static int run;

//    static CEDAR* cedar = &CEDAR::Instance();

	int NVrtx = e.NVertex();
	int Nevent = e.UniqueEvNum();
	int N_tracks;			// number of incoming tracks
	int N_hits_FI;
	int N_hits_SI;
	

    static bool first(true);
    if(first){ // histograms and Ntupes booking block
    Phast::Ref().HistFileDir("UserEvent108");

    h1[1]  = new TH1D("mom","outgoing particle momentum distribution 1; GeV/c", 200, 0, 200);
	h1[2]  = new TH1D("z_vertex","Z primary vertex; cm", 400, -200, 200);
	h1[3]  = new TH1D("y_vertex","Y primary vertex; cm", 50, -4, 4);
	h1[4]  = new TH1D("x_vertex","X primary vertex; cm", 50, -4, 4);
	h1[5]  = new TH1D("N_tracks", "Number of beam tracks", 15, 0, 15);

	h2[1]  = new TH2D("XvsY","XY primary vertex; X [cm]; Y [cm]", 100, -5, 5, 100, -5, 5);
	h2[1] -> SetOption("colz");
	h2[2]  = new TH2D("XvsZ","XZ primary vertex; Z [cm]; X [cm]", 200, -100, 100, 80, -5, 5);
	h2[2] -> SetOption("colz");
	h2[3]  = new TH2D("YvsZ","YZ primary vertex; Z [cm]; Y [cm]", 200, -100, 100, 80, -5, 5);
	h2[3] -> SetOption("colz");
	h2[4]  = new TH2D("SciFivsSi","Beam track hits; Silicon Hits; SciFi Hits", 12, 0, 12, 8, 0, 8);
	h2[4] -> SetOption("colz");

    tree = new TTree("USR108","User Ntuple example");
	
	tree->Branch("mom", &mom, "mom/D");
	tree->Branch("Zvertex", &Zvertex, "Zvertex/D");
	tree->Branch("Yvertex", &Yvertex, "Yvertex/D");
	tree->Branch("Xvertex", &Xvertex, "Xvertex/D");
	tree->Branch("N_tracks", &N_tracks, "N_tracks/I");


    first=false;
    }

    //pbar 				= 0;
	int Bvtx 		    = -1;
	bool extrap 		= false;
	double Xmin 		= 100;
	int ivok 			= e.iBestPrimaryVertex();

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
	//cout << "indice vertice scelto " << ivok << endl; 

	Nout = v.NOutParticles();
	N_tracks = v.InParticle();


	for(int in_p = 0; in_p < N_tracks; in_p++){				// loop over incoming particles in target
	 	
	 	int inc_index = v.InParticle();
		const PaParticle& beam_p = e.vParticle(inc_index);			// creates a beam particle for every incoming track

		if(beam_p.iTrack() == -1) {
			cout << " Evento " << (Nevent - 197248024)/4 << " incoming particle's index -1" << endl;
			continue;
		}
		const PaTrack& beam_tr = e.vTrack(beam_p.iTrack());

		N_hits_FI = beam_tr.NHitsFoundInDetect("FI"); 			// should return n of hits in SciFi
		N_hits_SI = beam_tr.NHitsFoundInDetect("SI");
	//	cout << " Evento " << (Nevent - 197248024)/4 << " n hits in SciFi "  << N_hits_FI << " n hits in Silicons " << N_hits_SI << endl;

	}

	run = e.RunNum();
//  n_rich = 1 + PaMetaDB::Ref().NminusOne(run)/1e6; //refraction index
//  if(PaMetaDB::Ref().NminusOne(run) < 0) n_rich = 1 + PaSetup::Ref().NminusOne()/1e6;
	
	
	Zvertex = v.Z();
	Yvertex = v.Y();
	Xvertex = v.X();
	
	
	//incoming particle identification
/*	cedar->Search(e);
	
	// if you want to identify a proton (if at least one CEDAR was set on protons)
	 if (((cedar->PID(1)==2) || (cedar->PID(2)==2)) && (cedar->PID(1)!=0) && (cedar->PID(2)!=0)){
		 countproton++;
	    }
       	  else return;
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


		   mom = p.ParInVtx(ivok).Mom();
		   q = p.Q();
/*		  
		  if(Nout < 2 && theta < 0.008 && mom > 188){
					Rqe++;
					continue;
		  }		  

		  if(zlast < 363) continue; // zfirst < zSM1
		  zfirstOK++;
		  h1[16] -> Fill(2.);

		  if(Mtime > 1e4) continue;
		  MtimeOK++;
		  h1[16] -> Fill(3.);

		  if(Zvertex < -68.52 || Zvertex > -28.52) continue;  //PV outside target region
		  zvertexOK++;
		  h1[16] -> Fill(4.);

		  if(Xrid > 10) continue;
		  XridOK++;
		  h1[16] -> Fill(5.);
			
		  if(abs(Yvertex) > 1.75) continue;

	      if(abs(Xvertex) > 1.75) continue;

		  if(radL > 10) continue; //no muons

		  if( r < 5 ) continue;   //scarto beam pipe rich

		  if (mom < 10 || mom > 45) continue;
		  
		  if (theta < 0.01 || theta > 0.18) continue;
		  
		  if(tm & 0 << 0 ) continue;	

		  if(mom > 190) continue;

 		 int count = 0;
  		//int N = 0;
		  //cout << "likelihood ";
		  for(int i = 0; i < 6; i++) {
			  like[i] = pid.GetLike(i,pt);
		      count += like[i];
			//cout << like[i] << " " ;
		     }
		  //cout << endl;
*/

	// ******************** filling histos **********************
		  h1[1] -> Fill(mom);
		  h1[2] -> Fill(Zvertex);
		  h1[3] -> Fill(Yvertex);
		  h1[4] -> Fill(Xvertex);
		  h1[5] -> Fill(N_tracks + 1);

		  h2[1] -> Fill(Xvertex, Yvertex);
		  h2[2] -> Fill(Zvertex, Xvertex);
		  h2[3] -> Fill(Zvertex, Yvertex); 
		  h2[4] -> Fill(N_hits_SI, N_hits_FI); 
		 
		  
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