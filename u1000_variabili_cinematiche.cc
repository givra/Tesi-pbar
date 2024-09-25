#include <iostream>
#include <cmath>
#include <stdio.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"
#include "PaParticle.h"
#include "PaAlgo.h"
#include "G3part.h"
#include "PaMCvertex.h"
#include "PaMCtrack.h"
#include "PaPid.h"


 void UserEvent1000(PaEvent& e){
	 
	 static TH1D* h1[50];
     static TTree* tree0(NULL);
	 static TTree* tree1(NULL);
	 static int q;
	 static int mc_q;
	 static double mom;
	 static double MCmom;
	 static double zfirst;
	 static double zlast;
	 static double MC_mass;
	 static double radL;
	 static double x;
	 static double y;
	 static double z;
	 static double like[7];
	 static int N0 = 0;
	 static int N1 = 0;
	 static int N2 = 0;
	 static int hits;
	 static double chi2;
	 static int ndf;
	 static double Mtime;
	 static int npar;
	 static double dip;
	 static double azi;
	 static double phi;
	 static double theta;
	 static int NOut;
	 PaPid pid;

    static bool first(true);
    if(first){ // histograms and Ntupes booking block
    Phast::Ref().HistFileDir("UserEvent1000");
	
	h1[1]  = new TH1D("h01","particle charge", 10, -2, 2); 
	h1[2]  = new TH1D("h02","MC charge", 10,-2, 2);
	h1[3]  = new TH1D("h03","particle momentum distribution", 100, -10, 250);
	h1[4]  = new TH1D("h04","particle MC momentum distribution", 100, -10, 250);
	h1[5]  = new TH1D("h05","Z first", 100, 0, 2000);
	h1[6]  = new TH1D("h06","Z last", 100, 0, 7000);
	h1[8]  = new TH1D("h08","Outgoing particles", 30, 0.5, 30.5);
	h1[9]  = new TH1D("h09","MC mass distribution", 100, -1, 1);
	h1[10]  = new TH1D("h10","radiation length", 100, -2, 40); //40 for lhe and 400 for dli data
	h1[11]  = new TH1D("h11","wrong charge", 10, -2, 2);
	h1[12]  = new TH1D("h12","particles counter", 40, 0, 10);
	
	h1[20]  = new TH1D("h20","hits", 100, -0.5, 600);
	h1[21]  = new TH1D("h21","X^2", 100, -0.5, 600);
	h1[22]  = new TH1D("h22","mean time", 100, -0.5, 600);
	h1[23]  = new TH1D("h23","# parameters", 10, 0, 10);
	h1[24]  = new TH1D("h24","dip", 100, -0.5, 0.5);
	h1[25]  = new TH1D("h25","azi", 100, -0.5, 0.5);
	h1[26]  = new TH1D("h26","phi", 100, -5, 5);
	h1[27]  = new TH1D("h27","theta", 100, 0, 0.5);
	
	
   tree0 = new TTree("USR1000","User Ntuple example"); //for all particles
   tree1 = new TTree("USR1001","User Ntuple example 2"); //wrong charge
   
   tree0->Branch("q",    &q,    "q/I");
   tree0->Branch("mc_q",    &mc_q,    "mc_q/I");
   tree0->Branch("mom", &mom, "mom/D");
   tree0->Branch("Mcmom", &MCmom, "MCmom/D");
   tree0->Branch("zfirst", &zfirst, "zfirst/D");
   tree0->Branch("zlast", &zlast, "zlast/D");
   tree0->Branch("MC_mass", &MC_mass, "MC_mass/D");
   tree0->Branch("radL", &radL, "radL/D");
   tree0->Branch("x", &x, "x/D");
   tree0->Branch("y", &y, "y/D");
   tree0->Branch("z", &z, "z/D");
   tree0->Branch("NOut", &NOut, "NOut/I");
   tree1->Branch("hits", &hits, "hits/I");
   tree1->Branch("chi2", &chi2, "chi2/D");
   tree1->Branch("ndf", &ndf, "ndf/I");
   tree1->Branch("Mtime", &Mtime, "Mtime/D");
   tree1->Branch("npar", &npar, "npar/I");
   tree1->Branch("dip", &dip, "dip/D");
   tree1->Branch("azi", &azi, "azi/D");
   tree1->Branch("phi", &phi, "phi/D");
   tree1->Branch("theta", &theta, "theta/D");
   tree1->Branch("like", &like, "like[7]/D");
   
   first=false;
	}
	int NVrtx = e.NVertex();
	int Nevent = e.UniqueEvNum();
	//cout << endl << endl;
	cout << "# event " << Nevent ;
	//cout << "# particles" << e.NParticle() << endl;
	
	for(int iv = 0; iv < NVrtx ; iv++){ // loop over reconstructed vertices
    const PaVertex& v = e.vVertex(iv);
	int Bvtx  =  e.iBestPrimaryVertex();
    NOut = v.NOutParticles();
    if(Bvtx !=-1 && Bvtx != iv) continue; //per evitare che conti la stessa particella più volte	
	
	   x = v.X();
	   y = v.Y();
	   z = v.Z();
	   
	   for(int ip = 0; ip < NOut; ip++){ //loop over outgoing particles from each vertex
		   //cout << "at vertex " << iv << " there are " << NOut << " outgoing particles" << endl;
		   
		   int index = v.iOutParticle(ip);
		   const PaParticle& p = e.vParticle(index); 
		   
		   if(p.iTrack() == -1) continue;
		   
		   const PaTrack& pt = e.vTrack(p.iTrack());
		   const PaTPar& par = pt.vTPar(0);
		   
		   //int mc_index = pt.iMCtrack();
		   
		   //if(mc_index == -1 ) continue;
		   //const PaMCtrack& p_mc = e.vMCtrack(mc_index);
			
			    mom = p.ParInVtx(iv).Mom();
			    q = p.Q();
			    //MCmom = p_mc.P(); 
			    zfirst = pt.ZFirst();
			    zlast = pt.ZLast();
			    /*const string name = p_mc.Name();
			    MC_mass = p_mc.MCmass();
			    mc_q = p_mc.Q();
			    */radL = pt.XX0();
				npar = pt.NTPar();
				
				double count = 0;
				
			   //cout << name << " " << "likelihood ";
			   cout << " likelihood: ";
			  for(int i = 0; i < 6; i++) {like[i] = pid.GetLike(i,pt);
			    cout << like[i] << " ";
				count += like[i];
			   }
			   //cout << "somma comp " << count << endl;
               cout << endl;
			   
			   if(count == 0) N0++;
			   else if(count == -6) N1++;
			   else N2++;
				
			   
			   if(q != mc_q){
				   /*cout << name << " " << "likelihood ";
			   for(int i = 0; i < 6; i++) {like[i] = pid.GetLike(i,pt);
			     cout << like[i] << " ";
			     count += like[i];
			   }*/
			   //cout << endl;
			   //cout << " " << name << " charge " << q << " MC charge " << mc_q << endl;
			   h1[11]->Fill(q);
			   hits = pt.NHits();
			   chi2 = pt.Chi2tot();
			   ndf = pt.Ndf();
			   Mtime = pt.MeanTime();
			   //npar = pt.NTPar();
			   dip = par.Dip();
			   azi = par.Azi();
			   phi = par.Phi();
			   theta = par.Theta();
			   //cout << name << " mean time " << Mtime << endl;
			   h1[20]->Fill(hits);
			   h1[21]->Fill(chi2);
			   h1[22]->Fill(Mtime);
			   //h1[23]->Fill(npar);
			   h1[24]->Fill(dip);
			   h1[25]->Fill(azi);
			   h1[26]->Fill(phi);
			   h1[27]->Fill(theta);
			   
			   //------------------------------------------------------------------------------
			   //positron, electron, kaon, proton, muon, pion, other
			   
			   
			   /*const int& kind = p_mc.Pid();
			   switch (kind)
			   { case 2 :
					h1[12] -> Fill(0.); //positron
					break;
			     case 3 :
					h1[12] -> Fill(1.); //electron
					break;
			     case 10 :
				 case 11 : 
				 case 12 :
					h1[12] -> Fill(2.); //kaon
					break;
			     case 14 :
					h1[12] -> Fill(3.); //proton
					break;
				 case 5 :
				 case 6 :
				    h1[12] -> Fill(4.); //muon
				 case 7 :
				 case 8 :
				 case 9 :
					h1[12] -> Fill(5.); //pion
					break;
			    default: 
					h1[12]-> Fill(6.); //other
				}*/
				
			   tree1->Fill();
			   }
			   
			   h1[1]->Fill(q);
			   h1[2]->Fill(mc_q);
		       h1[3]->Fill(mom);
			   h1[4]->Fill(MCmom);
			   h1[5]->Fill(zfirst);
			   h1[6]->Fill(zlast);
			   h1[8]->Fill(NOut);
			   h1[9]->Fill(MC_mass);
			   h1[10]->Fill(radL);
			 h1[23]->Fill(npar);
            
		tree0->Fill();
	   }
	}
	//if(Nevent == 1048577 + 5147){ //+5147 for lhe, +1913 for dli data
		//cout << endl << "nulli: " << N0 << " tutti -1: " << N1 << " normali: " << N2 << " TOT: " << N0+N1+N2 << endl;
		//cout << endl << "concerning h13: 0 = positron | 1 = electron | 2 = kaon | 3 = proton | 4 = muon | 5 = pion | 6 = other ";  
 }
 
  
 

