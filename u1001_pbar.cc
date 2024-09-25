#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
//#include "TLorentzVector.h"
//#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaDetect.h"
#include "PaRich.h"
#include "PaRichDet.h"
#include "PaEvent.h"
#include "PaParticle.h"
#include "PaAlgo.h"
#include "G3part.h"
//#include "PaMCvertex.h"
//#include "PaMCtrack.h"
#include "PaPid.h"
#include "PaMetaDB.h"
//#include "CEDAR_Helper.h"
#include "PaMtx.h"
#include "TMath.h"


 // /cmps-data/2009_HData/RunsW29/mDST-77105-2-7.root

class PaPid1: public PaPid{
	
 public: 
int IsProton(double mom, double mom_p, double mom_k, double mom_pi, double like[6], double count, int charge){ 
	 //cout  << "mom p " << mom_p << " mom k " << mom_k << " mom pi " << mom_pi << endl;
	 
	 int delta = 3;

	if( mom > (mom_p - delta) ){   //SOPRA SOGLIA

		double pi_k = like[0]/like[1];
		double pi_p = like[0]/like[2];
		double pi_bkg = like[0]/like[5];
		if(pi_k > 1.01 && pi_p > 1.01 && pi_bkg > 2.02){
			 if(charge < 0) return 0; //pion-
			 else return 4; //pion+
		}
		
		double k_pi = like[1]/like[0];
		double k_p = like[1]/like[2];
		double k_bkg = like[1]/like[5];
		if(k_pi > 1.08 && k_p > 1.08 && k_bkg > 1.5){
			 if(charge < 0) return 1; //kaon-
			 else return 5; //kaon+
		}

		
		double p_pi = like[2]/like[0];
		double p_k = like[2]/like[1];
		double p_bkg = like[2]/like[5];
		if(p_pi > 1. && p_k > 1. && p_bkg > 1.){
			 if(charge < 0) return 2; //pbar
			 else return 6; //proton
		}
	
	}
	
	if(mom < (mom_p + delta)){   //SOTTO SOGLIA

			
			double pi_k = like[0]/like[1];
			double pi_p = like[0]/like[2];
			double pi_bkg = like[0]/like[5];
			if(pi_k > 1.01 && pi_p > 1.01 && pi_bkg > 2.02){
				 if(charge < 0) return 0; //pion-
				 else return 4; //pion+
			}
		
		
			
			double k_pi = like[1]/like[0];
			double k_p = like[1]/like[2];
			double k_bkg = like[1]/like[5];
			if(k_pi > 1.08 && k_p > 1.08 && k_bkg > 1.5){
				 if(charge < 0 ) return 1; //kaon-
				 else return 5; //kaon+
			}
		
		
			double valueK = like[1]/like[5];
			double valuePi = like[0]/like[5];
			if(charge == 1){
				if((valueK < 1.5 && valuePi < 2.) || count == 0 ){
					 return 6; //proton
				}
			}
			else{
				if((valueK < 1.5 && valuePi < 2.) || count == 0 ){
					 return 2;  //pbar
				}
			}
		
	
	}
	return 3;
   }
};

double efficienze [12][3][3] = {
					{{0.967826,0.0251999,0.0697725},{0.00102249,0.836424,0.0535878},{0.0300532,0.126551,0.862242}},    		//dati00
					{{0.979647,0.02426,0.0642111},{0.00162168,0.947048,0.084379},{0.0160502,0.0230714,0.79981}},			//dati01
					{{0.982884,0.014773,0.0457637},{0.00131405,0.966517,0.0382535},{0.00331241,0.00751702,0.785299}},		//dati02
					{{0.981186,7.24089e-08,0.0284578},{0.00141763,0.959957,0.00913336},{0.00263844,0.0081946,0.925603}},	//dati03
					{{0.958464,0.00565461,0.0147763},{0.00366605,0.948489,0.00667184},{0.00227172,0.000601679,0.954042}},	//dati04
					{{0.857817,2.2104e-07,0.00503792},{0.010988,0.818245,0.00619428},{0.0017985,0.00529184,0.965162}},		//dati05
					{{0.967812,0.00144269,0.0227235},{0.00248539,0.876341,0.0722492},{0.0286409,0.115437,0.891952}},		//dati10
					{{0.971345,0.0226833,0.0271211},{0.0030422,0.936809,0.0658102},{0.0230537,0.0377323,0.878434}},			//dati11
					{{0.974118,3.71147e-10,0.0192414},{0.00194095,0.951378,0.0232427},{0.00929441,0.0257791,0.87768}},		//dati12
					{{0.975525,5.29429e-11,0.0122957},{0.00157261,0.959423,0.00646673},{0.00667612,0.0216172,0.952654}},	//dati13
					{{0.953714,3.97383e-11,0.00822642},{0.00447053,0.981585,0.00508667},{0.00511445,0.0103265,0.96623}},	//dati14
					{{0.859777,1.20066e-08,6.85031e-10},{0.0170868,0.989092,0.00118705},{0.00192132,1.29595e-05,0.983289}}	//dati15
								};
								
double errori [12][3][3] = {
					{{0.000152164,0.0179319,0.00158208},{3.61073e-05,0.0163739,0.00100831},{0.00014552,0.0059659,0.00186016}},
					{{0.000138027,0.00975522,0.0011794},{5.02397e-05,0.0099065,0.000929932},{0.000120575,0.00271357,0.00152446}},
					{{0.00015713,0.00888223,0.00102449},{5.81096e-05,0.00915076,0.000704431},{8.06075e-05,0.00191275,0.00141642}},
					{{0.000206632,0.000873041,0.00116966},{7.32427e-05,0.00414791,0.000427565},{0.000100291,0.00235344,0.00145836}},
					{{0.000358885,0.0140383,0.00101044},{0.000118652,0.0146241,0.000407208},{0.00013157,0.00340817,0.00130197}},
					{{0.000604875,0.0375065,0.000780519},{0.000184397,0.0324712,0.00035922},{0.000130466,0.00270179,0.00111302}}, 
					{{0.000115092,0.00242289,0.000872943},{3.63036e-05,0.00540367,0.000809974},{0.000107905,0.00473567,0.00119535}},
					{{0.000204387,0.0259114,0.000973983},{7.66489e-05,0.0255019,0.00081961},{0.000182526,0.00452862,0.0013233}},
					{{0.000390634,0.00203627,0.00157913},{0.00013888,0.00867371,0.000817997},{0.000268089,0.00640309,0.00227446}},
					{{0.00075334,0.0065859,0.00256997},{0.000259211,0.0139097,0.000848292},{0.00047908,0.00955645,0.00294325}},
					{{0.00174887,0.0257759,0.00312922},{0.000615972,0.0361922,0.00114284},{0.000924548,0.0178792,0.00383248}},
					{{0.00427428,0.0528177,0.000864027},{0.00160424,0.306597,0.00193717},{0.00144956,0.171047,0.00444718}}
							};


    int countproton = 0;
	int countproton1 = 0;
	int countAproton = 0;	
	int NoutPV = 0;
	int trackindex = 0;
	int zfirstOK = 0;
	int MtimeOK = 0;
	int zvertexOK = 0;
	int XridOK = 0;
	int Ri = 0; //inelastic events 
	int Rs = 0; //events with 1 or more pbar
	int Rqe = 0; //quasi elastic events
	int Eventifile = 0;
	int Totevent;  //total number of events per file
	int Procevent; //# of processed events
	int contatorevertici = 0;
	double N[12][3];  //matrice in cui andranno le molteplicità corrette, ogni riga corrisponde ad un bin
	const double sigma_inel = 31.8682; //sigma inelastica totale
	const float err_sig = 0.73; //errore sigma inel

	cout << "inziamo dall'evento..." << endl;

void UserEvent1001(PaEvent& e){

	Totevent = Phast::Ref().NevTot;
	Procevent = Phast::Ref().NevProc;
			

	static double like[6];
	static TH1D* h1[50];
	static TH2D* h2[50];
	static TTree* tree(NULL);
    static int Nout = 0;
	static double mom;
	static int q;
	static double Zvertex;
	static double Yvertex;
	static double Xvertex;
	static double zfirst;
    static double zlast;
	static double radL;
	static int hits;
	static double chi2;  // X^2 for particle
	static double X2;   // X^2 for vertex
	static double Mtime;
	static int ndf;
	static double Xrid;
	static double theta;
	static double eta;
	static double Xtraj; 
	static double Ytraj; 
	PaPid1 pid;
	static int proton;
	static int run;
	static double n_rich;
	const double M_p = G3partMass[14]; //proton's mass
	const double M_k = G3partMass[10]; //kaon's mass
	const double M_pi = G3partMass[9]; //pion-'s mass
	static double mom_p; //threshold momentum protons
	static double mom_k; //threshold momentum kaons
	static double mom_pi; 
	static int tm;
	static int pbar;
	int numPi = 0;
	int numK = 0;
	int numPbar = 0;
	static float determinant[12]; //ogni componente è il det di una matrice in un dato bin (p,theta)
	static double err_det_sq[12];
	static double err_det[12];
	static double err_cofatt[12][3][3];
	static double err_inv[12][3][3];
	static float inversa[12][3][3];
	static double cofatt[12][3][3];
	static double err_N[12][3];
	

//	static CEDAR* cedar = &CEDAR::Instance();

	int NVrtx = e.NVertex();
	int Nevent = e.UniqueEvNum();
	
	double theta_bins[3] = {0.01,0.03,0.18};
	double p_bins[7] = {10,15,20,25,30,35,45};
	
	static bool first(true);
    if(first){ // histograms and Ntupes booking block
    Phast::Ref().HistFileDir("UserEvent1001");
	
	
	h1[1]  = new TH1D("p","outgoing particle momentum distribution 1; GeV/c", 100, -10, 220);
	h1[2]  = new TH1D("z_vertex","Z primary vertex; cm", 100, -144, 56);
	h1[3]  = new TH1D("y_vertex","Y primary vertex; cm", 100, -4, 4);
	h1[4]  = new TH1D("x_vertex","X primary vertex; cm", 100, -4, 4);
	h1[5]  = new TH1D("particles_counter","particles counter", 8, 0, 8);
	h1[6]  = new TH1D("z_first","Z first ; cm ", 100, 0, 600);
	h1[7]  = new TH1D("z_last","Z last ; cm", 100, 0, 4000);
	h1[8]  = new TH1D("XX0","radiation length ; cm", 100, 0, 3);
	h1[9]  = new TH1D("hits","hits ", 100, -0.5, 150);
	h1[12]  = new TH1D("X^2","X^2 ", 100, -0.5, 400);
	h1[13]  = new TH1D("x_traj","X trajectory", 100, -50, 50);
	h1[14]  = new TH1D("y_traj","Y trajectory", 100, -50, 50);
	h1[15]  = new TH1D("X^2_vertex","X^2 vertex", 100, -0.5, 100);
	h1[16]  = new TH1D("statistical_table","statistic table", 7, 0, 7);
	h1[17]  = new TH1D("#pi","molteplicità pi-", 6, 0, 6);
	h1[18]  = new TH1D("#k","molteplicità k-", 6, 0, 6);
	h1[19]  = new TH1D("#pbar","molteplicità pbar", 6, 0, 6);
	
	//Xsec projection on p axis
	h1[20]  = new TH1D("pi_Xsec", "pi X-section; GeV/c; b", 6, p_bins);
	h1[21]  = new TH1D("k_Xsec", "k X-section; GeV/c; b,", 6, p_bins);
	h1[22]  = new TH1D("pbar_Xsec", "pbar X-section; GeV/c; b", 6, p_bins);

	
	//2-D histogram
	h2[1]  = new TH2D("theta_p_pi","theta vs p |pi; GeV/c; rad", 6, p_bins, 2, theta_bins); //first 3 param are for p, the rest for theta
	h2[1] -> SetOption("lego");
	h2[2]  = new TH2D("theta_p_k","theta vs p |k; GeV/c; rad", 6, p_bins, 2, theta_bins); //first 3 param are for p, the rest for theta
	h2[2] -> SetOption("lego");
	h2[3]  = new TH2D("theta_p_pbar","theta vs p |pbar; GeV/c; rad", 6, p_bins, 2, theta_bins); //first 3 param are for p, the rest for theta
	h2[3] -> SetOption("lego");
	//the following 3 concern the corrected count
	h2[6]  = new TH2D("correzione_pi","theta vs p corr |pi; GeV/c; rad", 6, p_bins, 2, theta_bins); 
	h2[6] -> SetOption("lego");
	h2[7]  = new TH2D("correzione_k","theta vs p corr |k; GeV/c; rad", 6, p_bins, 2, theta_bins); 
	h2[7] -> SetOption("lego");
	h2[8]  = new TH2D("correzione_pbar","theta vs p corr |pbar; GeV/c; rad", 6, p_bins, 2, theta_bins); 
	h2[8] -> SetOption("lego");
	//Xsec histos
	h2[10]  = new TH2D("Xsec_pi","Xsec pi; GeV/c; rad", 6, p_bins, 2, theta_bins); 
	h2[10] -> SetOption("lego");
	h2[11]  = new TH2D("Xsec_k","Xsec k; GeV/c; rad", 6, p_bins, 2, theta_bins); 
	h2[11] -> SetOption("lego");
	h2[12]  = new TH2D("Xsec_pbar","Xsec pbar; GeV/c; rad", 6, p_bins, 2, theta_bins); 
	h2[12] -> SetOption("lego");

	h2[4]  = new TH2D("eta","eta vs p; GeV/c; eta", 100, 0, 50, 100, 0, 10);   
	h2[4] -> SetOption("colz");
	h2[5]  = new TH2D("XvsY","XY primary vertex; cm; cm", 100, -3, 3, 100, -3, 3);
	h2[5] -> SetOption("colz");

	tree = new TTree("USR1001","User Ntuple example");
	
	tree->Branch("mom", &mom, "mom/D");
	tree->Branch("Zvertex", &Zvertex, "Zvertex/D");
	tree->Branch("Yvertex", &Yvertex, "Yvertex/D");
	tree->Branch("Xvertex", &Xvertex, "Xvertex/D");
	tree->Branch("proton", &proton, "proton/I");
	tree->Branch("zfirst", &zfirst, "zfirst/D");
   	tree->Branch("zlast", &zlast, "zlast/D");
	tree->Branch("radL", &radL, "radL/D");
	tree->Branch("hits", &hits, "hits/I");
    tree->Branch("chi2", &chi2, "chi2/D");
	tree->Branch("theta", &theta, "theta/D");
	tree->Branch("Mtime", &Mtime, "Mtime/D");
	
	first=false;
	}
	
	//cout << "****************************\n" ;
	//cout << endl;
	cout << Nevent << endl;
	
	for(int i = 0; i < 12; i++){
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){
					cofatt[i][j][k] = 0.;
					err_cofatt[i][j][k] = 0.;
					inversa[i][j][k] = 0.;
					err_inv[i][j][k] = 0.;
				}
			
			err_N[i][j] = 0.;
			}
		determinant[i] = 0.;
		err_det_sq[i] = 0.;
		err_det[i] = 0.;
	}

		//calcolo determinanti (copiato brutalmente da internet)
	for(int j = 0; j < 12; j++){
		for(int i = 0; i < 3; i++){
			determinant[j] = determinant[j] + (efficienze[j][0][i] * (efficienze[j][1][(i+1)%3] * efficienze[j][2][(i+2)%3] - efficienze[j][1][(i+2)%3] * efficienze[j][2][(i+1)%3]));
		//errore determinante ^2
			err_det_sq[j] = err_det[j] + pow((errori[j][0][i] * (efficienze[j][1][(i+1)%3] * efficienze[j][2][(i+2)%3] - efficienze[j][1][(i+2)%3] * efficienze[j][2][(i+1)%3])),2);
		
		}
		err_det[j] = sqrt(err_det_sq[j]);
	}
	
	
	//calcolo matrici inverse (same as above)
		for(int i = 0; i < 12; i++){
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){
					inversa[i][j][k] = ((efficienze[i][(k+1)%3][(j+1)%3] * efficienze[i][(k+2)%3][(j+2)%3]) - (efficienze[i][(k+1)%3][(j+2)%3] * efficienze[i][(k+2)%3][(j+1)%3]))/ determinant[i];
				}
		}	
	}
	
	for(int i = 0; i < 12; i++){
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){


					//errore elementi cofatt ^2
					err_cofatt[i][j][k] = sqrt( pow((efficienze[i][(k+2)%3][(j+2)%3]* errori[i][k+1][j+1]),2.) + pow((efficienze[i][(k+2)%3][(j+1)%3]* errori[i][k+1][j+2]),2.) + pow((efficienze[i][(k+1)%3][(j+2)%3]*errori[i][k+2][j+1]),2.) + pow((efficienze[i][(k+1)%3][(j+1)%3]*errori[i][k+2][j+2]),2.) );
				}
			}
		}
		
	for(int i = 0; i < 12; i++){
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){
					
					//errore inversa
					err_inv[i][j][k] = sqrt( pow((err_cofatt[i][j][k]/determinant[i]),2.) + pow((cofatt[i][j][k]*err_det[i])/(pow(determinant[i],2)),2.) );
				}
			}	
		}
	
	
	pbar 				= 0;
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
	
	Nout = v.NOutParticles();
	contatorevertici ++;
	h1[16] -> Fill(0.);
	proton = -1;
	
		

	run = e.RunNum();
	n_rich = 1 + PaMetaDB::Ref().NminusOne(run)/1e6; //refraction index
    if(PaMetaDB::Ref().NminusOne(run) < 0) n_rich = 1 + PaSetup::Ref().NminusOne()/1e6;
	
	mom_p = M_p/sqrt((n_rich*n_rich) - 1);
	mom_k = M_k/sqrt((n_rich*n_rich) - 1);
	mom_pi = M_pi/sqrt((n_rich*n_rich) - 1);
	
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

		
		for(int ip = 0; ip < Nout; ip++){ //loop over outgoing particles from each vertex
		
		   int index = v.iOutParticle(ip);
		   const PaParticle& p = e.vParticle(index); 
		   if(p.iTrack() == -1) continue;
		   trackindex ++;
		   h1[16] -> Fill(1.);
		   
		   NoutPV += Nout;

		   const PaTrack& pt = e.vTrack(p.iTrack());
		   const PaTPar& par = pt.vTPar(0);

		   const PaSetup& setup = PaSetup::Ref();
		   const double z_rich = setup.Rich().DetPos(0).Z();
		   PaTPar Hout;
		   extrap = pt.Extrapolate(z_rich, Hout);  //extrapolates trajectory parameters at z = z_rich and Hout is the result


		   mom = p.ParInVtx(ivok).Mom();
		   zfirst = pt.ZFirst();
		   zlast = pt.ZLast();
		   radL = pt.XX0();
		   hits = pt.NHits();
		   chi2 = pt.Chi2tot();
		   q = p.Q();
		   Mtime = pt.MeanTime();
		   ndf = pt.Ndf();
		   Xrid = chi2/ndf;
		   theta = par.Theta();
		   eta = -log(tan(theta/2));
		   Xtraj = Hout.X();
		   Ytraj = Hout.Y();
		   double r = sqrt((Xtraj*Xtraj) + (Ytraj*Ytraj));
		   
		   tm = e.TrigMask();
		  
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


		  h1[1] -> Fill(mom);
		  h1[6]-> Fill(zfirst);
		  h1[7]-> Fill(zlast);
		  h1[8]-> Fill(radL);
		  h1[9]-> Fill(hits);
		  h1[12]-> Fill(chi2);
		  h1[13]-> Fill(Xtraj);
		  h1[14]-> Fill(Ytraj);
		  h1[2] -> Fill(Zvertex);
		  h1[3] -> Fill(Yvertex);
		  h1[4] -> Fill(Xvertex);
		  h1[15]-> Fill(X2);
		  h2[4] -> Fill(mom, eta);
		  h2[5] -> Fill(Xvertex, Yvertex);
		  
		  //***************CORREZIONE NUOVA EFFICIENZE********************
		  
		 
		  
		  proton = pid.IsProton(mom, mom_p, mom_k, mom_pi, like, count, q);


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
		    default: 
			cout << " caso default # event " << (Nevent - 1048580)/20 << " proton " << proton << endl;	
			break;}
			
		 tree -> Fill();

		if(proton >= 3) continue; 

		static int zone;
		int range;
		range = (theta > 0.03);  //if theta > 0.03, range = 1; = 0 otherwise
		
		if((mom >= 10 && mom <= 15) && range == 0) zone = 0; 		//bin00
		if((mom > 15 && mom <= 20) && range == 0) zone = 1;  		//bin01
		if((mom > 20 && mom <= 25) && range == 0) zone = 2;  		//bin02
		if((mom > 25 && mom <= 30) && range == 0) zone = 3;  		//bin03
		if((mom > 30 && mom <= 35) && range == 0) zone = 4;  		//bin04
		if((mom > 35 && mom <= 45) && range == 0) zone = 5;  		//bin05
		if((mom >= 10 && mom <= 15) && range == 1)zone = 6;     //bin10
		if((mom > 15 && mom <= 20) && range == 1)zone = 7;     //bin11
		if((mom > 20 && mom <= 25) && range == 1)zone = 8;     //bin12
		if((mom > 25 && mom <= 30) && range == 1)zone = 9;     //bin13
		if((mom > 30 && mom <= 35) && range == 1)zone = 10;    //bin14
		if((mom > 35 && mom <= 45) && range == 1)zone = 11;	   //bin15
		
		
		for(int i = 0; i < 3; i++){
			N[zone][i] = N[zone][i] + inversa[zone][i][proton];
			err_N[zone][i] = err_N[zone][i] + err_inv[zone][i][proton];
		}
			 
      //********************FINE CORREZIONE EFFICIENZE********************
		  


		}   //end particles loop
		
		h1[17] -> Fill(numPi);
		h1[18] -> Fill(numK);
		h1[19] -> Fill(numPbar);
	
    
		if(pbar > 0) Rs++;  //if at least 1 pbar is produced, count the event
	
	
	    if(Procevent + 1 == Totevent ){ 
		 
			float areaA = 0.02*5.; //area zona 1,2,3,4,5
			float areaB = 0.15*5.; //area zona 7,8,9,10,11

		  Ri = Totevent - Rqe;
		  
/*for(int i = 0; i < 12; i++){
		for(int j = 0; j < 3; j++){
			cout << " i:" << i << " j:" << j << endl << N[i][j] << endl;
		}
	cout << endl;
}*/

		//si filla l'isto per efficienze corrette
		
			for(int j = 0; j < 3; j++){
				
				//efficiency-corrected histos

				h2[6+j] -> SetBinContent(1,1,N[0][j]);
				h2[6+j] -> SetBinContent(2,1,N[1][j]);
				h2[6+j] -> SetBinContent(3,1,N[2][j]);
				h2[6+j] -> SetBinContent(4,1,N[3][j]);
				h2[6+j] -> SetBinContent(5,1,N[4][j]);
				h2[6+j] -> SetBinContent(6,1,N[5][j]);
				h2[6+j] -> SetBinContent(1,2,N[6][j]);
				h2[6+j] -> SetBinContent(2,2,N[7][j]);
				h2[6+j] -> SetBinContent(3,2,N[8][j]);
				h2[6+j] -> SetBinContent(4,2,N[9][j]);
				h2[6+j] -> SetBinContent(5,2,N[10][j]);
				h2[6+j] -> SetBinContent(6,2,N[11][j]);
		
				//normalized Xsection histos

				h2[10+j] -> SetBinContent(1,1,(N[0][j]/(Ri*areaA))*sigma_inel);
				h2[10+j] -> SetBinContent(2,1,(N[1][j]/(Ri*areaA))*sigma_inel);
				h2[10+j] -> SetBinContent(3,1,(N[2][j]/(Ri*areaA))*sigma_inel);
				h2[10+j] -> SetBinContent(4,1,(N[3][j]/(Ri*areaA))*sigma_inel);
				h2[10+j] -> SetBinContent(5,1,(N[4][j]/(Ri*areaA))*sigma_inel);
				h2[10+j] -> SetBinContent(6,1,(N[5][j]/(Ri*areaA*2))*sigma_inel);
				h2[10+j] -> SetBinContent(1,2,(N[6][j]/(Ri*areaB))*sigma_inel);
				h2[10+j] -> SetBinContent(2,2,(N[7][j]/(Ri*areaB))*sigma_inel);
				h2[10+j] -> SetBinContent(3,2,(N[8][j]/(Ri*areaB))*sigma_inel);
				h2[10+j] -> SetBinContent(4,2,(N[9][j]/(Ri*areaB))*sigma_inel);
				h2[10+j] -> SetBinContent(5,2,(N[10][j]/(Ri*areaB))*sigma_inel);
				h2[10+j] -> SetBinContent(6,2,(N[11][j]/(Ri*areaB*2))*sigma_inel);
				
				//Xsec projections on x axis

				h1[20+j] -> SetBinContent(1,((N[0][j]+N[6][j])/(Ri*areaA))*sigma_inel);
				h1[20+j] -> SetBinContent(2,((N[1][j]+N[7][j])/(Ri*areaA))*sigma_inel);
				h1[20+j] -> SetBinContent(3,((N[2][j]+N[8][j])/(Ri*areaA))*sigma_inel);
				h1[20+j] -> SetBinContent(4,((N[3][j]+N[9][j])/(Ri*areaA))*sigma_inel);
				h1[20+j] -> SetBinContent(5,((N[4][j]+N[10][j])/(Ri*areaA))*sigma_inel);
				h1[20+j] -> SetBinContent(6,((N[5][j]+N[11][j])/(Ri*areaA*4))*sigma_inel);

			//Xsec errors

				h1[20+j] -> SetBinError(1,(1/(Ri*areaA)*sqrt( pow(((N[0][j]+N[6][j])*err_sig),2) + pow((sigma_inel*(err_N[0][j]+err_N[6][j])),2) )));
				h1[20+j] -> SetBinError(2,(1/(Ri*areaA)*sqrt( pow(((N[1][j]+N[7][j])*err_sig),2) + pow((sigma_inel*(err_N[1][j]+err_N[7][j])),2) )));
				h1[20+j] -> SetBinError(3,(1/(Ri*areaA)*sqrt( pow(((N[2][j]+N[8][j])*err_sig),2) + pow((sigma_inel*(err_N[2][j]+err_N[8][j])),2) )));
				h1[20+j] -> SetBinError(4,(1/(Ri*areaA)*sqrt( pow(((N[3][j]+N[9][j])*err_sig),2) + pow((sigma_inel*(err_N[3][j]+err_N[9][j])),2) )));
				h1[20+j] -> SetBinError(5,(1/(Ri*areaA)*sqrt( pow(((N[4][j]+N[10][j])*err_sig),2) + pow((sigma_inel*(err_N[4][j]+err_N[10][j])),2) )));
				h1[20+j] -> SetBinError(6,(1/(Ri*areaA*4)*sqrt( pow(((N[5][j]+N[11][j])*err_sig),2) + pow((sigma_inel*(err_N[5][j]+err_N[11][j])),2) )));

		 }	
		

		 
		 cout << "Rs " << Rs << " Ri " << Ri << " Rqe " << Rqe << " RBK1 " << RBK1 << " eventi validi " << contatoreEventi << endl; 
		 //cout << endl << "contatore vertici " << contatorevertici << endl << " # out particles from PV " << NoutPV << endl << " # con indice traccia definito " << trackindex << endl << " z first < z SM1 " << zfirstOK << endl << " mean time definito " << MtimeOK << endl << " z vertex ok " << zvertexOK << endl << " chi ridotto < 10 " << XridOK << endl;
		 cout << endl << "p-pbar counter: 0 = pion | 1 = kaon | 2 = pbar | 3 = neither | 4 = proton " << endl;
		 cout << " incoming protons " << countproton << " # out A-protons " << countAproton << " # out protons " << countproton1 << endl;}


}		