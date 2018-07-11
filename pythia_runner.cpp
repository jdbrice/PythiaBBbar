
#include "TPythia6.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"

// STL
#include <iostream>
#include <string> 
#include <fstream>
#include <map>

// #include <cstdio>
// #include <iostream>
// #include <fstream>
// #include <vector>

// #include "TMath.h"
// #include "TString.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TFile.h"
// #include "TCanvas.h"
// #include "TStyle.h"
// #include "TVector3.h"
// #include "TTree.h"
// #include "TGraphErrors.h"
// #include "TF1.h"
// #include "TH2.h"
// #include "TLorentzVector.h"
// #include "TH1.h"
// #include "TList.h"
// #include "TMultiGraph.h"

// 
// #include "TParticle.h"
// #include "TParticlePDG.h"


#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1I.h"

using namespace std;


TPythia6 		*pythia = NULL;
TTree 			*tree = NULL;
TClonesArray 	*pyTracks = NULL;
TFile 			*tFile = nullptr;
TNtuple 		*ntuple = nullptr;
int 			nParticles = 0;

#define PID_b 5
#define PID_e 11
#define PID_mu 13
#define PID_string 92

#define MASS_mu 0.1056583745 // GeV/c^2

#define K_ID 2
#define K_PARENT 3


TH1 * hParentId = nullptr,
	*hPt = nullptr,
	*hEta = nullptr,
	*hPhi = nullptr,
	*hPairPt = nullptr,
	*hPairEta = nullptr,
	*hPairY = nullptr,
	*hPairPhi = nullptr; 

TH2 *hFullAcc_dNdM_pT = nullptr, // full phase space
	*hPairCut_dNdM_pT = nullptr, // |y| < 0.5
	*hAccCut0_dNdM_pT = nullptr, // +|\eta| < 0.5
	*hAccCut1_dNdM_pT = nullptr; // +p_T > 1.1

// Electrons
TH1 *heParentId = nullptr,
	*hePt = nullptr,
	*heEta = nullptr,
	*hePhi = nullptr,
	*hePairPt = nullptr,
	*hePairEta = nullptr,
	*hePairY = nullptr,
	*hePairPhi = nullptr; 

TH2 *heFullAcc_dNdM_pT = nullptr, // full phase space
	*hePairCut_dNdM_pT = nullptr, // |y| < 1.0
	*heAccCut0_dNdM_pT = nullptr, // +|\eta| < 1.0
	*heAccCut1_dNdM_pT = nullptr; // +p_T > 0.2


// NOTE
// I did not turn off any decay channels
// so I don't need to do the cross section weighting like in cc-bar
// I decided not to do this because with bb-bar you need to do 2 levels of 
// br scaling, and it is not worth it.


// map<size_t, double> branchingRatios = {
// 	{411, 0.176},
// 	{421,0.067},
// 	{431,0.065},
// 	{4122,0.045},
// 	{511, 0.105}, // B0->mu + anything, from pythia table
// 	{}
// };


unsigned long long int get_seed(){
	unsigned long long int random_value = 0; //Declare value to store data into
	size_t size = sizeof(random_value); //Declare size of data
	ifstream urandom("/dev/urandom", ios::in|ios::binary); //Open stream
	if(urandom) //Check if stream is open
	{
		urandom.read(reinterpret_cast<char*>(&random_value), size); //Read from urandom
		if(urandom) {
			return random_value;
		}
		else { //Read failed
			return 0;
		}
		urandom.close(); //close stream
	} else { //Open failed
		std::cerr << "Failed to open /dev/urandom" << std::endl;
	}
	return 0;
}


void setupPythia( int trigger, long int seed = 0 ){
	pythia = new TPythia6();

	const int MSEL_MINBIAS = 1;
	const int MSEL_BBBAR_TRIG = 5;
	if ( MSEL_MINBIAS != trigger && MSEL_BBBAR_TRIG != trigger ){
		cout << "WARNING: trigger must be (1=mb) or (5=bbbar), got " << trigger << endl;
	}

	pythia->SetMSEL(trigger);

	// Tune from JOEY BUTTERWORTH
	pythia->SetPARP(91,1.0); //<kT>
	pythia->SetPARP(67,1.0);  //mstp*4

	
	unsigned long long int gSeed = seed;
	if ( 0 == seed ){
		gSeed = get_seed();
		cout << "Generating random seed from /dev/urandom = " << gSeed << endl;
	} else {
		cout << "Using user seed = " << gSeed << endl;
	}

	pythia->SetMRPY(1, gSeed);

	// pythia->SetMRPY(1, iseed); 
	pythia->Initialize("CMS","p","p",200);

	pythia->Pylist(12);
	pythia->Pystat(1);
	pythia->Pystat(4);
	pythia->Pystat(5);


	string name = "./pythia_bbbar_TRIG_" + to_string( trigger ) + "_seed_" + to_string( gSeed ) +".root" ;
	tFile = new TFile( name.c_str(), "RECREATE" );
	tFile->cd();

	hParentId = new TH1I( "hParentId", "", 1000, 0, 1000 );
	hPt  = new TH1D( "hPt", "", 150, 0, 15 );
	hEta = new TH1D( "hEta", "", 500, -5, 5 );
	hPhi = new TH1D( "hPhi", "", 320, -3.2, 3.2 );
	hFullAcc_dNdM_pT = new TH2D( "hFullAcc_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );
	hPairCut_dNdM_pT = new TH2D( "hPairCut_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );
	hAccCut0_dNdM_pT = new TH2D( "hAccCut0_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );
	hAccCut1_dNdM_pT = new TH2D( "hAccCut1_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );

	hPairPt   = new TH1D( "hPairPt", "", 150, 0, 15 );
	hPairEta  = new TH1D( "hPairEta", "", 500, -5, 5 );
	hPairY    = new TH1D( "hPairY", "", 500, -5, 5 );
	hPairPhi  = new TH1D( "hPairPhi", "", 320, -3.2, 3.2 );


	heParentId = new TH1I( "heParentId", "", 1000, 0, 1000 );
	hePt  = new TH1D( "hePt", "", 150, 0, 15 );
	heEta = new TH1D( "heEta", "", 500, -5, 5 );
	hePhi = new TH1D( "hePhi", "", 320, -3.2, 3.2 );
	heFullAcc_dNdM_pT = new TH2D( "heFullAcc_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );
	hePairCut_dNdM_pT = new TH2D( "hePairCut_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );
	heAccCut0_dNdM_pT = new TH2D( "heAccCut0_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );
	heAccCut1_dNdM_pT = new TH2D( "heAccCut1_dNdM_pT", ";M (GeV/c^2); p_{T} (GeV/c)", 1000, 0, 10, 1500, 0, 15 );

	hePairPt   = new TH1D( "hePairPt", "", 150, 0, 15 );
	hePairEta  = new TH1D( "hePairEta", "", 500, -5, 5 );
	hePairY    = new TH1D( "hePairY", "", 500, -5, 5 );
	hePairPhi  = new TH1D( "hePairPhi", "", 320, -3.2, 3.2 );


	// ntuple = new TNtuple("bbbar", "bbbar to mumu","nPt:nEta:nPhi:nMass:nParentId:pPt:pEta:pPhi:pMass:pParentId:pairPt:pairEta:pairPhi:pairMass:pairY");
}

void printPlc( int i ){
	cout << "K(1)=" << pythia->GetK(i, 1) << ", K(2)=" << pythia->GetK(i, 2) << ", K(3)=" << pythia->GetK(i, 3) << ", K(4)=" << pythia->GetK(i, 4) << ", K(5)=" << pythia->GetK(i, 5) << endl;
}

int state( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 1 );
	return -1;
}
int plcId( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 2 );
	return -1;
}
int parentIndex( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 3 );
	return -1;
}
int posX( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 4 );
	return -1;
}
int posY( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 5 );
	return -1;
}
TLorentzVector lvec( int i ){
	TLorentzVector lv;
	if ( i>0 && i<= nParticles ){
		lv.SetPxPyPzE(  pythia->GetP(i,1),
						pythia->GetP(i,2),
						pythia->GetP(i,3),
						pythia->GetP(i,4));
	}
	return lv;
}


/* findStrings
 * counts the number of strings in the event
 * for ccbar we want exactly 2 strings to be consistent with dielectron papers
 */
int findStrings(){
	int nStrings = 0;
	for(Int_t i = 0; i < nParticles; i++){
		int id = abs(pythia->GetK(i+1, K_ID ) );
		if(id == PID_string ){  
			int parentIdIndex = pythia->GetK( i+1, K_PARENT);
			
			if(abs(pythia->GetK(parentIdIndex,2)) == PID_b){
				nStrings++;
			}//charm

		}//string
	}
	return nStrings;
}


/* isMuon
 * checks particle i to see if it is decayed muon from some bottom or charm meson
 */
bool isMuon( int i ){

	if ( abs( plcId( i ) ) != PID_mu )
		return false;
	

	if ( (posX(i) != 0 || posY(i) != 0) && state( i ) != 1 ){
		cout << "ERROR state pos mismatch" << endl;
		return false;
	}

	int pIndex = parentIndex( i );
	if ( pIndex <= 0 )
		return false;
	// {
	// 	// printPlc( pIndex );

		int parentId = abs(plcId( pIndex ));
		// reject secondaries from b decays, like phi, eta, rho, pi0 etc.
		if ( abs(parentId) < 400 ) return false;



	// 	if(parentId==411||parentId==421||parentId==431||parentId==4122){
	// 		// printPlc( i );
	// 		// cout << "\tPARENT: "; printPlc( pIndex );
	// 		return true;
	// 	}

	// } else {
	// 	cout << "\t reject no parent" << endl;
	// }
	



	return true;
}


/* isElectron
 * checks particle i to see if it is decayed electron from some bottom or charm meson
 */
bool isElectron( int i ){

	if ( abs( plcId( i ) ) != PID_e )
		return false;
	
	if ( (posX(i) != 0 || posY(i) != 0) && state( i ) != 1 ){
		cout << "ERROR state pos mismatch" << endl;
		return false;
	}

	int pIndex = parentIndex( i );
	if ( pIndex <= 0 )
		return false;

	int parentId = abs(plcId( pIndex ));
	// reject secondaries from b decays, like phi, eta, rho, pi0 etc.
	if ( abs(parentId) < 400 ) 
		return false;

	return true;
}

void findMuons(){
	// cout << " START EVENT =============" << endl;
	TLorentzVector plv, nlv, lv;
	TLorentzVector plve, nlve, lve;
	int pParentId = -1, nParentId = -1;
	int peParentId = -1, neParentId = -1;
	bool foundPos = false; 
	bool foundNeg = false;

	bool foundPosElec = false; 
	bool foundNegElec = false;

	for(Int_t i = 0; i < nParticles; i++){
		int iPlc = i+1;
		bool isMu = isMuon( iPlc );
		bool isE = isElectron( iPlc );
		int pId = plcId( iPlc );
		int ppId = plcId( parentIndex( iPlc ) );
		int pppId = plcId( parentIndex( parentIndex( iPlc ) ) );

		// REJECT decays of non-prompt J/Psi
		if ( 443 == abs(ppId) ) {
			continue;
		}


		if (  isMu && pId == PID_mu ){
			// set pos
			nlv = lvec( iPlc );
			nParentId = plcId( parentIndex( iPlc ) );
			foundNeg = true;
		} else if ( isMu && pId == -PID_mu ){
			plv = lvec( iPlc );
			pParentId = plcId( parentIndex( iPlc ) );
			foundPos = true;
		}

		if (  isE && pId == PID_e ){
			// set pos
			nlve = lvec( iPlc );
			neParentId = plcId( parentIndex( iPlc ) );
			foundNegElec = true;
		} else if ( isE && pId == -PID_e ){
			plve = lvec( iPlc );
			peParentId = plcId( parentIndex( iPlc ) );
			foundPosElec = true;
		}

		if ( foundPos && foundNeg && foundPosElec && foundNegElec) {
			break;
		}
	}


	// Muons
	if ( foundPos && foundNeg ){
		lv = plv + nlv;

		hParentId->Fill( abs(pParentId) );
		hParentId->Fill( abs(nParentId) );

		hPt->Fill( plv.Pt() );
		hPt->Fill( nlv.Pt() );

		hEta->Fill( plv.Eta() );
		hEta->Fill( nlv.Eta() );

		hPhi->Fill( plv.Phi() );
		hPhi->Fill( nlv.Phi() );


		hPairPt->Fill( lv.Pt() );
		hPairEta->Fill( lv.Eta() );
		hPairY->Fill( lv.Rapidity() );
		hPairPhi->Fill( lv.Phi() );

		hFullAcc_dNdM_pT->Fill( lv.M(), lv.Pt() );
		if ( abs(lv.Rapidity()) < 0.5 ){
			hPairCut_dNdM_pT->Fill( lv.M(), lv.Pt() );

			if ( abs(plv.Eta()) < 0.5 && abs(nlv.Eta()) < 0.5 ){
				hAccCut0_dNdM_pT->Fill( lv.M(), lv.Pt() );
				if ( plv.Pt() > 1.1 && nlv.Pt() > 1.1 ){
					hAccCut1_dNdM_pT->Fill( lv.M(), lv.Pt() );
				} // pT > 1.1
			} // |eta| < 0.5
		} // |y| < 0.5

		// ntuple->Fill( data );
	} else {
		// cout << "COULD NOT find mu+mu-" << endl;
	}


	// Electrons
	if ( foundPosElec && foundNegElec ){
		lve = plve + nlve;

		heParentId->Fill( abs(peParentId) );
		heParentId->Fill( abs(neParentId) );

		hePt->Fill( plve.Pt() );
		hePt->Fill( nlve.Pt() );

		heEta->Fill( plve.Eta() );
		heEta->Fill( nlve.Eta() );

		hePhi->Fill( plve.Phi() );
		hePhi->Fill( nlve.Phi() );


		hePairPt->Fill( lve.Pt() );
		hePairEta->Fill( lve.Eta() );
		hePairY->Fill( lve.Rapidity() );
		hePairPhi->Fill( lve.Phi() );

		heFullAcc_dNdM_pT->Fill( lve.M(), lve.Pt() );
		if ( abs(lve.Rapidity()) < 1.0 ){
			hePairCut_dNdM_pT->Fill( lve.M(), lve.Pt() );

			if ( abs(plve.Eta()) < 1.0 && abs(nlve.Eta()) < 1.0 ){
				heAccCut0_dNdM_pT->Fill( lve.M(), lve.Pt() );
				if ( plve.Pt() > 0.2 && nlve.Pt() > 0.2 ){
					heAccCut1_dNdM_pT->Fill( lve.M(), lve.Pt() );
				} // pT > 1.1
			} // |eta| < 0.5
		} // |y| < 0.5

		// ntuple->Fill( data );
	} else {
		// cout << "COULD NOT find mu+mu-" << endl;
	}

	// cout << " ============= STOP EVENT" << endl;
}


void genEvents( ULong_t _nEvents ){

	ULong_t iEvent = 0;
	while ( iEvent < _nEvents ){

		pythia->GenerateEvent();
		nParticles = pythia->GetNumberOfParticles();

		int nStrings = findStrings();
		// cout << "Found " << nStrings << " strings" << endl;
		if ( 2 == nStrings )
			findMuons( );

		if ( iEvent % 1000 == 0 )
			cout << "." << std::flush;

		iEvent++;
	}
	cout << endl;

}


// NOTE: at first i thought it looked wrong to assume only 1 ccbar -> mu mu but it is so rare that this is fine

Int_t main( Int_t argc, Char_t **argv){
	cout << "argc = " << argc << endl;
	if ( argc < 3 ) {
		cout  << "USAGE\n GENERATOR trig(1=mb,5=bbbar) nEvents rndSeed=0" << endl;
		return 1;
	}

	int trigger = atoi( argv[1] );
	int nEvents = atoi( argv[2] );
	
	long int seed = 0;
	
	if (argc >= 4)
		seed = atol( argv[3] );

	setupPythia( trigger, seed );
	genEvents( nEvents );

	cout<<"::::::::Printing out stats:::::::"<<endl;
	pythia->Pystat(1);
	pythia->Pystat(4);
	pythia->Pystat(5);

	cout << "SAVE THEN END" << endl;
	// ntuple->Write();
	tFile->Write();
	tFile->Close();
	delete tFile;
	cout << "END" << endl;



	return 0;
}
