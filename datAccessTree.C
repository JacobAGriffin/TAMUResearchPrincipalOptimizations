#define datAccessTree_cxx
#include "datAccessTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TH1F.h>
#include <TString.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <cmath>
#include <TEfficiency.h>
#include <TF1.h>

vector<int> MatchLevel(vector<int> array2){
	int array[6][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}};
	array[0][0] = array2[0];
	array[0][1] = 1;
	for (int i = 1; i < 6; ++i) {
		for (int j = 0; j < 5; ++j) {
			if(array2[i] == array[j][0]) {
				++array[j][1];
				break;
			} else if(array[j+1][1] == 0) {
				array[j+1][0] = array2[i];
				array[j+1][1] = 1;
				break;
			}
		}
	}
	int num = 0;
	int pos = 0;
	for (int i = 0; i < 6; ++i) {
		if (array[i][1] > num && array[i][0] != -1) {
			pos = i;
			num = array[i][1];
		}
	}
	vector<int> final;
	final.push_back(num);
	final.push_back(array[pos][0]);
	return final;
}

vector<vector<long double> > compare(vector<vector<long double> > array, vector<vector<long double> > array2) {
	vector<vector<long double> > outArray;
	if (array.size() == 0) {
		outArray = array2;
	} else {
		outArray = array;
		for (unsigned i = 0; i < array2.size(); ++i) {
			for (unsigned j = 0; j < outArray.size(); ++j) {
				if(outArray[j][0] == array2[i][0] && outArray[j][1] == array2[i][1]){
					break;
				} else if(j + 1 == outArray.size()) {
					outArray.push_back(array2.at(i));
					break;
				}
			}
		}
	}
	return outArray;
}

vector<vector<long double> > singleCompare(vector<vector<long double> > array, vector<vector<long double> > array2) {
	vector<vector<long double> > outArray;
	for (unsigned i = 0; i < array2.size(); ++i) {
		for (unsigned j = 0; j < array.size(); ++j) {
			if(array[j][0] == array2[i][0] && array[j][1] == array2[i][1]) {
				outArray.push_back(array2.at(i));
				break;
			}
		}
	}
	return outArray;
}

int boolArrayToInt(vector<bool> temp){
	int count = 0;
	for (int i = 0; i < 6; ++i) {
		if (temp[i])
			++count;
	}
	return count;
}

vector<bool> radiusStubTPidChecker(vector<float> *tempf,vector<int> *tempi) {
	vector<bool> tempb = {0,0,0,0,0,0};
	for (int i = 0; i < tempf->size(); ++i) {
		if (tempi->at(i) != 0)
			continue;
		if (tempf->at(i) < 20)
			continue;
		if (tempf->at(i) >= 20 && tempf->at(i) <= 25)
			tempb[0] = true;
		else if (tempf->at(i) >= 33 && tempf->at(i) <= 38)
			tempb[1] = true;
		else if (tempf->at(i) >= 48 && tempf->at(i) <= 53)
			tempb[2] = true;
		else if (tempf->at(i) >= 66 && tempf->at(i) <= 71)
			tempb[3] = true;
		else if (tempf->at(i) >= 86 && tempf->at(i) <= 91)
			tempb[4] = true;
		else if (tempf->at(i) >= 106 && tempf->at(i) <= 110)
			tempb[5] = true;
	}
	return tempb;
}

double genTrackTransverse(const double &pt, const double &phi0, const double &d0, const int charge, const double &R) {
	double rho = charge*pt/(3.8114*0.003);
	double phiGen = phi0 - asin((d0*d0 + 2*d0*rho + R*R)/(2*R*(rho+d0)));
	return phiGen;
}


double genTrackLongitudinal(const double &z0, const double &cotTheta, const double &pt, const double &d0, const int charge, const double &R) {
  //std::cout<<"Longitudinal input: "<<z0<<", "<<cotTheta<<", "<<pt<<", "<<d0<<", "<<charge<<", "<<B<<", "<<R<<", "<<z<<std::endl;
	double rho = charge*pt/(3.8114*0.003);
	double zGen = z0 + charge*rho*cotTheta*acos((pow(rho,2) + pow(d0+rho,2) - pow(R,2))/(2*rho*(rho+d0)));
	return zGen;
}

bool Containment_tt27Test(const double &pt, const double &phi0, const double &z0, const double &eta, const int charge, const double &vx, const double &vy) {
	const double d0 = vx*sin(phi0)-vy*cos(phi0);
	const double cotTheta = 1./tan(2.*std::atan(std::exp(-1*eta)));
	float Rmin[6]={20.,33.,48.,66.,86.,106.};
	float Rmax[6]={25.,38.,53.,71.,91.,110.};
	float PhiMin[6]={0.564430, 0.653554, 0.641981, 0.717273, 0.658179, 0.618448};
	float PhiMax[6]={1.791765, 1.710419, 1.756567, 1.638922, 1.673851, 1.778293};
	float zMin[6]={-6.7127, -6.7797, -5.2542, -9.5318, -9.5318, -9.5318};
	float zMax[6]={26.9799, 36.7048, 47.7511, 59.4103, 78.7372, 88.9935};
	bool test=true;
	for(unsigned i=0; i<6; ++i){
		const double phi1=genTrackTransverse(pt, phi0, d0, charge, Rmin[i]);
		const double phi2=genTrackTransverse(pt, phi0, d0, charge, Rmax[i]);
		const double z1=genTrackLongitudinal(z0, cotTheta, pt, d0, charge, Rmin[i]);
		const double z2=genTrackLongitudinal(z0, cotTheta, pt, d0, charge, Rmax[i]);
		if(phi1>PhiMax[i] || phi1<PhiMin[i] || phi2>PhiMax[i] || phi2<PhiMin[i] || z1>zMax[i] || z1<zMin[i] || z2>zMax[i] || z2<zMin[i]) test=false;
		if(!test) break;
	}
	return test;
}

void datAccessTree::loop() {
		//   In a ROOT session, you can do:
		//      root> .L datAccessTree.C
		//      root> datAccessTree t
		//      root> t.GetEntry(12); // Fill t data members with entry number 12
		//      root> t.Show();       // Show values of entry 12
		//      root> t.Show(16);     // Read and show values of entry 16
		//      root> t.loop();       // Loop on all entries
		//	
		//    This is the loop skeleton where:
		//    jentry is the global entry number in the chain
		//    ientry is the entry number in the current Tree
		//    Note that the argument to GetEntry must be:
		//    jentry for TChain::GetEntry
		//    ientry for TTree::GetEntry and TBranch::GetEntry
		//
		//    To read only selected branches, Insert statements like:
		//    METHOD1:
		//    fChain->SetBranchStatus("*",0);  // disable all branches
		//    fChain->SetBranchStatus("branchname",1);  // activate branchname
		//    METHOD2: replace line
		//    fChain->GetEntry(jentry);       //read all branches
		//    by  b_branchname->GetEntry(ientry); //read only this branch
	
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	Long64_t nbytes = 0, nb = 0;
	
	TString type = "Muon";
	//TString type = "Pion";
	//TString type = "Electron";
	//TString type = "PileUp";
	//TString type = "T~T";
	TString fsss="6x6";
	//TString fsss="5x6";
	//TString fsss="5x5";
	TString NPP="PU";
	//TString NPP="NPU";
	bool TakeEverything = false, bMM = false;
	if (NPP == "PU") {
		TakeEverything = true;
	}
	if (type == "T~T" || type == "PileUP") {
		bMM = true;
	}


	int num = 12;
	if (fsss == "6x6") {
		num = 12;
	}
	
	vector<vector<vector<long double> > > gsubPrincipal;
	vector<vector<vector<long double> > > bsubPrincipal;
	vector<vector<long double> > midstep;
	for(int i=0; i<8; ++i){
		gsubPrincipal.push_back(midstep);
		bsubPrincipal.push_back(midstep);
	}
	
		//histogram definition section
	vector<TH1F*> h1_Principals_Good, h1_Principals_Bad, h1_Principals_Gooda, h1_Principals_Bada, h1_Principals_VBada;
	TH1F *test=new TH1F("","",5000,0,1000);
	TH1F *testa=new TH1F("","",400,-20,20);
	TH1F *histoCounter = new TH1F("Counter","Counter",7,-0.5,6.5);
	TH1F *histoCounter2 = new TH1F("Counter2","Counter2",7,-0.5,6.5);
	TH1F *histoTrackCandidates = new TH1F("Track Candidates","Track Candidates",61,-0.5,60.5);
	for(unsigned i=0; i<num; ++i){
		h1_Principals_Good.push_back((TH1F*)test->Clone(TString::Format("goodPrincipal%1.0d",i)));
		h1_Principals_Bad.push_back((TH1F*)test->Clone(TString::Format("badPrincipal%1.0d",i)));
		h1_Principals_Gooda.push_back((TH1F*)testa->Clone(TString::Format("goodPrincipal%1.0d",i)));
		h1_Principals_Bada.push_back((TH1F*)testa->Clone(TString::Format("badPrincipal%1.0d",i)));
		h1_Principals_VBada.push_back((TH1F*)testa->Clone(TString::Format("vbadPrincipal%1.0d",i)));
	}
	double bins[21]={1.,2.,3.,4.,5.,6.,7.,8.,10.,14.,18.,25.,40.,60.,100.,150.,250.,400.,800.,1500.,5000.};
	TEfficiency *PTDifT = new TEfficiency("PTDifT","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDifTChi = new TEfficiency("PTDifTChi","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif = new TEfficiency("PTDif","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_6 = new TEfficiency("PTDif_6","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_5P10 = new TEfficiency("PTDif_5P10","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_5P12 = new TEfficiency("PTDif_5P12","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_4 = new TEfficiency("PTDif_4","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_G4T = new TEfficiency("PTDif_G4T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_G5T = new TEfficiency("PTDif_G5T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_G6T = new TEfficiency("PTDif_G6T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_GnTT = new TEfficiency("PTDif_GnTT","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_Gn1T = new TEfficiency("PTDif_Gn1T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_Gn2T = new TEfficiency("PTDif_Gn2T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_Gn3T = new TEfficiency("PTDif_Gn3T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_Gn4T = new TEfficiency("PTDif_Gn4T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_Gn5T = new TEfficiency("PTDif_Gn5T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_Gn6T = new TEfficiency("PTDif_Gn6T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_AMGn4T = new TEfficiency("PTDif_AMGn4T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_AMGn5T = new TEfficiency("PTDif_AMGn5T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_AMGn6T = new TEfficiency("PTDif_AMGn6T","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTDif_AMGnnTT = new TEfficiency("PTDif_AMGnnTT","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PTChiDif = new TEfficiency("PTChiDif","PT Differential;p_{T};#epsilon",20,bins);
	TEfficiency *PhiDif_Gn1T = new TEfficiency("PhiDif_Gn1T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_Gn2T = new TEfficiency("PhiDif_Gn2T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_Gn3T = new TEfficiency("PhiDif_Gn3T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_Gn4T = new TEfficiency("PhiDif_Gn4T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_Gn5T = new TEfficiency("PhiDif_Gn5T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_Gn6T = new TEfficiency("PhiDif_Gn6T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_AMGn4T = new TEfficiency("PhiDif_AMGn4T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_AMGn5T = new TEfficiency("PhiDif_AMGn5T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *PhiDif_AMGn6T = new TEfficiency("PhiDif_AMGn6T","Phi Differential;#phi;#epsilon",100,.4,1.8);
	TEfficiency *Mismatch = new TEfficiency("Mismatch","Mismatch;# of Mismatched Stubs;#epsilon",7,0,7);
	
		//Other Histograms
	double maxchi2 = 0.;
	TH1F *chi2_good = new TH1F("Chi2_Good","Chi2_Good",10000,0,1000);
	TH1F *chi2_notgood = new TH1F("Chi2_NotGood","Chi2_NotGood",10000,0,1000);
	TH1F *chi2_bad = new TH1F("Chi2_Bad","Chi2_Bad",10000,0,1000);
	TH1F *chi2_verybad = new TH1F("Chi2_VeryBad","Chi2_VeryBad",10000,0,1000);
	TH1F *pt_good = new TH1F("","",800,0,200);
	TH1F *pt_bad = new TH1F("","",800,0,200);
	TH1F *phi_good = new TH1F("","",200,0.7,1.6);
	TH1F *phi_bad = new TH1F("","",200,0.7,1.6);
	TH1F *eta_good = new TH1F("","",200,-.1,.9);
	TH1F *eta_bad = new TH1F("","",200,-.1,.9);
	TH1F *goodCounter = new TH1F("goodCounter","",100,0,100);
	TH1F *badCounter = new TH1F("badCounter","",100,0,100);
	TH1F *eventCounter = new TH1F("EventCounter","",100,0,100);
	TH1F *hMismatch = new TH1F("hMismatch","Mismatch;# of Mismatched Stubs;percent",7,0,7);
	unsigned matchn = nentries;
	int errorCounter = 0;
	int AMDenom = 0;
	
	int GoodCounter = 0;
	int BadCounter = 0;
	int TGTN = 0;
	int TBTN = 0;
	vector<vector<long double> > gmasterlist, bmasterlist;
	
	std::cout<<type<<" "<<NPP<<std::endl;
	
	double cp0 = 100000,
	cp1 = 33,
	cp2 = 22,
	cp3 = 12,
	cp6 = 100000,
	cp7 = 3,
	cp8 = 3,
	cp9 = 3;

	//nentries = 100000;

		//event loop
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		int GTN = 0;
		int BTN = 0;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) {
			break;
		}
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		bool contained = Containment_tt27Test(trkParts_pt->at(0), trkParts_phi->at(0), trkParts_vz->at(0), trkParts_eta->at(0), trkParts_charge->at(0), trkParts_vx->at(0), trkParts_vy->at(0));
		if(trkParts_phi->at(0)<.8 || trkParts_phi->at(0)>1.55)
			contained=false; //special provision because of bank limitations
		if (type != "PileUp" && type != "T~T")	
			if (!contained)
				continue;
		pair<int,bool> MaxMatch = {0,false};
		pair<int,bool> MaxMatch2 = {0,false};
		for (double i = 0; i < AMTTTracks_principals->size(); ++i) {
			histoTrackCandidates->Fill(AMTTTracks_ndof->at(i));
			int counter = 0;
			int counter2 = 0;
			int previous = -2;
			bool match = true;
			bool noDummy = true;
			bool cut = false;
			vector<int> array2;
			for (unsigned j = 0; j < AMTTTracks_stubRefs->at(i).size(); ++j) {
				const int stub = AMTTTracks_stubRefs->at(i)[j];
				if (stub != 999999999) {
					int tp = TTStubs_tpId->at(stub);
					if (previous == -2)
						previous=tp;
					else if (tp == -1 || previous != tp)
						match = false;
					array2.push_back(tp);
				} else {
					array2.push_back(-1);
					if (AMTTTracks_ndof->at(i) != 10)
						++errorCounter;
					noDummy=false;
				}
			}
			if(noDummy && AMTTTracks_ndof->at(i) != 12)
				++errorCounter;
			if (AMTTTracks_ndof->at(i) == 12) {
				fsss = "6x6";
			} else if (AMTTTracks_ndof->at(i) == 10) {
				fsss = "5x6";
			} else std::cout<<"ERROR"<<std::endl;
			histoCounter->Fill(counter);
			histoCounter2->Fill(counter2);
			vector<int> temptemp = MatchLevel(array2);
			/*std::cout<<temptemp[1]<<" : "<<trkParts_pt->size()<<" : "<<TTStubs_r->size()<<" : "<<AMTTTracks_principals->size()<<std::endl;
			if (TTStubs_r->size() == 7)
				for (int r = 0; r<TTStubs_r->size();++r)
					std::cout<<TTStubs_r->at(r)<<std::endl;*/
			/*if (temptemp[1] != -1)
				if (trkParts_pt->at(temptemp[1]) < 3)
					continue;*/
			if (match && (TakeEverything || previous == 0)) {
				for (double j = 0; j < AMTTTracks_principals->at(i).size(); j++) {
					h1_Principals_Good[j]->Fill(fabs(AMTTTracks_principals->at(i)[j]));
					h1_Principals_Gooda[j]->Fill(AMTTTracks_principals->at(i)[j]);
					vector<long double> temp = {static_cast<long double>(jentry), i, fabs(AMTTTracks_principals->at(i)[j])};
					if (fsss == "6x6"){
						if (j == 0 && fabs(AMTTTracks_principals->at(i)[j]) > cp0) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 1 && fabs(AMTTTracks_principals->at(i)[j]) > cp1) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 2 && fabs(AMTTTracks_principals->at(i)[j]) > cp2) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 3 && fabs(AMTTTracks_principals->at(i)[j]) > cp3) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 6 && fabs(AMTTTracks_principals->at(i)[j]) > cp6) {
							gsubPrincipal[j-2].push_back(temp);
							cut = true;
						} else if (j == 7 && fabs(AMTTTracks_principals->at(i)[j]) > cp7) {
							gsubPrincipal[j-2].push_back(temp);
							cut = true;
						} else if (j == 8 && fabs(AMTTTracks_principals->at(i)[j]) > cp8) {
							gsubPrincipal[j-2].push_back(temp);
							cut = true;
						} else if (j == 9 && fabs(AMTTTracks_principals->at(i)[j]) > cp9) {
							gsubPrincipal[j-2].push_back(temp);
							cut = true;
						}
					} else {
						if (j == 0 && fabs(AMTTTracks_principals->at(i)[j]) > cp1) {
							gsubPrincipal[j+1].push_back(temp);
							cut = true;
						} else if (j == 1 && fabs(AMTTTracks_principals->at(i)[j]) > cp2) {
							gsubPrincipal[j+1].push_back(temp);
							cut = true;
						} else if (j == 2 && fabs(AMTTTracks_principals->at(i)[j]) > cp3) {
							gsubPrincipal[j+1].push_back(temp);
							cut = true;
						} else if (j == 5 && fabs(AMTTTracks_principals->at(i)[j]) > cp7) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 6 && fabs(AMTTTracks_principals->at(i)[j]) > cp8) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 7 && fabs(AMTTTracks_principals->at(i)[j]) > cp9) {
							gsubPrincipal[j].push_back(temp);
							cut = true;
						}
					}
				}//principal loop
				vector<long double> tempp = {static_cast<long double>(jentry), i};
				if (cut) {
					gmasterlist.push_back(tempp);
					PTDif->Fill(0,trkParts_pt->at(previous));
				} else {
					PTDif->Fill(1,trkParts_pt->at(previous));
				}
				chi2_good->Fill(AMTTTracks_chi2->at(i)/AMTTTracks_ndof->at(i));
				pt_good->Fill(AMTTTracks_pt->at(i));
				phi_good->Fill(AMTTTracks_phi0->at(i));
				eta_good->Fill(AMTTTracks_eta->at(i));
				++GoodCounter;
				++GTN;
			} else if (!match) {
				for (double j = 0; j < AMTTTracks_principals->at(i).size(); j++) {
					h1_Principals_Bad[j]->Fill(fabs(AMTTTracks_principals->at(i)[j]));
					if(temptemp[0] > (num/2)-2)
						h1_Principals_Bada[j]->Fill(AMTTTracks_principals->at(i)[j]);
					else
						h1_Principals_VBada[j]->Fill(AMTTTracks_principals->at(i)[j]);
					vector<long double> temp = {static_cast<long double>(jentry), i, fabs(AMTTTracks_principals->at(i)[j])};
					if (fsss == "6x6"){
						if (j == 0 && fabs(AMTTTracks_principals->at(i)[j]) > cp0) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 1 && fabs(AMTTTracks_principals->at(i)[j]) > cp1) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 2 && fabs(AMTTTracks_principals->at(i)[j]) > cp2) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 3 && fabs(AMTTTracks_principals->at(i)[j]) > cp3) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 6 && fabs(AMTTTracks_principals->at(i)[j]) > cp6) {
							bsubPrincipal[j-2].push_back(temp);
							cut = true;
						} else if (j == 7 && fabs(AMTTTracks_principals->at(i)[j]) > cp7) {
							bsubPrincipal[j-2].push_back(temp);
							cut = true;
						} else if (j == 8 && fabs(AMTTTracks_principals->at(i)[j]) > cp8) {
							bsubPrincipal[j-2].push_back(temp);
							cut = true;
						} else if (j == 9 && fabs(AMTTTracks_principals->at(i)[j]) > cp9) {
							bsubPrincipal[j-2].push_back(temp);
							cut = true;
						}
					} else {
						if (j == 0 && fabs(AMTTTracks_principals->at(i)[j]) > cp1) {
							bsubPrincipal[j+1].push_back(temp);
							cut = true;
						} else if (j == 1 && fabs(AMTTTracks_principals->at(i)[j]) > cp2) {
							bsubPrincipal[j+1].push_back(temp);
							cut = true;
						} else if (j == 2 && fabs(AMTTTracks_principals->at(i)[j]) > cp3) {
							bsubPrincipal[j+1].push_back(temp);
							cut = true;
						} else if (j == 5 && fabs(AMTTTracks_principals->at(i)[j]) > cp7) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 6 && fabs(AMTTTracks_principals->at(i)[j]) > cp8) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						} else if (j == 7 && fabs(AMTTTracks_principals->at(i)[j]) > cp9) {
							bsubPrincipal[j].push_back(temp);
							cut = true;
						}
					}
				}//principal loop
				vector<long double> tempp = {static_cast<long double>(jentry), i};
				if (cut) {
					bmasterlist.push_back(tempp);
				}
				chi2_notgood->Fill(AMTTTracks_chi2->at(i)/AMTTTracks_ndof->at(i));
				if(temptemp[0] > (num/2)-2)
					chi2_bad->Fill(AMTTTracks_chi2->at(i)/AMTTTracks_ndof->at(i));
				else
					chi2_verybad->Fill(AMTTTracks_chi2->at(i)/AMTTTracks_ndof->at(i));
				if(AMTTTracks_chi2->at(i)/AMTTTracks_ndof->at(i) > maxchi2)
					maxchi2 = AMTTTracks_chi2->at(i)/AMTTTracks_ndof->at(i);
				pt_bad->Fill(AMTTTracks_pt->at(i));
				phi_bad->Fill(AMTTTracks_phi0->at(i));
				eta_bad->Fill(AMTTTracks_eta->at(i));
				++BadCounter;
				++BTN;
			}//principal loop
			vector<int> tempml = MatchLevel(array2);
			if (tempml[1] == 0) {
				if (tempml[0] == 6) {
					PTDif_6->Fill(!cut,trkParts_pt->at(0));
				} else if (tempml[0] == 5){
					if (AMTTTracks_ndof->at(i) == 12){
						PTDif_5P12->Fill(!cut,trkParts_pt->at(0));
					} else if (AMTTTracks_ndof->at(i) == 10){
						PTDif_5P10->Fill(!cut,trkParts_pt->at(0));
					}
				} else if (tempml[0] == 4){
					PTDif_4->Fill(!cut,trkParts_pt->at(0));
				}
				if (tempml[0] == 6 || (tempml[0] == 5 && AMTTTracks_ndof->at(i) == 10)) {
					PTChiDif->Fill(MaxMatch2.second,trkParts_pt->at(0));
				}
			}
			if (tempml[1] == 0 || bMM) {
				if(tempml[0] > MaxMatch.first) {
					MaxMatch.first = tempml[0];
					MaxMatch.second = !cut;
				} else if(tempml[0] == MaxMatch.first) {
					if(!cut)
						MaxMatch.second = true;
				}
				if(tempml[0] > MaxMatch2.first) {
					MaxMatch2.first = tempml[0];
					MaxMatch2.second = !cut;
				} else if(tempml[0] == MaxMatch2.first) {
					if(AMTTTracks_chi2->at(i) < 14.6)
						MaxMatch2.second = true;
				}
			}
			if (fsss == "6x6") {
				Mismatch->Fill(!cut,6-temptemp[0]);
				hMismatch->Fill(6-temptemp[0]);
			} else {
				Mismatch->Fill(!cut,5-temptemp[0]);
				hMismatch->Fill(5-temptemp[0]);
			}
		}//track loop
		goodCounter->Fill(GTN);
		badCounter->Fill(BTN);
		eventCounter->Fill(BTN + GTN);
		vector<bool> tempb = radiusStubTPidChecker(TTStubs_r,TTStubs_tpId);
		int tempCount = boolArrayToInt(tempb);
		//std::cout<<tempCount<<std::endl;
		PTDif_G4T->Fill((tempCount>=4),trkParts_pt->at(0));
		PTDif_G5T->Fill((tempCount>=5),trkParts_pt->at(0));
		PTDif_G6T->Fill((tempCount>=6),trkParts_pt->at(0));
		PTDif_Gn1T->Fill(tempb[0],trkParts_pt->at(0));
		PTDif_Gn2T->Fill(tempb[1],trkParts_pt->at(0));
		PTDif_Gn3T->Fill(tempb[2],trkParts_pt->at(0));
		PTDif_Gn4T->Fill(tempb[3],trkParts_pt->at(0));
		PTDif_Gn5T->Fill(tempb[4],trkParts_pt->at(0));
		PTDif_Gn6T->Fill(tempb[5],trkParts_pt->at(0));
		PhiDif_Gn1T->Fill(tempb[0],trkParts_phi->at(0));
		PhiDif_Gn2T->Fill(tempb[1],trkParts_phi->at(0));
		PhiDif_Gn3T->Fill(tempb[2],trkParts_phi->at(0));
		PhiDif_Gn4T->Fill(tempb[3],trkParts_phi->at(0));
		PhiDif_Gn5T->Fill(tempb[4],trkParts_phi->at(0));
		PhiDif_Gn6T->Fill(tempb[5],trkParts_phi->at(0));
		/*if(MaxMatch.first == 6){
			PTDif_6->Fill(MaxMatch.second,trkParts_pt->at(0));
		} else if(MaxMatch.first == 5 && AMTTTracks_ndof->at(i) == 10){
			PTDif_5P10->Fill(MaxMatch.second,trkParts_pt->at(0));
		} else if(MaxMatch.first == 4){
			PTDif_4->Fill(MaxMatch.second,trkParts_pt->at(0));
		}*/
		if (tempCount >= 5) {
			PTDif_AMGn4T->Fill(MaxMatch.first >= 4,trkParts_pt->at(0));
			PhiDif_AMGn4T->Fill(MaxMatch.first >= 4,trkParts_phi->at(0));
		}
		if (tempCount >= 5) {
			PhiDif_AMGn5T->Fill(MaxMatch.first >= 5,trkParts_phi->at(0));
			PTDif_AMGn5T->Fill(MaxMatch.first >= 5,trkParts_pt->at(0));
		}
		if (tempCount >= 5) {
			PTDif_AMGn6T->Fill(MaxMatch.first >= 6,trkParts_pt->at(0));
			PhiDif_AMGn6T->Fill(MaxMatch.first >= 6,trkParts_phi->at(0));
		}
		PTDif_GnTT->Fill((tempCount>=4),trkParts_pt->at(0));
		PTDif_AMGnnTT->Fill(MaxMatch.first >= 4,trkParts_pt->at(0));
		if (tempCount >= 4 && MaxMatch.first >= 4) {
			PTDifT->Fill(MaxMatch.second,trkParts_pt->at(0));
			PTDifTChi->Fill(MaxMatch2.second,trkParts_pt->at(0));
		} else {
			PTDifT->Fill(false,trkParts_pt->at(0));
			PTDifTChi->Fill(false,trkParts_pt->at(0));
		}
		++AMDenom;
		TGTN += GTN;
		TBTN += BTN;
	}//event loop

	std::cout<<"HERE"<<std::endl;
	std::cout<<nentries<<std::endl;
	
	//std::cout<<"HERE with "<<input1<<" and "<<input2<<std::endl;
	
	/*double goodRetention99 = .99, badRetention99 = 1;
	 double goodRetention98 = .98, badRetention98 = 1;
	 for (double a = input1; a > input2; a -= 1) {
	 for (double b = 33; b > 12; b -= 1) {
	 for (double c = 22; c > 8; c -= 1) {
	 for (double d = 14; d > 6; d -= 1) {
	 for (double e = 27; e > 2; e -= 1) {
	 for (double f = 8; f > 2; f -= 1) {
	 for (double g = 3; g > 2; g -= 1) {
	 for (double h = 5; h > 2; h -= 1) {
	 vector<vector<vector<long double> > > ngsubPrincipal;
	 vector<vector<vector<long double> > > nbsubPrincipal;
	 vector<vector<long double> > nmidstep;
	 for(int i=0; i<8; ++i){
	 ngsubPrincipal.push_back(nmidstep);
	 nbsubPrincipal.push_back(nmidstep);
	 }
	 
	 for (unsigned long j = 0; j < gsubPrincipal.size(); ++j) {
	 for (double i = 0; i < gsubPrincipal[j].size(); i++) {
	 vector<long double> temp = gsubPrincipal[j][i];
	 if (j == 0 && gsubPrincipal[j][i][2] > a)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 1 && gsubPrincipal[j][i][2] > b)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 2 && gsubPrincipal[j][i][2] > c)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 3 && gsubPrincipal[j][i][2] > d)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 4 && gsubPrincipal[j][i][2] > e)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 5 && gsubPrincipal[j][i][2] > f)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 6 && gsubPrincipal[j][i][2] > g)
	 ngsubPrincipal[j].push_back(temp);
	 if (j == 7 && gsubPrincipal[j][i][2] > h)
	 ngsubPrincipal[j].push_back(temp);
	 }
	 }
	 
	 for (unsigned long j = 0; j < bsubPrincipal.size(); ++j) {
	 for (double i = 0; i < bsubPrincipal[j].size(); i++) {
	 vector<long double> temp = bsubPrincipal[j][i];
	 if (j == 0 && bsubPrincipal[j][i][2] > a)
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 1 && bsubPrincipal[j][i][2] > b)
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 2 && bsubPrincipal[j][i][2] > c)
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 3 && bsubPrincipal[j][i][2] > d)
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 4 && bsubPrincipal[j][i][2] > e)
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 5 && bsubPrincipal[j][i][2] > f)d
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 6 && bsubPrincipal[j][i][2] > g)
	 nbsubPrincipal[j].push_back(temp);
	 if (j == 7 && bsubPrincipal[j][i][2] > h)
	 nbsubPrincipal[j].push_back(temp);
	 }
	 }
	 
	 vector<vector<long double> > gmaster, bmaster;
	 for (int i = 0; i < 8; ++i) {
	 gmaster = compare(gmaster,ngsubPrincipal[i]);
	 }
	 for (int i = 0; i < 8; ++i) {
	 bmaster = compare(bmaster,nbsubPrincipal[i]);
	 }
	 pair<double, double> returnedValue = {1-(gmaster.size()/GoodCounter), 1-(bmaster.size()/BadCounter)};
	 if (get<0>(returnedValue) >= .99 && get<1>(returnedValue) < badRetention99) {
	 goodRetention99 = get<0>(returnedValue);
	 badRetention99 = get<1>(returnedValue);
	 std::cout<<a<<" : "<<b<<" : "<<c<<" : "<<d<<" : "<<e<<" : "<<f<<" : "<<g<<" : "<<h<<" : "<<goodRetention99<<" : "<<badRetention99<<"========"<<std::endl;
	 }
	 if (get<0>(returnedValue) >= .98 && get<1>(returnedValue) < badRetention98) {
	 goodRetention98 = get<0>(returnedValue);
	 badRetention98 = get<1>(returnedValue);
	 std::cout<<a<<" : "<<b<<" : "<<c<<" : "<<d<<" : "<<e<<" : "<<f<<" : "<<g<<" : "<<h<<" : "<<goodRetention98<<" : "<<badRetention98<<"*******"<<std::endl;
	 }
	 }
	 }
	 }
	 }
	 }
	 }
		}
	 }*/
	
	
	TFile *savefile=new TFile(type + NPP + "histograms.root","RECREATE");
	
	/*double gocorrelation[num-4][num-4], bocorrelation[num-4][num-4];
	for (int i = 0; i < num-4; ++i) {
		for (int j = 0; j < num-4; ++j) {
			gocorrelation[i][j] = compare(gsubPrincipal[i],gsubPrincipal[j]).size();
			bocorrelation[i][j] = compare(bsubPrincipal[i],bsubPrincipal[j]).size();
		}
	}

	std::cout<<"Between OR and AND"<<std::endl;
	
	double gacorrelation[num-4][num-4], bacorrelation[num-4][num-4];
	for (int i = 0; i < num-4; ++i) {
		for (int j = 0; j < num-4; ++j) {
			gacorrelation[i][j] = singleCompare(gsubPrincipal[i],gsubPrincipal[j]).size();
			bacorrelation[i][j] = singleCompare(bsubPrincipal[i],bsubPrincipal[j]).size();
		}
	}
	
		//Phi Coefficient
	TH2F *tgcorrelation = new TH2F("Good Match Correlations","Good Match Correlations",8,-0.5,7.5,8,-0.5,7.5);
	TH2F *tbcorrelation = new TH2F("Bad Match Correlations","Bad Match Correlations",8,-0.5,7.5,8,-0.5,7.5);
	for (int i = 0; i < num-4; ++i) {
		for (int j = 0; j < num-4; ++j) {
			double ga = gacorrelation[i][j],
			gb = gsubPrincipal[i].size() - gacorrelation[i][j],
			gc = gsubPrincipal[j].size() - gacorrelation[i][j],
			gd = GoodCounter - gocorrelation[i][j],
			ba = bacorrelation[i][j],
			bb = bsubPrincipal[i].size() - bacorrelation[i][j],
			bc = bsubPrincipal[j].size() - bacorrelation[i][j],
			bd = BadCounter - bocorrelation[i][j];
			tgcorrelation->SetBinContent(tgcorrelation->FindBin(i,j),((ga*gd-gb*gc)/sqrt((ga+gb)*(ga+gc)*(gb+gd)*(gc+gd))));
			tbcorrelation->SetBinContent(tbcorrelation->FindBin(i,j),((ba*bd-bb*bc)/sqrt((ba+bb)*(ba+bc)*(bb+bd)*(bc+bd))));
		}
	}
	TH2F *subcorrelation=(TH2F*)tgcorrelation->Clone("subcorrelaton");
	subcorrelation->SetName("Subtraction Correlations");
	subcorrelation->SetTitle("Subtraction Correlations");
	subcorrelation->Add(tbcorrelation,-1.);
	TCanvas *tempg = new TCanvas("Good Match Correlations","Good Match Correlations",1000,1000);
	tgcorrelation->Draw("colz, text");
	tgcorrelation->SetStats(0);
	tgcorrelation->GetXaxis()->SetTickLength(0);
	tgcorrelation->GetYaxis()->SetTickLength(0);
	tempg->Write();
	TCanvas *tempb = new TCanvas("Bad Match Correlations","Bad Match Correlations",1000,1000);
	tbcorrelation->Draw("colz, text");
	tbcorrelation->SetStats(0);
	tbcorrelation->GetXaxis()->SetTickLength(0);
	tbcorrelation->GetYaxis()->SetTickLength(0);
	tempb->Write();
	TCanvas *temps = new TCanvas("Subtraction Correlations","Subtraction Correlations",1000,1000);
	subcorrelation->Draw("colz, text");
	subcorrelation->SetStats(0);
	subcorrelation->GetXaxis()->SetTickLength(0);
	subcorrelation->GetYaxis()->SetTickLength(0);
	temps->Write();*/
	
	for (int i = 0; i < num-4; ++i) {
		std::cout<<gsubPrincipal[i].size()<<" : "<<bsubPrincipal[i].size()<<" : "<<((gsubPrincipal[i].size()*100.)/bsubPrincipal[i].size())<<std::endl;
	}
	std::cout<<TGTN<<" : "<<GoodCounter<<" : "<<gmasterlist.size()<<" : "<<((gmasterlist.size()*100.)/GoodCounter)<<" : "<<100-((gmasterlist.size()*100.)/GoodCounter)<<std::endl;
	std::cout<<TBTN<<" : "<<BadCounter<<" : "<<bmasterlist.size()<<" : "<<(bmasterlist.size()*100.)/BadCounter<<" : "<<100-((bmasterlist.size()*100.)/BadCounter)<<std::endl;
	
	std::cout<<AMDenom<<" : "<<(GoodCounter + BadCounter)<<std::endl;
	
	//std::cout<<"Number of principal miss fits : "<<errorCounter<<std::endl;
	
	//std::cout<<matchn<<" : "<<nentries<<std::endl;
	
	//std::cout<<"I did things"<<std::endl;
	
		//Intregals
	std::vector<TGraph*> graphs;
	for(double i=0; i<num; ++i) {
		double gdiv = h1_Principals_Good[i]->Integral(0,-1);
		double bdiv = h1_Principals_Bad[i]->Integral(0,-1);
		float x[5000],y[5000];
		for(int j=h1_Principals_Good[i]->GetNbinsX(); j!=-1; --j) {
			x[j]=h1_Principals_Good[i]->Integral(0,j)/gdiv;
			y[j]=h1_Principals_Bad[i]->Integral(0,j)/bdiv;
			if (x[j] > 0.7 && x[j] < .9999) {
				//std::cout<<i<<" : "<<h1_Principals_Good[i]->GetBinLowEdge(j)<<" : "<<y[j]<<" : "<<x[j]<<std::endl;
			}
		}
		/*graphs.push_back(new TGraph(5000,x,y));
		TCanvas *temp=new TCanvas(TString::Format("%1.0f",i),TString::Format("principal%1.0f",i),1000,1000);
		graphs[i]->SetTitle(TString::Format("Efficiency of Principal %1.0f",i));
		TH2F *background=new TH2F("background",";#epsilon_{good};#epsilon_{bad}",1,0,1,1,0,1);
		background->Draw();
		background->SetStats(0);
		graphs[i]->Draw("same");
		TLine *line=new TLine(0,0,1,1);
		line->SetLineStyle(2);
		line->SetLineColor(2);
		line->Draw();
		temp->Write();*/
	}
	//std::cout<<"I did other things"<<std::endl;
	
		//Chi^2 Integral
	/*double gchidiv = chi2_good->Integral(0,-1);
	double bchidiv = chi2_notgood->Integral(0,-1);
	float x2[10000],y2[10000];
	float x3[1],y3[1];
	x3[0] = 1-((float)(gmasterlist.size())/((float)GoodCounter));
	y3[0] = 1-((float)(bmasterlist.size())/((float)BadCounter));
	for(int j=chi2_good->GetNbinsX(); j!=-1; --j) {
		x2[j]=chi2_good->Integral(0,j)/gchidiv;
		y2[j]=chi2_notgood->Integral(0,j)/bchidiv;
		if (chi2_good->GetBinLowEdge(j) > 14 && chi2_good->GetBinLowEdge(j) < 100) {
			std::cout<<chi2_good->GetBinLowEdge(j)<<" : "<<j<<" : "<<x2[j]<<" : "<<y2[j]<<std::endl;
		}
	}
	TGraph *graph = new TGraph(chi2_good->GetNbinsX(),x2,y2);
	TGraph *graph2 = new TGraph(1,x3,y3);
	graph2->SetMarkerStyle(20);
	TCanvas *tempchi = new TCanvas("Chi2Efficiency","Chi2Efficiency",1000,1000);
	graph->SetTitle("Chi2Efficiency");
	graph->SetName("Chi2EfficiencyGraph");
	TH2F *background=new TH2F("background",";#epsilon_{good};#epsilon_{bad}",1,0,1,1,0,1);
	background->Draw();
	background->SetStats(0);
	graph->Draw("same");
	graph2->Draw("same,p");
	TLine *line=new TLine(0,0,1,1);
	line->SetLineStyle(2);
	line->SetLineColor(2);
	line->Draw();
	graph->Write();
	tempchi->Write();*/
	//std::cout<<"I did other other things"<<std::endl;
	
		//save Histograms
	/*for(double i=0; i<num; ++i) {
		TCanvas *temp = new TCanvas(TString::Format("GoodAndBad%1.0f",i),TString::Format("GoodAndBad%1.0f",i),1000,1000);
		h1_Principals_Gooda[i]->DrawNormalized();
		h1_Principals_Bada[i]->SetLineColor(2);
		h1_Principals_Bada[i]->DrawNormalized("same");
		h1_Principals_VBada[i]->SetLineColor(3);
		h1_Principals_VBada[i]->DrawNormalized("same");
		temp->SetLogy();
		TLegend *leg = new TLegend(0.8,0.6,0.95,0.75,"Collections");
		leg->AddEntry(h1_Principals_Gooda[i],"Good","l");
		leg->AddEntry(h1_Principals_Bada[i],"Bad","l");
		leg->AddEntry(h1_Principals_VBada[i],"Very Bad","l");
		leg->SetFillColor(0);
		leg->Draw();
		temp->SetTitle(TString::Format("Principal %1.0f",i));
		h1_Principals_Gooda[i]->Write();
		h1_Principals_Bada[i]->Write();
		h1_Principals_VBada[i]->Write();
		temp->Write();
	}*/
	
	/*TCanvas *temp = new TCanvas("C2GoodAndBad","C2GoodAndBad",1000,1000);
	chi2_good->GetXaxis()->SetRangeUser(0,100);
	chi2_good->DrawNormalized();
	chi2_bad->SetLineColor(2);
	chi2_bad->DrawNormalized("same");
	chi2_verybad->SetLineColor(3);
	chi2_verybad->DrawNormalized("same");
	temp->SetLogy();
	TLegend *leg = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg->AddEntry(chi2_good,"Good","l");
	leg->AddEntry(chi2_bad,"Bad","l");
	leg->AddEntry(chi2_verybad,"VeryBad","l");
	leg->SetFillColor(0);
	leg->Draw();
	chi2_bad->Write();
	chi2_good->Write();
	chi2_verybad->Write();
	temp->Write();
	
	TCanvas *temp1 = new TCanvas("PtGoodAndBad","PtGoodAndBad",1000,1000);
	pt_bad->Draw();
	pt_bad->SetLineColor(2);
	pt_good->Draw("same");
	temp1->SetLogy();
	TLegend *leg1 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg1->AddEntry(pt_good,"Good","l");
	leg1->AddEntry(pt_bad,"Bad","l");
	leg1->SetFillColor(0);
	leg1->Draw();
	temp1->Write();
	
	TCanvas *temp2 = new TCanvas("PhiGoodAndBad","PhiGoodAndBad",1000,1000);
	phi_bad->Draw();
	phi_bad->SetLineColor(2);
	phi_good->Draw("same");
	temp2->SetLogy();
	TLegend *leg2 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg2->AddEntry(phi_good,"Good","l");
	leg2->AddEntry(phi_bad,"Bad","l");
	leg2->SetFillColor(0);
	leg2->Draw();
	temp2->Write();
	
	TCanvas *temp3 = new TCanvas("EtaGoodAndBad","EtaGoodAndBad",1000,1000);
	eta_bad->Draw();
	eta_bad->SetLineColor(2);
	eta_good->Draw("same");
	temp3->SetLogy();
	TLegend *leg3 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg3->AddEntry(eta_good,"Good","l");
	leg3->AddEntry(eta_bad,"Bad","l");
	leg3->SetFillColor(0);
	leg3->Draw();
	temp3->Write();
	
	TCanvas *temp4 = new TCanvas("Counter","Counter",1000,1000);
	histoCounter->Draw();
	temp4->SetLogy();
	temp4->Write();
	
	TCanvas *temp4a = new TCanvas("Counter2","Counter2",1000,1000);
	histoCounter2->Draw();
	temp4a->SetLogy();
	temp4a->Write();
	
	TCanvas *temp5 = new TCanvas("Track Candidates","Track Candidates",1000,1000);
	histoTrackCandidates->Draw();
	temp5->Write();
	
	TCanvas *temp6 = new TCanvas("Good Tracks per Event","Good Tracks per Event",1000,1000);
	goodCounter->Draw();
	temp6->SetLogy();
	temp6->Write();
	
	TCanvas *temp7 = new TCanvas("Bad Tracks per Event","Bad Tracks per Event",1000,1000);
	badCounter->Draw();
	temp7->SetLogy();
	temp7->Write();
	
	TCanvas *temp8 = new TCanvas("Tracks per Event","Bad Tracks per Event",1000,1000);
	eventCounter->Draw();
	temp8->SetLogy();
	temp8->Write();

	TCanvas *temp9 = new TCanvas("PTDiff","PT Differential;p_{T};#epsilon",1000,1000);
	PTDif_6->Draw();
	PTDif_5P12->SetLineColor(2);
	PTDif_5P12->Draw("same");
	PTDif_5P10->SetLineColor(6);
	PTDif_5P10->Draw("same");
	PTDif_4->SetLineColor(4);
	PTDif_4->Draw("same");
	TLegend *leg4 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg4->AddEntry(PTDif_6,"6 Match","l");
	leg4->AddEntry(PTDif_5P12,"5 Match P12","l");
	leg4->AddEntry(PTDif_5P10,"5 Match P10","l");
	leg4->AddEntry(PTDif_4,"4 Match","l");
	leg4->SetFillColor(0);
	leg4->Draw();
	PTDif_6->Write();
	PTDif_5P12->Write();
	PTDif_5P10->Write();
	PTDif_4->Write();
	temp9->SetLogx();
	temp9->Write();

	TCanvas *temp10 = new TCanvas("PTDiff_GnT","PT Differential;p_{T};#epsilon",1000,1000);
	PTDif_G6T->Draw();
	PTDif_G5T->SetLineColor(2);
	PTDif_G5T->Draw("same");
	PTDif_G4T->SetLineColor(4);
	PTDif_G4T->Draw("same");
	TLegend *leg5 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg5->AddEntry(PTDif_G6T,"=6 Match","l");
	leg5->AddEntry(PTDif_G5T,">=5 Match","l");
	leg5->AddEntry(PTDif_G4T,">=4 Match","l");
	leg5->SetFillColor(0);
	leg5->Draw();
	PTDif_G6T->Write();
	PTDif_G5T->Write();
	PTDif_G4T->Write();
	temp10->SetLogx();
	temp10->Write();

	TCanvas *temp11 = new TCanvas("PTDiff_Add","PT Differential;p_{T};#epsilon",1000,1000);
	TEfficiency *PTDif_Add = (TEfficiency*)PTDif_6->Clone("PTDif_Add");
	*PTDif_Add+=*PTDif_5P10;
	PTDif_Add->Draw();
	PTDif_Add->Write();
	PTChiDif->SetLineColor(2);
	PTChiDif->Draw("same");
	PTChiDif->Write();
	TLegend *leg6 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg6->AddEntry(PTDif_Add,"Principal Cut","l");
	leg6->AddEntry(PTChiDif,"Chi^2 Cut","l");
	leg6->Draw();
	temp11->SetLogx();
	temp11->Write();

	TCanvas *temp12 = new TCanvas("TeMismatch","Mismatch;# of Mismatched Stubs;#epsilon",1000,1000);
	Mismatch->Draw();
	Mismatch->Write();
	temp12->Write();

	TCanvas *temp13 = new TCanvas("hiMismatch","Mismatch;# of Mismatched Stubs;#epsilon",1000,1000);
	hMismatch->DrawNormalized();
	hMismatch->Write();
	temp13->Write();

	TCanvas *temp18 = new TCanvas("PhiDiff_GnnT","Phi Differential;#phi;#epsilon",1000,1000);
	PhiDif_Gn1T->Draw();
	PhiDif_Gn2T->SetLineColor(2);
	PhiDif_Gn2T->Draw("same");
	PhiDif_Gn3T->SetLineColor(3);
	PhiDif_Gn3T->Draw("same");
	PhiDif_Gn4T->SetLineColor(4);
	PhiDif_Gn4T->Draw("same");
	PhiDif_Gn5T->SetLineColor(6);
	PhiDif_Gn5T->Draw("same");
	PhiDif_Gn6T->SetLineColor(8);
	PhiDif_Gn6T->Draw("same");
	TLegend *leg10 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg10->AddEntry(PhiDif_Gn1T,"Layer 1 Missing","l");
	leg10->AddEntry(PhiDif_Gn2T,"Layer 2 Missing","l");
	leg10->AddEntry(PhiDif_Gn3T,"Layer 3 Missing","l");
	leg10->AddEntry(PhiDif_Gn4T,"Layer 4 Missing","l");
	leg10->AddEntry(PhiDif_Gn5T,"Layer 5 Missing","l");
	leg10->AddEntry(PhiDif_Gn6T,"Layer 6 Missing","l");
	leg10->SetFillColor(0);
	leg10->Draw();
	PhiDif_Gn1T->Write();
	PhiDif_Gn2T->Write();
	PhiDif_Gn3T->Write();
	PhiDif_Gn4T->Write();
	PhiDif_Gn5T->Write();
	PhiDif_Gn6T->Write();
	temp18->Write();
	
	TCanvas *temp19 = new TCanvas("PTDiff_GnnT","PT Differential;p_{T};#epsilon",1000,1000);
	PTDif_Gn1T->Draw();
	PTDif_Gn2T->SetLineColor(2);
	PTDif_Gn2T->Draw("same");
	PTDif_Gn3T->SetLineColor(3);
	PTDif_Gn3T->Draw("same");
	PTDif_Gn4T->SetLineColor(4);
	PTDif_Gn4T->Draw("same");
	PTDif_Gn5T->SetLineColor(6);
	PTDif_Gn5T->Draw("same");
	PTDif_Gn6T->SetLineColor(8);
	PTDif_Gn6T->Draw("same");
	TLegend *leg11 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg11->AddEntry(PTDif_Gn1T,"Layer 1 Missing","l");
	leg11->AddEntry(PTDif_Gn2T,"Layer 2 Missing","l");
	leg11->AddEntry(PTDif_Gn3T,"Layer 3 Missing","l");
	leg11->AddEntry(PTDif_Gn4T,"Layer 4 Missing","l");
	leg11->AddEntry(PTDif_Gn5T,"Layer 5 Missing","l");
	leg11->AddEntry(PTDif_Gn6T,"Layer 6 Missing","l");
	leg11->SetFillColor(0);
	leg11->Draw();
	PTDif_Gn1T->Write();
	PTDif_Gn2T->Write();
	PTDif_Gn3T->Write();
	PTDif_Gn4T->Write();
	PTDif_Gn5T->Write();
	PTDif_Gn6T->Write();
	temp19->SetLogx();
	temp19->Write();

	TCanvas *temp20 = new TCanvas("PhiDiff_AMGnnT","Phi Differential;#phi;#epsilon",1000,1000);
	PhiDif_AMGn4T->Draw();
	PhiDif_AMGn5T->SetLineColor(2);
	PhiDif_AMGn5T->Draw("same");
	PhiDif_AMGn6T->SetLineColor(3);
	PhiDif_AMGn6T->Draw("same");
	TLegend *leg12 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg12->AddEntry(PhiDif_AMGn4T,">=4 Layer Match","l");
	leg12->AddEntry(PhiDif_AMGn5T,">=5 Layer Match","l");
	leg12->AddEntry(PhiDif_AMGn6T,"=6 Layer Match","l");
	leg12->SetFillColor(0);
	leg12->Draw();
	PhiDif_Gn4T->Write();
	PhiDif_Gn5T->Write();
	PhiDif_Gn6T->Write();
	temp20->Write();
	
	TCanvas *temp21 = new TCanvas("PTDiff_AMGnnT","PT Differential;p_{T};#epsilon",1000,1000);
	PTDif_AMGn4T->Draw();
	PTDif_AMGn5T->SetLineColor(2);
	PTDif_AMGn5T->Draw("same");
	PTDif_AMGn6T->SetLineColor(3);
	PTDif_AMGn6T->Draw("same");
	TLegend *leg13 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg13->AddEntry(PTDif_AMGn4T,">=4 Layer Match","l");
	leg13->AddEntry(PTDif_AMGn5T,">=5 Layer Match","l");
	leg13->AddEntry(PTDif_AMGn6T,"=6 Layer Match","l");
	leg13->SetFillColor(0);
	leg13->Draw();
	PTDif_AMGn4T->Write();
	PTDif_AMGn5T->Write();
	PTDif_AMGn6T->Write();
	temp21->SetLogx();
	temp21->Write();

	TCanvas *temp22 = new TCanvas("PTDiffT","PT Differential;p_{T};#epsilon",1000,1000);
	PTDifT->Draw();
	PTDifT->Write();
	PTDifTChi->SetLineColor(2);
	PTDifTChi->Draw("same");
	PTDifTChi->Write();
	PTDif_AMGnnTT->SetLineColor(3);
	PTDif_AMGnnTT->SetLineStyle(2);
	PTDif_AMGnnTT->Draw("same");
	PTDif_AMGnnTT->Write();
	PTDif_GnTT->SetLineColor(4);
	PTDif_GnTT->Draw("same");
	PTDif_GnTT->Write();
	TLegend *leg14 = new TLegend(0.8,0.6,0.95,0.75,"Collections");
	leg14->AddEntry(PTDifT,"Total Efficiency with Principal Cut","l");
	leg14->AddEntry(PTDifTChi,"Total Efficiency with Chi^2 < 14.6 Cut","l");
	leg14->AddEntry(PTDif_AMGnnTT,"AM Efficiency","l");
	leg14->AddEntry(PTDif_GnTT,"Stub Efficiency","l");
	leg14->SetFillColor(0);
	leg14->Draw();
	temp22->SetLogx();
	temp22->Write();

	
	for (unsigned i = 0; i<num; ++i)
		h1_Principals_Bada[i]->Write();
	for (unsigned i = 0; i<num; ++i)
		h1_Principals_Gooda[i]->Write();
	for (unsigned i = 0; i<num; ++i)
		graphs[i]->Write();
	*/
	savefile->Close();
}
