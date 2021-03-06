#define CATEfficiency_cxx
#include "CATEfficiency.h"
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
	for (unsigned i = 0; i < tempf->size(); ++i) {
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

pair<vector<int>, vector<float> > pCut(vector<float> principal, vector<float> cuts) {
	vector<int> cutNumber;
	vector<float> cutValue;
	if (principal.size() == 10) {
		principal.insert(principal.begin(),0);
		principal.insert(principal.begin() + 6,0);
	}
	for (unsigned i = 0; i < principal.size(); ++i) {
		if (cuts[i] != -1)
			if (fabs(principal[i]) > cuts[i]) {
				cutNumber.push_back(i);
				cutValue.push_back(fabs(principal[i]));
			}
	}
	pair<vector<int>, vector<float> > pair = make_pair(cutNumber, cutValue);
	return pair;
}

void	CATEfficiency::point(vector<float> cuts, bool highPt) {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	
	Long64_t nbytes = 0, nb = 0;

	bool matchPileUp = true;
	//vector<pair< vector<int>, vector<float> > > goodTracks, badTracks;
	int goodCut = 0, badCut = 0, goodTracks = 0, badTracks = 0;

	float tpMRet = 0;
	float tpMatchesSize = 0;
	float primaryTpMRet = 0;
	float primaryTpMatchesSize = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		bool contained = Containment_tt27Test(trkParts_pt->at(0), trkParts_phi->at(0), trkParts_vz->at(0), trkParts_eta->at(0), trkParts_charge->at(0), trkParts_vx->at(0), trkParts_vy->at(0));
		if (trkParts_phi->at(0) < .8 || trkParts_phi->at(0) > 1.55)
			contained = false; //special provision because of bank limitations
		if (!contained)
			continue;
		vector<pair<int,unsigned> > tpMatches;
		for (unsigned i = 0; i < AMTTTracks_principals->size(); ++i) {
			if (AMTTTracks_preEstimatedPt->at(i) < 10 && highPt)
				continue;
			if (AMTTTracks_preEstimatedPt->at(i) >= 10 && !highPt)
				continue;
			//int previous = -2;
			//bool match = true;
			vector<int> candidateStubs;
			for (unsigned j = 0; j < AMTTTracks_stubRefs->at(i).size(); ++j) {
				const int stub = AMTTTracks_stubRefs->at(i)[j];
				if (stub != 999999999) {
					int tp = TTStubs_tpId->at(stub);
					if(tp >= 0)
						if(trkParts_pt->at(tp) <= 3)
							tp = -1;
					/*if (previous == -2)
						previous=tp;
					else if (tp == -1 || previous != tp)
						match = false;*/
					candidateStubs.push_back(tp);
				} else
					candidateStubs.push_back(-1);
			}
			vector<int> matchLevel = MatchLevel(candidateStubs);
			pair<vector<int>, vector<float> > cut;
			if (matchLevel[0] == ((AMTTTracks_ndof->at(i) / 2) + 2) && (matchPileUp || matchLevel[1] == 0) && matchLevel[1]>=0) {
				int tpMatchPos=-1;
				for(unsigned p=0; p<tpMatches.size(); ++p){
					if(tpMatches[p].first==matchLevel[1]){
						tpMatchPos=p;
						break;
					}
				}
				if(tpMatchPos==-1){
					tpMatches.push_back(make_pair(matchLevel[1],0));
					tpMatchPos=tpMatches.size()-1;
				}
				cut = pCut(AMTTTracks_principals->at(i),cuts);
				++goodTracks;//goodTracks.push_back(cut);
				if (cut.first.size() > 0)
					++goodCut;
				else
					tpMatches[tpMatchPos].second++;
			} else {
				cut = pCut(AMTTTracks_principals->at(i),cuts);
				++badTracks;//badTracks.push_back(cut);
				if (cut.first.size() > 0)
					++badCut;
			}
		}//track loop
		tpMatchesSize += tpMatches.size();
		for (int z = 0; z < tpMatches.size(); ++z) {
			if (tpMatches[z].second)
				++tpMRet;
			if (tpMatches[z].first == 0) {
				++primaryTpMatchesSize;
				if (tpMatches[z].second)
					++primaryTpMRet;
			}
		}
	}//event loop
	std::cout<<cuts[1]<<" "<<cuts[2]<<" "<<cuts[3]<<" "<<cuts[7]<<" "<<cuts[8]<<" "<<cuts[9]<<" "<<nentries<<" "<<100*primaryTpMRet/primaryTpMatchesSize<<" "<<100*tpMRet/tpMatchesSize<<" "<<100 - ((goodCut*100.)/goodTracks)<<" "<<goodCut<<" "<<goodTracks<<" "<<100 - ((badCut*100.)/badTracks)<<" "<<badTracks - badCut<<" "<<badTracks<<std::endl;
}

void CATEfficiency::loop(int in1, int in2, bool highPt) {
//   In a ROOT session, you can do:
//      Root > .L CATEfficiency.C
//      Root > CATEfficiency t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//    Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
//    METHOD1:
//    	fChain->SetBranchStatus("*",0);  // disable all branches
//    	fChain->SetBranchStatus("branchname",1);  // activate branchname
//    METHOD2: replace line
//    	    fChain->GetEntry(jentry);       //read all branches
//      by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	int count = 0;
	float array1[] = {45,50,55,60};
	float array2[] = {35,40,45};
	float array3[] = {25,30,40};
	float array4[] = {3,3.5,4,5};
	float array5[] = {3};
	float array6[] = {3};
	for (float a : array1) {
		for (float b  : array2) {
			for (float c  : array3) {
				for (float d  : array4) {
					for (float e  : array5) {
						for (float f  : array6) {
							std::cout << "HERE " << count << std::endl;
							if (count >= in1 && count < in2) {
								//vector<pair< vector<int>, vector<float> > > goodTracks, badTracks;
								int goodCut = 0, badCut = 0, goodTracks = 0, badTracks = 0;

								//vector<float> cuts = {-1,33,22,12,-1,-1,-1,3,3,3,-1,-1};
								vector<float> cuts = {-1,a,b,c,-1,-1,-1,d,e,f,-1,-1};

								point(cuts, highPt);
								
							} else if (count >= in2)
								break;
							++count;
						}
					}
				}
			}
		}
	}
	//std::cout<<goodCut<<" "<<goodTracks<<" "<<((goodCut*100.)/goodTracks)<<" "<<100 - ((goodCut*100.)/goodTracks)<<std::endl;
	//std::cout<<badCut<<" "<<badTracks<<" "<<((badCut*100.)/badTracks)<<" "<<100 - ((badCut*100.)/badTracks)<<std::endl;

	/*int goodCutByPrincipal[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	int badCutByPrincipal[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
	for (unsigned i = 0; i < goodTracks.size(); ++i) {
		for (unsigned j = 0; j < goodTracks[i].first.size(); ++j) {
			++goodCutByPrincipal[goodTracks[i].first[j]];
		}
	}

	for (unsigned i = 0; i < badTracks.size(); ++i) {
		for (unsigned j = 0; j < badTracks[i].first.size(); ++j) {
			++badCutByPrincipal[badTracks[i].first[j]];
		}
	}

	for (int i = 0; i < 12; ++i) {
		std::cout<<i<<" "<<cuts[i]<<" "<<goodCutByPrincipal[i]<<" "<<badCutByPrincipal[i]<<std::endl;
	}*/
}
