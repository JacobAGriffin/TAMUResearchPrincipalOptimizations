#define principalIntegral_cxx
#include "principalIntegral.h"
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

void principalIntegral::loop(bool highPt)
{
//   In a ROOT session, you can do:
//      Root > .L principalIntegral.C
//      Root > principalIntegral t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	
	vector<TH1F*> h1_Principals_Good, h1_Principals_Bad;
	TH1F *test=new TH1F("","",10000,0,500);
	for(unsigned i=0; i<12; ++i){
		h1_Principals_Good.push_back((TH1F*)test->Clone(TString::Format("goodPrincipal%1.0d",i)));
		h1_Principals_Bad.push_back((TH1F*)test->Clone(TString::Format("badPrincipal%1.0d",i)));
	}

	//TString type = "Muon";
	//TString type = "Pion";
	TString type = "Electron";
	//TString type = "PileUp";
	//TString type = "T~T";
	TString matchLogic = "6x6";
	//TString matchLogic = "5x6";
	//TString matchLogic = "5x5";
	TString NPP = "PU";
	//TString NPP = "NPU";
	bool matchPileUp = false;
	if (NPP == "PU"	|| type == "T~T" || type == "PileUP") {
		matchPileUp = true;
	}

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		bool contained = Containment_tt27Test(trkParts_pt->at(0), trkParts_phi->at(0), trkParts_vz->at(0), trkParts_eta->at(0), trkParts_charge->at(0), trkParts_vx->at(0), trkParts_vy->at(0));
		if (trkParts_phi->at(0) < .8 || trkParts_phi->at(0) > 1.55)
			contained = false; //special provision because of bank limitations
		if (type != "PileUp" && type != "T~T")
			if (!contained)
				continue;
		for (double i = 0; i < AMTTTracks_principals->size(); ++i) {
			if (AMTTTracks_preEstimatedPt->at(i) < 10 && highPt)
				continue;
			if (AMTTTracks_preEstimatedPt->at(i) >= 10 && !highPt)
				continue;
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
			//std::cout << ((AMTTTracks_ndof->at(i) / 2) + 2) << " : " << matchLevel[0] << " : " << matchLevel[1] << endl;
			if (matchLevel[0] == ((AMTTTracks_ndof->at(i) / 2) + 2) && (matchPileUp || matchLevel[1] == 0) && matchLevel[1]>=0) {
				for (double j = 0; j < AMTTTracks_principals->at(i).size(); ++j) {
					if (((AMTTTracks_ndof->at(i) / 2) + 2) == 5 && (j == 0 || j == 6))
						continue;
					h1_Principals_Good[j]->Fill(fabs(AMTTTracks_principals->at(i)[j]));
				}
			} else {
				for (double j = 0; j < AMTTTracks_principals->at(i).size(); ++j) {
					if (((AMTTTracks_ndof->at(i) / 2) + 2) == 5 && (j == 0 || j == 6))
						continue;
					h1_Principals_Bad[j]->Fill(fabs(AMTTTracks_principals->at(i)[j]));
				}
			}
		}
	}

	TString highLow = (highPt) ? "highPt" : "lowPt";
	TFile *savefile = new TFile(type + NPP + "histograms" + highLow + ".root","RECREATE");

	std::vector<TGraph*> graphs;
	for(double i=0; i<12; ++i) {
		double gdiv = h1_Principals_Good[i]->Integral(0,-1);
		double bdiv = h1_Principals_Bad[i]->Integral(0,-1);
		float x[10001],y[10001];
		for(int j=h1_Principals_Good[i]->GetNbinsX(); j!=-1; --j) {
			x[j]=h1_Principals_Good[i]->Integral(0,j)/gdiv;
			y[j]=h1_Principals_Bad[i]->Integral(0,j)/bdiv;
			if (x[j] > 0.98 && x[j] < .9999) {
				std::cout<<i<<" : "<<h1_Principals_Good[i]->GetBinLowEdge(j)<<" : "<<y[j]<<" : "<<x[j]<<std::endl;
			}
		}
		graphs.push_back(new TGraph(5000,x,y));
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
		temp->Write();
	}

	for(double i=0; i<12; ++i) {
		TCanvas *temp = new TCanvas(TString::Format("GoodAndBad%1.0f",i),TString::Format("GoodAndBad%1.0f",i),1000,1000);
		h1_Principals_Good[i]->DrawNormalized();
		h1_Principals_Bad[i]->SetLineColor(2);
		h1_Principals_Bad[i]->DrawNormalized("same");
		temp->SetLogy();
		TLegend *leg = new TLegend(0.8,0.6,0.95,0.75,"Collections");
		leg->AddEntry(h1_Principals_Good[i],"Good","l");
		leg->AddEntry(h1_Principals_Bad[i],"Bad","l");
		leg->SetFillColor(0);
		leg->Draw();
		temp->SetTitle(TString::Format("Principal %1.0f",i));
		h1_Principals_Good[i]->Write();
		h1_Principals_Bad[i]->Write();
		temp->Write();
	}

	savefile->Close();
}
