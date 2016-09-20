#include "CATEfficiency.C"

void runCATEfficiencyPoint(int highPt){
	vector<float> cuts = {-1,33,22,12,-1,-1,-1,3,3,3,-1,-1};
	CATEfficiency t;
	t.point(cuts, highPt);
}
