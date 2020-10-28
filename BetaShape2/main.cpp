
#include"CreatBetaShape.h"
#include"GammaShape.h"





int main() {

	string filename = "funkyCube";
	string format = "ply";
	int i;
	cin >> i;
	if (i == 1) {
	CreatBetaShape cbs(filename, format, 10, 10);
	cbs.creatBetaShape();
	cbs.showPointCloud();
	}
	else {
	GammaShape gp(10, filename, format, 10, 10,0.8);
	gp.creatGammaShape();
	}
}


