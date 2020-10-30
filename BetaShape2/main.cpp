
#include"CreatBetaShape.h"
#include"GammaShape.h"





int main() {

	string filename = "cow-2";
	//filename = "funkyCube";
	string format = "ply";
	int i=0;
	
	if (i == 1) {
	CreatBetaShape cbs(filename, format, 10, 10);
	cbs.creatBetaShape();
	cbs.showPointCloud();
	}
	else {
	GammaShape gp(10, filename, format, 10, 30,0.95);
	gp.creatGammaShape();
	}
}


