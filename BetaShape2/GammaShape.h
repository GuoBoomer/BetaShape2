#pragma once
#include "CreatBetaShape.h"
#include<algorithm>
#include<list>


class GammaShape :
	public CreatBetaShape{
 public:
	//存储块索引
	struct blocks {
		vector<int> indexs;
		blocks* next=NULL;
	};
	//三角面片索引
	struct tri_index {
		vector<PointXYZ> tir_points;
		tri_index* next = NULL;
	};
	//球心索引
	struct BallCenter {
		//将球心分为两侧
		tri_index set1; 
		tri_index set2;
		BallCenter *next = NULL;

	};


 private:

	float front_rates;	//前置比率
	float Threshold;	//停驻点判别阈值
	vector<PointXYZ> tris;
	vector<int> visited;
	vector<int> stop_points;
	int stride; //步长,应该是不需要的元素
	//通过种子点开始生长
	void regionGrow(int seed, vector<int> &noise, blocks *block,
		tri_index *tri,BallCenter *bc, vector<int> &stop_point, vector<int> &current_queue);
	//种子生长法完全体
	void regionGrow2(int seed, vector<int> &noise, blocks *block,
		tri_index *tri, BallCenter *bc, vector<int> &stop_point, vector<int> &current_queue);
	
	bool judgeExistenceOfPoints(int a, int b, int c,tri_index *tri);
	//判断圆心集合
	void judgeSetsOfBallCenter(BallCenter *bc ,double Circumball[][3], int seed,
		vector<PointXYZ> &tri);
	//判断停驻点
	bool judgeStopPoints(BallCenter *bc,double Circumball[][3],int ind);
	bool judgeStopPoints2(BallCenter *bc, double Circumball[][3],
		vector<int> ind);
	bool judgeOtherPoints(double Circumball[][3],vector<int> &indexs,double radius);

 public:

	GammaShape(int mystride,string filename, 
		string format, double nbeta, int nk,float threshold);
	//开始创建gamma shape
	void creatGammaShape();
	void showPointCloud();

};
