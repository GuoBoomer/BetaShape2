#pragma once
#include "CreatBetaShape.h"
#include<algorithm>
#include<list>


class GammaShape :
	public CreatBetaShape{
 public:
	//�洢������
	struct blocks {
		vector<int> indexs;
		blocks* next=NULL;
	};
	//������Ƭ����
	struct tri_index {
		vector<PointXYZ> tir_points;
		tri_index* next = NULL;
	};
	//��������
	struct BallCenter {
		//�����ķ�Ϊ����
		tri_index set1; 
		tri_index set2;
		BallCenter *next = NULL;

	};


 private:

	float front_rates;	//ǰ�ñ���
	float Threshold;	//ͣפ���б���ֵ
	vector<PointXYZ> tris;
	vector<int> visited;
	vector<int> stop_points;
	int stride; //����,Ӧ���ǲ���Ҫ��Ԫ��
	//ͨ�����ӵ㿪ʼ����
	void regionGrow(int seed, vector<int> &noise, blocks *block,
		tri_index *tri,BallCenter *bc, vector<int> &stop_point, vector<int> &current_queue);
	//������������ȫ��
	void regionGrow2(int seed, vector<int> &noise, blocks *block,
		tri_index *tri, BallCenter *bc, vector<int> &stop_point, vector<int> &current_queue);
	
	bool judgeExistenceOfPoints(int a, int b, int c,tri_index *tri);
	//�ж�Բ�ļ���
	void judgeSetsOfBallCenter(BallCenter *bc ,double Circumball[][3], int seed,
		vector<PointXYZ> &tri);
	//�ж�ͣפ��
	bool judgeStopPoints(BallCenter *bc,double Circumball[][3],int ind);
	bool judgeStopPoints2(BallCenter *bc, double Circumball[][3],
		vector<int> ind);
	bool judgeOtherPoints(double Circumball[][3],vector<int> &indexs,double radius);

 public:

	GammaShape(int mystride,string filename, 
		string format, double nbeta, int nk,float threshold);
	//��ʼ����gamma shape
	void creatGammaShape();
	void showPointCloud();

};
