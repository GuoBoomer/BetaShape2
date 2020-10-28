#pragma once
#include<iostream>
#include<string>
#include<sstream>
#include<pcl/io/pcd_io.h>
#include<pcl/io/ply_io.h>
#include<pcl/point_types.h>
#include<pcl/visualization/cloud_viewer.h>
#include <pcl/kdtree/kdtree_flann.h>
#include<vector>
#include <math.h>
#include<cmath>


using namespace std;
using namespace pcl;

class CreatBetaShape
{
private:
	double rates = 5;//����
	int kneighbor = 10; //���������k����
	vector<PointXYZRGB> tri_indexs; //�洢��Ƭ����
	//�ж��Ƿ��Ѿ�������Щ��
	bool judegExistenceOfPoints(int i1, int i2, int i3);
	
	
public:
	double beta = 10; //betaֵ�����뾶
	double upperbound;	//�Ͻ�
	
	PointCloud<PointXYZRGB>::Ptr cloudA = NULL; //�������
	string pcname;	//��������
	string format;	//��ʽ
	boost::shared_ptr<visualization::PCLVisualizer> viewer;	//���ӻ�����
	KdTreeFLANN<PointXYZRGB> kdtree;

	//���Ҿ�������Ϊindex�ĵ�������k����
	bool findKnearest(int index, vector<int>& indexs);
	//�ж�indexs���Ƿ��������index
	bool existOfIndexs(int index, vector<int>& indexs);
	//�ж�һ����ĸ�˹ӳ��
	bool judgewithGaussianMap(vector<int> indexs, PointXYZRGB p);
	//�ж��Ƿ��Ѿ�������Щ��
	//bool judegExistenceOfPoints(int i1, int i2, int i3);
	//Ѱ���������
	bool findCircumball(double ch[][3], vector<PointXYZRGB> &pd);
	
	bool judegeGuassian(vector<PointXYZRGB> &normals);

	CreatBetaShape(string filename, string nformat, double nbeta, int nk);
	void showPointCloud();	//���ӻ�cloudA
	void creatBetaShape();	//����beta shape



};