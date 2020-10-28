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
	double rates = 5;//倍率
	int kneighbor = 10; //距离最近的k个点
	vector<PointXYZRGB> tri_indexs; //存储面片索引
	//判断是否已经存在这些点
	bool judegExistenceOfPoints(int i1, int i2, int i3);
	
	
public:
	double beta = 10; //beta值，即半径
	double upperbound;	//上界
	
	PointCloud<PointXYZRGB>::Ptr cloudA = NULL; //读入点云
	string pcname;	//点云名称
	string format;	//格式
	boost::shared_ptr<visualization::PCLVisualizer> viewer;	//可视化窗体
	KdTreeFLANN<PointXYZRGB> kdtree;

	//查找距离索引为index的点的最近的k个点
	bool findKnearest(int index, vector<int>& indexs);
	//判断indexs里是否存在索引index
	bool existOfIndexs(int index, vector<int>& indexs);
	//判断一个点的高斯映射
	bool judgewithGaussianMap(vector<int> indexs, PointXYZRGB p);
	//判断是否已经存在这些点
	//bool judegExistenceOfPoints(int i1, int i2, int i3);
	//寻找外接球心
	bool findCircumball(double ch[][3], vector<PointXYZRGB> &pd);
	
	bool judegeGuassian(vector<PointXYZRGB> &normals);

	CreatBetaShape(string filename, string nformat, double nbeta, int nk);
	void showPointCloud();	//可视化cloudA
	void creatBetaShape();	//创建beta shape



};