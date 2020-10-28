#include "CreatBetaShape.h"

CreatBetaShape::CreatBetaShape(string filename, string nformat, double nbeta, int nk)
{	//文件读入
	beta = nbeta;
	kneighbor = nk + 1;
	pcname = filename;
	format = nformat;
	string path = "./data/";
	path = path + filename + '.' + nformat; //文件路径
	PointCloud<PointXYZRGB>::Ptr cloud2(new PointCloud<PointXYZRGB>);

	if (format == "ply") {
		if (io::loadPLYFile<PointXYZRGB>(path, *cloud2))
			PCL_ERROR("文件读入失败，请检查目录\n");

	}
	else if (format == "pcd") {
		if (io::loadPCDFile<PointXYZRGB>(path, *cloud2))
			PCL_ERROR("文件读入失败，请检查目录\n");

	}
	else {
		cout << "请输入正确的点云格式（pcd或者ply）" << endl;
		return;
	}
	cloudA = cloud2;
	boost::shared_ptr<visualization::PCLVisualizer> viewer1(new visualization::PCLVisualizer("Green"));
	viewer = viewer1;
	kdtree.setInputCloud(cloudA);
}






//判断lines里面是否有本条直线,没有就返回true
struct line {
	int a = -1;
	int b = -1;

};
bool existOfLines(vector<line> &lines, int a, int b) {
	if (a > b) {
		int tem = a;
		a = b;
		b = a;
	}

	line l;
	l.a = a;
	l.b = b;

	if (lines.size() == 0) {
		lines.push_back(l);
		return true;
	}

	for (int i = 0; i < lines.size(); i++)
		if (a == lines[i].a&&b == lines[i].b)
			return false;
	lines.push_back(l);
	return true;

}




void CreatBetaShape::showPointCloud()
{
	visualization::PointCloudColorHandlerCustom<PointXYZRGB> singal_color(cloudA, 255, 255, 255);

	//存储有哪些边

	vector<line> lines;
	viewer->addPointCloud<PointXYZRGB>(cloudA, singal_color, "green");
	int linenumber = 0;

	for (int i = 0; i < tri_indexs.size(); i++) {
		if (existOfLines(lines, tri_indexs[i].x, tri_indexs[i].y)) {
			viewer->addLine(cloudA->points[tri_indexs[i].x], cloudA->points[tri_indexs[i].y], to_string(linenumber));
			linenumber++;
		}
		if (existOfLines(lines, tri_indexs[i].x, tri_indexs[i].z)) {
			viewer->addLine(cloudA->points[tri_indexs[i].x], cloudA->points[tri_indexs[i].z], to_string(linenumber));
			linenumber++;
		}
		if (existOfLines(lines, tri_indexs[i].y, tri_indexs[i].z)) {
			viewer->addLine(cloudA->points[tri_indexs[i].y], cloudA->points[tri_indexs[i].z], to_string(linenumber));
			linenumber++;
		}
		if (i % 10000 == 0)
			cout << "已完成面片：" << i << endl;
	}
	cout << "线段总数：" << linenumber << endl;


	while (!viewer->wasStopped()) {

		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(1000));
	}
}


bool CreatBetaShape::judegExistenceOfPoints(int i1, int i2, int i3) {
	int tem[3];
	tem[0] = i1, tem[1] = i2, tem[2] = i3;
	sort(tem, tem + 3);
	PointXYZRGB p;
	p.x = tem[0];
	p.y = tem[1];
	p.z = tem[2];

	//如果是第一个 则直接加入
	if (tri_indexs.size() == 0) {
		tri_indexs.push_back(p);
		return true;
	}
	//如果已经存在，则不再添加
	for (int i = 0; i < tri_indexs.size(); i++)
		if (p.x == tri_indexs[i].x&&p.y == tri_indexs[i].y&&p.z == tri_indexs[i].z)
			return false;
	tri_indexs.push_back(p);

	return true;
}


void CreatBetaShape::creatBetaShape()
{
	vector<int> noise;	//噪点索引
	long int tri_num = 0;


	for (int i = 0; i < cloudA->size(); i++) {

		vector<int> indexs;		//存储索引

		//找最近的K个点的索引,并更新beta
		if (!findKnearest(i, indexs))
			continue;
		if (indexs.size() < 2)			//不够连接数
			continue;

		//判断高斯映射,不符合条件跳出循环
		if (judgewithGaussianMap(indexs, cloudA->points[i])) {
			noise.push_back(i);
			cout << "Noise Points(near the main body) Found！Its index is " << i << endl;
			continue;
		}

		//开始找点
		for (int a = 0; a < indexs.size() - 1; a++)
			for (int b = a; b < indexs.size(); b++) {
				double heart_of_Circumball[2][3]; //两个外接球球心
				vector<PointXYZRGB> points_of_tri;

				points_of_tri.push_back(cloudA->points[indexs[i]]);	//三个点的坐标
				points_of_tri.push_back(cloudA->points[indexs[a]]);
				points_of_tri.push_back(cloudA->points[indexs[b]]);

				//如果外接圆圆心大于beta,则返回true,进入下一次循环
				//findCircumball(heart_of_Circumball, points_of_tri);
				if (findCircumball(heart_of_Circumball, points_of_tri))
					continue;
				//计算距离
				bool flag_of_dis = true;
				double dis1;
				double dis2;
				for (int c = 0; c < indexs.size(); c++) {
					if (c == a || c == b|| c==i)
						continue;
					//到两个圆心的距离
					//到两个圆心的距离
					dis1 = pow(heart_of_Circumball[0][0] - cloudA->points[indexs[c]].x, 2) +
						pow(heart_of_Circumball[0][1] - cloudA->points[indexs[c]].y, 2) +
						pow(heart_of_Circumball[0][2] - cloudA->points[indexs[c]].z, 2);
					dis2 = pow(heart_of_Circumball[1][0] - cloudA->points[indexs[c]].x, 2) +
						pow(heart_of_Circumball[1][1] - cloudA->points[indexs[c]].y, 2) +
						pow(heart_of_Circumball[1][2] - cloudA->points[indexs[c]].z, 2);
					dis1 = pow(dis1, 0.5);
					dis2 = pow(dis2, 0.5);
					if (dis1 < beta || dis2 < beta) {
						flag_of_dis = false;
						break;

					}

				}

				//符合标准
				if (flag_of_dis && judegExistenceOfPoints(i, indexs[a], indexs[b]))
					tri_num++;

			}

	}
	cout << "tri number is:" << tri_num << endl;
	cout << "Noise points (near the main body) number is " << noise.size() << endl;

}



bool  CreatBetaShape::findCircumball(double ch[][3], vector<PointXYZRGB> &pd) {

	double a1, b1, c1, d1;
	double a2, b2, c2, d2;
	double a3, b3, c3, d3;

	double x1 = pd[0].x, y1 = pd[0].y, z1 = pd[0].z;
	double x2 = pd[1].x, y2 = pd[1].y, z2 = pd[1].z;
	double x3 = pd[2].x, y3 = pd[2].y, z3 = pd[2].z;

	a1 = (y1*z2 - y2 * z1 - y1 * z3 + y3 * z1 + y2 * z3 - y3 * z2);
	b1 = -(x1*z2 - x2 * z1 - x1 * z3 + x3 * z1 + x2 * z3 - x3 * z2);
	c1 = (x1*y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
	d1 = -(x1*y2*z3 - x1 * y3*z2 - x2 * y1*z3 + x2 * y3*z1 + x3 * y1*z2 - x3 * y2*z1);

	a2 = 2 * (x2 - x1);
	b2 = 2 * (y2 - y1);
	c2 = 2 * (z2 - z1);
	d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

	a3 = 2 * (x3 - x1);
	b3 = 2 * (y3 - y1);
	c3 = 2 * (z3 - z1);
	d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;

	//外接圆圆心
	double centerpoint[3];
	centerpoint[0] = -(b1*c2*d3 - b1 * c3*d2 - b2 * c1*d3 + b2 * c3*d1 + b3 * c1*d2 - b3 * c2*d1) / (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);
	centerpoint[1] = (a1*c2*d3 - a1 * c3*d2 - a2 * c1*d3 + a2 * c3*d1 + a3 * c1*d2 - a3 * c2*d1) / (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);
	centerpoint[2] = -(a1*b2*d3 - a1 * b3*d2 - a2 * b1*d3 + a2 * b3*d1 + a3 * b1*d2 - a3 * b2*d1) / (a1*b2*c3 - a1 * b3*c2 - a2 * b1*c3 + a2 * b3*c1 + a3 * b1*c2 - a3 * b2*c1);

	//半径的平方
	double radius = pow((centerpoint[0] - x1), 2) + pow((centerpoint[1] - y1), 2) + pow((centerpoint[2] - z1), 2);
	if (pow(radius,0.5) > beta)
		return true;
	//圆心到球心的距离
	double dis = pow(pow(beta, 2) - radius, 0.5);
	//求单位法向量
	double hx = ((y2 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1));
	double hy = ((z2 - z1)*(x3 - x1) - (x2 - x1)*(z3 - z1));
	double hz = ((x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1));
	double tem = pow(pow(hx, 2) + pow(hy, 2) + pow(hz, 2), 0.5);
	hx = hx / tem;
	hy = hy / tem;
	hz = hz / tem;

	//计算球心,出现很小的点
	ch[0][0] = centerpoint[0] + dis * hx; 
	ch[0][1] = centerpoint[1] + dis * hy; 
	ch[0][2] = centerpoint[2] + dis * hz;
	ch[1][0] = centerpoint[0] - dis * hx; 
	ch[1][1] = centerpoint[1] - dis * hy; 
	ch[1][2] = centerpoint[2] - dis * hz;
	return false;
}

bool CreatBetaShape::judgewithGaussianMap(vector<int> indexs, PointXYZRGB p)
{
	vector<PointXYZRGB> normals;	//存法向

	//返回ture 表明该点为噪点
	for (int i = 0; i < indexs.size(); i++) {
		PointXYZRGB a;
		//计算向量，原点为（0，0，0）
		a.x = p.x - cloudA->points[indexs[i]].x;
		a.y = p.y - cloudA->points[indexs[i]].y;
		a.z = p.z - cloudA->points[indexs[i]].z;

		//单位化
		double dis = pow(pow(a.x, 2) + pow(a.y, 2) + pow(a.z, 2), 0.5);
		if (dis == 0)
			continue;
		a.x = a.x / dis;
		a.y = a.y / dis;
		a.z = a.z / dis;

		normals.push_back(a);

	}

	return judegeGuassian(normals);
}


/*
	超过45°大小的夹角达到三次以上，则认为该点不为噪声点
*/
bool CreatBetaShape::judegeGuassian(vector<PointXYZRGB>& normals)
{

	double cos45 = cos(3.1415926535 * 135 / 180);	// cos45的值
	int count = 0;		//大于45°范围的数量

	for (int i = 0; i < normals.size() - 1; i++) {
		for (int a = i + 1; a < normals.size(); a++) {
			double apb = normals[i].x*normals[a].x + normals[i].y*normals[a].y + normals[i].z*normals[a].z;
			if (cos45 > apb) {
				count++;
				break;
			}

		}

	}
	//count大于3，表示有若干角度很大的点，证明该点为非噪声点
	if (count > 3)
		return false;

	return true;
}


bool CreatBetaShape::findKnearest(int index, vector<int>& indexs)
{
	if (cloudA == NULL)
		return false;

	int cloudsize = cloudA->size();//点云大小


	PointXYZRGB searchPoint = cloudA->points[index];

	std::vector<int> pointIdxNKNSearch(kneighbor);
	std::vector<float> pointNKNSquaredDistance(kneighbor);

	if (kdtree.nearestKSearch(searchPoint, kneighbor, pointIdxNKNSearch, pointNKNSquaredDistance) == 0)
		return false;

	float dis = 0;
	for (int i = 0; i < pointNKNSquaredDistance.size(); i++)
		dis += pow(pointNKNSquaredDistance[i], 0.5);

	//更新beta，设置为k邻域内的距离均值，-1为消除其自身的影响
	if (index == 0)
		upperbound = dis / (pointNKNSquaredDistance.size() - 1);

	beta = dis / (pointNKNSquaredDistance.size() - 1);

	if (beta < upperbound)
		upperbound = beta;
	if (beta > upperbound * rates)
		return false;

	indexs = pointIdxNKNSearch;
	vector<int>::iterator k = indexs.begin(); //移除其自身
	indexs.erase(k);
	return true;

}

bool CreatBetaShape::existOfIndexs(int index, vector<int>& indexs)
{
	if (indexs.size() == 0)
		return false;
	for (int i = 0; i < indexs.size(); i++)
		if (index == indexs[i])
			return true;
	return false;
}

