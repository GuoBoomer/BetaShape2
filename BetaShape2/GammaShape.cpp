#include "GammaShape.h"

GammaShape::GammaShape(int mystride, string filename, string format, double nbeta, int nk,float threshold)
	:CreatBetaShape(filename,format,nbeta,nk){
	stride = mystride;
	Threshold = threshold;
	int r = 10;
	for (int i = 0; i < cloudA->points.size(); i++) {
		cloudA->points[i].x = cloudA->points[i].x * r;
		cloudA->points[i].y = cloudA->points[i].y * r;
		cloudA->points[i].z = cloudA->points[i].z * r;

	}
}

void GammaShape::creatGammaShape()
{
	blocks *blocks_head = new blocks; //分块头指针
	tri_index *tris_head = new tri_index;	//每一块对应的面片
	BallCenter *bc_head = new BallCenter;	//球心索引


	blocks *tem_block=blocks_head; //暂存指针
	tri_index *tem_tris=tris_head;
	BallCenter *tem_bc = bc_head;	//球心索引

	vector<int> noise;	//存噪声

	//遍历点云所有的点
	for (int i = 0; i < cloudA->points.size(); i++) {
		if (existOfIndexs(i, visited))
			continue;

		vector<int> current_queue;	//广度优先遍历队列
		blocks *block = new blocks; 
		tri_index *tri = new tri_index;
		BallCenter *bc = new BallCenter;
		vector<int> stop_point;		//停驻点索引

		block->indexs.push_back(i); //?
		int current;
		current_queue.push_back(i);
		while(current_queue.size()!=0){
			current = current_queue[0];
			vector<int>::iterator k = current_queue.begin(); //队首移除
			current_queue.erase(k);
			//区域生长
			//regionGrow(current, noise,block,tri,bc,stop_point,current_queue);
			if (!existOfIndexs(current, visited)) {
				regionGrow2(current,noise,block,tri,bc,stop_point,current_queue);
			    visited.push_back(current);
			}
		
		}
		//停驻点存储
		if(stop_point.size()!=0)
			stop_points.insert(stop_points.end(), stop_point.begin(), stop_point.end());

		
		//链表连接
		if (block->indexs.size() != 0) {
		tem_block->next = block;
		tem_block = block;
		}
		if (tri->tir_points.size() != 0) {
		tem_tris->next = tri;
		tem_tris = tri;
		}
		if (bc->set1.tir_points.size() != 0) {
		tem_bc->next = bc;
		tem_bc = bc;
		}
	}


	tem_tris = tris_head->next;
	tem_block = blocks_head;
	

	int cnum = 0;
	int anum = 0;

	while (tem_block != NULL) {
		tem_block = tem_block->next;
		anum++;
	}

	while (tem_tris != NULL) {
		tem_tris = tem_tris->next;
		cnum++;
	}

	cout << "停驻点长度:" << stop_points.size() << endl;
	cout << "分片点区块长度：" << anum << endl;
	cout <<"三角面片区块长度："<< cnum << endl;
	cout << "遍历点的数量：" << visited.size() << endl;


	


	showPointCloud();
	return;

	tem_block = blocks_head->next;
	while (tem_block != NULL) {
		if (tem_block->indexs.size() > 1) {
			for (int i = 0; i < tem_block->indexs.size(); i++)
				//cout << tem_block->indexs[i]<<" "<< tem_block->indexs[i]<<" "<< tem_block->indexs[i]<<endl;
				cout << tem_block->indexs[i]<<" ";
			cout << endl<<"===================================="<<endl;
		}
		tem_block = tem_block->next;
	}


}

void GammaShape::showPointCloud()
{
	//visualization::PointCloudColorHandlerCustom<PointXYZRGB> singal_color(cloudA, 255, 255, 255);

	//存储有哪些边
	int i = 0;
	while(i< cloudA->points.size()){
		if (existOfIndexs(i, stop_points)) {
			cloudA->points[i].r = 255;
			cloudA->points[i].g = 0;
			cloudA->points[i].b = 0;
		}
		else {
			cloudA->points[i].r = 255;
			cloudA->points[i].g = 255;
			cloudA->points[i].b = 255;
		}
		i++;
	}
	
	viewer->addPointCloud<PointXYZRGB>(cloudA, "green");
	

	while (!viewer->wasStopped()) {

		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(1000));
	}
}


void GammaShape::regionGrow(int seed,vector<int> &noise,
	blocks *block,tri_index *tri,BallCenter *bc, vector<int> &stop_point,vector<int> &current_queue)
{	//判断是否已经访问过
	if (existOfIndexs(seed, visited))
		return;
	visited.push_back(seed);//还要判断是否为停驻点
	vector<int> indexs;
	
	
	if (!findKnearest(seed, indexs))
		return;
	
	//判断高斯映射,不符合条件跳出循环
	/*if (judgewithGaussianMap(indexs, cloudA->points[seed])) {
		noise.push_back(seed);
		return;
	}*/

	//开始找点
	for (int a = 0; a < indexs.size() - 1; a++)
		for (int b = a+1; b < indexs.size(); b++) {
			double heart_of_Circumball[2][3]; //两个外接球球心
			vector<PointXYZRGB> points_of_tri;

			points_of_tri.push_back(cloudA->points[seed]);	//三个点的坐标
			points_of_tri.push_back(cloudA->points[a]);
			points_of_tri.push_back(cloudA->points[b]);

			//如果外接圆圆心大于beta,则返回true,进入下一次循环
			if (findCircumball(heart_of_Circumball, points_of_tri))
				continue;
			//计算距离
			bool flag_of_dis = true;
			double dis1;
			double dis2;
			for (int c = 0; c < indexs.size(); c++) {
				if (c == a || c == b)
					continue;
				//到两个圆心的距离,圆心到同一点的距离
				dis1 = pow(heart_of_Circumball[0][0] - points_of_tri[1].x, 2) +
					pow(heart_of_Circumball[0][1] - points_of_tri[1].y, 2) +
					pow(heart_of_Circumball[0][2] - points_of_tri[1].z, 2);
				dis2 = pow(heart_of_Circumball[1][0] - points_of_tri[1].x, 2) +
					pow(heart_of_Circumball[1][1] - points_of_tri[1].y, 2) +
					pow(heart_of_Circumball[1][2] - points_of_tri[1].z, 2);
				dis1 = pow(dis1, 0.5);
				dis2 = pow(dis2, 0.5);
				if (dis1 > beta || dis2 > beta) {
					flag_of_dis = false;	
					break;
				}

			}

			//符合标准
			if (flag_of_dis && judgeExistenceOfPoints(seed, indexs[a], indexs[b], tri)) {
				if (!existOfIndexs(indexs[a], block->indexs)) {
					block->indexs.push_back(indexs[a]);
					current_queue.push_back(indexs[a]);
					if (!existOfIndexs(indexs[a], visited))
						visited.push_back(indexs[a]);
				}
				if (!existOfIndexs(indexs[b], block->indexs)) {
					block->indexs.push_back(indexs[b]);
					current_queue.push_back(indexs[b]);
					if (!existOfIndexs(indexs[b], visited))
						visited.push_back(indexs[b]);
				}
			}
			
		}

}



void GammaShape::regionGrow2(int seed, vector<int>& noise, blocks * block, tri_index * tri, BallCenter * bc,
	vector<int>& stop_point, vector<int>& current_queue)
{
	if (existOfIndexs(seed, visited))
		return;
	//成为了种子点才算已访问
	
	vector<int> indexs;


	if (!findKnearest(seed, indexs))
		return;

	//判断高斯映射,不符合条件跳出循环
	/*if (judgewithGaussianMap(indexs, cloudA->points[seed])) {
		noise.push_back(seed);
		return;
	}*/

	for (int a = 0; a < indexs.size() - 1; a++)
		for (int b = a + 1; b < indexs.size(); b++) {
			
			

			double heart_of_Circumball[2][3]; //两个外接球球心
			vector<PointXYZRGB> points_of_tri;

			points_of_tri.push_back(cloudA->points[seed]);	//三个点的坐标
			points_of_tri.push_back(cloudA->points[indexs[a]]);
			points_of_tri.push_back(cloudA->points[indexs[b]]);

			//如果外接圆圆心大于beta,则返回true,进入下一次循环，如果不做判断，将会以外接圆的半径作为外接球的半径
			if (findCircumball(heart_of_Circumball, points_of_tri))
				continue;
			//计算距离
			bool flag_of_dis = true;
			double dis1;
			double dis2;
			//判断其他点是否在球内
			for (int c = 0; c < indexs.size(); c++) {
				if (c == a || c == b|| c== seed)
					continue;
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
					//break; //为甚要break?
				}
			}
			if (!flag_of_dis)
				continue;
			//找哪个三角面片与其共点
			int length = -1;
			if(tri->tir_points.size()>0)
				for (int c = 0; c < tri->tir_points.size(); c++) {
					if (tri->tir_points[c].x == seed || tri->tir_points[c].y == seed || tri->tir_points[c].z == seed) {
						length = c;
						break;
					}
			
				}

			//判断是否有停驻点, 是否需要保持return
			if (judgeStopPoints(bc, heart_of_Circumball, length)) {
				if ((!existOfIndexs(seed, stop_points)) && (!existOfIndexs(seed, stop_point))) {
					stop_point.push_back(seed);
					return;
				}
			}
			

			//符合标准 是否已经存在相同的面片
			if (judgeExistenceOfPoints(seed, indexs[a], indexs[b], tri)) {
				//将不同的圆心加入不同的集合中
				judgeSetsOfBallCenter(bc, heart_of_Circumball, seed, tri->tir_points);
				if (!existOfIndexs(indexs[a], block->indexs)) {
					block->indexs.push_back(indexs[a]);
					current_queue.push_back(indexs[a]);
				}
				if (!existOfIndexs(indexs[b], block->indexs)) {
					block->indexs.push_back(indexs[b]);
					current_queue.push_back(indexs[b]);
				}
			}

		}
}




//ind 应该是在tri_index中的索引 不是在CloudA里的索引
bool GammaShape::judgeStopPoints(BallCenter *bc,  double Circumball[][3],int ind)
{	

	//第一个点 或者没找到公共点 返回false 表示不是Stop Points
	if (bc->set1.tir_points.size() == 0||ind==-1)
		return false;
	

	//发生了越界 返回true 作为Stop Points 处理
	if (bc->set1.tir_points.size() <= ind) {
		cout << "judgeStopPoints 越界:"<<ind<<" "<<bc->set1.tir_points.size()<<endl;
		return false;
	}


	//对应的领域点
	PointXYZ tem = bc->set1.tir_points[ind];
	PointXYZ tem2 = bc->set2.tir_points[ind];
	PointXYZ heart1,heart2,tem3;

	heart1.x = Circumball[0][0];
	heart2.x = Circumball[1][0];
	heart1.y = Circumball[0][1];
	heart2.y = Circumball[1][1];
	heart1.z = Circumball[0][2];
	heart2.z = Circumball[1][2];
	
	float d[4];
	tem = bc->set1.tir_points[ind];
	tem2 = bc->set2.tir_points[ind];
	//离同一个球心的距离
	d[0] = pow(heart1.x - tem.x, 2) +
		pow(heart1.y - tem.y, 2) +
		pow(heart1.z - tem.z, 2);
	d[1] = pow(heart2.x - tem.x, 2) +
		pow(heart2.y - tem.y, 2) +
		pow(heart2.z - tem.z, 2);

	d[2] = pow(heart1.x - tem2.x, 2) +
		pow(heart1.y - tem2.y, 2) +
		pow(heart1.z - tem2.z, 2);
	d[3] = pow(heart2.x - tem2.x, 2) +
		pow(heart2.y - tem2.y, 2) +
		pow(heart2.z - tem2.z, 2);

	float min1 = d[0] < d[1] ? d[0] : d[1];
	float min2 = d[2] < d[3] ? d[2] : d[3];
	float min = min1 < min2 ? min1 : min2;

	//点集分割
	if (min == d[1] || min == d[3]) {
		tem3 = heart2;
		heart2 = heart1;
		heart1 = tem3;
	}
	
	//两个尺度
	float scale2= pow(heart2.x - tem2.x, 2) +
		pow(heart2.y - tem2.y, 2) +
		pow(heart2.z - tem2.z, 2);
	float scale1= pow(heart1.x - tem.x, 2) +
		pow(heart1.y - tem.y, 2) +
		pow(heart1.z - tem.z, 2);

	scale1 = pow(scale1, 0.5);
	scale2 = pow(scale2, 0.5);
	if (scale1 > scale2) {
		float a = scale1;
		scale1 = scale2;
		scale2 = a;
	}


	//小于阈值，非停驻点,否则为停驻点
	if (scale1 / scale2 > Threshold) {
		return false;
	}
	else{
		//cout << "发现停驻点" << endl;
		return true;
	}
}

bool GammaShape::judgeStopPoints2(BallCenter *bc, double Circumball[][3], vector<int> ind) {
//
//	第一个点 返回false 表示不是stop points
//	if (bc->set1.tir_points.size() == 0 || ind == 0)
//		return false;
//
//
//	发生了越界 返回true 作为stop points 处理
//	if (bc->set1.tir_points.size() < ind) {
//		cout << "judgestoppoints 越界" << endl;
//		return true;
//	}
//
//
//	对应的领域点
//	pointxyz tem = bc->set1.tir_points[ind];
//	pointxyz tem2 = bc->set2.tir_points[ind];
//	pointxyz heart1, heart2, tem3;
//
//	heart1.x = circumball[0][0];
//	heart2.x = circumball[0][0];
//	heart1.y = circumball[0][1];
//	heart2.y = circumball[0][1];
//	heart1.z = circumball[0][2];
//	heart2.z = circumball[0][2];
//
//	float d[4];
//	tem = bc->set1.tir_points[ind];
//	tem2 = bc->set2.tir_points[ind];
//	离同一个球心的距离
//	d[0] = pow(heart1.x - tem.x, 2) +
//		pow(heart1.y - tem.y, 2) +
//		pow(heart1.z - tem.z, 2);
//	d[1] = pow(heart2.x - tem.x, 2) +
//		pow(heart2.y - tem.y, 2) +
//		pow(heart2.z - tem.z, 2);
//
//	d[2] = pow(heart1.x - tem2.x, 2) +
//		pow(heart1.y - tem2.y, 2) +
//		pow(heart1.z - tem2.z, 2);
//	d[3] = pow(heart2.x - tem2.x, 2) +
//		pow(heart2.y - tem2.y, 2) +
//		pow(heart2.z - tem2.z, 2);
//
//	float min1 = d[0] < d[1] ? d[0] : d[1];
//	float min2 = d[2] < d[3] ? d[0] : d[1];
//	float min = min1 < min2 ? min1 : min2;
//
//	点集分割
//	if (min == d[1] || min == d[3]) {
//		tem3 = heart2;
//		heart2 = heart1;
//		heart1 = heart2;
//	}
//
//	两个尺度
//	float scale2 = pow(heart2.x - tem2.x, 2) +
//		pow(heart2.y - tem2.y, 2) +
//		pow(heart2.z - tem2.z, 2);
//	float scale1 = pow(heart1.x - tem.x, 2) +
//		pow(heart1.y - tem.y, 2) +
//		pow(heart1.z - tem.z, 2);
//
//	scale1 = pow(scale1, 0.5);
//	scale2 = pow(scale2, 0.5);
//	if (scale1 > scale2) {
//		float a = scale1;
//		scale1 = scale2;
//		scale2 = a;
//	}
//
//
//	小于阈值，非停驻点,否则为停驻点
//	if (scale1 / scale2 < threshold) {
//		return false;
//	}
//	else {
//		cout << "发现停驻点" << endl;
//		return true;
//	}
//
//
//
}

bool GammaShape::judgeOtherPoints(double Circumball[][3], vector<int> &indexs,double radius)
{	
	PointXYZ heart;
	heart.x = Circumball[0][0];
	heart.y = Circumball[0][1];
	heart.z = Circumball[0][2];

	for (int i = 0; i < indexs.size(); i++) {
		double dis = pow(heart.x - cloudA->points[indexs[i]].x, 2) +
					pow(heart.y - cloudA->points[indexs[i]].y, 2) +
					pow(heart.z - cloudA->points[indexs[i]].z, 2);
						
	
	}

	return false;
}


bool GammaShape::judgeExistenceOfPoints(int a, int b, int c,tri_index *tri)
{
	int tem[3];
	tem[0] = a, tem[1] = b, tem[2] = c;
	sort(tem, tem + 3);
	PointXYZ p;
	p.x = tem[0];
	p.y = tem[1];
	p.z = tem[2];
	if (tri->tir_points.size() == 0) {
		tri->tir_points.push_back(p);
		return true;
	}
	for (int i = 0; i < tri->tir_points.size(); i++)
		if (p.x == tri->tir_points[i].x&&p.y == tri->tir_points[i].y&&p.z == tri->tir_points[i].z)
			return false;
	tri->tir_points.push_back(p);
	return true;
}



void GammaShape::judgeSetsOfBallCenter(BallCenter * bc, double Circumball[][3], int seed, vector<PointXYZ> &tri)
{	//还没有点集进入
	PointXYZ tem;
	PointXYZ tem2;
	PointXYZ heart1;
	PointXYZ heart2;

	heart1.x = Circumball[0][0];
	heart2.x = Circumball[1][0];
	heart1.y = Circumball[0][1];
	heart2.y = Circumball[1][1];
	heart1.z = Circumball[0][2];
	heart2.z = Circumball[1][2];

	//如果没有还没有球心加入，则随机加入到两个集合中
	if (bc->set1.tir_points.size() == 0) {
		bc->set1.tir_points.push_back(heart1);
		bc->set2.tir_points.push_back(heart2);
		return;

	}

	//目前只找一个面片
	vector<int> index;
	int ind;
	for (int i = 0; i < tri.size(); i++) {
		if (tri[i].x == seed || tri[i].y == seed || tri[i].z == seed) {
			index.push_back(i);
			ind = i;
			break;
		}
	}

	//两组球心的两两距离
	float d[4];
	tem = bc->set1.tir_points[ind];
	tem2 = bc->set2.tir_points[ind];
	//离同一个球心的距离
	d[0] = pow(heart1.x - tem.x, 2) +
		pow(heart1.y - tem.y, 2) +
		pow(heart1.z - tem.z, 2);
	d[1] = pow(heart2.x - tem.x, 2) +
		pow(heart2.y - tem.y, 2) +
		pow(heart2.z - tem.z, 2);
	
	d[2] = pow(heart1.x - tem2.x, 2) +
		pow(heart1.y - tem2.y, 2) +
		pow(heart1.z - tem2.z, 2);
	d[3] = pow(heart2.x - tem2.x, 2) +
		pow(heart2.y - tem2.y, 2) +
		pow(heart2.z - tem2.z, 2);

	float min1 = d[0] < d[1] ? d[0] : d[1];
	float min2 = d[2] < d[3] ? d[2] : d[3];
	float min = min1 < min2 ? min1 : min2;
	

	//点集分开
	if (min == d[0]||min==d[2]) {
		bc->set1.tir_points.push_back(heart1);
		bc->set2.tir_points.push_back(heart2);
	}
	else if (min == d[1]||min==d[3]) {
		bc->set1.tir_points.push_back(heart2);
		bc->set2.tir_points.push_back(heart1);
	}
	

}


