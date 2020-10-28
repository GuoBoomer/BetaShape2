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
	blocks *blocks_head = new blocks; //�ֿ�ͷָ��
	tri_index *tris_head = new tri_index;	//ÿһ���Ӧ����Ƭ
	BallCenter *bc_head = new BallCenter;	//��������


	blocks *tem_block=blocks_head; //�ݴ�ָ��
	tri_index *tem_tris=tris_head;
	BallCenter *tem_bc = bc_head;	//��������

	vector<int> noise;	//������

	//�����������еĵ�
	for (int i = 0; i < cloudA->points.size(); i++) {
		if (existOfIndexs(i, visited))
			continue;

		vector<int> current_queue;	//������ȱ�������
		blocks *block = new blocks; 
		tri_index *tri = new tri_index;
		BallCenter *bc = new BallCenter;
		vector<int> stop_point;		//ͣפ������

		block->indexs.push_back(i); //?
		int current;
		current_queue.push_back(i);
		while(current_queue.size()!=0){
			current = current_queue[0];
			vector<int>::iterator k = current_queue.begin(); //�����Ƴ�
			current_queue.erase(k);
			//��������
			//regionGrow(current, noise,block,tri,bc,stop_point,current_queue);
			if (!existOfIndexs(current, visited)) {
				regionGrow2(current,noise,block,tri,bc,stop_point,current_queue);
			    visited.push_back(current);
			}
		
		}
		//ͣפ��洢
		if(stop_point.size()!=0)
			stop_points.insert(stop_points.end(), stop_point.begin(), stop_point.end());

		
		//��������
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

	cout << "ͣפ�㳤��:" << stop_points.size() << endl;
	cout << "��Ƭ�����鳤�ȣ�" << anum << endl;
	cout <<"������Ƭ���鳤�ȣ�"<< cnum << endl;
	cout << "�������������" << visited.size() << endl;


	


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

	//�洢����Щ��
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
{	//�ж��Ƿ��Ѿ����ʹ�
	if (existOfIndexs(seed, visited))
		return;
	visited.push_back(seed);//��Ҫ�ж��Ƿ�Ϊͣפ��
	vector<int> indexs;
	
	
	if (!findKnearest(seed, indexs))
		return;
	
	//�жϸ�˹ӳ��,��������������ѭ��
	/*if (judgewithGaussianMap(indexs, cloudA->points[seed])) {
		noise.push_back(seed);
		return;
	}*/

	//��ʼ�ҵ�
	for (int a = 0; a < indexs.size() - 1; a++)
		for (int b = a+1; b < indexs.size(); b++) {
			double heart_of_Circumball[2][3]; //�������������
			vector<PointXYZRGB> points_of_tri;

			points_of_tri.push_back(cloudA->points[seed]);	//�����������
			points_of_tri.push_back(cloudA->points[a]);
			points_of_tri.push_back(cloudA->points[b]);

			//������ԲԲ�Ĵ���beta,�򷵻�true,������һ��ѭ��
			if (findCircumball(heart_of_Circumball, points_of_tri))
				continue;
			//�������
			bool flag_of_dis = true;
			double dis1;
			double dis2;
			for (int c = 0; c < indexs.size(); c++) {
				if (c == a || c == b)
					continue;
				//������Բ�ĵľ���,Բ�ĵ�ͬһ��ľ���
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

			//���ϱ�׼
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
	//��Ϊ�����ӵ�����ѷ���
	
	vector<int> indexs;


	if (!findKnearest(seed, indexs))
		return;

	//�жϸ�˹ӳ��,��������������ѭ��
	/*if (judgewithGaussianMap(indexs, cloudA->points[seed])) {
		noise.push_back(seed);
		return;
	}*/

	for (int a = 0; a < indexs.size() - 1; a++)
		for (int b = a + 1; b < indexs.size(); b++) {
			
			

			double heart_of_Circumball[2][3]; //�������������
			vector<PointXYZRGB> points_of_tri;

			points_of_tri.push_back(cloudA->points[seed]);	//�����������
			points_of_tri.push_back(cloudA->points[indexs[a]]);
			points_of_tri.push_back(cloudA->points[indexs[b]]);

			//������ԲԲ�Ĵ���beta,�򷵻�true,������һ��ѭ������������жϣ����������Բ�İ뾶��Ϊ�����İ뾶
			if (findCircumball(heart_of_Circumball, points_of_tri))
				continue;
			//�������
			bool flag_of_dis = true;
			double dis1;
			double dis2;
			//�ж��������Ƿ�������
			for (int c = 0; c < indexs.size(); c++) {
				if (c == a || c == b|| c== seed)
					continue;
				//������Բ�ĵľ���
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
					//break; //Ϊ��Ҫbreak?
				}
			}
			if (!flag_of_dis)
				continue;
			//���ĸ�������Ƭ���乲��
			int length = -1;
			if(tri->tir_points.size()>0)
				for (int c = 0; c < tri->tir_points.size(); c++) {
					if (tri->tir_points[c].x == seed || tri->tir_points[c].y == seed || tri->tir_points[c].z == seed) {
						length = c;
						break;
					}
			
				}

			//�ж��Ƿ���ͣפ��, �Ƿ���Ҫ����return
			if (judgeStopPoints(bc, heart_of_Circumball, length)) {
				if ((!existOfIndexs(seed, stop_points)) && (!existOfIndexs(seed, stop_point))) {
					stop_point.push_back(seed);
					return;
				}
			}
			

			//���ϱ�׼ �Ƿ��Ѿ�������ͬ����Ƭ
			if (judgeExistenceOfPoints(seed, indexs[a], indexs[b], tri)) {
				//����ͬ��Բ�ļ��벻ͬ�ļ�����
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




//ind Ӧ������tri_index�е����� ������CloudA�������
bool GammaShape::judgeStopPoints(BallCenter *bc,  double Circumball[][3],int ind)
{	

	//��һ���� ����û�ҵ������� ����false ��ʾ����Stop Points
	if (bc->set1.tir_points.size() == 0||ind==-1)
		return false;
	

	//������Խ�� ����true ��ΪStop Points ����
	if (bc->set1.tir_points.size() <= ind) {
		cout << "judgeStopPoints Խ��:"<<ind<<" "<<bc->set1.tir_points.size()<<endl;
		return false;
	}


	//��Ӧ�������
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
	//��ͬһ�����ĵľ���
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

	//�㼯�ָ�
	if (min == d[1] || min == d[3]) {
		tem3 = heart2;
		heart2 = heart1;
		heart1 = tem3;
	}
	
	//�����߶�
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


	//С����ֵ����ͣפ��,����Ϊͣפ��
	if (scale1 / scale2 > Threshold) {
		return false;
	}
	else{
		//cout << "����ͣפ��" << endl;
		return true;
	}
}

bool GammaShape::judgeStopPoints2(BallCenter *bc, double Circumball[][3], vector<int> ind) {
//
//	��һ���� ����false ��ʾ����stop points
//	if (bc->set1.tir_points.size() == 0 || ind == 0)
//		return false;
//
//
//	������Խ�� ����true ��Ϊstop points ����
//	if (bc->set1.tir_points.size() < ind) {
//		cout << "judgestoppoints Խ��" << endl;
//		return true;
//	}
//
//
//	��Ӧ�������
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
//	��ͬһ�����ĵľ���
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
//	�㼯�ָ�
//	if (min == d[1] || min == d[3]) {
//		tem3 = heart2;
//		heart2 = heart1;
//		heart1 = heart2;
//	}
//
//	�����߶�
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
//	С����ֵ����ͣפ��,����Ϊͣפ��
//	if (scale1 / scale2 < threshold) {
//		return false;
//	}
//	else {
//		cout << "����ͣפ��" << endl;
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
{	//��û�е㼯����
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

	//���û�л�û�����ļ��룬��������뵽����������
	if (bc->set1.tir_points.size() == 0) {
		bc->set1.tir_points.push_back(heart1);
		bc->set2.tir_points.push_back(heart2);
		return;

	}

	//Ŀǰֻ��һ����Ƭ
	vector<int> index;
	int ind;
	for (int i = 0; i < tri.size(); i++) {
		if (tri[i].x == seed || tri[i].y == seed || tri[i].z == seed) {
			index.push_back(i);
			ind = i;
			break;
		}
	}

	//�������ĵ���������
	float d[4];
	tem = bc->set1.tir_points[ind];
	tem2 = bc->set2.tir_points[ind];
	//��ͬһ�����ĵľ���
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
	

	//�㼯�ֿ�
	if (min == d[0]||min==d[2]) {
		bc->set1.tir_points.push_back(heart1);
		bc->set2.tir_points.push_back(heart2);
	}
	else if (min == d[1]||min==d[3]) {
		bc->set1.tir_points.push_back(heart2);
		bc->set2.tir_points.push_back(heart1);
	}
	

}


