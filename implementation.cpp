#include <math.h>
#include <limits.h>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <ilcplex/ilocplex.h>
#include "predefine.h"
#include <ctime>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "float.h"


using namespace std;

/// Point 类的实现
Point::Point() :T(type::POINT) {
	loc[0] = 0;
	loc[1] = 0;
}

Point::Point(double* l, type type) : T(type) {
	for (int i = 0; i < 2; i++)
		loc[i] = l[i];
}

const double* Point::get_loc() const
{
	return this->loc;
}

const int Point::get_index()
{
	return this->index;
}

double Point::cal_distance(const Point& p)
{
	const double* loc1 = this->loc;
	const double* loc2 = p.get_loc();
	double dis = sqrt(pow(loc1[0] - loc2[0], 2) + pow(loc1[1] - loc2[1], 2));
	return dis;
}

void Point::set_location(double lx, double ly)
{
	loc[0] = lx;
	loc[1] = ly;
}

bool Point::operator<(const Point& p) const {
	const double* a = get_loc();
	const double* b = p.get_loc();

	double al = sqrt(pow(a[0], 2) + pow(a[1], 2));
	double bl = sqrt(pow(b[0], 2) + pow(b[1], 2));

	double cosa = 0;
	double cosb = 0;

	if (al != 0) cosa = a[1] / al;
	if (bl != 0) cosb = b[1] / bl;


	if (cosa < cosb)
		return true;
	else if (cosa == cosb && al < bl)
		return true;
	return false;
}


/// <summary>
/// Device 类的实现
/// </summary>
Device::Device(double * location, int i, int da, const Resource& t, double eb, double em)
{
	index = i;
	loc[0] = location[0];
	loc[1] = location[1];
	data = da;
	taskdemand.copy(t);
	T = type::DEVICE;
	e_unitbit = eb;
	e_unitm = em;
}

void Device::print_device()
{
	std::cout << "x=" << loc[0] << " y=" << loc[1] << " data="<<data <<" e_unitbit="<<e_unitbit<<" e_unitm" << ' ';
	taskdemand.print_();
	std::cout << "\n";
}

/// <summary>
/// Disk 类的实现
/// </summary>

double Disk::theta = 2;
Disk::Disk() {
	device = -1;
	base = -1;
	power = 0;
}

Disk::Disk(int base_, int device_)
{
	device = device_;
	base = base_;
	power = 0;
}

Disk::Disk(int base_, int device_, const double& radii, const Resource& capacity_)
{
	device = device_;
	base = base_;
	power = cal_power(radii);
	capacity.copy(capacity_);
}

double Disk::cal_power(const double& radii)
{
	double power = C * pow(radii, theta);		
	return power;
}


void Disk::cover_devices(const std::vector<int>& db, const int& i)
{
	devices_covered_by_disk = db;
	cover_devices(i);
}

void Disk::print_disk() const		
{
	std::cout << std::setprecision(2) << "base station" << base << " device" << device << " power" << power;
	capacity.print_();
	for (int i = 0; i < devices_covered_by_disk.size(); i++)
		std::cout << devices_covered_by_disk[i]<<" ";
	std::cout << "\n";
}

void Disk::remove_covered_devices(std::vector<int>& devs)
{
	for (auto dev : devs)
		devices_covered_by_disk.erase(remove_if(devices_covered_by_disk.begin(), devices_covered_by_disk.end(), [dev](int device_covered_by_disk) {return device_covered_by_disk == dev; }), devices_covered_by_disk.end());
}

/// <summary>
/// BaseStation 类的实现
/// </summary>

BaseStation::BaseStation(double* location, Resource& c, int b, double f_, double pf_)
{
	loc[0] = location[0];
	loc[1] = location[1];
	capacity.copy(c);
	T = type::BASESTATION;
	weight = 0;  
	za = -1;     
	index = b;
	basef = f_;
	basepf = pf_;
}

BaseStation::BaseStation()
{
	loc[0] = 0;
	loc[1] = 0;
	T = type::BASESTATION;
	weight = 0;
	za = -1;
	basef = 0;
	basepf = 0;
}

void BaseStation::print_basestation()
{
	std::cout << "x=" << loc[0] << " y=" << loc[1] << " basef=" << basef << " basepf ="<< basepf<< ' ';
	capacity.print_();
	std::cout << "\n";
}

/// <summary>
/// Cloud 类的实现
/// </summary>

void Cloud::set_cloud(Resource& c, double fvm_, double pfvm_)
{
	capacity.copy(c);
	cloudfvm = fvm_;
	cloudpfvm = pfvm_;
}

Cloud::Cloud()
{
	cloudfvm = 0;
	cloudpfvm = 0;
}

/// <summary>
/// Cover 类的实现
/// </summary>
Cover::Cover() 
{
	m = 0;
	n = 0;

	E_base_ex = nullptr;  
	E_base_tr = nullptr;
	E_cloud_ex = nullptr;
	E_cloud_tr = nullptr;
	E_cloud_sum_tr = nullptr;

	distance = nullptr;
	dist_index = nullptr;
	device_order = nullptr;

	alpha_i = nullptr; 
	beta_id = nullptr;   
	gamma_id = nullptr; 
	delta_d = nullptr; 
	epsilon_d = nullptr;
	mu_b = nullptr;
}

Cover::Cover(int m_, int n_)
{
	m = m_;
	n = n_;
	all_disk.resize(m);  

	E_base_ex = new double* [m];
	E_base_tr = new double* [m];
	E_cloud_ex = new double [n];
	E_cloud_tr = new double [n];
	E_cloud_sum_tr = new double* [m];

	distance = new double* [m];
	dist_index = new int* [m];
	device_order = new int* [m];

	for (int b = 0; b < m; b++) {
		E_base_ex[b] = new double[n];
		E_base_tr[b] = new double[n];
		E_cloud_sum_tr[b] = new double[n];

		distance[b] = new double[n];
		dist_index[b] = new int[n];
		device_order[b] = new int[n];
	}

	alpha_i = new double[n];  
	beta_id = new double** [n];  
	gamma_id = new double** [n];  
	delta_d = new double* [m]; 
	epsilon_d = new double* [m];
	mu_b = new double[m]; 

	for (int i = 0; i < n; i++) {
		alpha_i[i] = 0;

		beta_id[i] = new double* [m];
		gamma_id[i] = new double* [m];
		for (int b = 0; b < m; b++) {
			beta_id[i][b] = new double[n];
			gamma_id[i][b] = new double[n];
			for (int j = 0; j < n; j++) {
				beta_id[i][b][j] = 0;
				gamma_id[i][b][j] = 0;
			}
		}
	}

	for (int b = 0; b < m; b++) {
		delta_d[b] = new double[n];
		epsilon_d[b] = new double[n];

		mu_b[b] = 0;

		for (int j = 0; j < n; j++) {
			delta_d[b][j] = 0;
			epsilon_d[b][j] = 0;
		}
	}
}

Cover::Cover(std::string path)
{

	read_data(path);
	all_disk.resize(m);

	E_base_ex = new double* [m];
	E_base_tr = new double* [m];
	E_cloud_ex = new double[n];
	E_cloud_tr = new double[n];
	E_cloud_sum_tr = new double* [m];

	distance = new double* [m];
	dist_index = new int* [m];
	device_order = new int* [m];

	for (int b = 0; b < m; b++) {
		E_base_ex[b] = new double[n];
		E_base_tr[b] = new double[n];
		E_cloud_sum_tr[b] = new double[n];

		distance[b] = new double[n];
		dist_index[b] = new int[n];
		device_order[b] = new int[n];
	}

	alpha_i = new double[n];  
	beta_id = new double** [n];  
	gamma_id = new double** [n];  
	delta_d = new double* [m]; 
	epsilon_d = new double* [m];
	mu_b = new double[m]; 


	for (int i = 0; i < n; i++) {
		alpha_i[i] = 0;

		beta_id[i] = new double* [m];
		gamma_id[i] = new double* [m];
		for (int b = 0; b < m; b++) {
			beta_id[i][b] = new double[n];
			gamma_id[i][b] = new double[n];
			for (int j = 0; j < n; j++) {
				beta_id[i][b][j] = 0;
				gamma_id[i][b][j] = 0;
			}
		}
	}

	for (int b = 0; b < m; b++) {
		delta_d[b] = new double[n];
		epsilon_d[b] = new double[n];

		mu_b[b] = 0;

		for (int j = 0; j < n; j++) {
			delta_d[b][j] = 0;
			epsilon_d[b][j] = 0;
		}
	}

	cal_all_distance();		
	order_dist();
	set_all_disk();  
	set_devices_order();   //由dist_index推出device_order

	cal_all_offbase_energy();
	cal_all_offcloud_energy();
}

void Cover::read_data(std::string fname)
{

	using namespace std;
	string data;
	ifstream infpoint;

	infpoint.open(fname);
	// 读取m, n
	getline(infpoint, data);
	istringstream istr1(data);
	istr1 >> m;
	istr1 >> n;

	int CPU;				// CPU
	int BW;				// 带宽
	double f;              //CPU频率
	double p;              //功耗
	double x;    //横坐标
	double y;  //纵坐标
	int da;   //任务数据量
	double eb;
	double em;

	// 云端
	getline(infpoint, data);
	istringstream istr2(data);
	istr2 >>CPU >> BW >> f >> p;
	Resource c = Resource( CPU, BW);
	cloud.set_cloud(c, f, p);


	// 基站
	for (int b = 0; b < m; b++) {
		getline(infpoint, data);
		istringstream istr3(data);
		istr3 >> x >> y >> CPU >> BW >> f >> p;
		double l[2] = { x, y };
		Resource capacity = Resource( CPU, BW);
		BaseStation basestation(l, capacity, b, f,p);
		basestations.push_back(basestation);
	}
	  
	// devices
	for (int i = 0; i < n; i++) {
		getline(infpoint, data);
		istringstream istr4(data);
		istr4 >> x >> y  >>da>> CPU >> BW >>eb>>em;
		Resource taskdemand = Resource(CPU, BW);
		double l[2] = { x, y };
		Device device(l, i, da, taskdemand,eb,em);
		devices.push_back(device);
	}
}

void Cover::set_all_disk()  //构造圆盘
{
	for (int b = 0; b < m; b++)
		for (size_t i = 0; i < n; i++)
		{
			Disk d = Disk(b, dist_index[b][i], distance[b][dist_index[b][i]], basestations[b].get_capacity());
			if (i == 0) d.cover_devices(dist_index[b][i]);
			else d.cover_devices(all_disk[b][i - 1].get_devices(), dist_index[b][i]);
			all_disk[b].push_back(d);
		}
}

void Cover::update_all_disk()
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			all_disk[i][j].set_power(all_disk[i][j].cal_power(distance[i][dist_index[i][j]]));
		}
}

void Cover::set_devices_order()
{
	for (int b = 0; b < m; b++)
		for (int i = 0; i < n; i++)
		{
			int device_index = dist_index[b][i];
			device_order[b][device_index] = i;
		}
}


void Cover::set_theta(double theta)
{
	Disk::set_theta(theta);
	update_all_disk();
}

void Cover::cal_all_distance()
{
	for (int b = 0; b < m; b++)
		for (int i = 0; i < n; i++)
		{
			distance[b][i] = basestations[b].cal_distance(devices[i]);
			dist_index[b][i] = i; 
		}
}

void Cover::cal_all_offbase_energy()
{
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			E_base_ex[b][i] = basestations[b].get_basepf() * (devices[i].get_taskdemand().RCPU / basestations[b].get_basef());
			E_base_tr[b][i] = devices[i].get_e_unitbit() * (1e-9) * devices[i].get_data()*1024*8 + devices[i].get_e_unitm() * (1e-9) * devices[i].get_data() * 1024 * 8* pow(distance[b][i], K);  //数据量以KB为单位
		}
	}
}

void Cover::cal_all_offcloud_energy()
{
	for (int i = 0; i < n; i++) {
		E_cloud_ex[i] = cloud.get_cloudpfvm() * (devices[i].get_taskdemand().RCPU / cloud.get_cloudfvm());
		E_cloud_tr[i] = devices[i].get_data() * ELECTRICITY_INTENSITY;
		for (int b = 0; b < m; b++) {
			E_cloud_sum_tr[b][i] = E_base_tr[b][i] + E_cloud_tr[i];
		}
	}
}


void Cover::order_dist()
{
	// 排序算法
	int min = 0;		

	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n - 1; i++) {
			min = i;
			for (int k = i + 1; k < n; k++)
				if (distance[b][dist_index[b][k]] < distance[b][dist_index[b][min]])
					min = k;
				else if (distance[b][dist_index[b][k]] == distance[b][dist_index[b][min]])
					if (devices[k] < devices[min])
						min = k;
			int temp = dist_index[b][i];
			dist_index[b][i] = dist_index[b][min];
			dist_index[b][min] = temp;
		}
	}
}

const double Cover::get_distance(int b_index, int d_index)
{
	return distance[b_index][d_index];
}

int Cover::get_order(int b_index, int order)
{
	return dist_index[b_index][order];
}

void Cover::print_test()
{
	using namespace std;
	cout << "basestation:" << '\n';
	for (int b = 0; b < m; b++)
		basestations[b].print_basestation();
	cout << '\n';

	cout << "devices:" << '\n';
	for (int i = 0; i < n; i++)
		devices[i].print_device();
	cout << '\n';

	cout << "distance:" << '\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			cout << distance[b][i] << "\t";
		}
		cout << '\n';
	}
	cout << '\n';

	cout << "dist_index:" << '\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			cout << dist_index[b][i] << "\t";
		}
		cout << '\n';
	}
	cout << '\n';

	cout << "device's order:" << '\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			cout << device_order[b][i] << "\t";
		}
		cout << '\n';
	}
	cout << '\n';

	cout << "ordered distance:" << '\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			cout << distance[b][dist_index[b][i]] << "\t";
		}
		cout << '\n';
	}
	cout << '\n';

	cout << "ordered power:" << '\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			double power = all_disk[b][i].get_power();
			cout << setprecision(4) << power << "\t";
		}
		cout << '\n';
	}
	cout << '\n';

	cout << "E_base: "<<'\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			cout << E_base_ex[b][i] + E_base_tr[b][i]<<"\t";
		}
		cout << '\n';
	}
	cout << '\n';

	cout << "E_cloud: " << '\n';
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			cout << E_cloud_sum_tr[b][i] << "\t";
		}
		cout << '\n';
	}
	cout << '\n';
}

void Cover::print_allDisk(const std::vector<std::vector<Disk>>& disks)
{
	for (int b = 0; b < disks.size(); b++)
	{
		for (int i = 0; i < disks[b].size(); i++)
			disks[b][i].print_disk();
		std::cout << '\n';
	}
}

template <typename type>
void Cover::print_xyz(type** y_d0, type*** x_id0, type*** z_id0)
{
	cout << "y_d:\n";
	for (int b = 0; b < m; b++)
	{
		cout << "basestation" << b << ": ";
		for (int i = 0; i < n; i++)
			cout << y_d0[b][i] << ' ';
		cout << '\n';
	}
	cout << "x_id:\n";
	for (int i = 0; i < n; i++)
	{
		cout << "device" << i << ": \n";
		for (int b = 0; b < m; b++) {
			for (int ii = 0; ii < n; ii++)
				cout << x_id0[i][b][ii] << ' ';
			cout << '\n';
		}
		cout << '\n';
	}
	cout << "z_id:\n";
	for (int i = 0; i < n; i++)
	{

		cout << "device" << i << ": \n";
		for (int b = 0; b < m; b++) {
			for (int ii = 0; ii < n; ii++)
				cout << z_id0[i][b][ii] << ' ';
			cout << '\n';
		}
		cout << '\n';
	}
}

void Cover::print_result(Result r)
{
	/// <summary>
	/// </summary>
	/// <param name="r"></param>
	using namespace std;
	cout << "basestation:\n";
	for (int b = 0; b < m; b++)
	{
		const vector<int>& sd = basestations[b].get_direct_served_devices();
		cout << "b" << b << " is direct_serving: ";
		for (int j = 0; j < sd.size(); j++) cout << sd[j] << "\t"<<"\t";

		cout << "disk" << basestations[b].get_chosed_disk() << " is chosed.";
		cout << '\n';
	}
	cout << '\n';

	for (int b = 0; b < m; b++)
	{
		const vector<int>& sd = basestations[b].get_indirect_served_devices();
		cout << "b" << b << " is indirect_serving: ";
		for (int j = 0; j < sd.size(); j++) cout << sd[j] << "\t" << "\t";
		cout << '\n';
	}
	cout << '\n';

	cout << "devices direct_served by:";
	for (int i = 0; i < n; i++) {
		int a = devices[i].get_served_by();
		if(a>=0) 
		    cout << "d" << i << ": " << a << ',' << basestations[a].get_chosed_disk() << ' '<<"\n";
	}
	cout << "------------------------------------\n";
	cout << "total consumpution: ";


	cout << setprecision(5) << r.get_eng_cspt() << '\n';
	cout << "------------------------------------\n";
}

bool sortdevice_byCPU(Device& d1, Device& d2) {   //对devices按CPU需求降序排列
	return d1.get_taskdemand().RCPU > d2.get_taskdemand().RCPU;
}



//本文方法
Result Cover::PrimalDual_Result()
{
	using namespace std;

	clock_t start_, end_;
	start_ = clock();

	Result r;

	int** y_d0 = new int* [m];		// m * n  y_d0[b][i]，原始变量
	int** y_d0_temp = new int* [m];  
	for (int b = 0; b < m; b++) {
		y_d0[b] = new int[n];
		y_d0_temp[b] = new int[n];
		memset(y_d0[b], 0, sizeof(int) * n);
		memset(y_d0_temp[b], 0, sizeof(int) * n);
	}

	int*** x_id0 = new int** [n]; //n*m*n x_id[n][m][n]，原始变量
	int*** x_id0_temp = new int** [n];  
	for (int i = 0; i < n; i++) {
		x_id0[i] = new int* [m];
		x_id0_temp[i] = new int* [m];
		for (int b = 0; b < m; b++) {
			x_id0[i][b] = new int[n];
			x_id0_temp[i][b] = new int[n];
			memset(x_id0[i][b], 0, sizeof(int) * n);
			memset(x_id0_temp[i][b], 0, sizeof(int) * n);
		}
	}

	int*** z_id0 = new int** [n];	// n * m * n z_id0[n][m][n]，原始变量
	int*** z_id0_temp = new int** [n]; 
	for (int i = 0; i < n; i++) {
		z_id0[i] = new int* [m];
		z_id0_temp[i] = new int* [m];
		for (int b = 0; b < m; b++) {
			z_id0[i][b] = new int[n];
			z_id0_temp[i][b] = new int[n];
			memset(z_id0[i][b], 0, sizeof(int) * n);
			memset(z_id0_temp[i][b], 0, sizeof(int) * n);
		}
	}


	double* alpha_i_temp = new double[n];  //对偶变量，
	double*** beta_id_temp = new double** [n];  //对偶变量
	double*** gamma_id_temp = new double** [n];  //对偶变量
	double** delta_d_temp = new double* [m]; //对偶变量
	double** epsilon_d_temp = new double* [m];//对偶变量
	double* mu_b_temp = new double[m]; //对偶变量

	memset(alpha_i_temp, 0, sizeof(double) * n);
	memset(mu_b_temp, 0, sizeof(double)*m);

	for (int i = 0; i < n; i++) {
		beta_id_temp[i] = new double* [m];
		gamma_id_temp[i] = new double* [m];
		for (int b = 0; b < m; b++) {
			beta_id_temp[i][b] = new double[n];
			gamma_id_temp[i][b] = new double[n];
			memset(beta_id_temp[i][b],0,sizeof(double)*n);
			memset(gamma_id_temp[i][b], 0, sizeof(double) * n);
		}
	}

	for (int b = 0; b < m; b++) {
		delta_d_temp[b] = new double[n];
		epsilon_d_temp[b] = new double[n];
		memset(delta_d_temp[b], 0, sizeof(double) * n);
		memset(epsilon_d_temp[b], 0, sizeof(double) * n);
	}

	vector<Device> device_covered_by_maxdisk;  //被max_disk覆盖的devices
	vector<int> devices_served_by_maxdisk;  //被max_disk直接服务和间接服务的devices
	int base;
	int device;


	//猜测一个半径最大的圆盘，猜测m*n次
	Disk max_final;
	for (auto samebase_disks : all_disk) {
		for (int i = samebase_disks.size()-1; i >= 0; i--) {
	        Disk max_disk = samebase_disks[i];
			base = max_disk.get_base();
			device = max_disk.get_device();
			//if (distance[base][device] <= 120) {
				//输出max_disk的相关信息
				std::cout << "\n";
				max_disk.print_disk();
				std::cout << "\n";
				device_covered_by_maxdisk.clear();
				//初始化device_covered_by_maxdisk
				for (int j = 0; j < max_disk.get_devices().size(); j++) {
					for (int jj = 0; jj < n; jj++)
						if (devices[jj].get_index() == max_disk.get_devices()[j]) {
							device_covered_by_maxdisk.push_back(devices[jj]);
							break;
						}
				}

				//重置每一轮猜测的对偶变量，对偶变量初始化
				for (int i = 0; i < n; i++) {
					alpha_i_temp[i] = 0;
					for (int b = 0; b < m; b++)
						for (int ii = 0; ii < n; ii++) {
							beta_id_temp[i][b][ii] = 0;
							gamma_id_temp[i][b][ii] = 0;
							x_id0_temp[i][b][ii] = 0;
							z_id0_temp[i][b][ii] = 0;
						}
				}
				for (int b = 0; b < m; b++) {
					mu_b_temp[b] = 0;
					for (int i = 0; i < n; i++) {
						delta_d_temp[b][i] = 0;
						epsilon_d_temp[b][i] = 0;
						y_d0_temp[b][i] = 0;
					}
				}

				std::sort(device_covered_by_maxdisk.begin(), device_covered_by_maxdisk.end(), sortdevice_byCPU);


				//确定被max_disk覆盖且迁移到basestation或cloud的devices
				int maxdisk_CPU;
				int maxdisk_BW;
				maxdisk_CPU = max_disk.get_capacity().RCPU;
				maxdisk_BW = max_disk.get_capacity().RBW;
	
				for (auto dev : device_covered_by_maxdisk) {
					if (maxdisk_BW >= dev.get_taskdemand().RBW) {
						if (maxdisk_CPU >= dev.get_taskdemand().RCPU) {
							max_disk.cover_devicestobase(dev.get_index());
							maxdisk_CPU -= dev.get_taskdemand().RCPU;
							maxdisk_BW -= dev.get_taskdemand().RBW;
						}
						else {
							max_disk.cover_devicestocloud(dev.get_index());
							maxdisk_BW -= dev.get_taskdemand().RBW;
						}
					}
				}

				std::vector<int> a1 = max_disk.get_devicestobase();
				std::vector<int> a2 = max_disk.get_devicestocloud();


				devices_served_by_maxdisk.clear();
				devices_served_by_maxdisk.insert(devices_served_by_maxdisk.end(), a1.begin(), a1.end());
				devices_served_by_maxdisk.insert(devices_served_by_maxdisk.end(), a2.begin(), a2.end());


				//确定max_disk对应的basestation直接服务和间接服务的devices
				basestations[base].clear_direct_served_devices();
				basestations[base].clear_indirect_served_devices();
				for (int i = 0; i < max_disk.get_devicestobase().size(); i++)
					basestations[base].direct_served_device(max_disk.get_devicestobase()[i]);
				for (int i = 0; i < max_disk.get_devicestocloud().size(); i++)
					basestations[base].indirect_served_device(max_disk.get_devicestocloud()[i]);


				std::vector<std::vector<Disk> > instanceDisk = construct_instanceDisk(max_disk, devices_served_by_maxdisk);
				std::vector<int> unserved_devices = construct_unserved_devices(max_disk, devices_served_by_maxdisk);

				if (is_instanceDisk_feasible(instanceDisk, unserved_devices))
				{

					//设置max_disk对应的y_d，以及max_disk直接服务和间接服务的devices对应的x_id,z_id为1
					y_d0_temp[base][device] = 1;
					for (auto direct_de : max_disk.get_devicestobase())
						x_id0_temp[direct_de][base][device] = 1;
					for (auto indirect_de : max_disk.get_devicestocloud())
						z_id0_temp[indirect_de][base][device] = 1;


					instancePD(unserved_devices, instanceDisk, x_id0_temp, y_d0_temp, z_id0_temp, alpha_i_temp, beta_id_temp, gamma_id_temp, delta_d_temp, epsilon_d_temp, mu_b_temp);

					//移除同心的多个圆盘
					remove_unneccessary_disks(instanceDisk, x_id0_temp, y_d0_temp, z_id0_temp);

					if (is_capacity_satisfy(x_id0_temp, y_d0_temp, z_id0_temp)) {
						Result r_temp = cal_result(y_d0_temp, x_id0_temp, z_id0_temp);

						if (r_temp.get_eng_cspt() < r.get_eng_cspt()) {   
							r = r_temp;

							for (int i = 0; i < n; i++) {
								alpha_i[i] = alpha_i_temp[i];
								for (int b = 0; b < m; b++)
									for (int ii = 0; ii < n; ii++) {
										x_id0[i][b][ii] = x_id0_temp[i][b][ii];
										z_id0[i][b][ii] = z_id0_temp[i][b][ii];
										beta_id[i][b][ii] = beta_id_temp[i][b][ii];
										gamma_id[i][b][ii] = gamma_id_temp[i][b][ii];
									}
							}

							for (int b = 0; b < m; b++) {
								mu_b[b] = mu_b_temp[b];
								for (int i = 0; i < n; i++) {
									y_d0[b][i] = y_d0_temp[b][i];
									delta_d[b][i] = delta_d_temp[b][i];
									epsilon_d[b][i] = epsilon_d_temp[b][i];
								}
							}

							max_final = max_disk;   //能耗最小的那一轮的max_disk

						}

					}
					else
						break;
				}
				else
					break;
			//}
		}
	}	  

	cout << "\n max_final: " << max_final.get_base()<<","<<max_final.get_device()<< '\n';
   // print_xyz(y_d0, x_id0, z_id0);


	//Result r = Result();

	r = cal_result(y_d0, x_id0, z_id0);



			for (int i = 0; i < n; i++) {
				for (int b = 0; b < m; b++) {
					delete[] x_id0_temp[i][b];
					delete[] z_id0_temp[i][b];
				}
				delete[] x_id0_temp[i];
				delete[] z_id0_temp[i];
			}
			delete[] x_id0_temp;
			delete[] z_id0_temp;

			for (int b = 0; b < m; b++) {
				delete[] y_d0_temp[b];
			}
			delete[] y_d0_temp;


	end_ = clock();
	r.set_time((double)end_ - start_);

	return r;
}

std::vector<std::vector<Disk>> Cover::construct_instanceDisk(Disk& max_disk, std::vector<int>& devices_served_by_maxdisk)
{
	std::vector<std::vector<Disk> > insDisk;   
	insDisk.resize(static_cast<int>(m) - static_cast<int>(1));

	int ba = max_disk.get_base();
	int de = max_disk.get_device();

	//将与max_disk同心的全部圆盘去除，包括max_disk本身
	if (ba == 0) {
		for (int b = 1; b < m; b++) 
			for (int i = 0; i < all_disk[b].size(); i++) 
				insDisk[b-1].push_back(all_disk[b][i]);
	}
	else if (ba == m - 1) {
		for (int b = 0; b < m - 1; b++)
			for (int i = 0; i < all_disk[b].size(); i++)
				insDisk[b].push_back(all_disk[b][i]);
	}
	else {
		for (int b = 0; b < ba; b++)
			for (int i = 0; i < all_disk[b].size(); i++)
				insDisk[b].push_back(all_disk[b][i]);
		for (int b = ba + 1; b < m; b++)
			for (int j = 0; j < all_disk[b].size(); j++)
				insDisk[b-1].push_back(all_disk[b][j]);
	}


	const double maxdisk_radii = get_distance(max_disk.get_base(),max_disk.get_device());

	int index;  //存放半径大于max_disk半径的第一个圆盘的第二维下标

	vector<int> getdevices;


	//去除与max_disk非同心，半径大于max_disk半径的圆盘，并更新它们覆盖的devices
	for (int b = 0; b < m - 1; b++) {
		index = insDisk[b].size()-1;
		for (int i = 0; i < insDisk[b].size(); i++){
			if (get_distance(insDisk[b][i].get_base(), insDisk[b][i].get_device()) > maxdisk_radii) {
				index = i;
				break;
			}
			getdevices.clear();                                             //
			getdevices = insDisk[b][i].get_devices();
			for (int j = 0; j < getdevices.size(); j++)
				for (int ii = 0; ii < devices_served_by_maxdisk.size(); ii++) {
					if (getdevices[j] == devices_served_by_maxdisk[ii])
						getdevices.erase(getdevices.begin() + j);
				}
			insDisk[b][i].get_devices() = getdevices;                      //
		}
		
		insDisk[b].erase(insDisk[b].begin() + index, insDisk[b].end());
	}

	basestations_temp.clear();
	basestations_temp = basestations;
	basestations_temp.erase(basestations_temp.begin()+ba);


	return insDisk;
}

std::vector<int> Cover::construct_unserved_devices(Disk& max_disk, std::vector<int>& devices_served_by_maxdisk)
{

	std::vector<int> unserved_des;
	for (int i = 0; i < n; i++) {
		unserved_des.push_back(i);
	}
                                              

	for (int j = 0; j < devices_served_by_maxdisk.size(); j++) {
		vector<int>::iterator it;
		it = find(unserved_des.begin(), unserved_des.end(), devices_served_by_maxdisk[j]);
		if(it!=unserved_des.end())
		   unserved_des.erase(it);
	}

	return unserved_des;
}

bool Cover::is_instanceDisk_feasible(std::vector<std::vector<Disk>>& instanceDisk, std::vector<int> unserved_devices)
{ //判断一个实例是否有可行解
	
	int bBW = 0;  
	int dBW = 0;  

	for (int b = 0; b < instanceDisk.size(); b++)
		bBW += basestations[instanceDisk[b][0].get_base()].get_capacity().RBW;
	cout << "\n the total BW of basestations in instanceDisk:" << bBW << '\n';
	for (int i = 0; i < unserved_devices.size(); i++)
		dBW += devices[unserved_devices[i]].get_taskdemand().RBW;
	cout << "\n the total BW of devices in unserved_devices:" << dBW << '\n';
	
	vector<int> covered_devices;
	int i = 0;
	
	//unserved_devices被全部覆盖或者遍历完instanceDisk跳出循环
	while (unserved_devices.size()>0 && i<instanceDisk.size()){
		for (int j = 0; j < instanceDisk[i].size(); j++) {
			covered_devices.clear();
			covered_devices = instanceDisk[i][j].get_devices();

			for (int j = 0; j < covered_devices.size(); j++) {
				vector<int>::iterator it;
				it = find(unserved_devices.begin(), unserved_devices.end(), covered_devices[j]);
				if(it!=unserved_devices.end())
					unserved_devices.erase(it);
			}
		}
		i++;
	}


	if ((unserved_devices.size() >0) || (bBW<dBW))
		return false;
	else 
		return true;
}

//构造对偶可行解
void Cover::instancePD(std::vector<int> unserved_devices, std::vector<std::vector<Disk>>& instanceDisk, int*** x_id0_temp, int** y_d0_temp, int*** z_id0_temp, double* alpha_i_temp, double*** beta_id_temp, double*** gamma_id_temp, double** delta_d_temp, double** epsilon_d_temp, double* mu_b_temp)
{
	vector<vector<Disk>> disks_selected_temp; //临时被选中的圆盘
	disks_selected_temp.resize(static_cast<int>(m) - static_cast<int>(1));

	int sum_dCPU1;   
	int sum_dBW1,sum_dBW2;   
	Disk min_disk_coveri1;    //事件1
	Disk min_disk_coveri2;    //事件2


	std::vector<int> baseCPU; //基站当前剩余的CPU资源
	std::vector<int> baseBW;  //基站当前剩余的BW资源
	for (int b = 0; b < basestations_temp.size(); b++) {
		baseCPU.push_back(basestations_temp[b].get_capacity().RCPU);
		baseBW.push_back(basestations_temp[b].get_capacity().RBW);
	}


    //是否冻结的标记，用于判断是否需要增长对偶变量
	bool* is_frozen_alpha_i = new bool[n];
	bool*** is_frozen_beta_id = new bool** [n];
	bool*** is_frozen_gamma_id = new bool** [n];  
	bool** is_frozen_delta_d = new bool* [m];
	bool** is_frozen_epsilon_d= new bool* [m];
	bool* is_frozen_mu_b = new bool[m];

	memset(is_frozen_alpha_i, true, sizeof(bool) * n);  //初始化时is_frozen_alpha_i为true，其他为false
	memset(is_frozen_mu_b, false, sizeof(bool) * m);
	

	for (int i = 0; i < n; i++) {
		is_frozen_beta_id[i] = new bool* [m];
		is_frozen_gamma_id[i] = new bool* [m];
		for (int b = 0; b < m; b++) {
			is_frozen_beta_id[i][b] = new bool[n];
			is_frozen_gamma_id[i][b] = new bool[n];
			memset(is_frozen_beta_id[i][b], false, sizeof(bool) * n); //第三维按devices下标排序
			memset(is_frozen_gamma_id[i][b], false, sizeof(bool) * n);
		}
	}
				

	for (int b = 0; b < m; b++) {
		is_frozen_delta_d[b] = new bool[n];
		is_frozen_epsilon_d[b] = new bool[n];
		memset(is_frozen_delta_d[b], false, sizeof(bool) * n);
		memset(is_frozen_epsilon_d[b], false, sizeof(bool) * n);
	}

	//cout << "test instanceDisk's devices_covered: " << '\n';
	//print_allDisk(instanceDisk);

	while (unserved_devices.size() > 0) {

		//先判断是否要更新对偶变量
		for (int i = 0; i < unserved_devices.size(); i++)   //未被冻住的devices都要增长alpha_i
			if (is_frozen_alpha_i[unserved_devices[i]])
			alpha_i_temp[unserved_devices[i]] += INCREASE;


		for (int i = 0; i < unserved_devices.size(); i++)
			for (int b = 0; b < basestations_temp.size(); b++) {
				for (int ii = 0; ii < instanceDisk[b].size(); ii++) {  
					if (is_frozen_beta_id[unserved_devices[i]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]])    //第三维是devices的index
						beta_id_temp[unserved_devices[i]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]] += INCREASE;
					if (is_frozen_gamma_id[unserved_devices[i]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]])
						gamma_id_temp[unserved_devices[i]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]] += INCREASE;
				}
				for (int ii = 0; ii <n; ii++) {   
					if (is_frozen_delta_d[basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]])
						delta_d_temp[basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]] += INCREASE;  //第二维是devices的index
					if (is_frozen_epsilon_d[basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]])
						epsilon_d_temp[basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][ii]] += INCREASE;
				}
			}

		//mu_b更新的增量与其他对偶变量不一样
	    for (int b = 0; b < basestations_temp.size(); b++) 
			if (is_frozen_mu_b[basestations_temp[b].get_index()])
				mu_b_temp[basestations_temp[b].get_index()] += INCREASE*(static_cast<double>(basestations_temp[b].get_capacity().RCPU) + static_cast<double>(basestations_temp[b].get_capacity().RBW));


		//事件1
		vector<int> record_unserved_devices1;    //记录事件1unserved_devices中需要删除的devices
		record_unserved_devices1.clear();
		for (int undevice = 0; undevice < unserved_devices.size(); undevice++) {
			for (int basesta = 0; basesta < basestations_temp.size(); basesta++) {
				if (alpha_i_temp[unserved_devices[undevice]] <= (E_base_ex[basestations_temp[basesta].get_index()][unserved_devices[undevice]] + E_base_tr[basestations_temp[basesta].get_index()][unserved_devices[undevice]]) && fabs(alpha_i_temp[unserved_devices[undevice]]- (E_base_ex[basestations_temp[basesta].get_index()][unserved_devices[undevice]] + E_base_tr[basestations_temp[basesta].get_index()][unserved_devices[undevice]])) < ERROR) {
					cout << "\n E_base: " <<setprecision(6)<< E_base_ex[basestations_temp[basesta].get_index()][unserved_devices[undevice]] + E_base_tr[basestations_temp[basesta].get_index()][unserved_devices[undevice]] <<", unserved_devices[undevice]:"<< unserved_devices[undevice]<<", basestation:"<< basestations_temp[basesta].get_index()<<'\n';
					
					//找到以b为中心、覆盖i的圆盘中半径最小的圆盘：all_disk[basestations_temp[basesta].get_index()][device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]]
					if (device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]] < instanceDisk[basesta].size()) {    //找到以b为中心、覆盖i的圆盘中半径最小的圆盘，这样的圆盘不一定存在于instanceDisk中，所以判断
						min_disk_coveri1 = instanceDisk[basesta][device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]]];
						cout << "\n min_disk_coveri1: " << min_disk_coveri1.get_base() << " " << min_disk_coveri1.get_device() << '\n';
						sum_dCPU1 = 0;
						sum_dBW1 = 0;
						vector<Disk> unselectdisk_coveri1;
						vector<Disk> selectdisk_coveri1;
						bool selectdisk_flag1 = false;

						for (int i_order = device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]]; i_order < instanceDisk[basesta].size(); i_order++) {
							if (!instanceDisk[basesta][i_order].get_state())    //覆盖i的圆盘中有被选中的
								selectdisk_coveri1.push_back(instanceDisk[basesta][i_order]);
							if (instanceDisk[basesta][i_order].get_state())    //覆盖i的圆盘中有未被选中
								unselectdisk_coveri1.push_back(instanceDisk[basesta][i_order]);
						}

						//事件1.1
						if (selectdisk_coveri1.size() > 0) {   //覆盖i的圆盘中有被选中的
							//计算覆盖i、已被选中（开放）的圆盘中半径最小的那个的CPU和BW带宽剩余资源是否能满足i
							if (baseCPU[basesta] - devices[unserved_devices[undevice]].get_taskdemand().RCPU >= 0 && baseBW[basesta] - devices[unserved_devices[undevice]].get_taskdemand().RBW >= 0) {
								for (int i = 0; i < selectdisk_coveri1.size(); i++) {
									instanceDisk[basesta][device_order[basestations_temp[basesta].get_index()][selectdisk_coveri1[i].get_device()]].cover_devicestobase(unserved_devices[undevice]);
								}

								is_frozen_alpha_i[unserved_devices[undevice]] = false;   //冻住alpha_i_temp[unserved_devices[undevice]]

								for (int b = 0; b < basestations_temp.size(); b++)
									for (int i = 0; i < instanceDisk[b].size(); i++) {
										if (i >= device_order[basestations_temp[b].get_index()][unserved_devices[undevice]]) {
											is_frozen_beta_id[unserved_devices[undevice]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
											is_frozen_gamma_id[unserved_devices[undevice]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
										}
									}
								//记录要移除的devices
								record_unserved_devices1.push_back(unserved_devices[undevice]);
								//更新对应基站的资源容量
								baseCPU[basesta] -= devices[unserved_devices[undevice]].get_taskdemand().RCPU;
								baseBW[basesta] -= devices[unserved_devices[undevice]].get_taskdemand().RBW;
								selectdisk_flag1 = true;   //标记device被已选中的某个圆盘满足，不再执行事件1.2
							}
							else if (baseCPU[basesta] - devices[unserved_devices[undevice]].get_taskdemand().RCPU < 0)
								for (int select = 0; select < selectdisk_coveri1.size(); select++)
									is_frozen_delta_d[basestations_temp[basesta].get_index()][selectdisk_coveri1[select].get_device()] = true;
							else if (baseBW[basesta] - devices[unserved_devices[undevice]].get_taskdemand().RBW < 0)
								for (int select = 0; select < selectdisk_coveri1.size(); select++)
									is_frozen_epsilon_d[basestations_temp[basesta].get_index()][selectdisk_coveri1[select].get_device()] = true;
						}

						//事件1.2
						if (unselectdisk_coveri1.size() > 0 && !selectdisk_flag1) {
							for (int i_min = 0; i_min < unselectdisk_coveri1[0].get_devices().size(); i_min++) {
								for(int unde=0;unde<unserved_devices.size();unde++)
									if (unselectdisk_coveri1[0].get_devices()[i_min] == unserved_devices[unde]) {
										sum_dCPU1 += devices[unserved_devices[unde]].get_taskdemand().RCPU;
										sum_dBW1 += devices[unserved_devices[unde]].get_taskdemand().RBW;
									}
							}

							//如果CPU和BW都满足，增长同心、未被选中且半径大于等于unselectdisk_coveri[0]的圆盘的beta_id
							if (sum_dCPU1 <= baseCPU[basesta] && sum_dBW1 <= baseBW[basesta]) {
								for (int unselect = 0; unselect < unselectdisk_coveri1.size(); unselect++)
									is_frozen_beta_id[unserved_devices[undevice]][basestations_temp[basesta].get_index()][unselectdisk_coveri1[unselect].get_device()] = true;
							}
							if (sum_dCPU1 > baseCPU[basesta])
								for (int unselect = 0; unselect < unselectdisk_coveri1.size(); unselect++)
									is_frozen_delta_d[basestations_temp[basesta].get_index()][unselectdisk_coveri1[unselect].get_device()] = true;
							if(sum_dBW1 > baseBW[basesta])
								for(int unselect = 0; unselect < unselectdisk_coveri1.size(); unselect++)
									is_frozen_epsilon_d[basestations_temp[basesta].get_index()][unselectdisk_coveri1[unselect].get_device()] = true;
						}
					}
				}
			}
		}

		vector<int>::iterator iter1;
		for (int i1 = 0; i1 < record_unserved_devices1.size(); i1++) {
			iter1 = find(unserved_devices.begin(), unserved_devices.end(), record_unserved_devices1[i1]);
				if (iter1 != unserved_devices.end())  //找到
					unserved_devices.erase(iter1);
		}

		//事件2
		vector<int> record_unserved_devices2;  //记录事件2unserved_devices中需要删除的devices
		record_unserved_devices2.clear();
		for (int undevice = 0; undevice < unserved_devices.size(); undevice++) {
			for (int basesta = 0; basesta < basestations_temp.size(); basesta++) {
				if (alpha_i_temp[unserved_devices[undevice]] <= (E_cloud_ex[unserved_devices[undevice]] + E_cloud_sum_tr[basestations_temp[basesta].get_index()][unserved_devices[undevice]]) && fabs(alpha_i_temp[unserved_devices[undevice]] - (E_cloud_ex[unserved_devices[undevice]] + E_cloud_sum_tr[basestations_temp[basesta].get_index()][unserved_devices[undevice]])) < ERROR) {
					cout << "\n E_cloud: " << E_cloud_ex[unserved_devices[undevice]] + E_cloud_sum_tr[basestations_temp[basesta].get_index()][unserved_devices[undevice]] << ", unserved_devices[undevice]: " << unserved_devices[undevice]<<", basestation: "<< basestations_temp[basesta].get_index()<<'\n';
				
					//找到以b为中心、覆盖i的圆盘中半径最小的圆盘
					if (device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]] < instanceDisk[basesta].size()) {    //找到以b为中心、覆盖i的圆盘中半径最小的圆盘
						min_disk_coveri2 = instanceDisk[basesta][device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]]];
						cout << "\n min_disk_coveri2: " << min_disk_coveri2.get_base() << " " << min_disk_coveri2.get_device() << '\n';
						sum_dBW2 = 0;
						vector<Disk> unselectdisk_coveri2;
						vector<Disk> selectdisk_coveri2;
						bool selectdisk_flag2 = false;

						for (int i_order = device_order[basestations_temp[basesta].get_index()][unserved_devices[undevice]]; i_order < instanceDisk[basesta].size(); i_order++) {
							if (!instanceDisk[basesta][i_order].get_state())    //覆盖i的圆盘中有被选中的
								selectdisk_coveri2.push_back(instanceDisk[basesta][i_order]);
							if (instanceDisk[basesta][i_order].get_state())    //覆盖i的圆盘中有未被选中
								unselectdisk_coveri2.push_back(instanceDisk[basesta][i_order]);
						}

						//事件2.1
						if (selectdisk_coveri2.size() > 0) {   //覆盖i的圆盘中有被选中的
							//计算覆盖i、已被选中（开放）的圆盘中半径最小的那个的BW带宽剩余资源是否能满足i
							if ( baseBW[basesta] - devices[unserved_devices[undevice]].get_taskdemand().RBW >= 0) {
								for (int i = 0; i < selectdisk_coveri2.size(); i++) {
									instanceDisk[basesta][device_order[basestations_temp[basesta].get_index()][selectdisk_coveri2[i].get_device()]].cover_devicestocloud(unserved_devices[undevice]);
								}
								
								is_frozen_alpha_i[unserved_devices[undevice]] = false;   //冻住alpha_i_temp[unserved_devices[undevice]]
								
								for (int b = 0; b < basestations_temp.size(); b++)
									for (int i = 0; i < instanceDisk[b].size(); i++) {
										if (i >= device_order[basestations_temp[b].get_index()][unserved_devices[undevice]]) {
											is_frozen_beta_id[unserved_devices[undevice]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
											is_frozen_gamma_id[unserved_devices[undevice]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
										}
									}
								//记录要移除的devices，用于后面更新未被服务的列表unserved_devices
								record_unserved_devices2.push_back(unserved_devices[undevice]);
								//更新对应基站的资源容量
								baseBW[basesta] -= devices[unserved_devices[undevice]].get_taskdemand().RBW;
								selectdisk_flag2 = true;   //标记device被已选中的某个圆盘满足，不再执行事件2.2
							}
							else if (baseBW[basesta] - devices[unserved_devices[undevice]].get_taskdemand().RBW < 0)
								for (int select = 0; select < selectdisk_coveri2.size(); select++)
									is_frozen_epsilon_d[basestations_temp[basesta].get_index()][selectdisk_coveri2[select].get_device()] = true;
						}

						//事件2.2
						if (unselectdisk_coveri2.size() > 0 && !selectdisk_flag2) {
							for (int i_min = 0; i_min < unselectdisk_coveri2[0].get_devices().size(); i_min++) {
								for (int unde = 0; unde < unserved_devices.size(); unde++)
									if (unselectdisk_coveri2[0].get_devices()[i_min] == unserved_devices[unde]) {
										sum_dBW2 += devices[unserved_devices[unde]].get_taskdemand().RBW;
									}
							}

							//如果BW满足，增长同心、未被选中且半径大于等于unselectdisk_coveri[0]的圆盘的gamma_id
							if (sum_dBW2 <= baseBW[basesta]) {
								for (int unselect = 0; unselect < unselectdisk_coveri2.size(); unselect++)
									is_frozen_gamma_id[unserved_devices[undevice]][basestations_temp[basesta].get_index()][unselectdisk_coveri2[unselect].get_device()] = true;
							}
							if (sum_dBW2 > baseBW[basesta])
								for (int unselect = 0; unselect < unselectdisk_coveri2.size(); unselect++)
									is_frozen_epsilon_d[basestations_temp[basesta].get_index()][unselectdisk_coveri2[unselect].get_device()] = true;
						}
					}
				}
			}
		}

		vector<int>::iterator iter2;
		for (int i2 = 0; i2 < record_unserved_devices2.size(); i2++) {
			iter2 = find(unserved_devices.begin(), unserved_devices.end(), record_unserved_devices2[i2]);
			if (iter2 != unserved_devices.end())  //找到
				unserved_devices.erase(iter2);
		}


		//事件3
		for (int basesta = 0; basesta < basestations_temp.size(); basesta++) {
			for (int disk = 0; disk < instanceDisk[basesta].size(); disk++) {
				double sum_beta_id = 0;
				double sum_gamma_id = 0;
				double left_sum = 0;

				double sum_beta_id1 = 0;
				double sum_gamma_id1 = 0;
				double left_sum1 = 0;

				//对于未选中的圆盘执行
				if (instanceDisk[basesta][disk].get_state()) {

					for (int dev = 0; dev < instanceDisk[basesta][disk].get_devices().size(); dev++) {
						sum_beta_id += beta_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()];
						sum_gamma_id += gamma_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()];

						if(!(beta_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]<INCREASE))
							sum_beta_id1 += (beta_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]-INCREASE);
						if(!(gamma_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]<INCREASE))
						    sum_gamma_id1 += (gamma_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]-INCREASE);
					}
					left_sum = sum_beta_id + sum_gamma_id + basestations_temp[basesta].get_capacity().RCPU * delta_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()] + basestations_temp[basesta].get_capacity().RBW* epsilon_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()];
					left_sum1 = sum_beta_id1 + sum_gamma_id1;
					if(!(delta_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]<INCREASE))
						left_sum1 += basestations_temp[basesta].get_capacity().RCPU * (delta_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]-INCREASE);
					if (!(epsilon_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]<INCREASE))
						left_sum1 += basestations_temp[basesta].get_capacity().RBW * (epsilon_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()] - INCREASE);

					//圆盘收紧
					if (left_sum >= instanceDisk[basesta][disk].get_power() * COVER_TIME && left_sum1 < instanceDisk[basesta][disk].get_power() * COVER_TIME) {
					
						cout << "\n instanceDisk[basesta][disk]: " << instanceDisk[basesta][disk].get_base() << "," << instanceDisk[basesta][disk].get_device()<<"  disk_energy:"<< instanceDisk[basesta][disk].get_power();
						cout << "\n instanceDisk[basesta][disk].get_devices(): ";
						for (int dev = 0; dev < instanceDisk[basesta][disk].get_devices().size(); dev++)
							cout << instanceDisk[basesta][disk].get_devices()[dev] << " ";

						cout << "\n beta_id_temp: ";
						for (int dev = 0; dev < instanceDisk[basesta][disk].get_devices().size(); dev++) {
							cout << beta_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()] << " ";
						}
						cout << "\n gamma_id_temp: ";
						for (int dev = 0; dev < instanceDisk[basesta][disk].get_devices().size(); dev++) {
							cout << gamma_id_temp[instanceDisk[basesta][disk].get_devices()[dev]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()] << " ";
						}
						cout << "\n delta_d_temp: " << delta_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()];
						cout << "\n epsilon_d_temp: " << epsilon_d_temp[instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()]<<'\n';
				
						instanceDisk[basesta][disk].open_disk();  //state=false
						is_frozen_mu_b[instanceDisk[basesta][disk].get_base()] = true;


						//更新圆盘直接服务的devices
						for (int i = 0; i < unserved_devices.size(); i++)
							
							if (beta_id_temp[unserved_devices[i]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()] > 0 && (baseCPU[basesta] - devices[unserved_devices[i]].get_taskdemand().RCPU)>=0 && (baseBW[basesta] - devices[unserved_devices[i]].get_taskdemand().RBW)>=0) {
								instanceDisk[basesta][disk].cover_devicestobase(unserved_devices[i]);
								
								baseCPU[basesta] -= devices[unserved_devices[i]].get_taskdemand().RCPU;
								baseBW[basesta] -= devices[unserved_devices[i]].get_taskdemand().RBW;
							}

						//更新圆盘间接服务的devices
						for (int j = 0; j < unserved_devices.size(); j++) {
							int count = 0;
							vector<int>::iterator it1;
							
							if (gamma_id_temp[unserved_devices[j]][instanceDisk[basesta][disk].get_base()][instanceDisk[basesta][disk].get_device()] > 0 && (baseBW[basesta] - devices[unserved_devices[j]].get_taskdemand().RBW) >= 0) {
								
								for (int ii = 0; ii < instanceDisk[basesta][disk].get_devicestobase().size(); ii++) {
									if (instanceDisk[basesta][disk].get_devicestobase()[ii] == unserved_devices[j])
										break;
									else
										count++;
									if (count == instanceDisk[basesta][disk].get_devicestobase().size()) {
										instanceDisk[basesta][disk].cover_devicestocloud(unserved_devices[j]);
										baseBW[basesta] -= devices[unserved_devices[j]].get_taskdemand().RBW;
									}
								}
								//}
							}
						}


						cout << "instanceDisk[basesta][disk].cover_devicestobase: ";
						for (int i = 0; i < instanceDisk[basesta][disk].get_devicestobase().size(); i++) {
							cout << instanceDisk[basesta][disk].get_devicestobase()[i] << " ";
						}

						cout << "instanceDisk[basesta][disk].cover_devicestocloud: ";
						for (int i = 0; i < instanceDisk[basesta][disk].get_devicestocloud().size(); i++) {
							cout << instanceDisk[basesta][disk].get_devicestocloud()[i] << " ";
						}
						cout << '\n';

						//冻住instanceDisk[basesta][disk].get_devicestobase()中的devices
						for (int debase = 0; debase < instanceDisk[basesta][disk].get_devicestobase().size(); debase++) {
							is_frozen_alpha_i[instanceDisk[basesta][disk].get_devicestobase()[debase]] = false;
							//instanceDisk[basesta][disk].get_devicestobase()[debase]对应的所有beta_id_temp、gamma_id_temp都停止增长
							for (int b = 0; b < basestations_temp.size(); b++)
								for (int i = 0; i < instanceDisk[b].size(); i++) {
									if (i >= device_order[basestations_temp[b].get_index()][instanceDisk[basesta][disk].get_devicestobase()[debase]]) {
										is_frozen_beta_id[instanceDisk[basesta][disk].get_devicestobase()[debase]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
										is_frozen_gamma_id[instanceDisk[basesta][disk].get_devicestobase()[debase]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
									}
								}
						}
						//冻住instanceDisk[basesta][disk].get_devicestocloud()中的devices
						for (int decloud = 0; decloud < instanceDisk[basesta][disk].get_devicestocloud().size(); decloud++) {
							is_frozen_alpha_i[instanceDisk[basesta][disk].get_devicestocloud()[decloud]] = false;
							//instanceDisk[basesta][disk].get_devicestocloud()[decloud]对应的所有beta_id_temp、gamma_id_temp都停止增长
							for (int b = 0; b < basestations_temp.size(); b++)
								for (int i = 0; i < instanceDisk[b].size(); i++) {
									if (i >= device_order[basestations_temp[b].get_index()][instanceDisk[basesta][disk].get_devicestocloud()[decloud]]) {
										is_frozen_beta_id[instanceDisk[basesta][disk].get_devicestocloud()[decloud]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
										is_frozen_gamma_id[instanceDisk[basesta][disk].get_devicestocloud()[decloud]][basestations_temp[b].get_index()][dist_index[basestations_temp[b].get_index()][i]] = false;
									}
								}
						}

						vector<int>::iterator it2,it3;
						//更新unserved_devices
							for (int debase = 0; debase < instanceDisk[basesta][disk].get_devicestobase().size(); debase++) {
								it2 = find(unserved_devices.begin(), unserved_devices.end(), instanceDisk[basesta][disk].get_devicestobase()[debase]);
								if (it2 != unserved_devices.end())  //找到了
									unserved_devices.erase(it2);
							}
							for (int decloud = 0; decloud < instanceDisk[basesta][disk].get_devicestocloud().size(); decloud++) {
								it3 = find(unserved_devices.begin(), unserved_devices.end(), instanceDisk[basesta][disk].get_devicestocloud()[decloud]);
								if (it3 != unserved_devices.end())  //找到了
									unserved_devices.erase(it3);
							}
							
					}
				}
				
			}
		}
	}

}

void Cover::remove_unneccessary_disks(std::vector<std::vector<Disk>>& instanceDisk, int*** x_id0_temp, int** y_d0_temp, int*** z_id0_temp)
{
	//依次考查每一个圆盘
	for (int base = 0; base < basestations_temp.size(); base++) {
		int disk = instanceDisk[base].size() - 1;
		int flag = 0;  
		int maxdisk = -1;   //找到每个基站对应的圆盘中、被选中的、半径最大的那个
		while(disk >= 0) {
			if (flag == 0) {
				if (!instanceDisk[base][disk].get_state()) {   //instanceDisk[base]中第一个被选中的圆盘
					maxdisk = disk;
					flag++;
					disk--;
					y_d0_temp[instanceDisk[base][maxdisk].get_base()][instanceDisk[base][maxdisk].get_device()] = 1;
					continue;
				}
			}
			else 
			{
				if (!instanceDisk[base][disk].get_state()) {
					vector<int>::iterator it1,it2;
					if(instanceDisk[base][disk].get_devicestobase().size()>0)
						for (int i = 0; i < instanceDisk[base][disk].get_devicestobase().size(); i++) {
						//if (it1 == instanceDisk[base][maxdisk].get_devicestobase().end())  //没有找到
							instanceDisk[base][maxdisk].cover_devicestobase(instanceDisk[base][disk].get_devicestobase()[i]);
						}
					if(instanceDisk[base][disk].get_devicestocloud().size()>0)
						for (int j = 0; j < instanceDisk[base][disk].get_devicestocloud().size(); j++) {
					   // if(it2 == instanceDisk[base][maxdisk].get_devicestocloud().end())   //没有找到
							instanceDisk[base][maxdisk].cover_devicestocloud(instanceDisk[base][disk].get_devicestocloud()[j]);
						}
				}
			}
			disk--;
		}

		if (flag != 0) {

			for (int i = 0; i < instanceDisk[base][maxdisk].get_devicestobase().size(); i++)
				x_id0_temp[instanceDisk[base][maxdisk].get_devicestobase()[i]][instanceDisk[base][maxdisk].get_base()][instanceDisk[base][maxdisk].get_device()] = 1;
			for (int j = 0; j < instanceDisk[base][maxdisk].get_devicestocloud().size(); j++)
				z_id0_temp[instanceDisk[base][maxdisk].get_devicestocloud()[j]][instanceDisk[base][maxdisk].get_base()][instanceDisk[base][maxdisk].get_device()] = 1;

		}

	}
}

bool Cover::is_capacity_satisfy(int*** x_id0_temp, int** y_d0_temp, int*** z_id0_temp)
{
	int count = 0;
	vector<int> sum_CPU, sum_BW;
	sum_CPU.resize(m);
	sum_BW.resize(m);

	for (int b = 0; b < m; b++) {
		sum_CPU[b] = 0;
		for (int i = 0; i < n; i++)
			for (int ii = 0; ii < n; ii++)
				sum_CPU[b] += devices[i].get_taskdemand().RCPU * x_id0_temp[i][b][ii];
	}

	for (int b = 0; b < m; b++) {
		sum_BW[b] = 0;
		for (int i = 0; i < n; i++)
			for (int ii = 0; ii < n; ii++)
				sum_BW[b] += devices[i].get_taskdemand().RBW * x_id0_temp[i][b][ii] + devices[i].get_taskdemand().RBW * z_id0_temp[i][b][ii];
	}

	for (int b = 0; b < m; b++)
		if (sum_CPU[b] <= basestations[b].get_capacity().RCPU && sum_BW[b] <= basestations[b].get_capacity().RBW)
			count++;

	if (count < m)
		return false;
	else
		return true;
}


typedef IloArray <IloNumVarArray> IloNumVarArray2;
typedef IloArray <IloNumVarArray2> IloNumVarArray3;


ILOSTLBEGIN

Result Cover::OPT_Result()
{
	clock_t start_, end_;
	start_ = clock();

	int** y_d0 = new int* [m];		// m * n  y_d0[b][i]
	for (int b = 0; b < m; b++)
		y_d0[b] = new int[n];

	int*** x_id0 = new int** [n];	// n * m * n  x_id0[i][b][ii]
	for (int i = 0; i < n; i++) {
		x_id0[i] = new int* [m];
		for (int b = 0; b < m; b++)
			x_id0[i][b] = new int[n];
	}

	int*** z_id0 = new int** [n];	// n * m * n z_id0[i][b][ii]
	for (int i = 0; i < n; i++) {
		z_id0[i] = new int* [m];
		for (int b = 0; b < m; b++)
			z_id0[i][b] = new int[n];
	}

	IloEnv env;
	double solution_value = 0;
	try {
		IloModel model(env);

		IloNumVarArray2 y_d(env, m);
		IloNumVarArray3 x_id(env, n);
		IloNumVarArray3 z_id(env, n);

		// 圆盘 b,i 对应原始变量 y_d[b][i]
		for (IloInt b = 0; b < m; b++)
		{
			y_d[b] = IloNumVarArray(env, n);
			for (IloInt i = 0; i < n; i++)
			{
				y_d[b][i] = IloNumVar(env, 0, 1, ILOINT);
				model.add(y_d[b][i]);
			}
		}
		// 用户 i, 圆盘b,ii 对应原始变量 x_id[i][b][ii]
		for (IloInt i = 0; i < n; i++)
		{
			// 第一维：n个用户
			x_id[i] = IloNumVarArray2(env, m);
			for (IloInt b = 0; b < m; b++) {
				// 二三维：mn 个disk，其中第三维是按照device下标排序
				// 注意与all_disk区分，all_disk第二维是按照半径大小排序
				x_id[i][b] = IloNumVarArray(env, n);
				for (IloInt ii = 0; ii < n; ii++)
				{
					x_id[i][b][ii] = IloNumVar(env, 0, 1, ILOINT);
					model.add(x_id[i][b][ii]);
				}
			}
		}

		// 用户 i, 圆盘b,ii 对应原始变量 z_id[i][b][ii]
		for (IloInt i = 0; i < n; i++)
		{
			// 第一维：n个用户
			z_id[i] = IloNumVarArray2(env, m);
			for (IloInt b = 0; b < m; b++) {
				// 二三维：mn 个disk，其中第三维是按照device下标排序
				// 注意与all_disk区分，all_disk第二维是按照半径大小排序
				z_id[i][b] = IloNumVarArray(env, n);
				for (IloInt ii = 0; ii < n; ii++)
				{
					z_id[i][b][ii] = IloNumVar(env, 0, 1, ILOINT);
					model.add(z_id[i][b][ii]);
				}
			}
		}


		// 约束1
		int dd = 0;   //disk的下标
		IloExpr cover_to_base(env);
		IloExpr cover_to_cloud(env);
		for (IloInt i = 0; i < n; i++) {
			cover_to_base.clear();
			cover_to_cloud.clear();
			for (IloInt b = 0; b < m; b++) {
				dd = device_order[b][i];
				for (IloInt ii = dd; ii < n; ii++) {
					IloInt iii = all_disk[b][ii].get_device();
					cover_to_base += x_id[i][b][iii];
					cover_to_cloud += z_id[i][b][iii];
				}
			}
			model.add(cover_to_base + cover_to_cloud >= 1);
		}

		// 约束2
		IloExpr devicecovered_tobase_chose_disk(env);
		vector<int> dcovered;
		for (IloInt b = 0; b < m; b++) {
			for (IloInt i = 0; i < n; i++) {
				dd = device_order[b][i];
				dcovered = all_disk[b][dd].get_devices();
				for (IloInt ii = 0; ii < dcovered.size(); ii++) {
					IloInt iii = dcovered[ii];
					devicecovered_tobase_chose_disk.clear();
					devicecovered_tobase_chose_disk = y_d[b][i] - x_id[iii][b][i];
					model.add(devicecovered_tobase_chose_disk >= 0);
				}
			}
		}

		//约束3
		IloExpr devicecovered_tocloud_chose_disk(env);
		for (IloInt b = 0; b < m; b++) {
			for (IloInt i = 0; i < n; i++) {
				dd = device_order[b][i];
				dcovered = all_disk[b][dd].get_devices();
				for (IloInt ii = 0; ii < dcovered.size(); ii++) {
					IloInt iii = dcovered[ii];
					devicecovered_tocloud_chose_disk.clear();
					devicecovered_tocloud_chose_disk = y_d[b][i] - z_id[iii][b][i];
					model.add(devicecovered_tocloud_chose_disk >= 0);
				}
			}
		}

		//约束4
		IloExpr diskCPU_capacity(env);
		IloExpr dcpu(env);
		for (IloInt b = 0; b < m; b++) {
			for (IloInt i = 0; i < n; i++) {
				diskCPU_capacity.clear();
				dcpu.clear();
				dd = device_order[b][i];
				dcovered = all_disk[b][dd].get_devices();
				for (IloInt ii = 0; ii < dcovered.size(); ii++)
				{
					IloInt iii = dcovered[ii];
					dcpu += devices[iii].get_taskdemand().RCPU * x_id[iii][b][i];
				}

				diskCPU_capacity = all_disk[b][dd].get_capacity().RCPU * y_d[b][i]- dcpu;
				model.add(diskCPU_capacity >= 0);
			}
		}

        //约束5
		IloExpr diskBW_capacity(env);
		IloExpr dbw_tobase(env);
		IloExpr dbw_tocloud(env);
		for (IloInt b = 0; b < m; b++) {
			for (IloInt i = 0; i < n; i++) {
				diskBW_capacity.clear();
				dbw_tobase.clear();
				dbw_tocloud.clear();
				dd = device_order[b][i];
				dcovered = all_disk[b][dd].get_devices();
				for (IloInt ii = 0; ii < dcovered.size(); ii++)
				{
					IloInt iii = dcovered[ii];
					dbw_tobase += devices[iii].get_taskdemand().RBW * x_id[iii][b][i];
					dbw_tocloud += devices[iii].get_taskdemand().RBW * z_id[iii][b][i];
				}
				diskBW_capacity = all_disk[b][dd].get_capacity().RBW * y_d[b][i] - dbw_tobase - dbw_tocloud;
				model.add(diskBW_capacity >=0);
			}
		}

		//约束6
		IloExpr onedisk_eachbase(env);
		for (IloInt b = 0; b < m; b++) {
			onedisk_eachbase.clear();
			for (IloInt i = 0; i < n; i++) {
				onedisk_eachbase += y_d[b][i];
			}
			model.add(onedisk_eachbase <= 1);
		}
		
		// 目标函数
		IloExpr obj(env);
	    for (IloInt b = 0; b < m; b++)
			for (IloInt i = 0; i < n; i++){
				dd = device_order[b][i];
				obj += all_disk[b][dd].get_power()* COVER_TIME * y_d[b][i];
			}	
			   
		for (IloInt i = 0; i < n; i++) {
			for (IloInt b = 0; b < m; b++) {
			   dd = device_order[b][i];
			   for (IloInt ii = dd; ii < n; ii++) {
					IloInt iii = all_disk[b][ii].get_device();
					obj += (E_base_ex[b][i]+E_base_tr[b][i]) * x_id[i][b][iii];
			   }
		    }
		}
		
		for (IloInt i = 0; i < n; i++){
			for (IloInt b = 0; b < m; b++) {
					dd = device_order[b][i];
					for (IloInt ii = dd; ii < n; ii++) {
						IloInt iii = all_disk[b][ii].get_device();
						obj += (E_cloud_ex[i] + E_cloud_sum_tr[b][i]) * z_id[i][b][iii];
					}
			}
		}

		model.add(IloMinimize(env, obj));			
				
		// 创建求解对象
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		solution_value = cplex.getObjValue();
		printf("Solution value = %lf\n", solution_value);	

		for (IloInt b = 0; b < m; b++)
		{
			IloNumArray yy_d(env);
			cplex.getValues(yy_d, y_d[b]);
			for (int i = 0; i < n; i++) {
				y_d0[b][i] = yy_d[i];
			}
		}
		for (IloInt i = 0; i < n; i++)
		{
			
			for (IloInt b = 0; b < m; b++) {
				IloNumArray xx_id(env);
				cplex.getValues(xx_id, x_id[i][b]);
				for (IloInt ii = 0; ii < n; ii++)
					x_id0[i][b][ii] = xx_id[ii];
			}
		}
		for (IloInt i = 0; i < n; i++)
		{
			
			for (IloInt b = 0; b < m; b++) {
				IloNumArray zz_id(env);
				cplex.getValues(zz_id, z_id[i][b]);
				for (IloInt ii = 0; ii < n; ii++)
					z_id0[i][b][ii] = zz_id[ii];
			}
		}

	}
	catch (IloException& e) { cerr << "Concert exception caught:" << e << endl; }
	catch (...) { cerr << "Unknuwn exception caught" << endl; }
	env.end();

	//print_xyz(y_d0, x_id0,z_id0);

	Result r = Result();

	r = cal_result(y_d0, x_id0, z_id0);

	
	for (int b = 0; b < m; b++)
		delete[] y_d0[b];
	delete[] y_d0;

	for (int i = 0; i < n; i++) {
		for (int b = 0; b < m; b++)
			delete[] x_id0[i][b];
		delete[] x_id0[i];
	}
	delete[] x_id0;

	for (int i = 0; i < n; i++) {
		for (int b = 0; b < m; b++)
			delete[] z_id0[i][b];
		delete[] z_id0[i];
	}
	delete[] z_id0;

	end_ = clock();
	r.set_time((double)end_ - start_);
	return r;
}


Result Cover::Greedy_Result()
{
	using namespace std;

	clock_t start_, end_;
	start_ = clock();

	int** y_d0 = new int* [m];		// m * n  y_d0[b][i]
	for (int b = 0; b < m; b++) {
		y_d0[b] = new int[n];
		memset(y_d0[b], 0, sizeof(int) * n); 
	}

	int*** x_id0 = new int** [n];	// n * m * n  x_id0[i][b][ii]
	for (int i = 0; i < n; i++) {
		x_id0[i] = new int* [m];
		for (int b = 0; b < m; b++) {
			x_id0[i][b] = new int[n];
			memset(x_id0[i][b], 0, sizeof(int) * n);
		}
	}

	int*** z_id0 = new int** [n];	// n * m * n z_id0[i][b][ii]
	for (int i = 0; i < n; i++) {
		z_id0[i] = new int* [m];
		for (int b = 0; b < m; b++) {
			z_id0[i][b] = new int[n];
			memset(z_id0[i][b], 0, sizeof(int) * n);
		}
	}

	std::vector<Device> unserved_devices=devices;  //初始化未被服务的devices

	std::vector<std::vector<double>> power_fulldisk;   //每个disk的覆盖功耗
	power_fulldisk.resize(m);
	for (int b = 0; b < m; b++) {
		power_fulldisk[b].resize(n);
	}

	int** flag=new int* [m];  //标记disk是否开放，并初始化
	for (int b = 0; b < m; b++)
		flag[b] = new int[n];

	//清空前一种方法的结果
	for (int b = 0; b < m; b++) {
		basestations[b].clear_direct_served_devices();
		basestations[b].clear_indirect_served_devices();
		for (int i = 0; i < n; i++) {
			all_disk[b][i].delete_devicetobase();
			all_disk[b][i].delete_devicetocloud();
		}
	}


	Resource R_temp=Resource();  //用于更新basestation[b_min]的资源容量

	//按CPU需求降序排列devices
	std::sort(unserved_devices.begin(), unserved_devices.end(), sortdevice_byCPU);


	std::vector<Disk> chosed_disk;  //记录每一轮循环most energy effectiveness、且同心的disks中半径最大的disk

	//测试基站的BW资源是否足够所有devices
	int total_capacityCPU = 0;
	int total_capacityBW = 0;
	int total_demandCPU = 0;
	int total_demandBW = 0;
	for (int b = 0; b < m; b++) {
		total_capacityCPU += basestations[b].get_capacity().RCPU;
		total_capacityBW += basestations[b].get_capacity().RBW;
	}
	cout <<'\n'<< "total_capacityCPU:" << total_capacityCPU << '\t' << "total_capacityBW:" << total_capacityBW;
	for (int i = 0; i < n; i++) {
		total_demandCPU += devices[i].get_taskdemand().RCPU;
		total_demandBW += devices[i].get_taskdemand().RBW;
	}
	cout << '\n' << "total_demandCPU:" << total_demandCPU << '\t' << "total_demandBW:" << total_demandBW;

	if (total_capacityBW >= total_demandBW) {
		cout << '\n' << "BW is enough!" << '\n';

		int cycle = 0;  
		while (unserved_devices.size() > 0) {
			cycle++;

			vector<std::vector<Disk> > full_disk = construct_fulldisk(unserved_devices);  //构造fulldisk

			//每一个full_disk的总能耗
			vector<std::vector<double> > energy_fulldisk;
			energy_fulldisk.resize(full_disk.size());  

			//每一个full_disk的energy effectiveness
			vector<std::vector<double> > energy_effectiveness;
			energy_effectiveness.resize(full_disk.size());  

			//初始化每个圆盘的覆盖功耗
			if (chosed_disk.size() == 0) {
				for (int b = 0; b < full_disk.size(); b++)
					for (int i = 0; i < full_disk[b].size(); i++)
						power_fulldisk[b][i] = full_disk[b][i].get_power();

			}

			//计算每一个full_disk的energy effectiveness
			for (int b = 0; b < full_disk.size(); b++) {
				for (int i = 0; i < full_disk[b].size(); i++) {
					//计算圆盘覆盖能耗
					energy_fulldisk[b].push_back(power_fulldisk[b][i] * COVER_TIME);
					//计算迁移到basestation的devices能耗
					for (int j1 = 0; j1 < full_disk[b][i].get_devicestobase().size(); j1++) {
						energy_fulldisk[b][i] += E_base_ex[b][full_disk[b][i].get_devicestobase()[j1]] + E_base_tr[b][full_disk[b][i].get_devicestobase()[j1]];
					}
					//计算迁移到cloud的devices能耗
					for (int j2 = 0; j2 < full_disk[b][i].get_devicestocloud().size(); j2++) {
						energy_fulldisk[b][i] += E_cloud_ex[full_disk[b][i].get_devicestocloud()[j2]] + E_cloud_sum_tr[b][full_disk[b][i].get_devicestocloud()[j2]];
					}
					//计算能效
					if ((full_disk[b][i].get_devicestobase().size() + full_disk[b][i].get_devicestocloud().size()) != 0)
						energy_effectiveness[b].push_back(energy_fulldisk[b][i] / (full_disk[b][i].get_devicestobase().size() + full_disk[b][i].get_devicestocloud().size()));
					else energy_effectiveness[b].push_back(DBL_MAX);
				}
			}

			//求full_disk中能效最小的disk
			int b_min = 0;
			int i_min = 0;
			double energy_min = energy_effectiveness[0][0];
			for (int b = 0; b < full_disk.size(); b++) {
				for (int i = 0; i < full_disk[b].size(); i++) {
					if (energy_effectiveness[b][i] < energy_min) {
						energy_min = energy_effectiveness[b][i];
						b_min = b;
						i_min = i;
					}
				}
			}

			cout << "most energy-effectiveness is: "<<b_min<<" "<<i_min<<'\n';

			//如果disk未开放，则开放它，并判断是否将它放入chosed_disk
			int count = 0;
			if (full_disk[b_min][i_min].get_state()) {
				full_disk[b_min][i_min].open_disk();
				flag[b_min][i_min] = 1;              //标记full_disk[b_min][i_min]开放
				count = count_1(flag[b_min]);  

				if (count == 1) {    //说明这是以b_min为中心的圆盘中第一个开放的圆盘
					chosed_disk.push_back(full_disk[b_min][i_min]);
					basestations[b_min].set_last_disk(dist_index[b_min][i_min]);  //设置最后选择的圆盘所对应device的index

					//更新同心圆盘的覆盖功耗
					for (int i = 0; i < n; i++) {
						if (i <= i_min) {   //full_disk是all_disk返回的，同心的圆盘是按半径从小到大排序的
							power_fulldisk[b_min][i] = 0;  //与full_disk[b_min][i_min]同心，半径小于等于full_disk[b_min][i_min]半径的圆盘的power置为0
						}
						else power_fulldisk[b_min][i] = full_disk[b_min][i].get_power() - full_disk[b_min][i_min].get_power();
					}
				}

				if (count > 1) {  //说明以b_min为中心的圆盘已经开放多个
					if (device_order[b_min][basestations[b_min].get_chosed_disk()] < i_min) { //当前的圆盘比basestation最后选中的圆盘半径大
						for (std::vector<Disk>::iterator pos = chosed_disk.begin(); pos != chosed_disk.end(); pos++)
							if ((*pos).get_base() == b_min) {
								chosed_disk.erase(pos);
								break;
							}

						//重新设置选中的圆盘
						basestations[b_min].set_last_disk(dist_index[b_min][i_min]);
						chosed_disk.push_back(full_disk[b_min][i_min]);

						//更新同心圆盘的覆盖功耗
						for (int i = 0; i < n; i++) {
							if (i <= i_min) {  
								power_fulldisk[b_min][i] = 0;  
							}
							else power_fulldisk[b_min][i] = full_disk[b_min][i].get_power() - full_disk[b_min][i_min].get_power();
						}
					}

				}
			}


			//更新basestations[b_min]对应的直接、间接服务的devices
			for (auto detobase : full_disk[b_min][i_min].get_devicestobase()) {
				basestations[b_min].direct_served_device(detobase);
			}
			for (auto detocloud : full_disk[b_min][i_min].get_devicestocloud())
				basestations[b_min].indirect_served_device(detocloud);


			//更新未服务的devices
			for (int j = 0; j < full_disk[b_min][i_min].get_devicestobase().size(); j++) {
				for (vector<Device>::iterator iter = unserved_devices.begin(); iter != unserved_devices.end(); iter++) {
					if ((*iter).get_index() == full_disk[b_min][i_min].get_devicestobase()[j]) {
						unserved_devices.erase(iter);
						break;
					}
				}
			}
			for (int j = 0; j < full_disk[b_min][i_min].get_devicestocloud().size(); j++) {
				for (vector<Device>::iterator iter = unserved_devices.begin(); iter != unserved_devices.end(); iter++) {
					if ((*iter).get_index() == full_disk[b_min][i_min].get_devicestocloud()[j]) {
						unserved_devices.erase(iter);
						break;
					}
				}
			}


			//更新basestations[b_min]的CPU和BW资源容量
			R_temp = Resource(basestations[b_min].get_capacity().RCPU, basestations[b_min].get_capacity().RBW);
			for (int j = 0; j < full_disk[b_min][i_min].get_devicestobase().size(); j++) {
				R_temp.RCPU -= devices[full_disk[b_min][i_min].get_devicestobase()[j]].get_taskdemand().RCPU;
				R_temp.RBW -= devices[full_disk[b_min][i_min].get_devicestobase()[j]].get_taskdemand().RBW;
			}
			for (int j = 0; j < full_disk[b_min][i_min].get_devicestocloud().size(); j++)
				R_temp.RBW -= devices[full_disk[b_min][i_min].get_devicestocloud()[j]].get_taskdemand().RBW;
			basestations[b_min].set_capacity(R_temp);

		
			vector<std::vector<double> >().swap(energy_fulldisk);
			vector<std::vector<double> >().swap(energy_effectiveness);
			
		}
		
	}

	
	
	//根据chosed_disk确定y_d0的值
	for (int j = 0; j < chosed_disk.size(); j++) {
		y_d0[chosed_disk[j].get_base()][chosed_disk[j].get_device()] = 1;
	}
	//根据basestation的direct_served_devices和direct_served_devices确定x_id0和z_id0的值
	for (int b = 0; b < m; b++) {
		for (int i = 0; i < basestations[b].get_direct_served_devices_num(); i++) {
			for (int j = 0; j < chosed_disk.size(); j++) {
				if (chosed_disk[j].get_base() == b)
					x_id0[basestations[b].get_direct_served_devices()[i]][b][chosed_disk[j].get_device()] = 1;
			}
		}
		for (int i = 0; i < basestations[b].get_indirect_served_devices_num(); i++) {
			for (int j = 0; j < chosed_disk.size(); j++) {
				if(chosed_disk[j].get_base()==b)
					z_id0[basestations[b].get_indirect_served_devices()[i]][b][chosed_disk[j].get_device()] = 1;
			}
		}
	}
	
	//print_xyz(y_d0, x_id0, z_id0);

	Result r=Result();

	r = cal_result(y_d0, x_id0, z_id0);


	for (int b = 0; b < m; b++)
		delete[] y_d0[b];
	delete[] y_d0;

	for (int i = 0; i < n; i++) {
		for (int b = 0; b < m; b++)
			delete[] x_id0[i][b];
		delete[] x_id0[i];
	}
	delete[] x_id0;

	for (int i = 0; i < n; i++) {
		for (int b = 0; b < m; b++)
			delete[] z_id0[i][b];
		delete[] z_id0[i];
	}
	delete[] z_id0;

	for (int b = 0; b < m; b++)
		delete[] flag[b];
	delete[] flag;

	end_ = clock();
	r.set_time((double)end_ - start_);
	return r;
}

//构造满圆盘
std::vector<std::vector<Disk>> Cover::construct_fulldisk(std::vector<Device>& ordered_unserved_devices)
{
	vector<std::vector<Disk>> full_D;
	full_D.resize(m);

	//初始化每一个basestation的本轮直接服务devices和间接服务devices
	std::vector<std::vector<int>> current_direct_served_devices;
	std::vector<std::vector<int>> current_indirect_served_devices;

	current_direct_served_devices.resize(m);
	current_indirect_served_devices.resize(m);

	//初始化每一个disk的本轮覆盖的devices
	for(int b=0;b<m;b++){
		for (int i = 0; i < n; i++) {
			all_disk[b][i].delete_devicetobase();
			all_disk[b][i].delete_devicetocloud();
		}
	}

	int* bCPU = new int[m];
	int* bBW = new int[m];
	for (int b = 0; b < m; b++) {
		bCPU[b] = basestations[b].get_capacity().RCPU;
		bBW[b] = basestations[b].get_capacity().RBW;
	}
	
	//遍历ordered_unserved_devices
	for (int b = 0; b < m; b++) {
		for (auto dev : ordered_unserved_devices) {
			if (bBW[b] >= dev.get_taskdemand().RBW) {
				if (bCPU[b] >= dev.get_taskdemand().RCPU) {  
					current_direct_served_devices[b].push_back(dev.get_index());
					bCPU[b] -= dev.get_taskdemand().RCPU;
					bBW[b] -= dev.get_taskdemand().RBW;
				}
				else {    
					current_indirect_served_devices[b].push_back(dev.get_index());
					bBW[b] -= dev.get_taskdemand().RBW;
				}
			}
		}
	}


	//求被每一个disk覆盖、迁移到云或者basestation的devices
	for (int b = 0; b < m; b++) 
		for (int i = 0; i < n; i++) {
			for (auto j1 : current_direct_served_devices[b]) {
				for(auto iter1:all_disk[b][i].get_devices())
					if(iter1==j1)
					all_disk[b][i].cover_devicestobase(j1);
			}
			for(auto j2:current_indirect_served_devices[b]){
				for (auto iter2 : all_disk[b][i].get_devices())
					if(iter2==j2)
					all_disk[b][i].cover_devicestocloud(j2);
		}
	}



	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++) {
			if (all_disk[b][i].get_devicestobase().size() != 0 || all_disk[b][i].get_devicestocloud().size() != 0)
				full_D[b].push_back(all_disk[b][i]);
		}
	}


	delete[] bCPU;
	delete[] bBW;

	vector<std::vector<int> >().swap(current_direct_served_devices);
	vector<std::vector<int> >().swap(current_indirect_served_devices);

	return all_disk;
}

Result Cover::LP_Result()
{
	Result r;
	return r;
}

Result Cover::DP_Result()
{
	Result r;
	return r;
}


Result Cover::cal_result(int** y_d, int*** x_id, int*** z_id) {
	using namespace std;
	Result r = Result();

	double eng_cspt = 0;   
	double ratio_base = 0; 
	double cpu_utilization = 0;   
	double bw_utilization = 0;   
	double maxdevices_bybase = 0;  
	double avedevices_bybase = 0; 
	double averadius = 0;  
	double maxadius = 0; 
	double numselecteddisk = 0;

	//根据实际需求计算所需的结果
	int* d_direct = new int[m];			
	int* d_indirect = new int[m];     
	for (int b = 0; b < m; b++) {
		d_direct[b] = 0;  
		d_indirect[b] = 0; 
	}

	int num_base_dirser=0;
	int num_diskselected = 0;  
		
	
	int* sum_allocpu = new int[m];    
	int* sum_allobw = new int[m];   
	for (int b = 0; b < m; b++) {
		sum_allocpu[b] = 0;
		sum_allobw[b] = 0;
	}

	for (int b = 0; b < m; b++) {
		basestations[b].clear_direct_served_devices();
		basestations[b].clear_indirect_served_devices();
	}

	for (int b = 0; b < m; b++) {
		for (int i = 0; i < n; i++)
			if (y_d[b][i] == 1) {
				// 圆盘bi被选中
				basestations[b].set_last_disk(i);
				eng_cspt += all_disk[b][device_order[b][i]].get_power()* COVER_TIME;
				num_diskselected++;

				averadius += get_distance(b,i);
				if (maxadius < get_distance(b, i))
					maxadius = get_distance(b,i);

				for (int ii = 0; ii < n; ii++) {
					if (x_id[ii][b][i] == 1) {
						basestations[b].direct_served_device(ii);
                        d_direct[b]++;
						eng_cspt += E_base_ex[b][ii] + E_base_tr[b][ii];
						devices[ii].set_served_by(b);
						sum_allocpu[b] += devices[ii].get_taskdemand().RCPU;
						sum_allobw[b] += devices[ii].get_taskdemand().RBW;
					}
					if (z_id[ii][b][i] == 1) {
						basestations[b].indirect_served_device(ii);
						d_indirect[b]++;
						eng_cspt += E_cloud_ex[ii] + E_cloud_sum_tr[b][ii];
						devices[ii].set_served_by(-1);
						sum_allobw[b] += devices[ii].get_taskdemand().RBW;
					}	
				}	
			}
		
		cpu_utilization += (double)sum_allocpu[b] / basestations[b].get_capacity().RCPU;
		bw_utilization += (double)sum_allobw[b] / basestations[b].get_capacity().RBW;
	}
	
	cpu_utilization /= num_diskselected;
	bw_utilization /= num_diskselected;

	for (int b = 0; b < m; b++) {
		ratio_base += d_direct[b];

		if (d_direct[b] > maxdevices_bybase)
			maxdevices_bybase = d_direct[b];

		avedevices_bybase += d_direct[b];
		if (basestations[b].get_direct_served_devices_num() > 0)
			num_base_dirser++;

	}

	ratio_base /= n;
	avedevices_bybase /= num_base_dirser;
	averadius /= num_diskselected;
	numselecteddisk = num_diskselected;

	r.set_eng_cspt(eng_cspt);
	r.set_ratiobase(ratio_base);
	r.set_cpu_uti(cpu_utilization);
	r.set_bw_uti(bw_utilization);
	r.set_maxdevices_bybase(maxdevices_bybase);
	r.set_avedevices_bybase(avedevices_bybase);
	r.set_averadius(averadius);
	r.set_maxadius(maxadius);
	r.set_numselecteddisk(numselecteddisk);

	delete[] d_direct;
	delete[] d_indirect;
	delete[] sum_allocpu;
	delete[] sum_allobw;
	return r;
}

bool Cover::is_assined_all()
{
	return false;
}

int Cover::find_device_index_in_all_disk(int base, int device) {
	for (int i = 0; i < n; i++)
		if (dist_index[base][i] == device)
			return i;
	return -1;
}

void Result::print_titled() {
	std::cout << "energy\t\ttime\t\tratio\t\tcpuutilization\t\tbwutilization\n";
	std::cout << setprecision(5) << eng_cspt << "\t\t" << time << "\t\t";
	std::cout << setprecision(5) << ratio_base << "\t\t" <<  cpu_utilization << "\t\t"<<bw_utilization<<"\t\t"<<"\n";
}

std::string Result::print_untitled() {
	string r;
	r += to_string(eng_cspt) + ' ';
	r += to_string(time) + ' ';
	r += to_string(ratio_base) + ' ';
	//r += to_string(variance) + ' ';
	r += to_string(cpu_utilization) + ' ';
	r += to_string(bw_utilization) + '\n';
	cout << r;
	return r;
}
