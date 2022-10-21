

#pragma once
#pragma once
#ifndef PREDEFINE_H
#define PREDEFINE_H
#include<vector>
#include<string>
#include <iostream>

#define C 1  
#define K 2   
#define ELECTRICITY_INTENSITY (0.06*3.6*1e6)/(1024*1024)  
#define COVER_TIME 1  
#define INCREASE 1    
#define ERROR 1      


enum class type { POINT, DEVICE, BASESTATION, CLOUD, DISK };
class Point;
class Device;
class BaseStation;
class Cloud;
class Disk;
class Result;


// 点类
class Point {
protected:
	double loc[2] = {};                         //坐标
	int index = -1;								
	type T;
public:
	Point(double* l, type type);
	Point();           
	~Point() {}        
	double cal_distance(const Point& p);		
	void set_location(double lx, double ly);	
	const double* get_loc() const;				
	const int get_index();						
	bool operator<(const Point& p) const;		
};

// Resource类 
class Resource {
public:
	int RCPU;				// CPU
	int RBW;				// 带宽
public:
	Resource() { RCPU = 0; RBW = 0; }
	Resource(int c, int bw) { RCPU = c; RBW = bw; }
	void copy(const Resource& r) { RCPU = r.RCPU; RBW = r.RBW; }
	void print_() const { std::cout << " CPU=" << RCPU << " BW=" << RBW << ' '; }
};

// Device类
class Device : public Point {
private:
	int y = 0;					
	int served_by = -2;			
	int data;  
	Resource taskdemand;				
	double e_unitbit=0;  
	double e_unitm=0;   

public:

	const int get_data() { return data; }  
	Resource& get_taskdemand() { return taskdemand; }  
	const double get_e_unitbit() { return e_unitbit; }  
	const double get_e_unitm() { return e_unitm; } 

	Device(double* location, int i, int da, const Resource& t, double eb, double em);		

	void set_served_by(int base) { served_by = base; }		
			
	const int get_served_by() { return served_by; }			
	void print_device();									
};

// Disk类
class Disk {
private:
	int device;												
	int base;													
	Resource capacity;										
	double e_d = DBL_MAX;					
	std::vector<int> devices_covered_by_disk;
	std::vector<int> devices_tobase_covered_by_disk;			
	std::vector<int> devices_tocloud_covered_by_disk;              
	double power;																		
	bool state = true;										
	static double theta;									

public:
	Disk();
	Disk(int base_, int device_);
	Disk(int base_, int device_, const double& radii, const Resource& capacity_);		
	~Disk() {}

	// 以下函数都是获取私有变量的函数
	const int get_device() { return device; }		
	const int get_base() { return base; }			
	const Resource &get_capacity() { return capacity; }	
	const double get_power() { return power; }
	std::vector<int> get_devices() { return devices_covered_by_disk; }
	std::vector<int> get_devicestobase() { return devices_tobase_covered_by_disk; }
	std::vector<int> get_devicestocloud() { return devices_tocloud_covered_by_disk; }
	const double get_ed() { return e_d; }
	const bool get_state() { return state; }
	static double get_theta() { return theta; }

	static void set_theta(double a) { theta = a; }

	// 设置相关变量的值
	void set_capacity(Resource &c) { capacity = c; }
	void set_power(double p) { power = p; }
	void set_ed(double ed) { e_d = ed; }
	void open_disk() { state = false; }   

	
	double cal_power(const double& radii);	// 根据半径计算覆盖功耗
	void cover_devices(const int i) { devices_covered_by_disk.push_back(i); };
	void cover_devices(const std::vector<int>& d, const int& i);
	void delete_device(const int i) { devices_covered_by_disk.erase(devices_covered_by_disk.begin() + i); }
	void cover_devicestobase(const int i) { devices_tobase_covered_by_disk.push_back(i); };										
	void cover_devicestocloud(const int i) { devices_tocloud_covered_by_disk.push_back(i); };   
	void delete_devicetobase() { devices_tobase_covered_by_disk.clear(); }
	void delete_devicetocloud(){ devices_tocloud_covered_by_disk.clear(); }

	void print_disk() const;	// 控制台输出圆盘相关参数
	
	//重定义操作符==
	bool operator ==(const Disk &di) const {
		return base == di.base;
	}

	void remove_covered_devices(std::vector<int>& devs);  
};


class BaseStation :public Point {	
private:
	Resource capacity;			
	double weight;			
	int za;					
	std::vector<int> direct_served_devices;
	std::vector<int> indirect_served_devices; 
	int chosed_disk_index = -1;	
	double basef;												
	double basepf;											


public:
	BaseStation(double* location, Resource& c, int b, double f_, double pf_);			
	BaseStation();

	Resource get_capacity() { return capacity; }
	int get_chosed_disk() { return chosed_disk_index; }	
	size_t get_direct_served_devices_num() { return direct_served_devices.size(); }
	size_t get_indirect_served_devices_num() { return indirect_served_devices.size(); }   
	const std::vector<int>& get_direct_served_devices() { return direct_served_devices; }
	const std::vector<int>& get_indirect_served_devices() { return indirect_served_devices; }

	const double get_basef() { return basef; }  
	const double get_basepf() { return basepf; }  

	void set_last_disk(int index) { chosed_disk_index = index; }			
	void direct_served_device(int device) { direct_served_devices.push_back(device); }
	void indirect_served_device(int device) { indirect_served_devices.push_back(device); }
	void clear_direct_served_devices() { direct_served_devices.clear(); }  
	void clear_indirect_served_devices() { indirect_served_devices.clear(); }

	void set_capacity(Resource r) { capacity.RCPU = r.RCPU; capacity.RBW = r.RBW; }
	
	void print_basestation();

};

class Cloud {
private:
	Resource capacity;			
	std::vector<int> served_devices;	
	double cloudfvm;												
	double cloudpfvm;												

public:
				
	Cloud();

	void set_cloud(Resource& c, double fvm_, double pfvm_);
	Resource get_capacity() { return capacity; }
	double get_cloudfvm() { return cloudfvm; }
	double get_cloudpfvm() { return cloudpfvm; }

	size_t get_served_devices_num() { return served_devices.size(); }

	void serve_device(int device) { served_devices.push_back(device); }
	void clear_served_devices() { served_devices.clear(); }

	const std::vector<int>& get_served_devices() { return served_devices; }

};

// Cover类
class Cover {
private:
	int m;  
	int n;   
	double result = 0;
	Cloud cloud;
	std::vector<Device> devices;
	std::vector<BaseStation> basestations;

	std::vector<BaseStation> basestations_temp;


	double** E_base_ex = nullptr;					
	double** E_base_tr = nullptr;					
	double* E_cloud_ex = nullptr;					
	double* E_cloud_tr = nullptr;					
	double** E_cloud_sum_tr = nullptr;				

	std::vector<std::vector<Disk>> all_disk;			
	double** distance = nullptr;			
	int** dist_index = nullptr;		
	int** device_order = nullptr;			

	double* alpha_i = nullptr;  //对偶变量
	double*** beta_id=nullptr;   //对偶变量
	double*** gamma_id = nullptr;  //对偶变量
	double** delta_d = nullptr; //对偶变量
	double** epsilon_d = nullptr;//对偶变量
	double* mu_b = nullptr; //对偶变量

public:

	Cover();
	Cover(int m_, int n_);
	Cover(std::string fname);                                                                                                                                                                                                                                                                                                                                                
	~Cover() {
		for (int i = 0; i < n; i++) {
			for (int b = 0; b < m; b++) {
	            delete[] beta_id[i][b];
				delete[] gamma_id[i][b];
			}
			delete[] beta_id[i];
			delete[] gamma_id[i];
		}
		for (int b = 0; b < m; b++) {
			delete[] E_base_ex[b];
			delete[] E_base_tr[b];
			delete[] E_cloud_sum_tr[b];
			delete[] distance[b];
			delete[] dist_index[b];
			delete[] device_order[b];
			delete[] delta_d[b];
			delete[] epsilon_d[b];
		}
		delete[] E_base_ex;
		delete[] E_base_tr;
		delete[] E_cloud_ex;
		delete[] E_cloud_tr;
		delete[] E_cloud_sum_tr;
		delete[] distance;
		delete[] dist_index;
		delete[] device_order;
		delete[] alpha_i;
		delete[] beta_id;
		delete[] gamma_id;
		delete[] delta_d;
		delete[] epsilon_d;
		delete[] mu_b;
	}

	void read_data(std::string fname);

	int get_m() { return m; }
	int get_n() { return n; }
	void set_all_disk(); //构造圆盘
	void update_all_disk();		
	void set_devices_order();	
	void set_theta(double theta);
	void print_theta() {
		for (int i = 0; i < 5; i++)
			std::cout << all_disk[0][i].get_theta() << '\n';
	}
	void cal_all_distance(); // 计算基站与用户两两之间的距离
	void cal_all_offbase_energy();  //计算用户迁移到基站的能耗
	void cal_all_offcloud_energy();   //计算用户迁移到云的能耗

	void order_dist();  //排序算法，设置dist_index

	const double get_distance(int b_index, int d_index);
	int get_order(int b_index, int order);  //返回距离基站第order远的设备index
	const std::vector<std::vector<Disk>>& get_alldisk() { return all_disk; }

	void print_test();
	void print_allDisk(const std::vector<std::vector<Disk>>& disks);
	template <typename type>
	void print_xyz(type** y_d0, type*** x_id0, type*** z_id0);
	void print_result(Result r);

	// 所有算法的实现
	Result PrimalDual_Result();	
	std::vector<std::vector<Disk>> construct_instanceDisk(Disk& max_disk,std::vector<int>& devices_served_by_maxdisk); //构造新实例
	std::vector<int> construct_unserved_devices(Disk& max_disk, std::vector<int>& devices_served_by_maxdisk);//构造未被服务的devices
	bool is_instanceDisk_feasible(std::vector<std::vector<Disk>>& instanceDisk,std::vector<int> unserved_devices);
	void instancePD(std::vector<int> unserved_devices, std::vector<std::vector<Disk>>& instanceDisk, int*** x_id0_temp, int** y_d0_temp, int*** z_id0_temp, double* alpha_i_temp, double*** beta_id_temp, double*** gamma_id_temp, double** delta_d_temp,double** epsilon_d_temp,double* mu_b_temp);
	void remove_unneccessary_disks(std::vector<std::vector<Disk>>& instanceDisk, int*** x_id0_temp, int** y_d0_temp, int*** z_id0_temp);
	bool is_capacity_satisfy(int*** x_id0_temp, int** y_d0_temp, int*** z_id0_temp);

	Result OPT_Result();	// 原整数规划
	Result LP_Result();		
	Result DP_Result();		

	Result Greedy_Result();	
	std::vector<std::vector<Disk>> construct_fulldisk(std::vector<Device>& ordered_unserved_devices);//构造fulldisk set

	Result cal_result(int** y_d, int*** x_id, int*** z_id);	
	bool is_assined_all();
	int find_device_index_in_all_disk(int base, int device);	
	int count_1(int* arr) {
		int count = 0;
		for (int j = 0; j < n; j++)
			if (arr[j] == 1)
				count++;
		return count;
	}
};

// Result类
class Result {
private:
	double eng_cspt = DBL_MAX;	
	double time = 0;		
	double ratio = 0;		
	double ratio_base = 0;		
	double cpu_utilization = 0;  
	double bw_utilization = 0;  
	double averageradius = 0;
	double maxdevices_bybase = 0;    
	double avedevices_bybase = 0;      
	double averadius = 0; 
	double maxadius = 0;
	double numselecteddisk = 0;

public:
	void set_eng_cspt(double eng) { eng_cspt = eng; }
	void set_time(double t) { time = t; }
	void set_ratio(double r) { ratio = r; }
	void set_ratiobase(double rbase) { ratio_base = rbase; }
	void set_cpu_uti(double c) { cpu_utilization= c; }
	void set_bw_uti(double b) { bw_utilization = b; }
	void set_maxdevices_bybase(double mdn) { maxdevices_bybase = mdn; }
	void set_avedevices_bybase(double adn) { avedevices_bybase = adn; }
	void set_averadius(double aa) { averadius = aa; }
	void set_maxadius(double ma) { maxadius = ma; }
	void set_numselecteddisk(double ns) { numselecteddisk = ns; }

	double get_eng_cspt() { return eng_cspt; }
	double get_time() { return time; }
	double get_ratio() { return ratio; }
	double get_ratiobase() { return ratio_base; }
	double get_cpu_uti() { return cpu_utilization; }
	double get_bw_uti() { return bw_utilization; }
	double get_maxdevices_bybase() { return maxdevices_bybase; }
	double get_avedevices_bybase() { return avedevices_bybase; }
	double get_averadius() { return averadius; }
	double get_maxadius() { return maxadius; }
	double get_numselecteddisk() { return numselecteddisk;  }

	void print_titled();
	std::string print_untitled();
	void average_result(int inst_num) {
		eng_cspt /= inst_num;
		time /= inst_num;
		ratio /= inst_num;
		ratio_base /= inst_num;
		cpu_utilization /= inst_num;
		bw_utilization /= inst_num;
		maxdevices_bybase /= inst_num;
		avedevices_bybase /= inst_num;
		averadius /= inst_num;
		maxadius /= inst_num;
		numselecteddisk/=inst_num;
	}

	void add_result(Result& add) {
		eng_cspt += add.get_eng_cspt();
		time += add.get_time();
		ratio += add.get_ratio();
		ratio_base += add.get_ratiobase();
		cpu_utilization += add.get_cpu_uti();
		bw_utilization += add.get_bw_uti();
		maxdevices_bybase += add.get_maxdevices_bybase();
		avedevices_bybase += add.get_avedevices_bybase();
		averadius += add.get_averadius();
		maxadius += add.get_maxadius();
		numselecteddisk += add.get_numselecteddisk();
	}

	void cal_ratio(Result& opt) {
		std::cout << eng_cspt << ' ' << opt.get_eng_cspt() << "\n========================\n";
		ratio = eng_cspt / opt.get_eng_cspt();
	}
};


#endif // PREDEFINE_H
