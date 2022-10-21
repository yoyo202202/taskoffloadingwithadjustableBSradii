#include<iostream>
#include"predefine.h"
#include<string>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace std;

#define CASE_NUM 5


//三组实验
string case1[CASE_NUM] = { "d20","d50","d100","d200","d500" };    //随devices数变化
string case2[CASE_NUM] = { "b4","b8","b15","b25","b40" };   //随basestation数变化
string case4[CASE_NUM] = { "r200","r500","r800","r1500","r2500" };    //随region变化

//string root = "C:\\Users\\SQ\\Desktop\\Alidataset0606\\data\\numberofdevices\\";
string root = "C:\\Users\\SQ\\Desktop\\Alidataset0606\\data\\numberofbasestations\\";
//string root = "P:\\Alidataset\\data\\region\\";

Result OPT[CASE_NUM];
Result Greedy[CASE_NUM];
Result PrimalDual[CASE_NUM];

void an_example();
void case_n();
void print_to_csv(Result* OPT,Result* Greedy,Result* PrimalDual);

int main()
{
	case_n();

}


int inst_num =5;    //每组实验的实例数
string case_name = "case2";
int order =3;//每组实验的第几种情况


void case_n() {
	for (int i = 0; i < CASE_NUM; i++) {   
		OPT[i] = Result();
		OPT[i].set_eng_cspt(0);
		Greedy[i] = Result();
		Greedy[i].set_eng_cspt(0);
		PrimalDual[i] = Result();
		PrimalDual[i].set_eng_cspt(0);
	}

	//for (int i = 0; i < CASE_NUM; i++) {   
		string path = root + case2[order];
		int nu = inst_num;

		// inst_num 指的是每组实验的实例数
		for (int j = 0; j < inst_num; j++) { 
			// 每次循环计算一个实例
			Result opt;
			Result gre;
			Result pd;
			string file_name = path + "\\instance" + to_string(j) + ".txt";
			cout << file_name << '\n';
			Cover cover = Cover(file_name);
			//cover.print_test();

			//pd = cover.PrimalDual_Result();
			//std::cout << "PrimalDual_Result:\n";
			//pd.print_titled();
			//cover.print_result(pd);

			//gre = cover.Greedy_Result();
			//std::cout << "Greedy_Result:\n";
		   // gre.print_titled();
			//cover.print_result(gre);

			//if (i <= 3)
			//{
				opt = cover.OPT_Result();
				std::cout << "OPT_Result:\n";
				opt.print_titled();
				cover.print_result(opt);

				//if (opt.get_eng_cspt() != 0)
				//{
					//pd.cal_ratio(opt);
				//	gre.cal_ratio(opt);
				//}

			//}

			OPT[order].add_result(opt);
		     //Greedy[order].add_result(gre);
			//if (pd.get_eng_cspt() == 0)
				//nu -= 1;
			//PrimalDual[order].add_result(pd);
		}

		OPT[order].average_result(inst_num);
		//Greedy[order].average_result(inst_num);
		//PrimalDual[order].average_result(nu);
	//}
	//cout << "print_to_csv"<<'\n';
	print_to_csv(OPT,Greedy,PrimalDual);
}


// 结果写入csv文件
void print_to_csv(Result* OPT, Result* Greedy, Result* PrimalDual) {
	ofstream output;
	string ename = "energy.csv";
	string tname = "time.csv";
	string rname = "ratiobase.csv";
	string cpuname = "cpuutilization.csv";
	string bwname = "bwutilization.csv";
	string mdname = "maxdevicesbybase.csv";
	string adname = "avedevicesbybase.csv";
	string arname = "averadius.csv";
	string maname = "maxadius.csv";
	string ndname = "numselecteddisk.csv";
	string path = "C:\\Users\\SQ\\Desktop\\CoverProblemsu20220606\\experiment\\";
	path += case_name + "\\";

	//for (int i = 0; i < 2; i++) {
		output.open(path + ename, ios::app);
		output << setprecision(5) << OPT[order].get_eng_cspt() << ',' << Greedy[order].get_eng_cspt() << ',' << PrimalDual[order].get_eng_cspt() << '\n';
		output.close();

		output.open(path + tname, ios::app);
		output << setprecision(5) << OPT[order].get_time() << ',' << Greedy[order].get_time() << ',' << PrimalDual[order].get_time() << '\n';
		output.close();

		output.open(path + rname, ios::app);
		output << setprecision(5) << OPT[order].get_ratiobase() << ',' << Greedy[order].get_ratiobase() << ',' << PrimalDual[order].get_ratiobase() << '\n';
		output.close();

		output.open(path + cpuname, ios::app);
		output << setprecision(5) << OPT[order].get_cpu_uti() << ',' << Greedy[order].get_cpu_uti() << ',' << PrimalDual[order].get_cpu_uti() << '\n';
		output.close();

		output.open(path + bwname, ios::app);
		output << setprecision(5) << OPT[order].get_bw_uti() << ',' << Greedy[order].get_bw_uti() << ',' << PrimalDual[order].get_bw_uti() << '\n';
		output.close();

		output.open(path + mdname, ios::app);
		output << setprecision(5) << OPT[order].get_maxdevices_bybase() << ',' << Greedy[order].get_maxdevices_bybase() << ',' << PrimalDual[order].get_maxdevices_bybase() << '\n';
		output.close();

		output.open(path + adname, ios::app);
		output << setprecision(5) << OPT[order].get_avedevices_bybase() << ',' << Greedy[order].get_avedevices_bybase() << ',' << PrimalDual[order].get_avedevices_bybase() << '\n';
		output.close();

		output.open(path + arname, ios::app);
		output << setprecision(5) << OPT[order].get_averadius() << ',' << Greedy[order].get_averadius() << ',' << PrimalDual[order].get_averadius() << '\n';
		output.close();

		output.open(path + maname, ios::app);
		output << setprecision(5) << OPT[order].get_maxadius() << ',' << Greedy[order].get_maxadius() << ',' << PrimalDual[order].get_maxadius() << '\n';
		output.close();

		output.open(path + ndname, ios::app);
		output << setprecision(5) << OPT[order].get_numselecteddisk() << ',' << Greedy[order].get_numselecteddisk() << ',' << PrimalDual[order].get_numselecteddisk() << '\n';
		output.close();
	//}
}


