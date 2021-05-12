#include <ilcplex/ilocplex.h>
#include <fstream>
#include "stdio.h"
#include "stdlib.h"
#include <cstring>
#include <iostream>
#include <sstream>

ILOSTLBEGIN

//typedef IloArray<IloIntArray>    IntMatrix;          //整数型数组
typedef IloArray<IloNumArray>    NumMatrix;          //数值数组
typedef IloArray<NumMatrix>      NumArray3;
typedef IloArray<IloNumVarArray> NumVarMatrix;      //变量矩阵
typedef IloArray<NumVarMatrix>   NumVarArray3;
typedef IloArray<NumVarArray3>  NumVarArray4;

//const IloInt T = 550;              //时间线
//const IloInt Horizon = T + 2;
ofstream out1put("M_0.txt");
ofstream out2put("Result.txt");
ofstream out3put("b.txt");
ofstream out4put("c.txt");
ofstream out5put("m.txt");

int main(int argc, char **argv)            //main函数包含两个参数，一个整形，一个指针的指针
{
	IloEnv env;                          //设定环境
										 //define_data(env);
										 //IloEnv env2 = env;
	try
	{                                //用try检验
	//	const char* filename = "H:/郜振鑫/example/example3_1.dat";
		const char* filename = "H:/郜振鑫/example/Data needed for computing Minimum LIM.dat";
		//const char* filename = "G:/job_shop/job_shop.dat";             //先是定义了一个字符型指针，指向文件
		if (argc > 1)                     //如果main函数中的第一个参数值大于1
		filename = argv[1];            //则将第二个参数赋值给filename，main的第二个参数存储了file的地址，也即第二参数指向了文件
		ifstream file(filename);    //打开文件的方式是在类ios中定义，其中ios::in表示文件以输入方式打开
		if (!file)
		{
			cerr << "ERROR: could not open file '" << filename
				<< "' for reading" << endl;
			cerr << "usage:   " << argv[0] << " <file>" << endl;
			throw(-1);
		}
		NumMatrix     Lamda(env);
		NumMatrix     Large(env);
		IloNumArray   Qi(env);
		NumArray3     Pro(env); 
		int R_num;
		int number;
		int Circle;
		int big_num ;
		int Chang ;
		file >> R_num;
		file >> number;
		file >> Circle;
		file >> big_num;
		file >> Chang;
		file >> Lamda;
		file >> Large;
		file >> Qi;
		Pro = NumArray3(env, Circle);
		for (IloInt k = 0; k < Circle; k++)
		{
			Pro[k] = NumMatrix(env, R_num);
			file >> Pro[k];
		}
	//	cout << "firstJ: " << Lamda << endl;
	//	cout << "firstJ: " << Large << endl;
		cout << " Qi: " << Qi << endl;
	//	cout << "Processtime1: " << Processtime1 << endl;
		cout << 1;
	
		//建立模型!!!!!!!!!!!!!!
		IloModel model(env);
		//创建而维数组变量Y
		NumVarMatrix a(env, Circle);
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			a[ii] = IloNumVarArray(env, R_num);
			for (IloInt t = 0; t < R_num; t++)
			{
				std::stringstream f_name;
				f_name << "a[" << ii << t << "]";
			//a[ii][t].setName(f_name.str().c_str());
			    if (Lamda[ii][t] > 0)
			 	{
					a[ii][t] = IloNumVar(env, 0, 1, ILOINT, f_name.str().c_str());
			 	}
				f_name.flush();
				//Y[i].add(Y[i][m]);
			}
		}
	//	NumVarMatrix b(env, Circle);
		NumVarArray3 x(env, Circle);
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			x[ii] = NumVarMatrix(env, R_num);
		//	x[ii] = IloNumVarArray(env, R_num);
			for (IloInt t = 0; t < R_num; t++)
			{
				x[ii][t] = IloNumVarArray(env, number);
				for (IloInt tt = 0; tt < number; tt++)
				{
					std::stringstream f_name;
					f_name << "x[" << ii << t << tt << "]";
					if (Pro[ii][t][tt] > 0)
					{
						x[ii][t][tt] = IloNumVar(env, 0, 1, ILOINT, f_name.str().c_str());
					}

					//	if (ii == 1 && t == 2)
					//	{
					//		b[ii][t] = IloNumVar(env, 0, 0, ILOINT, f_name.str().c_str());
					//	}
						//Y[i].add(Y[i][m]);
					f_name.flush();
				}
			}
		}
	
		NumVarMatrix c(env, Circle);
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			c[ii] = IloNumVarArray(env, R_num);
			for (IloInt t = 0; t < R_num; t++)
			{
				std::stringstream f_name;
				f_name << "c[" << ii << t << "]";
			 	if (Lamda[ii][t] > 0)
			 	{
					c[ii][t] = IloNumVar(env, 0, 1, ILOINT, f_name.str().c_str());
			 	}
				f_name.flush();
			}
		}
		NumVarMatrix m(env, Circle);
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			m[ii] = IloNumVarArray(env, R_num);
			for (IloInt t = 0; t < R_num; t++)
			{
				std::stringstream f_name;
				f_name << "m[" << ii << t << "]";
			 	if (Lamda[ii][t] > 0)
			 	{
					m[ii][t] = IloNumVar(env, 0, 1, ILOINT, f_name.str().c_str());
					f_name.flush();
			    }
			}
		}
		IloNumVarArray   M_0(env, R_num);
		IloInt R_max;
		for (IloInt t = 0; t < R_num; t++)
		{
			R_max = 0;
			for (IloInt r = 0; r < Circle; r++)
			{
			//	cout << Large[r][t]<< "-->";
				if (R_max < Large[r][t])
				{
					R_max = Large[r][t];
				}
			}
		//	cout << "最大值：" << R_max;
		//	cout << endl;
			std::stringstream f_name;
			f_name << "M_0[" << t << "]";
			M_0[t] = IloNumVar(env, Qi[t], R_max, ILOINT, f_name.str().c_str());
			f_name.flush();
		//	cout <<"下限"<< M_0[t].getLB()<< "  上限" << M_0[t].getUB();
		//	cout << endl;
		}
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
					model.add(M_0[t] + big_num*(a[ii][t] - 1)<= Lamda[ii][t] - 1);
				}
			}
		}
	//	IloInt flag;
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
				//	flag = 0;
						//	for (IloInt kk = 0; kk < number; kk++)
					for (IloInt kk = 0; kk < number; kk++)
					{
						if (Pro[ii][t][kk] > 0)
						{
						//	flag = 1;
							model.add(M_0[t] + big_num*(x[ii][t][kk] - 1) <= Lamda[ii][t] + Pro[ii][t][kk]);
						}
						
					}
				
				}
			}
		}
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
				//	if (ii == 0)
				//	{
						for (IloInt kk = 0; kk < number; kk++)
						{
							if (Pro[ii][t][kk] > 0)
							{
								model.add(M_0[t] - big_num*(x[ii][t][kk] - 1) >= Lamda[ii][t] + Pro[ii][t][kk]);
							}
							
						}
			
				}
			}
		}
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
					model.add(M_0[t] - big_num*(c[ii][t] - 1) >= Large[ii][t]);
				}
			}
		}
		IloInt flag;
		IloExpr costSum1(env);
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
					flag = 0;
					costSum1.clear();
					for (IloInt kk = 0; kk < number; kk++)
					{
						if (Pro[ii][t][kk] > 0)
						{
							flag = 1;
							costSum1 = costSum1 + x[ii][t][kk];
						}
					}
					if (flag == 1)
					{
						model.add((a[ii][t] + costSum1 + c[ii][t]) <= m[ii][t] * Chang);
					}
					else
					{
						 model.add((a[ii][t] + c[ii][t]) <= m[ii][t] * Chang);;
					}
				}
				costSum1.clear();
			}
		}
		IloInt flag1;
		IloExpr costSum2(env);
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
					flag1 = 0;
					costSum2.clear();
					for (IloInt kk = 0; kk < number; kk++)
					{
						if (Pro[ii][t][kk] > 0)
						{
							flag1 = 1;
							costSum2 = costSum2 + x[ii][t][kk];
						}
					}
					if (flag1 == 1)
					{
						model.add(m[ii][t] - big_num*(a[ii][t] + costSum2 + c[ii][t]) <= 0);
					}
					else
					{
						model.add(m[ii][t] - big_num*(a[ii][t] + c[ii][t]) <= 0);;
					}
				}
				costSum2.clear();
			}
		}
		for (IloInt ii = 0; ii< Circle; ii++)
		{
			IloExpr cost1(env);
			for (IloInt t = 0; t < R_num; t++)
			{
				if (Lamda[ii][t] > 0)
				{
					cost1 += m[ii][t];
				}
			}
			model.add(cost1 >= 1);
			cost1.clear();
		}
		// objective function 目标函数
		IloExpr costSum(env);
		for (IloInt t = 0; t < R_num; t++)
		{
			costSum = costSum + M_0[t];
		}
		IloObjective obj = IloMinimize(env, costSum);
	//	IloObjective obj = IloMaximize(env, costSum);
		model.add(obj);
		costSum.end();
		//模型求解
		IloCplex cplex(env);
		cplex.extract(model);
		cplex.exportModel("job shop.lp");
		//IloCplex cplex(model);

		//IloNum tolerance = cplex.getParam(IloCplex::EpInt);
		//cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
		if (!cplex.solve())
		{
			env.error() << "Failed to optimize this problem." << endl;
			throw(-1);
		}
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution Value = " << cplex.getObjValue() << endl;
		env.out() << "number of finished units of order i :" << endl;
		for (IloInt i = 0; i < R_num; i++)
		{
			env.out() << cplex.getValue(M_0[i]) << endl;
		}
		
	//	stringstream ss;
	//	ss << T;
		//string str = ss.str() + "output.txt";
	//	string str = "FMS0" + fms + "-ORDER" + order + ".txt";
	//	ofstream fout(str);
		for (IloInt i = 0; i < R_num; i++)
		{
			out2put << "M0[" << i << "]="<< cplex.getValue(M_0[i]) << " ";
		}
		out2put << "\n";
		out2put << "\n";
		for (IloInt i = 0; i != Circle; ++i) {
			for (IloInt mm = 0; mm != R_num; ++mm) {
					if (Lamda[i][mm] > 0)
					{
						//cout << "U[" << i << "][" << machine << "][" << t << "]=" << cplex.getValue(U[i][machine][t]) << "\t"<<endl;
						out2put << "a[" << i << "][" << mm << "]=" << cplex.getValue(a[i][mm]) << "\t";
				}
			//	out2put << "\n";
			}
		}
		out2put << "\n";
		out2put << "\n";
		for (IloInt i = 0; i != Circle; ++i) {
			for (IloInt mm = 0; mm != R_num; ++mm) {
				if (Lamda[i][mm] > 0)
				{
					for (IloInt kk = 0; kk != number; ++kk) {
						if(Pro[i][mm][kk] > 0)
						//cout << "U[" << i << "][" << machine << "][" << t << "]=" << cplex.getValue(U[i][machine][t]) << "\t"<<endl;
						out2put << "x[" << i << "][" << mm << "][" << kk <<"]="<< cplex.getValue(x[i][mm][kk]) << "\t";
					}
				}
			//	out2put << "\n";
			}
		}
		out2put << "\n";
		out2put << "\n";
		for (IloInt i = 0; i != Circle; ++i) {
			for (IloInt mm = 0; mm != R_num; ++mm) {
				if (Lamda[i][mm] > 0)
				{
					//cout << "U[" << i << "][" << machine << "][" << t << "]=" << cplex.getValue(U[i][machine][t]) << "\t"<<endl;
					out2put << "c[" << i << "][" << mm << "]=" << cplex.getValue(c[i][mm]) << "\t";
				}
			//	out2put << "\n";
			}
		}
		out2put << "\n";
		out2put << "\n";
		for (IloInt i = 0; i != Circle; ++i) {
			for (IloInt mm = 0; mm != R_num; ++mm) {
				if (Lamda[i][mm] > 0)
				{
					//cout << "U[" << i << "][" << machine << "][" << t << "]=" << cplex.getValue(U[i][machine][t]) << "\t"<<endl;
					out2put << "m[" << i << "][" << mm << "]=" << cplex.getValue(m[i][mm]) << "\t";
				}
			//	out2put << "\n";
			}
		}


	
	}
	catch (IloException& e)
	{
		cerr << "Concert Exception: " << e << endl;
	}
	catch (...)
	{
		cerr << "Other Exception" << endl;
	}
	env.end();
	return 0;
}