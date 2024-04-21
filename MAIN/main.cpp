#include <iostream>
#include <fstream>
#include <functional>
#include <io.h>
#include <random>
#include <cmath>
#include <omp.h>
#include <vector>
#include <BGAL/Optimization/LinearSystem/LinearSystem.h>
#include <BGAL/Optimization/ALGLIB/optimization.h>
#include <BGAL/Optimization/LBFGS/LBFGS.h>
#include <BGAL/BaseShape/Point.h>
#include <BGAL/BaseShape/Polygon.h>
#include <BGAL/Integral/Integral.h>
#include <BGAL/Model/ManifoldModel.h>
#include <BGAL/Model/Model_Iterator.h>
#include <BGAL/Optimization/GradientDescent/GradientDescent.h>
#include <BGAL/BaseShape/KDTree.h>
#include "BGAL/Optimization/ALGLIB/dataanalysis.h"
#include "MyRPD.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include <Eigen/Sparse>

using namespace alglib;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point1;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point1, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron1;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;



pair<double, double>  V3toV2(Eigen::Vector3d nor)
{
	double u = 0.0, v = 0.0;
	u = acos(nor.z());
	if (u == 0)
	{
		v = 0;
	}
	else
	{
		double tmp1 = abs(acos(nor.x() / sin(acos(nor.z()))));
		double tmp2 = abs(asin(nor.y() / sin(acos(nor.z()))));
		if (isnan(tmp1))
		{
			if (isnan(tmp2))
			{
				v = 0.0;
				pair<double, double> p(u, v);
				return p;
			}
			else
			{
				v = tmp2;
				pair<double, double> p(u, v);
				return p;
			}

		}
		if (isnan(tmp2))
		{
			if (isnan(tmp1))
			{
				v = 0.0;
				pair<double, double> p(u, v);
				return p;
			}
			else
			{
				v = tmp1;
				pair<double, double> p(u, v);
				return p;
			}
		}
		Eigen::Vector3d n11(sin(u) * cos(tmp1), sin(u) * sin(tmp1), cos(u));
		Eigen::Vector3d n22(sin(u) * cos(tmp2), sin(u) * sin(tmp2), cos(u));
		double tot1, tot2;
		tot1 = (n11 - nor).x() * (n11 - nor).x() + (n11 - nor).y() * (n11 - nor).y() + (n11 - nor).z() * (n11 - nor).z();
		tot2 = (n22 - nor).x() * (n22 - nor).x() + (n22 - nor).y() * (n22 - nor).y() + (n22 - nor).z() * (n22 - nor).z();
		if (tot1 < tot2)
		{
			v = tmp1;
		}
		else
		{
			v = tmp2;
		}
		if (abs(sin(u) * cos(v) - nor.x()) > 0.1)
		{
			u = -1.0 * u;
		}
		if (abs(sin(u) * sin(v) - nor.y()) > 0.1)
		{
			v = -1.0 * v;
		}

	}
	pair<double, double> p(u, v);
	return p;
}




double PI = 3.1415926535;
void WindingNumLBFGSTest(string modelpath, string model, vector<vector<Eigen::Vector3d>> VDPs, pair<vector<double>, map<MyPoint, set<MyPoint>>> areasNeibor,bool Ifdoublelyer)
{


	
	clock_t start, end;
	ifstream inNewPs(modelpath + model + ".xyz");
	vector<Eigen::Vector3d> Vall;
	vector<Eigen::Vector3d> Nall;
	map<MyPoint, pair<int, int>> VDPsMap;
	map<pair<int, int>,int> ID2ID;
	map<pair<int, int>,int> ID2insideID;
	map<pair<int, int>,int> ID2outsideID;
	map<pair<int, int>, pair<int, int>> VDPsIDMap;
	map<pair<int, int>, pair<int, int>> VDPs_DIDMap;
	int n = 0;
	double x11, y11, z11;
	while (inNewPs >> x11 >> y11 >> z11)
	{
		n++;
		Eigen::Vector3d NewP(x11, y11, z11);
		Vall.push_back(NewP);

		
		inNewPs >> x11 >> y11 >> z11;
		Eigen::Vector3d normal(x11, y11, z11);
		normal.normalize();
		Nall.push_back(normal);
	}

	/*vector<Eigen::Vector3d> VDPs_P = VDPs.first;
	vector<double> VDPs_D = VDPs.second;*/
	vector<double> Areas = areasNeibor.first;
	auto Neibors = areasNeibor.second;
	vector<vector<Eigen::Vector3d>> VDPs_QC;
	vector<vector<double>> WdNum;
	vector<vector<double>> DisT;
	vector<pair<int, int>> InsidePs, OutsidePs;
	map<pair<int, int>, bool> Isinside;
	double Totn = 0, Qcn = 0, maxDis=-99999,Insn=0,Outsn=0;
	for (int i = 0; i < n; i++)
	{
		//Areas[i] = 9.445902 / 4000.0;
		vector<Eigen::Vector3d> VDPs_QC_tmp;
		vector<double> WdNum_tmp;vector<double> Dist_tmp;
		maxDis = -99999;
		for (int qpid = 0; qpid < VDPs[i].size(); qpid++)
		{
			
			Totn++;
			Dist_tmp.push_back((VDPs[i][qpid] - Vall[i]).norm());
			maxDis = max(maxDis, Dist_tmp[Dist_tmp.size()-1]);
			MyPoint mp(VDPs[i][qpid]);
			if (VDPsMap.find(mp) == VDPsMap.end())
			{
				ID2ID[make_pair(i, VDPs_QC_tmp.size())] = Qcn;
				VDPsMap[mp] = make_pair(i, VDPs_QC_tmp.size());
				VDPsIDMap[make_pair(i, qpid)] = make_pair(i, VDPs_QC_tmp.size());
				//VDPs_DIDMap[make_pair(i, VDPs_QC_tmp.size())] = make_pair(i, qpid);
				
				auto queryP = VDPs[i][qpid];
				double wwd = 0;
				for (int j = 0; j < n; j++)
				{
					auto YY = (Vall[j] - queryP);
					wwd += Areas[j] * ((YY.dot(Nall[j])) / (4 * PI * pow((YY.squaredNorm()), 1.5)));
				}
				WdNum_tmp.push_back(wwd);
				if (wwd > 0.5)
				{
					InsidePs.push_back(make_pair(i, VDPs_QC_tmp.size()));
					Isinside[make_pair(i, VDPs_QC_tmp.size())] = 1;
					ID2insideID[make_pair(i, VDPs_QC_tmp.size())] = Insn;
					Insn++;
				}
				else
				{
					OutsidePs.push_back(make_pair(i, VDPs_QC_tmp.size()));
					ID2outsideID[make_pair(i, VDPs_QC_tmp.size())] = Outsn;
					Isinside[make_pair(i, VDPs_QC_tmp.size())] = 0;
					Outsn++;
				}
				VDPs_QC_tmp.push_back(VDPs[i][qpid]);
				Qcn++;
			}
			else
			{
				//VDPs_DIDMap[make_pair(i, qpid)] = make_pair(i, qpid);
				//WdNum_tmp.push_back(-99999);
				VDPsIDMap[make_pair(i, qpid)] = VDPsMap[mp];
			}
		}
		VDPs_QC.push_back(VDPs_QC_tmp);
		WdNum.push_back(WdNum_tmp);
		DisT.push_back(Dist_tmp);
	}
	cout << "Total Query points : " << Totn << " After QuChong n: " << Qcn << endl;





	int fgNum = 0;
	std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
		= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
	{
		fgNum++;
		double func = 0;
		double wd = 0;

		
		double lambda_x = 10;
		double lambda_nor = 1000;
		double lambda_half = 50.0;
		if (Ifdoublelyer)
		{
			lambda_x = 1;
			lambda_nor = 10;
			lambda_half = 50.0;
		}


		double lambda_smooth = 0.5;

		double sigma = 0.2;
		//vector<double> WdNum;
		double lossNor = 0, lossFC = 0, losswdn = 0, loss_laplaceIn = 0, loss_laplaceOut = 0;
		vector<Eigen::Vector3d> Nall_new;
		ofstream outEnd;
		if (Ifdoublelyer)
		{
			outEnd.open(modelpath + "Out\\" + model +to_string(fgNum) + "_Double_End.xyz");
		}
		else
		{
			outEnd.open(modelpath + "Out\\" + model  + to_string(fgNum) + "_End.xyz");
		}
		
		vector<double> ff;

		vector<int> II;
		for (int j = 0; j < n; j++)
		{
			II.push_back(j);
			g(j * 2) = 0;
			g(j * 2 + 1) = 0;
			double u = X(j * 2);
			double v = X(j * 2 + 1);

			Eigen::Vector3d N(sin(u) * cos(v), sin(u) * sin(v), cos(u));
			Nall_new.push_back(N);
			outEnd << Vall[j].transpose() << " " << (1.0*Nall_new[j]).transpose() << endl;

		}
		outEnd.close();
		
		double maxWD = -99999, minWD = 99999;
		int ComputeN = 0;
		int allN = 0;
		Eigen::VectorXd WdInsideVec, WdOutsideVec;
		WdInsideVec.resize(Insn);
		WdOutsideVec.resize(Outsn);








		
		
#pragma omp parallel for //schedule(dynamic, 5)	
		for (int i = 0; i < n; ++i)  // compute winding number
		{

			for (int qpid = 0; qpid < VDPs_QC[i].size(); qpid++)
			{
				//ComputeN++;
				auto queryP = VDPs_QC[i][qpid];
				double wwd = 0;
				for (int j = 0; j < n; j++)
				{
					double u = X(j * 2);
					double v = X(j * 2 + 1);
					auto YY = (Vall[j] - queryP);

					wwd += Areas[j] * ((YY.dot(Nall_new[j])) / (4 * PI * pow((YY.squaredNorm()), 1.5)));
				}
				//wwd = totarea * wwd / (n);
				WdNum[i][qpid] = wwd;

			}
		}
		

		
		if (1)
		{
			for (int i = 0; i < n; i++)
			{
				double TotWd = 0;

				for (int qpid = 0; qpid < VDPs_QC[i].size(); qpid++)
				{
					auto queryP = VDPs_QC[i][qpid];

					double rate = 3;
					double move = 1;
					double Oneside = 1;
					wd = WdNum[i][qpid];

					
					if (Ifdoublelyer)
					{

						func = func + lambda_x * (pow((wd - 0.5) / (sqrt(0.5)), 4) - pow((wd - 0.5) / (sqrt(0.5)), 2) - wd / 4);
						losswdn += lambda_x * (pow((wd - 0.5) / (sqrt(0.5)), 4) - pow((wd - 0.5) / (sqrt(0.5)), 2) - wd / 4);
					}
					else
					{
						func = func + lambda_x * (-wd);
						losswdn += lambda_x * (-wd);
					}

#pragma omp parallel for //schedule(dynamic, 5)	
					for (int j = 0; j < n; j++)
					{
						double u = X(j * 2);
						double v = X(j * 2 + 1);
						auto YY = (Vall[j] - queryP);
						// (rate*exp(rate*(move - wd)))/(exp(rate*(move - wd)) + 1)^2 - (rate*exp(-rate*(wd - move + 2)))/(exp(-rate*(wd - move + 2)) + 1)^2
						if (Ifdoublelyer)
						{
							g(j * 2) = g(j * 2) + lambda_x * (4 * pow((wd - 0.5) / (sqrt(0.5)), 3) - 2 * ((wd - 0.5) / (sqrt(0.5))) - 0.25)
								* (YY.x() * cos(v) * cos(u) + YY.y() * sin(v) * cos(u) - YY.z() * sin(u)) * (Areas[j])
								/ (4 * PI * pow((YY.squaredNorm()), 1.5));
							g(j * 2 + 1) = g(j * 2 + 1) + lambda_x * (4 * pow((wd - 0.5) / (sqrt(0.5)), 3) - 2 * ((wd - 0.5) / (sqrt(0.5)))- 0.25)
								* (-1.0 * YY.x() * sin(v) * sin(u) + YY.y() * cos(v) * sin(u)) * (Areas[j])
								/ (4 * PI * pow((YY.squaredNorm()), 1.5));
						}
						else
						{
							g(j * 2) = g(j * 2) + lambda_x * (-1)
								* (YY.x() * cos(v) * cos(u) + YY.y() * sin(v) * cos(u) - YY.z() * sin(u)) * (Areas[j])
								/ (4 * PI * pow((YY.squaredNorm()), 1.5));
							g(j * 2 + 1) = g(j * 2 + 1) + lambda_x * (-1)
								* (-1.0 * YY.x() * sin(v) * sin(u) + YY.y() * cos(v) * sin(u)) * (Areas[j])
								/ (4 * PI * pow((YY.squaredNorm()), 1.5));
						}
						
					}
				}
			}
		}
		

		if (1) //normal cons
		{

			for (int i = 0; i < n; i++)
			{
				double kk = VDPs[i].size();
				double TotWd = 0;
				vector<double> Gavg(2 * n, 0);
				for (int qpid = 0; qpid < VDPs[i].size(); qpid++)
				{
					//allN++;
					pair<int, int> Nid = VDPsIDMap[make_pair(i, qpid)];
					TotWd += WdNum[Nid.first][Nid.second];
				
					auto queryP = VDPs[i][qpid];

#pragma omp parallel for //schedule(dynamic, 5)	
					for (int j = 0; j < n; j++)
					{
						double u = X(j * 2);
						double v = X(j * 2 + 1);
						auto YY = (Vall[j] - queryP);
						Gavg[j * 2] = Gavg[j * 2] + (1 / kk)
							* (YY.x() * cos(v) * cos(u) + YY.y() * sin(v) * cos(u) - YY.z() * sin(u)) * (Areas[j])
							/ (4 * PI * pow((YY.squaredNorm()), 1.5));
						Gavg[j * 2 + 1] = Gavg[j * 2 + 1] + (1 / kk)
							* (-1.0 * YY.x() * sin(v) * sin(u) + YY.y() * cos(v) * sin(u)) * (Areas[j])
							/ (4 * PI * pow((YY.squaredNorm()), 1.5));
					}

				}
				double avg = TotWd / kk;

				for (int qpid = 0; qpid < VDPs[i].size(); qpid++)
				{
					auto queryP = VDPs[i][qpid];
					auto Nj = queryP - Vall[i];
					//auto Nj = Nall[i];
					pair<int, int> Nid = VDPsIDMap[make_pair(i, qpid)];
					double wdd = WdNum[Nid.first][Nid.second];
					double u = X(i * 2);
					double v = X(i * 2 + 1);
					Eigen::Vector3d Ni(sin(u) * cos(v), sin(u) * sin(v), cos(u));
					func = func + lambda_nor * (Ni.dot(Nj)) * wdd / kk;
					lossNor += lambda_nor * (Ni.dot(Nj)) * wdd / kk;
					func += -lambda_half / kk * (wdd - avg) * (wdd - avg);
					lossFC += -lambda_half / kk * (wdd - avg) * (wdd - avg);
					g(i * 2) = g(i * 2) + lambda_nor * ((Nj.x() * cos(u) * cos(v) - Nj.z() * sin(u) + Nj.y() * cos(u) * sin(v)) * wdd)/kk;
					g(i * 2 + 1) = g(i * 2 + 1) + lambda_nor * (Nj.y() * cos(v) * sin(u) - Nj.x() * sin(u) * sin(v)) * wdd /kk;
#pragma omp parallel for //schedule(dynamic, 5)	
					for (int j = 0; j < n; j++)
					{
						double u = X(j * 2);
						double v = X(j * 2 + 1);
						auto YY = (Vall[j] - queryP);
						g(j * 2) = g(j * 2) + lambda_nor * ((Ni.dot(Nj)) * ((YY.x() * cos(v) * cos(u) + YY.y() * sin(v) * cos(u) - YY.z() * sin(u)) * (Areas[j]) 
							/ (4 * PI  * pow((YY.squaredNorm()), 1.5))))/kk;
						g(j * 2 + 1) = g(j * 2 + 1) + lambda_nor * ((Ni.dot(Nj)) * ((-1.0 * YY.x() * sin(v) * sin(u) + YY.y() * cos(v) * sin(u)) * (Areas[j]) 
							/ (4 * PI  * pow((YY.squaredNorm()), 1.5))))/kk;
						g(j * 2) = g(j * 2) - lambda_half * (2 * (wdd - avg))
							* (((YY.x() * cos(v) * cos(u) + YY.y() * sin(v) * cos(u) - YY.z() * sin(u)) * (Areas[j])
								/ (4 * PI * pow(((Vall[j] - queryP).squaredNorm()), 1.5)))
								- Gavg[j * 2]) / kk;
						g(j * 2 + 1) = g(j * 2 + 1) - lambda_half * (2 * (wdd - avg))
							* (((-1.0 * YY.x() * sin(v) * sin(u) + YY.y() * cos(v) * sin(u)) * (Areas[j])
								/ (4 * PI * pow(((Vall[j] - queryP).squaredNorm()), 1.5)))
								- Gavg[j * 2 + 1]) / kk;
					}

				}
			}
		}

		maxWD = 1.1;
		minWD = -0.1;
		
		ofstream outQp;
		ofstream outwd;

		//if (Ifdoublelyer)
		//{
		//	outQp.open(modelpath + "Out\\" + model + to_string(fgNum) + "_Double_QueryPoints.txt");
		//	outwd.open(modelpath + "Out\\" + model + to_string(fgNum) + "_Double_WindingNums.txt");
		//}
		//else
		//{
		//	outQp.open(modelpath + "Out\\" + model + "_QueryPoints.txt");
		//	outwd.open(modelpath + "Out\\" + model + "_WindingNums.txt");
		//}
		//for (int i = 0; i < n; i++)
		//{
		//	for (int qpid = 0; qpid < VDPs_QC[i].size(); qpid++)
		//	{
		//		auto queryP = VDPs_QC[i][qpid];
		//		pair<int, int> Nid(i,qpid);
		//		if (WdNum[Nid.first][Nid.second] > maxWD )
		//		{
		//			outQp << queryP.transpose() << " " << 0 << " 0 1" << endl;
		//		}
		//		else
		//		{
		//			if (WdNum[Nid.first][Nid.second] < minWD)
		//			{
		//				outQp << queryP.transpose() << " " << 0 << " 1 0" << endl;
		//			}
		//			else
		//			{
		//				outQp << queryP.transpose() << " " << (WdNum[Nid.first][Nid.second] - minWD) / (maxWD - minWD) << " 0 0" << endl;
		//				
		//				
		//			}
		//		}
		//		outwd << WdNum[Nid.first][Nid.second] << "\n";
		//		//outwd << WdNum[Nid.first][Nid.second] <<"    "<<i<<" "<<qpid << "\n";
		//		continue;
		//	}
		//	
		//}
		//outQp.close();
		//outwd.close();
		
		//cout << "Func: " << func << " LossWindingNum: " << losswdn << " LossNormal: " << lossNor << " LossVariance: "<<lossFC<< " NormG: "<< g.norm()  << endl;
		return func;
	};
	BGAL::_LBFGS::_Parameter para;
	para.is_show = true;
	para.epsilon = 1e-8;
	para.max_iteration = 150;
	//para.max_time = 1000*1000;
		
	para.max_linearsearch = 20;
	
	BGAL::_LBFGS lbfgs(para);
	Eigen::VectorXd iterX(n * 2);
	Eigen::VectorXd RightX(n * 2);
	ofstream outStart(modelpath + "Out\\"+model+"_Start.xyz");
	for (int i = 0; i < n; ++i)
	{
		auto Q = V3toV2(Nall[i]);
		Eigen::Vector3d N_ball=Vall[i];
		N_ball.normalize();
		auto Q1 = V3toV2(N_ball);

		iterX(i * 2) = rand();
		iterX(i * 2 + 1) = rand(); // Here, we use random normal vectors to start the optimization.
		

		RightX(i * 2) = Q.first;
		RightX(i * 2 + 1) = Q.second;

		
		Eigen::Vector3d N(sin(iterX(i * 2)) * cos(iterX(i * 2 + 1)), sin(iterX(i * 2)) * sin(iterX(i * 2 + 1)), cos(iterX(i * 2)));
		outStart << Vall[i].transpose() << " " << N.transpose() << endl;
	}
	outStart.close();
	Eigen::VectorXd TempG(n * 2);
	// open this to check the gt normal loss, note that if using this, the first iter output will be GT
	// cout << "Init loss: " << fg(RightX, TempG) << endl;

	//return;
	cout << "Start opt...\n";
	lbfgs.minimize(fg, iterX);
	vector<Eigen::Vector3d> Nall_new;
	ofstream outEnd(modelpath + "Out\\" + model + "_End.xyz");
	for (int j = 0; j < n; j++)
	{
		double u = iterX(j * 2);
		double v = iterX(j * 2 + 1);
		Eigen::Vector3d N(sin(u) * cos(v), sin(u) * sin(v), cos(u));
		Nall_new.push_back(N);
		outEnd << Vall[j].transpose() << " " << Nall_new[j].transpose() << endl;
	}
	outEnd.close();
	end = clock();
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;

	
}


int main()
{
	// IMPORTANT NOTE: This code is not optimized for speed, but for clarity. 
	// Please open Openmp and AVX2 in Visual Studio to speed up the code.
	// Please set the floating point model to fast in Visual Studio to speed up the code.
	// The default number of Openmp parallel threads is 28, set according to an AMD Ryzen 5950x CPU, 
	// please set different number of threads according to the cpu you use to get the best running effect.
	// 
	// In order to allow you to view the optimization process in more detail, 
	// we have not set the optimization stop condition, you can manually stop the optimization, 
	// and view all iteration results in the data\out folder
	
	// A noise-free point cloud generally requires about 50 iterations, and a noisy point cloud may require more.


	string modelpath = "..\\..\\data\\";
	string modelname;
	omp_set_num_threads(28);
	{
		modelname = "BS_1000_torus"; // 

		cout << modelname << endl;
		bool Ifdoublelyer = 1; //change this to 0 if you want to try f_{01} = \sum_j^M -w_j


		auto VDPs = Comput_RPD(modelpath, modelname, 1); //Note that this is not a parallel version of Voronoi diagram computation. We will release the parallel version later.
		cout << "Compute 3D Voronoi DONE>> \n";


		std::cout << "====================WindingNumLBFGSTest" << std::endl;
		WindingNumLBFGSTest(modelpath, modelname, VDPs.first, VDPs.second, Ifdoublelyer); 

		std::cout << "successful!" << std::endl;
	
		//break;
	}
	return 0;
}
