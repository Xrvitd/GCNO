#include <iostream>
#include <fstream>
#include <functional>
#include <io.h>
#include <random>
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
#include "nanoflann.hpp"
#include "nanoflann/examples/utils.h"
#include "BGAL/Optimization/ALGLIB/dataanalysis.h"
#include "MyRPD.hpp"
#include "MyRPD_rnn.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>

using namespace nanoflann;
using namespace alglib;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point1;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point1, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron1;

int rnnnum = 80;

void Poisson(string modelpath, string model)
{
	clock_t start, end;
	start = clock();


	string outputPath = modelpath;

	//ifstream in("E:\\Dropbox\\MyProjects\\SIG-2022-Feature-preserving-recon\\data\\famous_pts_normal\\timing\\" + model + ".xyz");

	std::vector<Pwn> points;
	if (!CGAL::IO::read_points(CGAL::data_file_path(modelpath + "Denoise_Final_"+model+".xyz"), std::back_inserter(points),
		CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
		.normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
	{
		std::cerr << "Error: cannot read input file!" << std::endl;
		return;
	}
	Polyhedron1 output_mesh;
	double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
		(points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));
	if (CGAL::poisson_surface_reconstruction_delaunay
	(points.begin(), points.end(),
		CGAL::First_of_pair_property_map<Pwn>(),
		CGAL::Second_of_pair_property_map<Pwn>(),
		output_mesh, average_spacing))
	{
		std::ofstream out(outputPath + "\\model_poisson_" + model + ".off");
		out << output_mesh;
	}
	ifstream in(outputPath + "\\model_poisson_" + model + ".off");
	vector<Eigen::Vector3d> pts;
	vector<Eigen::Vector3i> facs;
	
	string line;
	in >> line;
	int q, w, e;
	in >> q >> w >> e;
	for (int i = 0; i < q; i++)
	{
		double x, y, z;
		in >> x >> y >> z;
		pts.push_back(Eigen::Vector3d(x, y, z));
	}
	for (int i = 0; i < w; i++)
	{
		int a, x, y, z;
		in >>a>> x >> y >> z;
		facs.push_back(Eigen::Vector3i(x, y, z));
	}
	std::ofstream outobj(outputPath + "\\model_poisson_" + model + ".obj");
	
	for (int i = 0; i < q; i++)
	{
		outobj << "v " << pts[i].transpose() << endl;
	}
	for (int i = 0; i < w; i++)
	{
		outobj << "f " << (facs[i]+ Eigen::Vector3i(1,1,1)).transpose() << endl;
	}
	outobj.close();


	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Poisson Running Time: " << endtime << endl;
}

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




void RFEPSTest(string modelpath,string model,bool ifdenoise)
{
	clock_t start, end;
	vector<Eigen::Vector3d> Vall, Nall;
	vector<vector<int>> neighboor;

	bool debugOutput = 0, IfoutputFile =1;

	ofstream out;
	//alglib::setglobalthreading(alglib::parallel);

	//string model = "angleWithNor";
	//string model = "0.005_50000_00040123_8fc7d06caf264003a242597a_trimesh_000";
	//string modelnn = "DenoisePoints";

	
	string outputPath = modelpath;

	//ifstream in("E:\\Dropbox\\MyProjects\\SIG-2022-Feature-preserving-recon\\data\\"+ outputFile +"\\" +model+"\\"+model+".xyz");
	//ifstream in("E:\\Dropbox\\MyProjects\\SIG-2022-Feature-preserving-recon\\data\\" + outputFile + "\\" + model + "\\01_"+model+".xyz");
	ifstream in;
	if (ifdenoise)
	{
		in.open(outputPath + "\\Denoise_" + model + ".xyz");
	}
	else
	{
		in.open(outputPath + "\\" + model + ".xyz");
	}

	while (!in.eof())
	{
		Eigen::Vector3d p, n;
		in >> p.x() >> p.y() >> p.z() >> n.x() >> n.y() >> n.z();
		Vall.push_back(p);
		Nall.push_back(n);
	}
	cout << "Read PointCloud. xyz File. \n";
	int r = Vall.size();


	// 0-1
	if (0)
	{
		Eigen::Vector3d maxp(-99999, -99999, -99999);
	Eigen::Vector3d minp(99999, 99999, 99999);

	for (int i = 0; i < r; i++)
	{
		maxp.x() = max(maxp.x(), Vall[i].x());
		maxp.y() = max(maxp.y(), Vall[i].y());
		maxp.z() = max(maxp.z(), Vall[i].z());
		minp.x() = min(minp.x(), Vall[i].x());
		minp.y() = min(minp.y(), Vall[i].y());
		minp.z() = min(minp.z(), Vall[i].z());
		Nall[i].normalize();
	}
	double minn = min(minp.x(), min(minp.y(), minp.z()));
	double maxn = max(maxp.x(), max(maxp.y(), maxp.z()));
	double maxl = max(maxp.x() - minp.x(), max(maxp.y() - minp.y(), maxp.z() - minp.z()));
	for (int i = 0; i < r; i++)
	{
		Vall[i].x() = (Vall[i].x() - minp.x()) / (maxl);
		Vall[i].y() = (Vall[i].y() - minp.y()) / (maxl);
		Vall[i].z() = (Vall[i].z() - minp.z()) / (maxl);
	} //
	out.open(outputPath + "\\01_" + model + ".xyz");
	for (size_t i = 0; i < r; i++)
	{
		out << Vall[i].transpose() << " " << Nall[i].transpose() << endl;
	}
	out.close();
	return;
	}
	
	

	double radis = 0.002;
	double lambda = 0.05;

	// rnn

	PointCloud<double> cloud;
	int mink = 9999999, maxk = -9999999;

	cloud.pts.resize(r);
	for (int i = 0; i < r; i++)
	{
		cloud.pts[i].x = Vall[i].x();
		cloud.pts[i].y = Vall[i].y();
		cloud.pts[i].z = Vall[i].z();
	}
	typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<double, PointCloud<double> >,
		PointCloud<double>,
		3 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();

	do
	{
		radis += 0.001;
		neighboor.clear();
		mink = 99999; maxk = -99999;
		for (int i = 0; i < r; i++)
		{
			vector<int> tmp;
			neighboor.push_back(tmp);
			double query_pt[3] = { Vall[i].x(), Vall[i].y(), Vall[i].z() };
			const double search_radius = static_cast<double>((radis) * (radis));
			std::vector<std::pair<uint32_t, double> >   ret_matches;
			nanoflann::SearchParams params;
			const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
			for (size_t j = 0; j < nMatches; j++)
			{
				neighboor[i].push_back(ret_matches[j].first);
			}
			maxk = max(maxk, int(nMatches));
			mink = min(mink, int(nMatches));
		}
		mink = 99999; maxk = -99999;
		auto neighboor_denoise = neighboor;
		for (int i = 0; i < r; i++)
		{
			int k = neighboor[i].size();
			neighboor_denoise[i].clear();
			for (int j = 0; j < k; j++)
			{
				int pid = neighboor[i][j];
				double dis = (Nall[i] - Nall[pid]).norm();
				//if (dis < 1.5)
				{
					neighboor_denoise[i].push_back(pid);
				}
			}
			maxk = max(maxk, int(neighboor_denoise[i].size()));
			mink = min(mink, int(neighboor_denoise[i].size()));
		}
		neighboor = neighboor_denoise;
		
		cout << "maxk: " << maxk << "   mink:" << mink << endl;
	}  while (maxk < rnnnum );

	// add par later

	start = clock();

	map<int, double> R2, R3, S;
	for (int iter = 0; iter < r; iter++)
	{
		R2[iter] = 0.0; R3[iter] = 0; S[iter] = 0;
	}

	int omp_cnt = 0;
	omp_set_num_threads(24);
#pragma omp parallel for schedule(dynamic, 20)  //part 1
	for (int iter = 0; iter < r; iter++) // part 1
	{
		//cout << iter << endl;
		omp_cnt++;
		std::function<void(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)> fop_lambda
			= [&](const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr) -> void
		{
			int k = neighboor[iter].size();
			int r = Vall.size();
			vector<double > w;
			for (int i = 0; i < k; i++)
			{
				w.push_back(0);
			}
			double maxw = -99999999.0;
			for (int i = 0; i < k; i++)
			{
				double dis = (Vall[iter] - Vall[neighboor[iter][i]]).norm();
				if (dis < 0.0001)
				{
					w[i] = 0.0;
				}
				else
				{
					w[i] = 1.0 / (dis * dis);
					maxw = max(maxw, w[i]);
				}
			}
			for (int i = 0; i < k; i++)
			{
				w[i] = 1.0;
				//w[i] = w[i] / maxw;
			}

			double u1 = x[k], v1 = x[k + 1], u2 = x[k + 2], v2 = x[k + 3];
			Eigen::Vector3d n1(sin(u1) * cos(v1), sin(u1) * sin(v1), cos(u1));
			Eigen::Vector3d n2(sin(u2) * cos(v2), sin(u2) * sin(v2), cos(u2));
			n1.normalize();
			n2.normalize();
			func = 0;
			for (int i = 0; i < k; i++)
			{
				auto Nj = Nall[neighboor[iter][i]];
				func += x[i] * w[i] * (pow((Nj.x() - sin(u1) * cos(v1)), 2) + pow(Nj.y() - sin(u1) * sin(v1), 2) + pow(Nj.z() - cos(u1), 2))
					+ (1.0 - x[i]) * w[i] * (pow((Nj.x() - sin(u2) * cos(v2)), 2) + pow(Nj.y() - sin(u2) * sin(v2), 2) + pow(Nj.z() - cos(u2), 2));
				//func += x[i] * w[i] * (Nall[neighboor[iter][i]] - n1).squaredNorm()  + (1 - x[i]) * w[i] * (Nall[neighboor[iter][i]] - n2).squaredNorm();
			}

			for (int i = 0; i < k; i++)
			{
				//double a = abs(cos(x[k]))* abs(cos(x[k]));
				auto Nj = Nall[neighboor[iter][i]];
				grad[i] = w[i] * (pow((Nj.x() - sin(u1) * cos(v1)), 2) + pow(Nj.y() - sin(u1) * sin(v1), 2) + pow(Nj.z() - cos(u1), 2))
					- w[i] * (pow((Nj.x() - sin(u2) * cos(v2)), 2) + pow(Nj.y() - sin(u2) * sin(v2), 2) + pow(Nj.z() - cos(u2), 2));
			}

			grad[k] = 0; grad[k + 1] = 0; grad[k + 2] = 0; grad[k + 3] = 0;

			for (int i = 0; i < k; i++)
			{
				auto Nj = Nall[neighboor[iter][i]];
				grad[k] += 2 * x[i] * sin(u1) * w[i] * (Nj.z() - cos(u1))
					- (2 * x[i] * cos(v1) * cos(u1) * w[i] * (Nj.x() - sin(u1) * cos(v1))
						+ 2 * x[i] * cos(u1) * w[i] * sin(v1) * (Nj.y() - sin(u1) * sin(v1)));

				grad[k + 1] += 2 * x[i] * sin(u1) * sin(v1) * w[i] * (Nj.x() - sin(u1) * cos(v1))
					- 2 * x[i] * sin(u1) * cos(v1) * w[i] * (Nj.y() - sin(u1) * sin(v1));

				grad[k + 2] += 2 * w[i] * sin(u2) * (1 - x[i]) * (Nj.z() - cos(u2))
					- (2 * w[i] * cos(u2) * cos(v2) * (1 - x[i]) * (Nj.x() - sin(u2) * cos(v2))
						+ 2 * w[i] * cos(u2) * sin(v2) * (1 - x[i]) * (Nj.y() - sin(u2) * sin(v2)));

				grad[k + 3] += 2 * sin(u2) * sin(v2) * (1 - x[i]) * w[i] * (Nj.x() - sin(u2) * cos(v2))
					- 2 * sin(u2) * cos(v2) * (1 - x[i]) * w[i] * (Nj.y() - sin(u2) * sin(v2));
			}
		};


		if (omp_cnt % 10000 == 0)
		{
			cout << "Part 1 --- Point iter: " << omp_cnt << " \n";
		}


		int k = neighboor[iter].size();
		if (k < 5)
		{
			continue;
		}
		// kmeans

		clusterizerstate s;
		kmeansreport rep;
		//real_2d_array xy = "[[1,1],[1,2],[4,1],[2,3],[4,1.5]]";
		real_2d_array xy;
		xy.setlength(k, 3);
		for (int i = 0; i < k; i++)
		{
			xy[i][0] = Nall[neighboor[iter][i]].x();
			xy[i][1] = Nall[neighboor[iter][i]].y();
			xy[i][2] = Nall[neighboor[iter][i]].z();
			if (debugOutput)
				cout << xy[i][0] << " " << xy[i][1] << " " << xy[i][2] << endl;
		}

		alglib::clusterizercreate(s);
		alglib::clusterizersetpoints(s, xy, 2);
		alglib::clusterizersetkmeanslimits(s, 10, 0);
		alglib::clusterizerrunkmeans(s, 2, rep);// ?
		if (debugOutput)
			printf("%d\n", int(rep.terminationtype)); // EXPECTED: 1
		Eigen::Vector3d C1, C2;
		if (int(rep.terminationtype) == -3)
		{
			C1 = Nall[neighboor[iter][0]];
			C2 = Nall[neighboor[iter][0]];
			//其实可以直接结束
		}
		else
		{
			C1.x() = rep.c[0][0];
			C1.y() = rep.c[0][1];
			C1.z() = rep.c[0][2];
			C2.x() = rep.c[1][0];
			C2.y() = rep.c[1][1];
			C2.z() = rep.c[1][2];

		}
		if (debugOutput)
		{
			cout << rep.c[0][0] << " " << rep.c[0][1] << " " << rep.c[0][2] << endl;
			cout << rep.c[1][0] << " " << rep.c[1][1] << " " << rep.c[1][2] << endl;
		}


		real_1d_array x0;
		x0.setlength(k + 4);
		real_1d_array s0;
		s0.setlength(k + 4);
		for (int i = 0; i < k; i++)
		{
			x0[i] = 0.5;
			s0[i] = 1;
		}
		for (int i = k; i < k + 4; i++)
		{
			x0[i] = 0;
			s0[i] = 1;
		}

		C2.normalize();
		C1.normalize();

		auto Q = V3toV2(C1);
		x0[k] = Q.first;
		x0[k + 1] = Q.second;
		Q = V3toV2(C2);
		x0[k + 2] = Q.first;
		x0[k + 3] = Q.second;
		//kmeans[iter] = make_pair(make_pair(x0[k], x0[k + 1]), make_pair(x0[k + 2], x0[k + 3]));

		real_1d_array bndl;
		real_1d_array bndu;
		bndl.setlength(k + 4);
		bndu.setlength(k + 4);

		real_2d_array c;
		c.setlength(1, neighboor[iter].size() + 5);
		integer_1d_array ct = "[0]";
		for (int i = 0; i < k; i++)
		{
			bndl[i] = 0;
			bndu[i] = 1;
			c[0][i] = 1;
		}
		for (int i = k; i < k + 4; i++)
		{
			bndl[i] = -99999999999;
			bndu[i] = 99999999999;
			c[0][i] = 0;
		}
		c[0][k + 4] = double(k) / 2.0;
		minbleicstate state;
		double epsg = 0;
		double epsf = 0;
		double epsx = 0;
		ae_int_t maxits = 0.000001;
		alglib::minbleiccreate(x0, state);
		alglib::minbleicsetlc(state, c, ct);
		alglib::minbleicsetbc(state, bndl, bndu);
		alglib::minbleicsetscale(state, s0);
		alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);

		alglib::minbleicoptguardsmoothness(state);
		alglib::minbleicoptguardgradient(state, 0.0001);

		if (debugOutput)
		{
			for (int i = 0; i < k + 4; i++)
			{
				cout << x0[i] << " ";
			}cout << endl;
		}



		minbleicreport rep2;
		//minbleicoptimize(state, fop, nullptr, nullptr, alglib::parallel);
		if (1)
		{
			alglib::minbleicoptimize(state, fop_lambda);
			alglib::minbleicresults(state, x0, rep2);
			//cout << rep2.debugff << endl;
			double mn = 0;
			real_1d_array G_tmp;
			G_tmp.setlength(k + 4);
			fop_lambda(x0, mn, G_tmp, nullptr);
			if (debugOutput)
				cout << mn << endl;
			if (debugOutput)
			{
				printf("%d\n", int(rep2.terminationtype)); // EXPECTED: 4
				printf("%s\n", x0.tostring(2).c_str()); // EXPECTED: [2,4]

				optguardreport ogrep;
				minbleicoptguardresults(state, ogrep);
				printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
				printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
				printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
			}
		}
		

		double u1 = x0[k], v1 = x0[k + 1], u2 = x0[k + 2], v2 = x0[k + 3];
		Eigen::Vector3d n1(sin(u1) * cos(v1), sin(u1) * sin(v1), cos(u1));
		Eigen::Vector3d n2(sin(u2) * cos(v2), sin(u2) * sin(v2), cos(u2));
		/*n1.normalize();
		n2.normalize();*/
		if (debugOutput)
		{
			cout << n1 << endl;
			cout << n2 << endl;
		}
		//n1 = C1;
		//n2 = C2;

		double dis = sqrt((n1.x() - n2.x()) * (n1.x() - n2.x()) + (n1.y() - n2.y()) * (n1.y() - n2.y()) + (n1.z() - n2.z()) * (n1.z() - n2.z()));
		if (debugOutput)
			cout << "dis " << dis << endl;

		//compute angle between n1 n2
		double angle = acos(n1.dot(n2));
		// angel to drgee
		angle = angle * 180 / 3.1415926535;
		if (debugOutput)
			cout << "angle " << angle << endl;

		if (dis < 30 / 100.0) // importent
		{
			R3[iter] = 1; // normal point
		}
		else
		{
			R3[iter] = 0;
		}
		//if (mn < 0.25)
		{
			R2[iter] = 1;
		}
	}

	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "T2 Running Time: " << endtime << endl;
	cout << "T2 Running Time: " << endtime * 1000 << " ms " << endl;
	

	if (IfoutputFile)
	{
		out.open(outputPath  + "\\ShowColorR2_"+model+".txt");
		for (int i = 0; i < r; i++)
		{
			out << Vall[i].x() << " " << Vall[i].y() << " " << Vall[i].z() << " " << R2[i] << " 0.1 0.1 \n";
		}
		out.close();
		out.open(outputPath + "\\ShowColorR3_" + model + ".txt");
		for (int i = 0; i < r; i++)
		{
			out << Vall[i].x() << " " << Vall[i].y() << " " << Vall[i].z() << " " << R3[i] << " 0.1 0.1 \n";
		}
		out.close();
	}

	//return;
	// part2  normal
	//KNNSC
	auto neighboor_M = neighboor;
	cout << "ASDA\n";
	for (int i = 0; i < r; i++)
	{
		int k = neighboor[i].size();
		if (k < 5)
		{
			continue;
		}

		bool flag = 0;
		if (R3[i] == 1)
		{
			flag = 1;
		}
		/*for (auto p : neighboor_M[i])
		{
			if (R2[p] == 1)
			{
				flag = 1;
				break;
			}
		}*/
		if (flag)
		{
			
			continue;
		}
		int pid = 1;

		
		k = neighboor_M[i].size();
		neighboor_M[i].clear();

		double query_pt[3] = { Vall[i].x(), Vall[i].y(), Vall[i].z() };
		const double search_radius = static_cast<double>((radis * 3) * (radis * 3));
		std::vector<std::pair<uint32_t, double> > ret_matches;
		nanoflann::SearchParams params;
		const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

		neighboor_M[i].push_back(i);
		for (size_t j = 0; j < nMatches; j++)
		{
			if (R3[ret_matches[j].first] == 1)
			{
				neighboor_M[i].push_back(ret_matches[j].first);
				pid++;
				if (pid == k)
				{
					break;
				}
			}
		}
		double tmpr = radis * 3;
		while (pid != k)
		{
			tmpr *= 4;
			pid = 1;
			neighboor_M[i].clear();
			const double search_radius = static_cast<double>((tmpr) * (tmpr));
			std::vector<std::pair<uint32_t, double> > ret_matches;
			nanoflann::SearchParams params;
			const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

			neighboor_M[i].push_back(i);
			for (size_t j = 0; j < nMatches; j++)
			{
				if (R3[ret_matches[j].first] == 1)
				{
					neighboor_M[i].push_back(ret_matches[j].first);
					pid++;
					if (pid == k)
					{
						break;
					}
				}
			}
		}
	}
	cout << "ASDA\n";

	omp_cnt = 0;
	map<int, Eigen::Vector3d> Nall_new;
	for (int iter = 0; iter < r; iter++)
	{
		Nall_new[iter] = Nall[iter];
	}
	
	start = clock();
	debugOutput = 0;
#pragma omp parallel for schedule(dynamic, 20)  //part 2
	for (int iter = 0; iter < r; iter++) // part 2
	{
		omp_cnt++;

		if (omp_cnt % 10000 == 0)
			cout << "Part 2 --- Point iter: " << omp_cnt << " \n";

		int k = neighboor[iter].size();
		if (k < 5)
		{
			continue;
		}
		k = neighboor_M[iter].size();
		/*if (R2[iter] == 1 && R3[iter] == 1)
		{
			S[iter] = 1;
			continue;
		}*/

		if (1)
		{
			vector<int> v1, cnt;
			for (int i = 0; i < k; i++)
			{
				v1.push_back(i);
				cnt.push_back(0);
			}
			for (int i = 0; i < k - 1; i++)
			{
				for (int j = i + 1; j < k; j++)
				{
					double dij = (Nall[neighboor_M[iter][i]] - Nall[neighboor_M[iter][j]]).norm();
					if (dij < lambda)
					{
						v1[j] = v1[i];
					}
				}
			}
			for (int i = 0; i < k; i++)
			{
				cnt[v1[i]] = cnt[v1[i]] + 1;
			}
			int mx = -99999;
			int center = 0;
			for (int i = 0; i < k; i++)
			{
				if (cnt[i] > mx)
				{
					mx = cnt[i];
					center = i;
				}
			}
			Eigen::Vector3d nor(0, 0, 0);
			for (int i = 0; i < k; i++)
			{
				if (v1[i] == center)
				{
					nor = nor + Nall[neighboor_M[iter][i]];
				}
			}
			nor.x() /= double(mx);
			nor.y() /= double(mx);
			nor.z() /= double(mx);
			nor.normalize();
			S[iter] = 1;

			Nall_new[iter] = nor;
			continue;
		}
		// kmeans test
		if (0)
		{
			// kmeans

			clusterizerstate s;
			kmeansreport rep;
			//real_2d_array xy = "[[1,1],[1,2],[4,1],[2,3],[4,1.5]]";
			real_2d_array xy;
			xy.setlength(k, 3);
			for (int i = 0; i < k; i++)
			{
				xy[i][0] = Nall[neighboor_M[iter][i]].x();
				xy[i][1] = Nall[neighboor_M[iter][i]].y();
				xy[i][2] = Nall[neighboor_M[iter][i]].z();
				if (debugOutput)
					cout << xy[i][0] << " " << xy[i][1] << " " << xy[i][2] << endl;
			}

			alglib::clusterizercreate(s);
			alglib::clusterizersetpoints(s, xy, 2);
			alglib::clusterizersetkmeanslimits(s, 10, 0);
			alglib::clusterizerrunkmeans(s, 3, rep);// ?
			if (debugOutput)
				printf("%d\n", int(rep.terminationtype)); // EXPECTED: 1
			Eigen::Vector3d C1, C2 ,C3;
			if (int(rep.terminationtype) == -3)
			{
				C1 = Nall[neighboor_M[iter][0]];
				C2 = Nall[neighboor_M[iter][0]];
				C3 = Nall[neighboor_M[iter][0]];
				cout << "Error\n";
				//其实可以直接结束
			}
			else
			{
				C1.x() = rep.c[0][0];
				C1.y() = rep.c[0][1];
				C1.z() = rep.c[0][2];
				C2.x() = rep.c[1][0];
				C2.y() = rep.c[1][1];
				C2.z() = rep.c[1][2];
			
			}
			
			if (debugOutput)
			{
				cout << "C1: " << C1.x() << " " << C1.y() << " " << C1.z() << endl;
				cout << "C2: " << C2.x() << " " << C2.y() << " " << C2.z() << endl;
				cout << "C3: " << C3.x() << " " << C3.y() << " " << C3.z() << endl;
			}
			vector<int> cntk;
			cntk.resize(k);
			for (int i = 0; i < k; i++)
			{
				cntk[i] = 0;
			}
			for (int i = 0; i < k; i++)
			{
				cntk[rep.cidx[i]]++;
			}
			//find max cntk
			int max = -1;
			int center = 0;
			for (int i = 0; i < k; i++)
			{
				if (cntk[i] > max)
				{
					max = cntk[i];
					center = i;
				}
			}
			Nall_new[iter] = C1;
			continue;
		}
		// opt .  no use.


		clusterizerstate s;
		kmeansreport rep;
		//real_2d_array xy = "[[1,1],[1,2],[4,1],[2,3],[4,1.5]]";
		real_2d_array xy;
		xy.setlength(k, 3);
		for (int i = 0; i < k; i++)
		{
			xy[i][0] = Nall[neighboor_M[iter][i]].x();
			xy[i][1] = Nall[neighboor_M[iter][i]].y();
			xy[i][2] = Nall[neighboor_M[iter][i]].z();
			if (debugOutput)
				cout << xy[i][0] << " " << xy[i][1] << " " << xy[i][2] << endl;
		}

		alglib::clusterizercreate(s);
		alglib::clusterizersetpoints(s, xy, 2);
		alglib::clusterizersetkmeanslimits(s, 10, 0);
		alglib::clusterizerrunkmeans(s, 2, rep);// ?
		if (debugOutput)
			printf("%d\n", int(rep.terminationtype)); // EXPECTED: 1
		Eigen::Vector3d C1, C2;
		if (int(rep.terminationtype) == -3)
		{
			C1 = Nall[neighboor[iter][0]];
			C2 = Nall[neighboor[iter][0]];
			//其实可以直接结束
		}
		else
		{
			C1.x() = rep.c[0][0];
			C1.y() = rep.c[0][1];
			C1.z() = rep.c[0][2];
			C2.x() = rep.c[1][0];
			C2.y() = rep.c[1][1];
			C2.z() = rep.c[1][2];

		}
		if (debugOutput)
		{
			cout << rep.c[0][0] << " " << rep.c[0][1] << " " << rep.c[0][2] << endl;
			cout << rep.c[1][0] << " " << rep.c[1][1] << " " << rep.c[1][2] << endl;
		}


		real_1d_array x0;
		x0.setlength(k + 4);
		real_1d_array s0;
		s0.setlength(k + 4);
		for (int i = 0; i < k; i++)
		{
			x0[i] = 0.5;
			s0[i] = 1;
		}
		for (int i = k; i < k + 4; i++)
		{
			x0[i] = 0;
			s0[i] = 1;
		}

		C2.normalize();
		C1.normalize();

		auto Q = V3toV2(C1);
		x0[k] = Q.first;
		x0[k + 1] = Q.second;
		Q = V3toV2(C2);
		x0[k + 2] = Q.first;
		x0[k + 3] = Q.second;



		real_1d_array bndl;
		real_1d_array bndu;
		bndl.setlength(k + 4);
		bndu.setlength(k + 4);

		/*real_2d_array c;
		c.setlength(1, neighboor[iter].size() + 5);
		integer_1d_array ct = "[0]";*/
		for (int i = 0; i < k; i++)
		{
			bndl[i] = 0;
			bndu[i] = 1;
			//c[0][i] = 1;
		}
		for (int i = k; i < k + 4; i++)
		{
			bndl[i] = -9999999999;
			bndu[i] = 9999999999;
			//c[0][i] = 0;
		}
		//c[0][k + 4] = double(k) / 2.0;
		minbleicstate state;
		double epsg = 0;
		double epsf = 0;
		double epsx = 0.000001;
		ae_int_t maxits = 0;
		alglib::minbleiccreate(x0, state);
		//alglib::minbleicsetlc(state, c, ct);
		alglib::minbleicsetbc(state, bndl, bndu);
		alglib::minbleicsetscale(state, s0);
		alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);

		alglib::minbleicoptguardsmoothness(state);
		alglib::minbleicoptguardgradient(state, 0.0001);

		/*	for (int i = 0; i < k+4; i++)
			{
				cout << x0[i] << " ";
			}cout << endl;*/



		std::function<void(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)> fop_w_lambda
			= [&](const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr) -> void
		{
			int k = neighboor_M[iter].size();
			int r = Vall.size();
			vector<double > w;
			for (int i = 0; i < k; i++)
			{
				w.push_back(R3[neighboor_M[iter][i]]);
			}
			double maxw = -99999999.0;
			for (int i = 0; i < k; i++)
			{
				double dis = (Vall[iter] - Vall[neighboor_M[iter][i]]).norm();
				if (dis < 0.0001)
				{
					w[i] = 0.0;
				}
				else
				{
					w[i] = 1.0 / (dis * dis);
					maxw = max(maxw, w[i]);
				}
			}
			for (int i = 0; i < k; i++)
			{
				//w[i] = 1.0;
				w[i] = w[i] / maxw;
			}

			double u1 = x[k], v1 = x[k + 1], u2 = x[k + 2], v2 = x[k + 3];
			Eigen::Vector3d n1(sin(u1) * cos(v1), sin(u1) * sin(v1), cos(u1));
			Eigen::Vector3d n2(sin(u2) * cos(v2), sin(u2) * sin(v2), cos(u2));
			n1.normalize();
			n2.normalize();
			func = 0;
			for (int i = 0; i < k; i++)
			{
				auto Nj = Nall[neighboor_M[iter][i]];
				func += x[i] * w[i] * (pow((Nj.x() - sin(u1) * cos(v1)), 2) + pow(Nj.y() - sin(u1) * sin(v1), 2) + pow(Nj.z() - cos(u1), 2))
					+ (1.0 - x[i]) * w[i] * (pow((Nj.x() - sin(u2) * cos(v2)), 2) + pow(Nj.y() - sin(u2) * sin(v2), 2) + pow(Nj.z() - cos(u2), 2));
				//func += x[i] * w[i] * (Nall[neighboor[iter][i]] - n1).squaredNorm()  + (1 - x[i]) * w[i] * (Nall[neighboor[iter][i]] - n2).squaredNorm();
			}

			for (int i = 0; i < k; i++)
			{
				//double a = abs(cos(x[k]))* abs(cos(x[k]));
				auto Nj = Nall[neighboor_M[iter][i]];
				grad[i] = w[i] * (pow((Nj.x() - sin(u1) * cos(v1)), 2) + pow(Nj.y() - sin(u1) * sin(v1), 2) + pow(Nj.z() - cos(u1), 2))
					- w[i] * (pow((Nj.x() - sin(u2) * cos(v2)), 2) + pow(Nj.y() - sin(u2) * sin(v2), 2) + pow(Nj.z() - cos(u2), 2));
			}

			grad[k] = 0; grad[k + 1] = 0; grad[k + 2] = 0; grad[k + 3] = 0;

			for (int i = 0; i < k; i++)
			{
				auto Nj = Nall[neighboor_M[iter][i]];
				grad[k] += 2 * x[i] * sin(u1) * w[i] * (Nj.z() - cos(u1))
					- (2 * x[i] * cos(v1) * cos(u1) * w[i] * (Nj.x() - sin(u1) * cos(v1))
						+ 2 * x[i] * cos(u1) * w[i] * sin(v1) * (Nj.y() - sin(u1) * sin(v1)));

				grad[k + 1] += 2 * x[i] * sin(u1) * sin(v1) * w[i] * (Nj.x() - sin(u1) * cos(v1))
					- 2 * x[i] * sin(u1) * cos(v1) * w[i] * (Nj.y() - sin(u1) * sin(v1));

				grad[k + 2] += 2 * w[i] * sin(u2) * (1 - x[i]) * (Nj.z() - cos(u2))
					- (2 * w[i] * cos(u2) * cos(v2) * (1 - x[i]) * (Nj.x() - sin(u2) * cos(v2))
						+ 2 * w[i] * cos(u2) * sin(v2) * (1 - x[i]) * (Nj.y() - sin(u2) * sin(v2)));

				grad[k + 3] += 2 * sin(u2) * sin(v2) * (1 - x[i]) * w[i] * (Nj.x() - sin(u2) * cos(v2))
					- 2 * sin(u2) * cos(v2) * (1 - x[i]) * w[i] * (Nj.y() - sin(u2) * sin(v2));
			}
		};


		minbleicreport rep2;
		//minbleicoptimize(state, fop, nullptr, nullptr, alglib::parallel);
		alglib::minbleicoptimize(state, fop_w_lambda);
		alglib::minbleicresults(state, x0, rep2);
		//cout << rep2.debugff << endl;
		//double mn = 0;
		//real_1d_array G_tmp;
		//G_tmp.setlength(k + 4);
		//fop_w_lambda(x0, mn, G_tmp, nullptr);
		//cout << mn << endl;
		//printf("%d\n", int(rep2.terminationtype)); // EXPECTED: 4
		//printf("%s\n", x0.tostring(2).c_str()); // EXPECTED: [2,4]

		//optguardreport ogrep;
		//minbleicoptguardresults(state, ogrep);
		//printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
		//printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
		//printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false

		double u1 = x0[k], v1 = x0[k + 1], u2 = x0[k + 2], v2 = x0[k + 3];
		Eigen::Vector3d n1(sin(u1) * cos(v1), sin(u1) * sin(v1), cos(u1));
		Eigen::Vector3d n2(sin(u2) * cos(v2), sin(u2) * sin(v2), cos(u2));
		n1.normalize();
		n2.normalize();
		//cout << n1 << endl;
		//cout << n2 << endl;

		double sum = 0.0;
		for (int i = 0; i < k; i++)
		{
			sum += x0[i];
		}
		if (sum < double(k) / 2.0)
		{
			Nall_new[iter] = n2;
		}
		else
		{
			Nall_new[iter] = n1;
		}
	}

	

	// noisy Nall_new
	/*for (int i = 0; i < r; i++)
	{
		
		Eigen::Vector3d noise = Eigen::Vector3d(rand_num(), rand_num(), rand_num());
		noise.normalize();
		Nall_new[i] = Nall_new[i] + 0.2 * noise;
		Nall_new[i].normalize();
	}*/
	
	
	if (IfoutputFile)
	{
		out.open(outputPath + "\\ShowNormal_" + model + ".xyz");
		for (int i = 0; i < r; i++)
		{
			int k = neighboor[i].size();
			if (k < 5)
			{
				continue;
			}
			out << Vall[i].x() << " " << Vall[i].y() << " " << Vall[i].z() << " " << Nall_new[i].x() << " " << Nall_new[i].y() << " " << Nall_new[i].z() << "\n";
		}
		out.close();
	}

	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "T3 Running Time: " << endtime << endl;
	cout << "T3 Running Time: " << endtime * 1000 << " ms " << endl;
	
	// part 2.5 denoise 

	
	auto neighboor_denoise = neighboor;

	std::function<void(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)> fop_denoise
		= [&](const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr) -> void
	{
		func = 0.0;


		for (int i = 0; i < r; i++)
		{
			grad[i] = 0;
		}
		for (int i = 0; i < r; i++)
		{

			int k = neighboor_denoise[i].size();
			if (k < 5)
			{
				continue;
			}

			for (int j = 0; j < k; j++)
			{
				int pj = neighboor_denoise[i][j];
				Eigen::Vector3d Pj = Vall[pj] + x[pj] * Nall_new[pj];
				Eigen::Vector3d Pi = Vall[i] + x[i] * Nall_new[i];
				Eigen::Vector3d Vij = Pj - Pi;

				Eigen::Matrix3d a;
				a = Vij * Vij.transpose();
				Eigen::Vector3d ans = a * Nall_new[i];
				func += ans.norm();

				auto t0 = Vij;

				double t1 = t0.transpose() * Nall_new[i];
				double t2 = Nall_new[i].transpose() * t0;
				double t3 = t0.transpose() * t0;
				double t4 = Nall_new[i].transpose() * Nall_new[i];


				grad[i] += -1.0 * (2.0 * t1 * t1 * t2 + 2.0 * t1 * t3 * t4);
				double t5 = t0.transpose() * Nall_new[pj];
				double t6 = Nall_new[i].transpose() * Nall_new[pj];

				grad[pj] += 2.0 * t1 * t5 * t2 + 2.0 * t1 * t3 * t6;
			}
		}


		//return loss;
	};

	std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
		= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
	{
		double func = 0.0;


		for (int i = 0; i < r; i++)
		{
			g(i) = 0;
		}
		for (int i = 0; i < r; i++)
		{

			int k = neighboor_denoise[i].size();
			if (k < 5)
			{
				continue;
			}

			for (int j = 0; j < k; j++)
			{
				int pj = neighboor_denoise[i][j];
				Eigen::Vector3d Pj = Vall[pj] + X(pj) * Nall_new[pj];
				Eigen::Vector3d Pi = Vall[i] + X(i) * Nall_new[i];
				Eigen::Vector3d Vij = Pj - Pi;

				Eigen::Matrix3d a;
				a = Vij * Vij.transpose();
				Eigen::Vector3d ans = a * Nall_new[i];
				func += ans.norm();

				auto t0 = Vij;

				double t1 = t0.transpose() * Nall_new[i];
				double t2 = Nall_new[i].transpose() * t0;
				double t3 = t0.transpose() * t0;
				double t4 = Nall_new[i].transpose() * Nall_new[i];


				g[i] += -1.0 * (2.0 * t1 * t1 * t2 + 2.0 * t1 * t3 * t4);
				double t5 = t0.transpose() * Nall_new[pj];
				double t6 = Nall_new[i].transpose() * Nall_new[pj];

				g[pj] += 2.0 * t1 * t5 * t2 + 2.0 * t1 * t3 * t6;
			}
		}
		return func;
	};

	for (int i = 0; i < r; i++)
	{
		int k = neighboor[i].size();
		if (k < 5)
		{
			continue;
		}
		neighboor_denoise[i].clear();
		for (int j = 0; j < k; j++)
		{
			int pid = neighboor[i][j];
			double dis = (Nall_new[i] - Nall_new[pid]).norm();
			if (dis < 0.2)
			{
				neighboor_denoise[i].push_back(pid);
			}
		}
	}

	
	

	

	//cout << "loss = " << fop_denoise() << endl;
	// opt.

	vector<Eigen::Vector3d> Vall_new = Vall;
	start = clock();
	// lbfgs test
	if(1)
	{
		
		BGAL::_LBFGS::_Parameter param = BGAL::_LBFGS::_Parameter();
		param.epsilon = 1e-4;
		param.is_show = true;
		param.max_iteration = 35;
		param.max_linearsearch = 5;
		BGAL::_LBFGS lbfgs(param);


		Eigen::VectorXd iterX(r);
		for (int i = 0; i < r; i++)
		{
			iterX(i) = 0;
		}
		int n = lbfgs.minimize(fg, iterX);
		//int a = 53;
		//int n = lbfgs.test(a);
		//std::cout << iterX << std::endl;
		std::cout << "n: " << n << std::endl;


		for (int i = 0; i < r; i++)
		{
			Vall_new[i] = Vall[i] + iterX(i) * Nall_new[i];
		}


	}
	
	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "T4 Running Time: " << endtime << endl;
	cout << "T4 Running Time: " << endtime * 1000 << " ms " << endl;

	// rnn

	PointCloud<double> cloud1;
	mink = 9999999; maxk = -9999999;

	cloud1.pts.resize(r);
	for (int i = 0; i < r; i++)
	{
		cloud1.pts[i].x = Vall_new[i].x();
		cloud1.pts[i].y = Vall_new[i].y();
		cloud1.pts[i].z = Vall_new[i].z();
	}

	my_kd_tree_t   index1(3 /*dim*/, cloud1, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index1.buildIndex();


	auto neighboor_qc = neighboor;
	neighboor_qc.clear();
		mink = 99999; maxk = -99999;
		for (int i = 0; i < r; i++)
		{
			
			
			vector<int> tmp;
			neighboor_qc.push_back(tmp);
			double query_pt[3] = { Vall_new[i].x(), Vall_new[i].y(), Vall_new[i].z() };
			const double search_radius = static_cast<double>((radis) * (radis));
			std::vector<std::pair<uint32_t, double> >   ret_matches;
			nanoflann::SearchParams params;
			const size_t nMatches = index1.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
			for (size_t j = 0; j < nMatches; j++)
			{
				neighboor_qc[i].push_back(ret_matches[j].first);
			}
			maxk = max(maxk, int(nMatches));
			mink = min(mink, int(nMatches));
		}
		cout << "maxk: " << maxk << "   mink:" << mink << endl;
	

	ofstream fout(outputPath +"\\Denoise_Final_" + model + ".xyz");
	for (size_t i = 0; i < r; i++)
	{
		int k = neighboor_denoise[i].size();
		if (k < 5)
		{
			continue;
		}
	/*	k = neighboor_qc[i].size();
		if (k < 55555)
		{
			continue;
		}*/
		
		
		fout << Vall_new[i].transpose() << " " << Nall_new[i].transpose() << endl;
	}

	Vall = Vall_new;

	
	



	// part 3
	map<int, bool> flag3;

	for (int iter = 0; iter < r; iter++)
	{
		flag3[iter] = 0;
		if (R3[iter] == 1)
		{
			flag3[iter] = 1;
			//continue;
		}
	}
	cout << "Part3 pre compute\n";
	for (int iter = 0; iter < r; iter++)
	{
		if (flag3[iter])
		{
			continue;
		}
		int k = neighboor[iter].size();
		if (k < 5)
		{
			continue;
		}
		
		if (iter % 100000 == 0)
		{
			cout << "Part3 pre Compute iter: " << iter << endl;
		}
		double maxAng = -9999.0;
		/*for (int i = 0; i < k - 1; i++)
		{
			for (int j = i + 1; j < k; j++)
			{
				double sigma = acos(Nall_new[neighboor[iter][i]].dot(Nall_new[neighboor[iter][j]]) / (Nall_new[neighboor[iter][i]].norm() * Nall_new[neighboor[iter][j]].norm()));
				sigma = sigma / EIGEN_PI * 180.0;
				maxAng = max(maxAng, sigma);
			}
		}*/
		for (int i = 1; i < k; i++)
		{
				double sigma = acos(Nall_new[neighboor[iter][i]].dot(Nall_new[neighboor[iter][0]]) / (Nall_new[neighboor[iter][i]].norm() * Nall_new[neighboor[iter][0]].norm()));
				sigma = sigma / EIGEN_PI * 180.0;
				maxAng = max(maxAng, sigma);
		}

		

		if (maxAng < 60)
		{
			flag3[iter] = 1;
			continue;
		}

	}
	bool iw = 0;
	/*while (!iw)
	{
		iw = 1;
		for (int iter = 0; iter < r; iter++)
		{
			if (flag3[iter])
			{
				continue;
			}
			int cnt = 0;
			for (auto p : neighboor[iter])
			{
				if (flag3[p] == 0)
				{
					cnt++;
				}
			}
			if (cnt <= 3)
			{
				flag3[iter] = 1;
				iw = 0;
			}
		}
	}*/
	omp_cnt = 0;
	start = clock();
	//debugOutput = 1;
	map<int, Eigen::Vector3d> NewPoints;
//#pragma omp parallel for schedule(dynamic, 20) //part 3
	cout << "Begin Part3 feature points...\n";
	for (int iter = 0; iter < r; iter++)
	{
		if (flag3[iter])
		{
			continue;
		}
		int k = neighboor[iter].size();
		if (k < 5)
		{
			continue;
		}
	


		omp_cnt++;

		if (omp_cnt % 1000 == 0)
			cout << "Part 3 --- Point iter: " << omp_cnt << " \n";

		real_1d_array x0;
		x0.setlength(3);
		x0[0] = Vall[iter].x();
		x0[1] = Vall[iter].y();
		x0[2] = Vall[iter].z();
		real_1d_array s0 = "[1,1,1]";



		std::function<void(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)> fop_z_lambda
			= [&](const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr) -> void
		{
			int k = neighboor[iter].size();
			int r = Vall.size();


			Eigen::Vector3d z(x[0], x[1], x[2]);
			double mu = 0.01;
			func = 0;
			for (int i = 0; i < k; i++)
			{
				auto zp = Vall[neighboor[iter][i]] - z;
				func += pow(zp.dot(Nall_new[neighboor[iter][i]]), 2);
				//func += x[i] * w[i] * (Nall[neighboor[iter][i]] - n1).squaredNorm()  + (1 - x[i]) * w[i] * (Nall[neighboor[iter][i]] - n2).squaredNorm();
			}
			func = func + mu * (Vall[iter] - z).squaredNorm();

			Eigen::Vector3d g(0, 0, 0);
			g = 2 * mu * (z - Vall[iter]);
			for (int i = 0; i < k; i++)
			{
				//double a = abs(cos(x[k]))* abs(cos(x[k]));
				auto Pj = Vall[neighboor[iter][i]];
				auto nj = Nall_new[neighboor[iter][i]];
				g = g + 2 * (z - Pj).dot(nj) * nj;
			}

			grad[0] = g.x();
			grad[1] = g.y();
			grad[2] = g.z();


		};

		minbleicstate state;
		double epsg = 0;
		double epsf = 0;
		double epsx = 0;
		ae_int_t maxits = 0;
		alglib::minbleiccreate(x0, state);
		//alglib::minbleicsetlc(state, c, ct);
		//alglib::minbleicsetbc(state, bndl, bndu);
		alglib::minbleicsetscale(state, s0);
		alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);

		alglib::minbleicoptguardsmoothness(state);
		alglib::minbleicoptguardgradient(state, 0.000001);

		if (debugOutput)
		{
			for (int i = 0; i < 3; i++)
			{
				cout << x0[i] << " ";
			}cout << endl;
		}



		minbleicreport rep2;
		//minbleicoptimize(state, fop, nullptr, nullptr, alglib::parallel);
		alglib::minbleicoptimize(state, fop_z_lambda);
		alglib::minbleicresults(state, x0, rep2);
		//cout << rep2.debugff << endl;
		double mn = 0;
		real_1d_array G_tmp;
		G_tmp.setlength(3);
		fop_z_lambda(x0, mn, G_tmp, nullptr);
		if (debugOutput)
			cout << mn << endl;
		if (debugOutput)
		{
			printf("%d\n", int(rep2.terminationtype)); // EXPECTED: 4
			printf("%s\n", x0.tostring(2).c_str()); // EXPECTED: [2,4]

			optguardreport ogrep;
			minbleicoptguardresults(state, ogrep);
			printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
			printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
			printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
		}
		double dis = (Vall[iter] - Eigen::Vector3d(x0[0], x0[1], x0[2])).norm();
		if (dis > 0.0001)
		{
			NewPoints[iter] = Eigen::Vector3d(x0[0], x0[1], x0[2]);
			
		}
		else
		{
			flag3[iter] = true;
		}
		
		
		//NewPoints[iter] = Eigen::Vector3d(x0[0], x0[1], x0[2]);
	}



	for (auto np : NewPoints)
	{
		if (flag3[np.first])
		{
			continue;
		}
		double dis = (Vall[np.first] - np.second).norm();
		if (dis > radis*2 || dis < 0.00001)
		{
			flag3[np.first] = 1;
		}

	}

	if (IfoutputFile) {
		/*out.open("FinalPointCloud.xyz");
		for (int i = 0; i < r; i++)
		{
			out << Vall[i].transpose() << endl;
		}

		for (auto np : NewPoints)
		{
			if (flag3[np.first])
			{
				continue;
			}
			out << np.second.transpose() << endl;
		}
		out.close();*/

		out.open(outputPath + "\\FeaturePoints_" + model + ".xyz");
		for (auto np : NewPoints)
		{
			if (flag3[np.first])
			{
				continue;
			}
			out << np.second.transpose() << endl;
		}
		out.close();

		out.open(outputPath + "\\FinalPointCloud1_" + model + ".xyz");
		vector<bool> flag2(r, 0);
		for (auto np : NewPoints)
		{
			if (flag3[np.first])
			{
				continue;
			}
			out << np.second.transpose() <<" "<< Nall_new[np.first].transpose() << endl;
			flag2[np.first] = 1;
		}
		for (size_t i = 0; i < r; i++)
		{
			int k = neighboor_denoise[i].size();
			if (k < 5)
			{
				continue;
			}
			
			if(!flag2[i])
			out << Vall[i].transpose() << " " << Nall_new[i].transpose() << endl;
		}

		/*for (size_t i = 0; i < r; i++)
		{
			int k = neighboor_denoise[i].size();
			if (k < 5)
			{
				continue;
			}
			out << Vall[i].transpose() << " " << Nall_new[i].transpose() << endl;
		}*/
		out.close();

		out.open(outputPath + "\\FinalPointCloud2_" + model + ".xyz");
		vector<bool> flag22(r, 0);
		for (auto np : NewPoints)
		{
			if (flag3[np.first])
			{
				continue;
			}
			out << np.second.transpose() << " " << Nall_new[np.first].transpose() << endl;
			flag22[np.first] = 1;
		}
		for (size_t i = 0; i < r; i++)
		{
			int k = neighboor_denoise[i].size();
			if (k < 5)
			{
				continue;
			}

			//if (!flag2[i])
			out << Vall[i].transpose() << " " << Nall_new[i].transpose() << endl;
		}

		out.close();


		out.open(outputPath + "\\FinalPointCloud_WithoutFeature_" + model + ".xyz");
		
		for (auto np : NewPoints)
		{
			if (flag3[np.first])
			{
				continue;
			}
			//out << np.second.transpose() << " " << Nall_new[np.first].transpose() << endl;
			//flag2[np.first] = 1;
		}
		for (size_t i = 0; i < r; i++)
		{
			int k = neighboor_denoise[i].size();
			if (k < 5)
			{
				continue;
			}

			if(!flag2[i])
			out << Vall[i].transpose() << " " << Nall_new[i].transpose() << endl;
		}
		out.close();
		
		out.open(outputPath + "\\FeaturePointsNum_" + model + ".txt");
		int cnt = 0;
		for (auto np : NewPoints)
		{
			if (flag3[np.first])
			{
				continue;
			}
			cnt++;
			
		}
		out << cnt << endl;
		out.close();

		out.open(outputPath + "\\radis_" + model + ".txt");
		
		out << radis << endl;
		out.close();
		
	}


	end = clock();
	endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "T5 Running Time: " << endtime << endl;
	cout << "T5 Running Time: " << endtime * 1000 << " ms " << endl;

}


void DenoiseTest(string modelpath,string model)
{
	clock_t start, end;

	bool debugOutput = 0, IfoutputFile = 1;
	vector<Eigen::Vector3d> Vall, Nall;
	vector<vector<int>> neighboor;

	ofstream out;
	//alglib::setglobalthreading(alglib::parallel);

	//string model = "angleWithNor";
	//string model = "0.005_50000_00040123_8fc7d06caf264003a242597a_trimesh_000";

	
	string outputPath = modelpath;

	//ifstream in("E:\\Dropbox\\MyProjects\\SIG-2022-Feature-preserving-recon\\data\\Noise\\abc_chunk4\\0.0025\\" + model + ".xyz");

	ifstream in(outputPath + model + ".xyz");
	int ppid = 0;
	while (!in.eof())
	{
		Eigen::Vector3d p, n;
		in >> p.x() >> p.y() >> p.z() >> n.x() >> n.y() >> n.z();
		ppid++;
		//if (ppid % 10 == 0)
		{
			Vall.push_back(p);
			Nall.push_back(n);
		}
		
	}

	cout << "Read PointCloud. xyz File. \n";
	int r = Vall.size();
	Eigen::Vector3d maxp(-99999, -99999, -99999);
	Eigen::Vector3d minp(99999, 99999, 99999);

	for (int i = 0; i < r; i++)
	{
		maxp.x() = max(maxp.x(), Vall[i].x());
		maxp.y() = max(maxp.y(), Vall[i].y());
		maxp.z() = max(maxp.z(), Vall[i].z());
		minp.x() = min(minp.x(), Vall[i].x());
		minp.y() = min(minp.y(), Vall[i].y());
		minp.z() = min(minp.z(), Vall[i].z());
		Nall[i].normalize();
	}
	double minn = min(minp.x(), min(minp.y(), minp.z()));
	double maxn = max(maxp.x(), max(maxp.y(), maxp.z()));
	double maxl = max(maxp.x() - minp.x(), max(maxp.y() - minp.y(), maxp.z() - minp.z()));
	//for (int i = 0; i < r; i++)
	//{
	//	Vall[i].x() = (Vall[i].x() - minp.x()) / (maxl);
	//	Vall[i].y() = (Vall[i].y() - minp.y()) / (maxl);
	//	Vall[i].z() = (Vall[i].z() - minp.z()) / (maxl);
	//} //
	//out.open("01_" + model + ".xyz");
	//for (size_t i = 0; i < r; i++)
	//{
	//	out << Vall[i].transpose() << " " << Nall[i].transpose() << endl;
	//}
	//out.close();

	double radis = 0.001;
	double lambda = 0.05;

	// rnn

	PointCloud<double> cloud;
	int mink = 9999999, maxk = -9999999;

	cloud.pts.resize(r);
	for (int i = 0; i < r; i++)
	{
		cloud.pts[i].x = Vall[i].x();
		cloud.pts[i].y = Vall[i].y();
		cloud.pts[i].z = Vall[i].z();
	}
	typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<double, PointCloud<double> >,
		PointCloud<double>,
		3 /* dim */
	> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();

	do
	{
		radis += 0.001;
		neighboor.clear();
		mink = 99999; maxk = -99999;
		for (int i = 0; i < r; i++)
		{
			vector<int> tmp;
			neighboor.push_back(tmp);
			double query_pt[3] = { Vall[i].x(), Vall[i].y(), Vall[i].z() };
			const double search_radius = static_cast<double>((radis) * (radis));
			std::vector<std::pair<uint32_t, double> >   ret_matches;
			nanoflann::SearchParams params;
			const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
			for (size_t j = 0; j < nMatches; j++)
			{
				neighboor[i].push_back(ret_matches[j].first);
			}
			maxk = max(maxk, int(nMatches));
			mink = min(mink, int(nMatches));
		}
		cout << "maxk: " << maxk << "   mink:" << mink << endl;

	}  while (maxk < rnnnum );

	// smooth 
	//for (int i = 0; i < r; i++)
	//{
	//	int k = neighboor[i].size();
	//	if (k < 5)
	//	{
	//		continue;
	//	}
	//	
	//	Eigen::Vector3d sum(0, 0, 0);
	//	for (int j = 0; j < k; j++)
	//	{
	//		sum += Nall[neighboor[i][j]];
	//	}
	//	sum /= k;
	//	Nall[i] = sum;
	//}			


	
	// add par later

	//auto neighboor_denoise = neighboor;
	//for (int i = 0; i < r; i++)
	//{
	//	int k = neighboor[i].size();
	//	if (k < 5)
	//	{
	//		continue;
	//	}
	//	neighboor_denoise[i].clear();
	//	for (int j = 0; j < k; j++)
	//	{
	//		int pid = neighboor[i][j];
	//		double dis = (Nall[i] - Nall[pid]).norm();
	//		if (dis < 0.3)
	//		{
	//			neighboor_denoise[i].push_back(pid);
	//		}
	//	}
	//}
	//neighboor = neighboor_denoise;


	start = clock();


	auto Vall_new = Vall;
	auto Nall_new = Nall;

	std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
		= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
	{
		double func = 0.0;
		double lambda = 0.01;

		for (int i = 0; i < r*3; i++)
		{
			g(i) = 0;
		}
		auto Nall_tmp = Nall_new;
		for (int i = 0; i < r; i++)
		{
			double u = X(r + i * 2);
			double v = X(r + i * 2 + 1);
			Nall_tmp[i].x() = sin(u) * cos(v);
			Nall_tmp[i].y() = sin(u) * sin(v);
			Nall_tmp[i].z() = cos(u);
		}
		
		for (int i = 0; i < r; i++)
		{
			

			int k = neighboor[i].size();
			if (k < 5)
			{
				continue;
			}

			for (int j = 0; j < k; j++)
			{
				int pj = neighboor[i][j];
				Eigen::Vector3d Pj = Vall[pj] + X(pj) * Nall_tmp[pj];
				Eigen::Vector3d Pi = Vall[i] + X(i) * Nall_tmp[i];
				Eigen::Vector3d Vij = Pj - Pi;

				Eigen::Matrix3d a;
				a = Vij * Vij.transpose();
				Eigen::Vector3d ans = a * Nall_tmp[i];
				func += ans.norm();

				auto t0 = Vij;

				double t1 = t0.transpose() * Nall_tmp[i];
				double t2 = Nall_tmp[i].transpose() * t0;
				double t3 = t0.transpose() * t0;
				double t4 = Nall_tmp[i].transpose() * Nall_tmp[i];


				g[i] += -1.0 * (2.0 * t1 * t1 * t2 + 2.0 * t1 * t3 * t4);
				double t5 = t0.transpose() * Nall_tmp[pj];
				double t6 = Nall_tmp[i].transpose() * Nall_tmp[pj];

				g[pj] += 2.0 * t1 * t5 * t2 + 2.0 * t1 * t3 * t6;

				// g of ni
				{
					Eigen::Vector3d t0 = Vij;
					Eigen::Vector3d t1 = t0.transpose() * Nall_tmp[i] * t0;
					double t2 = 2 * X(i);
					double t3 = t0.transpose() * t0;
					double t4 = Nall_tmp[i].transpose() * t0;
					Eigen::Vector3d gni(0, 0, 0);
					
					gni = 2.0 * t3 * t1 -(t2 * t4 * t1 + t2 * t3 * t4 * Nall_tmp[i]);
					
					g[r + i * 2] += gni.x() * cos(X(r + i * 2)) * cos(X(r + i * 2 + 1))  + gni.y() * cos(X(r + i * 2)) * sin(X(r + i * 2 + 1)) + -1.0 * gni.z() * sin(X(r + i * 2));
					g[r + i * 2 + 1] += -1.0* gni.x() * sin(X(r + i * 2)) * sin(X(r + i * 2 + 1)) + gni.y() * sin(X(r + i * 2)) * cos(X(r + i * 2 + 1));
					
				}
				// g of nj
				{
					Eigen::Vector3d t0 = Vij;
					double t1 = 2*X(pj);
					double t2 = Nall_tmp[i].transpose() * t0;
					double t3 = t0.transpose() * Nall_tmp[i];
					double t4 = t0.transpose() * t0;
					Eigen::Vector3d gnj(0,0,0);
					gnj = t1 * t2 * t3 * t0 + t1 * t2 * t4*Nall_tmp[i];
					g[r + pj * 2] += gnj.x() * cos(X(r + pj * 2)) * cos(X(r + pj * 2 + 1)) + gnj.y() * cos(X(r + pj * 2)) * sin(X(r + pj * 2 + 1)) + -1.0 * gnj.z() * sin(X(r + pj * 2));
					g[r + pj * 2 + 1] += -1.0 * gnj.x() * sin(X(r + pj * 2)) * sin(X(r + pj * 2 + 1)) + gnj.y() * sin(X(r + pj * 2)) * cos(X(r + pj * 2 + 1));
				}

			}
		}
		for (int i = 0; i < r; i++)
		{
			Eigen::Vector3d Pi = Vall[i] + X(i) * Nall_tmp[i];
			func += lambda* (Pi - Vall[i]).squaredNorm();
			double tt = (Pi - Vall[i]).transpose() * Nall_tmp[i];
			g[i] += lambda * 2 * tt;
		}

		//double cnt = 0;
		//for (int i = 0; i < r; i++)
		//{
		//	cnt +=  X(i);
		//	//g[i] += lambda * 2 * X(i);
		//}
		//for (int i = 0; i < r; i++)
		//{
		//	//cnt += X(i);
		//	g[i] += lambda * 2 * cnt;
		//}
		//func += lambda * cnt * cnt;

		return func;
	};


	
	BGAL::_LBFGS::_Parameter param = BGAL::_LBFGS::_Parameter();
	param.epsilon = 1e-5;
	param.is_show = true;
	BGAL::_LBFGS lbfgs(param);


	Eigen::VectorXd iterX(r*3);
	for (int i = 0; i < r; i++)
	{
		iterX(i) = 0;
	}
	for (int i = 0; i < r; i++)
	{
		auto Q = V3toV2(Nall_new[i]);
		iterX(r+i*2) = Q.first;
		iterX(r+i*2+1) = Q.second;
	}
	int n = lbfgs.minimize(fg, iterX);
	//int a = 53;
	//int n = lbfgs.test(a);
	//std::cout << iterX << std::endl;
	std::cout << "n: " << n << std::endl;


	for (int i = 0; i < r; i++)
	{
		double u = iterX(r+i*2);
		double v = iterX(r+i*2+1);
		Nall_new[i].x() = sin(u) * cos(v);
		Nall_new[i].y() = sin(u) * sin(v);
		Nall_new[i].z() = cos(u);
		Vall_new[i] = Vall[i] + iterX(i) * Nall_new[i];
	}



	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "Running Time: " << endtime << endl;
	cout << "Running Time: " << endtime * 1000 << " ms " << endl;


	ofstream fout(outputPath+"\\Denoise_"+ model +".xyz");
	
	for (size_t i = 0; i < r; i++)
	{
		int k = neighboor[i].size();
		if (k < 5)
		{
			continue;
		}
		fout << Vall_new[i].transpose() << " " << Nall_new[i].transpose() << endl;
	}
	fout.close();

}



int main()
{
	string modelpath = "..//..//data//";
	
	string modelname = "01_82-block";


	std::cout << "====================RFEPSTest" << std::endl;

	rnnnum = 60;
	bool denoise = 0;
	//if your input point cloud is noisy, please open the following code
	if (denoise)
	{
		rnnnum *= 2.0;
		//DenoiseTest(modelpath,modelname);
	}

	RFEPSTest(modelpath,modelname, denoise);
	Poisson(modelpath, modelname);
	Comput_rnn(modelpath,modelname);
	Comput_RPD(modelpath,modelname);
	
	std::cout << "successful!" << std::endl;
	
	return 0;
}