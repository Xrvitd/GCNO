#include <iostream>
#include <fstream>
#include <functional>
#include <io.h>
#include <random>
#include <omp.h>
#include <BGAL/Optimization/LinearSystem/LinearSystem.h>
#include <BGAL/Optimization/ALGLIB/optimization.h>
#include <BGAL/Optimization/LBFGS/LBFGS.h>
#include <BGAL/BaseShape/Point.h>
#include <BGAL/BaseShape/Polygon.h>
#include <BGAL/Tessellation2D/Tessellation2D.h>
#include <BGAL/Draw/DrawPS.h>
#include <BGAL/Integral/Integral.h>
#include <BGAL/Model/ManifoldModel.h>
#include <BGAL/Model/Model_Iterator.h>
#include <BGAL/Optimization/GradientDescent/GradientDescent.h>
#include <BGAL/Tessellation3D/Tessellation3D.h>
#include <BGAL/BaseShape/KDTree.h>
#include <BGAL/PointCloudProcessing/Registration/ICP/ICP.h>
#include <BGAL/Reconstruction/MarchingTetrahedra/MarchingTetrahedra.h>
#include <BGAL/Geodesic/Dijkstra/Dijkstra.h>
#include <BGAL/CVTLike/CPD.h>
#include <BGAL/CVTLike/CVT.h>
//Test BOC sign
void BOCSignTest()
{
	BGAL::_BOC::_Sign res = BGAL::_BOC::sign_(1e-6);
	if (res == BGAL::_BOC::_Sign::ZerO)
	{
		std::cout << "0" << std::endl;
	}
	else if (res == BGAL::_BOC::_Sign::PositivE)
	{
		std::cout << "1" << std::endl;
	}
	else if (res == BGAL::_BOC::_Sign::NegativE)
	{
		std::cout << "-1" << std::endl;
	}
	BGAL::_BOC::set_precision_(1e-5);
	res = BGAL::_BOC::sign_(1e-6);
	if (res == BGAL::_BOC::_Sign::ZerO)
	{
		std::cout << "0" << std::endl;
	}
	else if (res == BGAL::_BOC::_Sign::PositivE)
	{
		std::cout << "1" << std::endl;
	}
	else if (res == BGAL::_BOC::_Sign::NegativE)
	{
		std::cout << "-1" << std::endl;
	}
}

//Test LinearSystem
void LinearSystemTest()
{
	Eigen::MatrixXd m(3, 3);
	m.setZero();
	m(0, 0) = 1;
	m(1, 1) = 2;
	m(2, 2) = 3;
	Eigen::SparseMatrix<double> h = m.sparseView();
	Eigen::VectorXd r(3);
	r(0) = 2;
	r(1) = 4;
	r(2) = 6;
	Eigen::VectorXd res = BGAL::_LinearSystem::solve_ldlt(h, r);
	std::cout << res << std::endl;
}
//****************************************

//Test ALGLIB
void function1_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
{
	//
	// this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
	// and its derivatives df/d0 and df/dx1
	//
	func = 100 * pow(x[0] + 3, 4) + pow(x[1] - 3, 4);
	grad[0] = 400 * pow(x[0] + 3, 3);
	grad[1] = 4 * pow(x[1] - 3, 3);
}
void ALGLIBTest()
{
	//
	// This example demonstrates minimization of
	//
	//     f(x,y) = 100*(x+3)^4+(y-3)^4
	//
	// using LBFGS method, with:
	// * initial point x=[0,0]
	// * unit scale being set for all variables (see minlbfgssetscale for more info)
	// * stopping criteria set to "terminate after short enough step"
	// * OptGuard integrity check being used to check problem statement
	//   for some common errors like nonsmoothness or bad analytic gradient
	//
	// First, we create optimizer object and tune its properties
	//
	alglib::real_1d_array x = "[0,0]";
	alglib::real_1d_array s = "[1,1]";
	double epsg = 0;
	double epsf = 0;
	double epsx = 0.0000000001;
	alglib::ae_int_t maxits = 0;
	alglib::minlbfgsstate state;
	alglib::minlbfgscreate(1, x, state);
	alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
	alglib::minlbfgssetscale(state, s);

	//
	// Activate OptGuard integrity checking.
	//
	// OptGuard monitor helps to catch common coding and problem statement
	// issues, like:
	// * discontinuity of the target function (C0 continuity violation)
	// * nonsmoothness of the target function (C1 continuity violation)
	// * erroneous analytic gradient, i.e. one inconsistent with actual
	//   change in the target/constraints
	//
	// OptGuard is essential for early prototyping stages because such
	// problems often result in premature termination of the optimizer
	// which is really hard to distinguish from the correct termination.
	//
	// IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
	//            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
	//
	//            Other OptGuard checks add moderate overhead, but anyway
	//            it is better to turn them off when they are not needed.
	//
	alglib::minlbfgsoptguardsmoothness(state);
	alglib::minlbfgsoptguardgradient(state, 0.001);

	//
	// Optimize and examine results.
	//
	alglib::minlbfgsreport rep;
	alglib::minlbfgsoptimize(state, function1_grad);
	alglib::minlbfgsresults(state, x, rep);
	printf("%s\n", x.tostring(2).c_str()); // EXPECTED: [-3,3]

	//
	// Check that OptGuard did not report errors
	//
	// NOTE: want to test OptGuard? Try breaking the gradient - say, add
	//       1.0 to some of its components.
	//
	alglib::optguardreport ogrep;
	alglib::minlbfgsoptguardresults(state, ogrep);
	printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
	printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
	printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false
}
//***************************************

//LBFGSTest
void LBFGSTest()
{
	BGAL::_LBFGS::_Parameter param = BGAL::_LBFGS::_Parameter();
	param.epsilon = 1e-10;
	param.is_show = true;
	BGAL::_LBFGS lbfgs(param);
	class problem
	{
	public:
		double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& g)
		{
			double fval = (x(0) - 1) * (x(0) - 1) + (x(1) - 1) * (x(1) - 1);
			g.setZero();
			g(0) = 2 * (x(0) - 1);
			g(1) = 2 * (x(1) - 1);
			return fval;
		}
	};
	problem fun = problem();
	Eigen::VectorXd iterX(2);
	iterX(0) = 53;
	iterX(1) = -68;
	int n = lbfgs.minimize(fun, iterX);
	int a = 53;
	//int n = lbfgs.test(a);
	std::cout << iterX << std::endl;
	std::cout << "n: " << n << std::endl;
}
//***********************************

//BaseShapeTest
void BaseShapeTest()
{
	BGAL::_Point2 p0(0.3, 0.25);
	BGAL::_Point3 p1(1, 2, 3);
	BGAL::_Point3 p2(2, 3, 4);
	BGAL::_Polygon boundary;
	boundary.start_();
	boundary.insert_(0, 0);
	boundary.insert_(1, 0);
	boundary.insert_(1, 1);
	boundary.insert_(0, 1);
	boundary.end_();
	if (boundary.is_in_(p0))
	{
		std::cout << "yes!" << std::endl;
	}
	std::cout << p1.dot_(p2) << std::endl;
}
//***********************************

//TessellationTest2D
void Tessellation2DTest()
{
	BGAL::_Polygon boundary;
	boundary.start_();
	boundary.insert_(BGAL::_Point2(0, 0));
	boundary.insert_(BGAL::_Point2(1, 0));
	boundary.insert_(BGAL::_Point2(1, 1));
	boundary.insert_(BGAL::_Point2(0, 1));
	boundary.end_();
	int num_sites = 2;
	std::vector<BGAL::_Point2> sites;
	//sites.push_back(BGAL::_Point2(0.5, 0.5));
	sites.push_back(BGAL::_Point2(0.5, 0));
	//sites.push_back(BGAL::_Point2(0.4, 0.3));
	//sites.push_back(BGAL::_Point2(0.3, 0.4));
	sites.push_back(BGAL::_Point2(0, 0.5));
	//for (int i = 0; i < num_sites; ++i)
	//{
	//    sites.push_back(BGAL::_Point2(BGAL::_BOC::rand_(), BGAL::_BOC::rand_()));
	//}
	BGAL::_Tessellation2D vor(boundary, sites);
	std::vector<BGAL::_Polygon> cells = vor.get_cell_polygons_();
	std::vector<std::vector<std::pair<int, int>>> edges = vor.get_cells_();
	std::ofstream out("data\\TessellationTest.ps");
	BGAL::_PS ps(out);
	ps.set_bbox_(boundary.bounding_box_());

	for (int i = 0; i < num_sites; ++i)
	{
		ps.draw_polygon_(cells[i], 0.001);
		ps.draw_point_(sites[i], 0.005, 1, 0, 0);
	}
	ps.end_();
	out.close();
}
//***********************************

//CVTLBFGSTest
void CVTLBFGSTest()
{
	BGAL::_Polygon boundary;
	boundary.start_();
	boundary.insert_(BGAL::_Point2(0, 0));
	boundary.insert_(BGAL::_Point2(1, 0));
	boundary.insert_(BGAL::_Point2(1, 1));
	boundary.insert_(BGAL::_Point2(0, 1));
	boundary.end_();
	std::vector<BGAL::_Point2> sites;
	int num = 100;
	for (int i = 0; i < num; ++i)
	{
		sites.push_back(BGAL::_Point2(BGAL::_BOC::rand_(), BGAL::_BOC::rand_()));
	}
	BGAL::_Tessellation2D voronoi(boundary, sites);
	std::vector<BGAL::_Polygon> cells = voronoi.get_cell_polygons_();
	std::function<double(BGAL::_Point2 p)> rho
		= [](BGAL::_Point2 p)
	{
		return 1.0;
	};
	std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
		= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
	{
		for (int i = 0; i < num; ++i)
		{
			BGAL::_Point2 p(X(i * 2), X(i * 2 + 1));
			if (!boundary.is_in_(p))
			{
				return std::numeric_limits<double>::max();
			}
			sites[i] = p;
		}
		voronoi.calculate_(boundary, sites);
		cells = voronoi.get_cell_polygons_();
		double energy = 0;
		g.setZero();
		for (int i = 0; i < num; ++i)
		{
			Eigen::VectorXd inte = BGAL::_Integral::integral_polygon_fast(
				[&](BGAL::_Point2 p)
				{
					Eigen::VectorXd r(3);
					r(0) = rho(p) * ((sites[i] - p).sqlength_());
					r(1) = 2 * rho(p) * (sites[i].x() - p.x());
					r(2) = 2 * rho(p) * (sites[i].y() - p.y());
					return r;
				}, cells[i]
					);
			energy += inte(0);
			g(i * 2) = inte(1);
			g(i * 2 + 1) = inte(2);
		}
		return energy;
	};
	BGAL::_LBFGS::_Parameter para;
	para.is_show = true;
	para.epsilon = 1e-4;
	BGAL::_LBFGS lbfgs(para);
	Eigen::VectorXd iterX(num * 2);
	for (int i = 0; i < num; ++i)
	{
		iterX(i * 2) = sites[i].x();
		iterX(i * 2 + 1) = sites[i].y();
	}
	lbfgs.minimize(fg, iterX);
	for (int i = 0; i < num; ++i)
	{
		sites[i] = BGAL::_Point2(iterX(i * 2), iterX(i * 2 + 1));
	}
	voronoi.calculate_(boundary, sites);
	cells = voronoi.get_cell_polygons_();
	std::ofstream out("data\\CVTLBFGSTest.ps");
	BGAL::_PS ps(out);
	ps.set_bbox_(boundary.bounding_box_());
	for (int i = 0; i < cells.size(); ++i)
	{
		ps.draw_point_(sites[i], 0.005, 1, 0, 0);
		ps.draw_polygon_(cells[i], 0.005);
	}
	ps.end_();
	out.close();
}
//

//DrawTest
void DrawTest()
{
	BGAL::_Polygon boundary;
	boundary.start_();
	boundary.insert_(BGAL::_Point2(0, 0));
	boundary.insert_(BGAL::_Point2(1, 0));
	boundary.insert_(BGAL::_Point2(1, 1));
	boundary.insert_(BGAL::_Point2(0, 1));
	boundary.end_();
	BGAL::_Polygon poly;
	poly.start_();
	poly.insert_(BGAL::_Point2(0.25, 0.25));
	poly.insert_(BGAL::_Point2(0.75, 0.25));
	poly.insert_(BGAL::_Point2(0.75, 0.75));
	poly.insert_(BGAL::_Point2(0.25, 0.75));
	poly.end_();
	std::ofstream out("data\\DrawTest.ps");
	BGAL::_PS ps(out);
	ps.set_bbox_(boundary.bounding_box_());
	ps.draw_polygon_(poly, 0.01);
	ps.end_();
	out.close();
}
//***********************************

//IntegralTest
void IntegralTest()
{
	BGAL::_Polygon boundary;
	boundary.start_();
	boundary.insert_(BGAL::_Point2(0, 0));
	boundary.insert_(BGAL::_Point2(1, 0));
	boundary.insert_(BGAL::_Point2(1, 1));
	boundary.insert_(BGAL::_Point2(0, 1));
	boundary.end_();
	Eigen::VectorXd r = BGAL::_Integral::integral_polygon_fast(
		[](BGAL::_Point2 p)
		{
			Eigen::VectorXd res(1);
			res(0) = p.x();
			return res;
		}, boundary
	);
	std::cout << r << std::endl;
}
//***********************************

//ModelTest
void ModelTest()
{
	BGAL::_ManifoldModel model("data\\sphere.obj");
	std::cout << "V number: " << model.number_vertices_() << std::endl;
	std::cout << "F number: " << model.number_faces_() << std::endl;
	model.initialization_PQP_();
	if (model.is_in_(BGAL::_Point3(1.25, 0.1, 0.05)))
	{
		std::cout << "in" << std::endl;
	}
}
//***********************************

//Tessellation3DTest
void Tessellation3DTest()
{
	BGAL::_ManifoldModel model("data\\sphere.obj");
	int num = 20;
	std::vector<BGAL::_Point3> sites;
	for (int i = 0; i < num; ++i)
	{
		double phi = BGAL::_BOC::PI() * 2.0 * BGAL::_BOC::rand_();
		double theta = BGAL::_BOC::PI() * BGAL::_BOC::rand_();
		sites.push_back(BGAL::_Point3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)));
	}
	std::vector<double> weights(num, 0);
	BGAL::_Restricted_Tessellation3D RVD(model, sites, weights);
	const std::vector<std::vector<std::tuple<int, int, int>>>& cells = RVD.get_cells_();
	std::ofstream out("data\\Tessellation3DTest.obj");
	out << "g 3D_Object\nmtllib BKLineColorBar.mtl\nusemtl BKLineColorBar" << std::endl;
	for (int i = 0; i < RVD.number_vertices_(); ++i)
	{
		out << "v " << RVD.vertex_(i) << std::endl;
	}
	for (int i = 0; i < cells.size(); ++i)
	{
		double color = (double)BGAL::_BOC::rand_();
		out << "vt " << color << " 0" << std::endl;
		for (int j = 0; j < cells[i].size(); ++j)
		{
			out << "f " << std::get<0>(cells[i][j]) + 1 << "/" << i + 1
				<< " " << std::get<1>(cells[i][j]) + 1 << "/" << i + 1
				<< " " << std::get<2>(cells[i][j]) + 1 << "/" << i + 1 << std::endl;
		}
	}
	//for (int i = 0; i < cells.size(); ++i)
	//{
	//	double color = (double)(i % 16) / (15.0);
	//	out << "vt " << color << " 0" << std::endl;
	//	for (int j = 0; j < cells[i].size(); ++j)
	//	{
	//		out << "f " << std::get<0>(cells[i][j]) + 1 << "/" << i + 1
	//			<< " " << std::get<1>(cells[i][j]) + 1 << "/" << i + 1
	//			<< " " << std::get<2>(cells[i][j]) + 1 << "/" << i + 1 << std::endl;
	//	}
	//}
	out.close();
}
//***********************************

////ReadFileTest
//void ReadFileTest()
//{
//	std::string path("data");
//	std::vector<string> files;
//	int num = BGAL::_BOC::search_files_(path, ".obj", files);
//	std::cout << "num: " << num << std::endl;
//	for (int i = 0; i < files.size(); ++i)
//	{
//		std::cout << files[i] << std::endl;
//	}
//}
////***********************************

//KDTreeTest
void KDTreeTest()
{
	std::vector<BGAL::_Point3> pts;
	int num = 100;
	std::ifstream ip("data\\KDTreeTest.txt");
	for (int i = 0; i < num; ++i)
	{
		double x, y, z;
		ip >> x >> y >> z;
		pts.push_back(BGAL::_Point3(x, y, z));
	}
	BGAL::_KDTree kdt(pts);
	double mind = 0;
	BGAL::_Point3 query(-0.03988281, -0.9964109, 0.02454429);
	int id = kdt.search_(query, mind);
	std::cout << id << " " << setprecision(20) << pts[id] << std::endl;
	std::cout << mind << "    " << (query - pts[id]).length_() << std::endl;
	double bmind = 10000;
	int bid = -1;
	for (int i = 0; i < num; ++i)
	{
		if ((pts[i] - query).length_() < bmind)
		{
			bmind = (pts[i] - query).length_();
			bid = i;
		}
	}
	std::cout << "my:  " << bid << "\t" << bmind << std::endl;
}
//***********************************

//ICPTest
void ICPTest()
{
	std::vector<BGAL::_Point3> pts;
	int num = 2895;
	std::ifstream ip("data\\ICPTest.txt");
	for (int i = 0; i < num; ++i)
	{
		double x, y, z;
		ip >> x >> y >> z;
		pts.push_back(BGAL::_Point3(x, y, z));
	}
	std::vector<BGAL::_Point3> dpts(num);
	Eigen::Matrix3d R;
	R.setIdentity();
	double theta = BGAL::_BOC::PI() * 0.5 * 0.125;
	R(0, 0) = 1 - 2 * sin(theta) * sin(theta);
	R(0, 1) = -2 * cos(theta) * sin(theta);
	R(1, 0) = 2 * cos(theta) * sin(theta);;
	R(1, 1) = 1 - 2 * sin(theta) * sin(theta);
	for (int i = 0; i < num; ++i)
	{
		dpts[i] = pts[i].rotate_(R) + BGAL::_Point3(0.8, 0.2, -0.15);
	}
	BGAL::_ICP icp(pts);
	Eigen::Matrix4d RTM = icp.registration_(dpts);
	std::cout << RTM << std::endl;
	Eigen::Matrix4d RRTM;
	RRTM.setIdentity();
	RRTM.block<3, 3>(0, 0) = R;
	RRTM(0, 3) = 0.8;
	RRTM(1, 3) = 0.2;
	RRTM(2, 3) = -0.15;
	std::cout << RRTM * RTM << std::endl;
}
//***********************************

//MarchingTetrahedraTest
void MarchingTetrahedraTest()
{
	double bboxl = 4.35;
	std::pair<BGAL::_Point3, BGAL::_Point3> bbox(
		BGAL::_Point3(-bboxl, -bboxl, -bboxl), BGAL::_Point3(bboxl, bboxl, bboxl));
	BGAL::_Marching_Tetrahedra MT(bbox, 6);
	MT.set_method_(0);
	function<std::pair<double, BGAL::_Point3>(BGAL::_Point3 p)> dis
		= [](BGAL::_Point3 p)
	{
		double x, y, z;
		x = p.x();
		y = p.y();
		z = p.z();
		double len = (x * x * x * x + y * y * y * y + z * z * z * z) / 16 - (x * x + y * y + z * z) / 4 + 0.4;
		//len = (1 - sqrt(x * x + y * y)) * (1 - sqrt(x * x + y * y)) + z * z - 0.16;
		double dx, dy, dz;
		dx = 4 * x * x * x / 16 - 2 * x / 4;
		dy = 4 * y * y * y / 16 - 2 * y / 4;
		dz = 4 * z * z * z / 16 - 2 * z / 4;
		//dx = -2 * (x / (sqrt(x * x + y * y)) - x);
		//dy = -2 * (y / (sqrt(x * x + y * y)) - y);
		//dz = 2 * z;
		double dl = (dx * dx + dy * dy + dz * dz);
		dx /= dl;
		dy /= dl;
		dz /= dl;
		BGAL::_Point3 rp(x - dx * len, y - dy * len, z - dz * len);
		std::pair<double, BGAL::_Point3> res(len, rp);
		return res;
	};
	BGAL::_ManifoldModel model = MT.reconstruction_(dis);
	model.save_obj_file_("data\\MarchingTetrahedraTest.obj");
}
//***********************************

void GeodesicDijkstraTest()
{
	BGAL::_ManifoldModel model("data\\sphere.obj");
	std::map<int, double> source;
	source[0] = 0;
	BGAL::Geodesic::_Dijkstra dijk(model, source);
	dijk.execute_();
	std::vector<double> distance = dijk.get_distances_();
	double maxd = *(std::max_element(distance.begin(), distance.end()));
	for (int i = 0; i < distance.size(); ++i)
	{
		distance[i] = distance[i] / maxd;
	}
	model.save_scalar_field_obj_file_("data\\GeodesicDijkstraTest.obj", distance);
}
//************************************

void CPDTest()
{
	BGAL::_ManifoldModel model("data\\sphere.obj");
	std::function<double(BGAL::_Point3& p)> rho = [](BGAL::_Point3& p)
	{
		return 1;
	};
	BGAL::_LBFGS::_Parameter para;
	para.is_show = true;
	para.epsilon = 5e-4;
	BGAL::_CPD3D cpd(model, rho, para);
	cpd._omt_eps = 5e-4;
	double sum_mass = 0;
	for (auto fit = model.face_begin(); fit != model.face_end(); ++fit)
	{
		auto tri = model.face_(fit.id());
		sum_mass += tri.area_();
	}
	int num = 200;
	std::vector<double> capacity(num, sum_mass / num);
	cpd.calculate_(capacity);
	const std::vector<BGAL::_Point3>& sites = cpd.get_sites();
	const std::vector<double>& weights = cpd.get_weights();
	const BGAL::_Restricted_Tessellation3D& RPD = cpd.get_RPD();
	const std::vector<std::vector<std::tuple<int, int, int>>>& cells = RPD.get_cells_();
	std::cout << "mass - capacity" << std::endl;
	for (int i = 0; i < num; ++i)
	{
		double cal_mass = 0;
		for (int j = 0; j < cells[i].size(); ++j)
		{
			BGAL::_Triangle3 tri(RPD.vertex_(std::get<0>(cells[i][j])), RPD.vertex_(std::get<1>(cells[i][j])), RPD.vertex_(std::get<2>(cells[i][j])));
			cal_mass += tri.area_();
		}
		std::cout << cal_mass << "\t" << capacity[i] << "\t" << cal_mass - capacity[i] << std::endl;
	}
	std::ofstream out("data\\CPD3DTest.obj");
	out << "g 3D_Object\nmtllib BKLineColorBar.mtl\nusemtl BKLineColorBar" << std::endl;
	for (int i = 0; i < RPD.number_vertices_(); ++i)
	{
		out << "v " << RPD.vertex_(i) << std::endl;
	}
	for (int i = 0; i < cells.size(); ++i)
	{
		double color = (double)BGAL::_BOC::rand_();
		out << "vt " << color << " 0" << std::endl;
		for (int j = 0; j < cells[i].size(); ++j)
		{
			out << "f " << std::get<0>(cells[i][j]) + 1 << "/" << i + 1
				<< " " << std::get<1>(cells[i][j]) + 1 << "/" << i + 1
				<< " " << std::get<2>(cells[i][j]) + 1 << "/" << i + 1 << std::endl;
		}
	}
	out.close();
}

void CVT3DTest()
{
	BGAL::_ManifoldModel model("data\\bunny.obj");
	std::function<double(BGAL::_Point3& p)> rho = [](BGAL::_Point3& p)
	{
		return 1;
	};
	BGAL::_LBFGS::_Parameter para;
	para.is_show = true;
	para.epsilon = 1e-4;
	BGAL::_CVT3D cvt(model, rho, para);
	int num = 300;
	cvt.calculate_(num);
	const std::vector<BGAL::_Point3>& sites = cvt.get_sites();
	const BGAL::_Restricted_Tessellation3D& RVD = cvt.get_RVD();
	const std::vector<std::vector<std::tuple<int, int, int>>>& cells = RVD.get_cells_();
	std::ofstream out("data\\CVT3DTest.obj");
	out << "g 3D_Object\nmtllib BKLineColorBar.mtl\nusemtl BKLineColorBar" << std::endl;
	for (int i = 0; i < RVD.number_vertices_(); ++i)
	{
		out << "v " << RVD.vertex_(i) << std::endl;
	}
	for (int i = 0; i < cells.size(); ++i)
	{
		double color = (double)BGAL::_BOC::rand_();
		out << "vt " << color << " 0" << std::endl;
		for (int j = 0; j < cells[i].size(); ++j)
		{
			out << "f " << std::get<0>(cells[i][j]) + 1 << "/" << i + 1
				<< " " << std::get<1>(cells[i][j]) + 1 << "/" << i + 1
				<< " " << std::get<2>(cells[i][j]) + 1 << "/" << i + 1 << std::endl;
		}
	}
	out.close();
}

/*************************************

Expect:
====================BOCSignTest
1
0
====================LinearSystemTest
2
2
2
====================ALGLIBTest
[-3.00,3.00]
false
false
false
====================LBFGSTest
0       0       0.001   172.8   7465
1       5       0.002   140.8   4956.19
2       6       0.003   3.17764e-14     2.52435e-28
reach the gradient tolerance
1
1
n: 2
====================BaseShapeTest
yes!
20
====================Tessellation2DTest
====================CVTLBFGSTest
0       0       0.007   0.00829834      0.00377827
1       2       0.014   0.00350175      0.0024125
2       3       0.022   0.00190224      0.00207575
3       4       0.03    0.00138204      0.00189172
4       5       0.037   0.001058        0.00179544
5       6       0.044   0.000630967     0.00174395
6       7       0.052   0.000453071     0.00171594
7       8       0.059   0.000407016     0.00169794
8       9       0.067   0.000382601     0.00168672
9       10      0.074   0.000236295     0.00167822
10      11      0.083   0.000188259     0.00167365
11      12      0.089   0.00018715      0.00166988
12      13      0.095   0.000191645     0.00166741
13      14      0.102   0.0001226       0.00166534
14      15      0.11    0.000109166     0.00166371
15      16      0.118   0.00012134      0.00166249
16      17      0.127   0.000142921     0.00166056
17      18      0.135   0.000219643     0.00165917
18      19      0.147   0.00015968      0.00165643
19      20      0.155   0.000126286     0.0016543
20      21      0.161   0.000124011     0.00165273
21      22      0.171   0.000113692     0.00165165
22      23      0.18    8.62573e-05     0.00165076
reach the gradient tolerance
====================DrawTest
====================IntegralTest
0.5
====================ModelTest
V number: 642
F number: 1280
====================Tessellation3DTest
====================KDTreeTest
61 0.0051334350000000004283 -0.95469340000000002533 0.18068470000000000364
0.16776960319175260317    0.16776960319175260317
my:  61 0.16776960319175260317
====================ICPTest
    0.92387953251128940302     0.38268343236508994831  6.7307270867900115263e-16     -0.8156403124820496009
   -0.38268343236509022587     0.92387953251128607235 -5.2041704279304212832e-17      0.1213708393898148552
 1.3444106938820254982e-15  1.0451708942760262744e-16      1.0000000000000013323       0.149999999999999023
                         0                          0                          0                          1
     1.0000000000000026645  4.4408920985006261617e-16   6.417535974601941462e-16 -2.2204460492503130808e-15
 5.5511151231257827021e-16     0.99999999999999933387  2.0949350896789406863e-16 -4.7184478546569152968e-16
 1.3444106938820254982e-15  1.0451708942760262744e-16      1.0000000000000013323 -9.7144514654701197287e-16
                         0                          0                          0                          1
====================MarchingTetrahedraTest
====================GeodesicDijkstraTest
====================CPDTest
0       0       2.8559999999999998721   0.97077838761131618472  1.5745422097489745195
1       1       3.7879999999999998117   0.50654419345871326552  0.83706428321869918996
2       2       4.3559999999999998721   0.10341526695700922756  0.54823767814940860266
3       3       4.9279999999999999361   0.064348966727981085634 0.53241666390643571649
4       4       5.4930000000000003268   0.055978077579559790133 0.52303719647387347802
5       5       6.0590000000000001634   0.027350443137352131728 0.51869437610611468514
6       6       6.6159999999999996589   0.025028562315720939702 0.51602315877088045237
7       7       7.1310000000000002274   0.041603689229046926512 0.51276448733047774731
8       8       7.6490000000000000213   0.028852734948012288135 0.51009318178425633317
9       9       8.1630000000000002558   0.018929838601576352147 0.50782103206131068429
10      10      8.6669999999999998153   0.016154404661631139445 0.50680290391269833261
11      11      9.1940000000000008384   0.01314021086265272989  0.50611520612709437472
12      12      9.7010000000000005116   0.012379146785242668358 0.50538449212462099869
13      13      10.180999999999999162   0.011617886072246796231 0.50490602085507330088
14      14      10.675000000000000711   0.0083309046909829185396        0.50453526747561872057
15      15      11.137999999999999901   0.0057183422454947386432        0.50432714568149383805
16      16      11.617000000000000881   0.0040848813100232556073        0.50421835653346636086
17      17      12.096000000000000085   0.0036973204459675658613        0.50412838219085553959
18      18      12.564999999999999503   0.0053717282361031475427        0.50404070920698462732
19      19      13.041999999999999815   0.0050809776275175227295        0.50392744977097758685
20      21      14.038000000000000256   0.0061808731807658189028        0.50381777140314776275
21      22      14.516000000000000014   0.0080891448900077683737        0.50370077706668570094
22      23      15.026999999999999247   0.0058684635027048412739        0.50356297222700674432
23      24      15.489000000000000767   0.0032095819943958178375        0.50348377130270893787
24      25      15.942999999999999616   0.0038988783924685118353        0.50342742924118755177
25      26      16.40899999999999892    0.0045883746557730644214        0.50337278805180996066
26      27      16.885999999999999233   0.0039180574454319282846        0.50331414136301411144
27      28      17.349000000000000199   0.0027244578038674569821        0.50326446135064273335
28      29      17.821000000000001506   0.0019379409723720758454        0.5032354190101470115
29      30      18.329000000000000625   0.0018087757414846634789        0.50321983470384168413
30      31      18.812000000000001165   0.0013028630258690508496        0.50320715636633683854
31      32      19.321000000000001506   0.0011783374032515757986        0.50319915304635931541
32      33      19.818999999999999062   0.0014676454032900041   0.50319086790881939475
33      34      20.353999999999999204   0.0015583810724837368535        0.50318094734710294702
34      35      20.833999999999999631   0.002362075227291897915 0.50316470947673608283
35      36      21.309999999999998721   0.0032787689329946251814        0.50314831577231877713
36      37      21.803999999999998494   0.0026986173532612236017        0.50312304385060069301
37      38      22.262000000000000455   0.0014219395958781934613        0.50310648530927826183
38      39      22.739000000000000767   0.001615130064985079595 0.50309781280288112804
39      40      23.193000000000001393   0.0010293988295718848724        0.50309054561560839769
40      41      23.653999999999999915   0.001150556612973379798 0.50308201595055312971
41      42      24.176999999999999602   0.0018894256762530183425        0.50306894262166490517
42      43      24.685999999999999943   0.0032389060252479831212        0.50303698963606058303
43      46      26.155999999999998806   0.0041089762752991101924        0.50301811622742298447
44      47      26.615999999999999659   0.0052529712037364022573        0.50298995415369451845
45      48      27.082000000000000739   0.0052228035541843978451        0.50291050212833288136
46      49      27.588000000000000966   0.0037558136639685708695        0.50283448424049492775
47      50      28.059000000000001052   0.0049192653062895708854        0.50274012469509332668
48      51      28.556000000000000938   0.010434469543566506078 0.5025991741619602049
49      52      29.042999999999999261   0.0094345597461962023983        0.50242554455932364466
50      53      29.55099999999999838    0.0054597914390715086147        0.50221842072561184711
51      54      30.047999999999998266   0.0059490895457608421876        0.50209167430241052887
52      55      30.553999999999998494   0.0063564344407932042366        0.50200891702364192071
53      56      31.030000000000001137   0.0043077053720921012689        0.50192435069084972987
54      57      31.545000000000001705   0.0027823923960697792557        0.50187037839039727594
55      58      32.029000000000003467   0.0019360024677122155638        0.50183171328020392821
56      59      32.529000000000003467   0.0033217704188521347645        0.50182773657329549089
57      60      33.045999999999999375   0.0010214169540488850663        0.50181491924833965257
58      61      33.567000000000000171   0.00064919895087481641321       0.50181214372791038691
59      62      34.08100000000000307    0.00055139179941391021095       0.50180937461070618255
60      63      34.624000000000002331   0.00078619052963909937578       0.50180825059119238407
61      64      35.082999999999998408   0.00027829333736625617426       0.5018074726296616328
62      65      35.61500000000000199    0.00016581815728628954663       0.50180731065158645787
63      66      36.06799999999999784    0.00012529229361561805794       0.50180715577438583797
64      67      36.533000000000001251   0.0001881051793069788672        0.50180708027468379218
65      68      37      7.1675672068046648237e-05       0.50180702933455167969
reach the gradient tolerance
mass - capacity
0.25012830819059872489  0.25012984804989535359  -1.5398592966286983597e-06
0.25013019810614184335  0.25012984804989535359  3.5005624648976052526e-07
0.25013058054167031097  0.25012984804989535359  7.3249177495737782806e-07
0.25012839111854828777  0.25012984804989535359  -1.4569313470658151743e-06
0.25012853615062691226  0.25012984804989535359  -1.3118992684413299799e-06
0.25013060255317565161  0.25012984804989535359  7.5450328029802449237e-07
0.25013039914674778386  0.25012984804989535359  5.5109685243026618195e-07
0.25012847467613846808  0.25012984804989535359  -1.3733737568855097777e-06
0.25013025518694881333  0.25012984804989535359  4.0713705345973849603e-07
0.25013061992165341874  0.25012984804989535359  7.7187175806514574106e-07
0.25013014735038452407  0.25012984804989535359  2.9930048917048424073e-07
0.25013043153240493988  0.25012984804989535359  5.8348250958628611329e-07
0.2501283195836258022   0.25012984804989535359  -1.5284662695513873132e-06
0.25012961148541623668  0.25012984804989535359  -2.3656447911690747787e-07
0.25013030405732905592  0.25012984804989535359  4.5600743370233232099e-07
0.25013015638657509765  0.25012984804989535359  3.0833667974405898349e-07
0.25012830901375759929  0.25012984804989535359  -1.5390361377543015919e-06
0.25013049676849052894  0.25012984804989535359  6.4871859517534602446e-07
0.25012961153793605851  0.25012984804989535359  -2.3651195929508261884e-07
0.25013055535693423659  0.25012984804989535359  7.0730703888299828463e-07
0.25013036070870031669  0.25012984804989535359  5.1265880496309534919e-07
0.25012834730090183211  0.25012984804989535359  -1.5007489935214834986e-06
0.25012812541727147408  0.25012984804989535359  -1.7226326238795053314e-06
0.25013048852504093933  0.25012984804989535359  6.4047514558573936938e-07
0.25013062244000516809  0.25012984804989535359  7.7439010981450451254e-07
0.25013009594294205451  0.25012984804989535359  2.4789304670091993898e-07
0.25013029708295431153  0.25012984804989535359  4.4903305895793721447e-07
0.25013016985467523279  0.25012984804989535359  3.2180477987919786642e-07
0.25012833494101649467  0.25012984804989535359  -1.5131088788589153182e-06
0.25013059658552866393  0.25012984804989535359  7.4853563331034322914e-07
0.25012839803272507444  0.25012984804989535359  -1.4500171702791497808e-06
0.25013017163012218891  0.25012984804989535359  3.2358022683531828534e-07
0.25013017696256534261  0.25012984804989535359  3.2891266998902324303e-07
0.25013038689742878029  0.25012984804989535359  5.3884753342670066445e-07
0.25013020554906495452  0.25012984804989535359  3.5749916960092775753e-07
0.25013012258386424502  0.25012984804989535359  2.7453396889143277804e-07
0.25013061224202304267  0.25012984804989535359  7.6419212768907840427e-07
0.25013031715688155421  0.25012984804989535359  4.6910698620061808128e-07
0.25013038684554850244  0.25012984804989535359  5.3879565314884914073e-07
0.2501304064087886414   0.25012984804989535359  5.5835889328781362906e-07
0.25013013896463942576  0.25012984804989535359  2.9091474407216821874e-07
0.25012840413401904449  0.25012984804989535359  -1.4439158763090986781e-06
0.25013006787565189581  0.25012984804989535359  2.1982575654222458184e-07
0.25012825838974489523  0.25012984804989535359  -1.5896601504583607323e-06
0.25013037638930035733  0.25012984804989535359  5.2833940500374154681e-07
0.25013031479326036655  0.25012984804989535359  4.6674336501295599078e-07
0.25013019806972219827  0.25012984804989535359  3.5001982684468302409e-07
0.25013060746784665511  0.25012984804989535359  7.5941795130152200954e-07
0.25013048499741385999  0.25012984804989535359  6.3694751850640329849e-07
0.25013061963995819603  0.25012984804989535359  7.7159006284244213703e-07
successful!

*******************************************/


int alltest()
{
	//std::cout << "====================GraphCutsTest" << std::endl;
	//GraphCutsTest();
	std::cout << "====================BOCSignTest" << std::endl;
	BOCSignTest();	
	std::cout << "====================LinearSystemTest" << std::endl;
	LinearSystemTest();
	std::cout << "====================ALGLIBTest" << std::endl;
	ALGLIBTest();
	std::cout << "====================LBFGSTest" << std::endl;
	LBFGSTest();
	std::cout << "====================BaseShapeTest" << std::endl;
	BaseShapeTest();
	std::cout << "====================Tessellation2DTest" << std::endl;
	Tessellation2DTest();
	std::cout << "====================CVTLBFGSTest" << std::endl;
	CVTLBFGSTest();
	std::cout << "====================DrawTest" << std::endl;
	DrawTest();
	std::cout << "====================IntegralTest" << std::endl;
	IntegralTest();
	std::cout << "====================ModelTest" << std::endl;
	ModelTest();
	std::cout << "====================Tessellation3DTest" << std::endl;
	Tessellation3DTest();
	//std::cout << "====================ReadFileTest" << std::endl;
	//ReadFileTest();
	std::cout << "====================KDTreeTest" << std::endl;
	KDTreeTest();
	std::cout << "====================ICPTest" << std::endl;
	ICPTest();
	std::cout << "====================MarchingTetrahedraTest" << std::endl;
	MarchingTetrahedraTest();
	std::cout << "====================GeodesicDijkstraTest" << std::endl;
	GeodesicDijkstraTest();
	std::cout << "====================CPDTest" << std::endl;
	CPDTest();
	std::cout << "====================CVT3DTest" << std::endl;
	CVT3DTest();
	std::cout << "successful!" << std::endl;
	return 0;
}

int main()
{
	alltest();
	return 0;
}