#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "BGAL/CVTLike/CVT.h"
#include "BGAL/Algorithm/BOC/BOC.h"
#include "BGAL/Integral/Integral.h"
#include "BGAL/Optimization/LinearSystem/LinearSystem.h"
#include <omp.h> 
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/IO/OBJ.h>
typedef CGAL::Simple_cartesian<double> K_T;
typedef K_T::FT FT;
typedef K_T::Point_3 Point_T;

typedef K_T::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K_T> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K_T, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef K_T::Plane_3                                     Plane;

namespace BGAL
{
	_CVT3D::_CVT3D(const _ManifoldModel& model) : _model(model), _RVD(model), _para()
	{
		_rho = [](BGAL::_Point3& p)
		{
			return 1;
		};
		_para.is_show = true;
		_para.epsilon = 5e-5;
	}
	_CVT3D::_CVT3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para) : _model(model), _RVD(model), _rho(rho), _para(para)
	{
		
	}
	_CVT3D::_CVT3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para, std::string modelname) : _model(model), _RVD(model), _rho(rho), _para(para)
	{
		_modelname = modelname;
	}
	void _CVT3D::calculate_(int num_sites)
	{		
		int num = num_sites;
		_sites.resize(num);
		for (int i = 0; i < num; ++i)
		{
			int fid = rand() % _model.number_faces_();
			double l0, l1, l2, sum;
			l0 = _BOC::rand_();
			l1 = _BOC::rand_();
			l2 = _BOC::rand_();
			sum = l0 + l1 + l2;
			l0 /= sum;
			l1 /= sum;
			l2 /= sum;
			_sites[i] = _model.face_(fid).point(0) * l0 + _model.face_(fid).point(1) * l1 + _model.face_(fid).point(2) * l2;
		}
		_RVD.calculate_(_sites);
		std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
			= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
		{
			for (int i = 0; i < num; ++i)
			{
				BGAL::_Point3 p(X(i * 3), X(i * 3 + 1), X(i * 3 + 2));
				// 这里是不是需要将点投影到mesh表面有待考虑
				_sites[i] = p;
			}
			_RVD.calculate_(_sites);
			const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RVD.get_cells_();
			double energy = 0;
			g.setZero();
			for (int i = 0; i < num; ++i)
			{
				for (int j = 0; j < cells[i].size(); ++j)
				{
					Eigen::VectorXd inte = BGAL::_Integral::integral_triangle3D(
						[&](BGAL::_Point3 p)
						{
							Eigen::VectorXd r(5);
							r(0) = _rho(p);
							r(1) = _rho(p) * ((_sites[i] - p).sqlength_());
							r(2) = 2 * _rho(p) * (_sites[i].x() - p.x());
							r(3) = 2 * _rho(p) * (_sites[i].y() - p.y());
							r(4) = 2 * _rho(p) * (_sites[i].z() - p.z());
							return r;
						}, _RVD.vertex_(std::get<0>(cells[i][j])), _RVD.vertex_(std::get<1>(cells[i][j])), _RVD.vertex_(std::get<2>(cells[i][j]))
							);
					energy += inte(1);
					g(i * 3) += inte(2);
					g(i * 3 + 1) += inte(3);
					g(i * 3 + 2) += inte(4);
				}
			}
			return energy;
		};
		BGAL::_LBFGS lbfgs(_para);
		Eigen::VectorXd iterX(num * 3);
		for (int i = 0; i < num; ++i)
		{
			iterX(i * 3) = _sites[i].x();
			iterX(i * 3 + 1) = _sites[i].y();
			iterX(i * 3 + 2) = _sites[i].z();
		}
		lbfgs.minimize(fg, iterX);
		for (int i = 0; i < num; ++i)
		{
			_sites[i] = BGAL::_Point3(iterX(i * 3), iterX(i * 3 + 1), iterX(i * 3 + 2));
		}
		_RVD.calculate_(_sites);
	}

	
	void _CVT3D::calculate_CapVT(std::vector<BGAL::_Point3>& sites)
	{
		Polyhedron polyhedron;
		string modelname = _modelname;
		std::ifstream input("..\\..\\data\\" + modelname + ".off");
		input >> polyhedron;
		input.close();
		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		
		int num = sites.size();
		_sites = sites;
		_RVD.calculate_(_sites);
		const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RVD.get_cells_();
		double TotArea = 0;
		for (int i = 0; i < cells.size(); ++i)
		{
			for (int j = 0; j < cells[i].size(); ++j)
			{
				double side[3];//存储三条边的长度;
				auto a = _RVD.vertex_(std::get<0>(cells[i][j]));
				auto b = _RVD.vertex_(std::get<1>(cells[i][j]));
				auto c = _RVD.vertex_(std::get<2>(cells[i][j]));

				side[0] = sqrt(pow(a.x() - b.x(), 2) + pow(a.y() - b.y(), 2) + pow(a.z() - b.z(), 2));
				side[1] = sqrt(pow(a.x() - c.x(), 2) + pow(a.y() - c.y(), 2) + pow(a.z() - c.z(), 2));
				side[2] = sqrt(pow(c.x() - b.x(), 2) + pow(c.y() - b.y(), 2) + pow(c.z() - b.z(), 2));
				double p = (side[0] + side[1] + side[2]) / 2;
				double area = sqrt(p * (p - side[0]) * (p - side[1]) * (p - side[2]));
				TotArea += area;
			}
		}
		

		std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
			= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
		{

			for (int i = 0; i < num; ++i)
			{
				
				// 这里是不是需要将点投影到mesh表面有待考虑
				//_model
				Point_T query(X(i * 3), X(i * 3 + 1), X(i * 3 + 2));
				Point_T closest = tree.closest_point(query);
				
				BGAL::_Point3 p(closest.x(), closest.y(), closest.z());
				//BGAL::_Point3 p(X(i * 3), X(i * 3 + 1), X(i * 3 + 2));
				_sites[i] = p;
			}
			_RVD.calculate_(_sites);
			const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RVD.get_cells_();
			double energy = 0;
			g.setZero();


			map<int,double> CellAreas;

			double AreaDiff = 0;

			for (int i = 0; i < cells.size(); ++i)
			{
				double CellsArea = 0;
				for (int j = 0; j < cells[i].size(); ++j)
				{
					double side[3];//存储三条边的长度;
					auto a = _RVD.vertex_(std::get<0>(cells[i][j]));
					auto b = _RVD.vertex_(std::get<1>(cells[i][j]));
					auto c = _RVD.vertex_(std::get<2>(cells[i][j]));
					side[0] = sqrt(pow(a.x() - b.x(), 2) + pow(a.y() - b.y(), 2) + pow(a.z() - b.z(), 2));
					side[1] = sqrt(pow(a.x() - c.x(), 2) + pow(a.y() - c.y(), 2) + pow(a.z() - c.z(), 2));
					side[2] = sqrt(pow(c.x() - b.x(), 2) + pow(c.y() - b.y(), 2) + pow(c.z() - b.z(), 2));
					double p = (side[0] + side[1] + side[2]) / 2;
					double area = sqrt(p * (p - side[0]) * (p - side[1]) * (p - side[2]));
					CellsArea += area;
				}
				CellAreas[i]=CellsArea;
				AreaDiff += (CellsArea - TotArea / num) * (CellsArea - TotArea / num);
			}
			
			energy = AreaDiff;
			
			auto Edges = _RVD.get_edges_();

			for (int i = 0; i < num; i++)
			{
				double sumx = 0.0, sumy = 0.0, sumz = 0.0;
				for (auto ee : Edges[i])
				{
					int j = ee.first;
					map<int,int> m;
					for (auto e : ee.second)
					{
						if (m.find(e.first) == m.end())
						{
							m[e.first] = 0;
						}
						if (m.find(e.second) == m.end())
						{
							m[e.second] = 0;
						}
						m[e.first]++;
						m[e.second]++;
					}
					vector<int> EdgPts;
					for (auto e : m)
					{
						if (e.second == 1)
						{
							EdgPts.push_back(e.first);
						}
					}
					if (EdgPts.size() != 2)
					{
						cout << "Error!!" << endl;
					}
					
					double tmp = (CellAreas[i]  - CellAreas[j]) / (_sites[j] - _sites[i]).length_();
					auto Epts1 = _RVD.vertex_(EdgPts[0]);
					auto Epts2 = _RVD.vertex_(EdgPts[1]);
					double l = (Epts1 - Epts2).length_();
					auto vecx = l * ((Epts1.x() + Epts2.x()) / 2.0 - _sites[i].x());
					auto vecy = l * ((Epts1.y() + Epts2.y()) / 2.0 - _sites[i].y());
					auto vecz = l * ((Epts1.z() + Epts2.z()) / 2.0 - _sites[i].z());
					sumx += tmp*vecx;
					sumy += tmp*vecy;
					sumz += tmp*vecz;
				}
				
				
				g(i * 3) += 2.0*sumx;
				g(i * 3 + 1) += 2.0*sumy;
				g(i * 3 + 2) += 2.0*sumz;
			}
			bool projection = 1;
			if (projection)
			{

				for (int i = 0; i < num; i++)
				{
					Point_T query(_sites[i].x(), _sites[i].y(), _sites[i].z());
					Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
					Point_T closest = pp.first;
					Polyhedron::Face_handle f = pp.second;

					Plane pl(f->halfedge()->vertex()->point(), f->halfedge()->next()->vertex()->point(), f->halfedge()->next()->next()->vertex()->point());
					auto nn = pl.orthogonal_vector();
					Eigen::Vector3d n(nn.x(), nn.y(), nn.z());
					Eigen::Vector3d u(g(i * 3), g(i * 3 + 1), g(i * 3 + 2));
					Eigen::Vector3d newG = u - n * ((u.dot(n)) / n.squaredNorm());
					g(i * 3) = newG.x();
					g(i * 3 + 1) = newG.y();
					g(i * 3 + 2) = newG.z();
				}
			}
			

			
			return energy;
		};
		auto _para2 = _para;
		//_para2.epsilon = 1e-4;
		_para2.max_iteration = 100;
		_para2.max_linearsearch = 10;
		//_para2.is_show = false;
		BGAL::_LBFGS lbfgs(_para2);
		Eigen::VectorXd iterX(num * 3);
		for (int i = 0; i < num; ++i)
		{
			iterX(i * 3) = _sites[i].x();
			iterX(i * 3 + 1) = _sites[i].y();
			iterX(i * 3 + 2) = _sites[i].z();
		}
		lbfgs.minimize(fg, iterX);
		for (int i = 0; i < num; ++i)
		{
			Point_T query(iterX(i * 3), iterX(i * 3 + 1), iterX(i * 3 + 2));
			Point_T closest = tree.closest_point(query);

			BGAL::_Point3 p(closest.x(), closest.y(), closest.z());
			//BGAL::_Point3 p(X(i * 3), X(i * 3 + 1), X(i * 3 + 2));
			_sites[i] = p;
			//_sites[i] = BGAL::_Point3(iterX(i * 3), iterX(i * 3 + 1), iterX(i * 3 + 2));


			
		}
		_RVD.calculate_(_sites);

		double AreaDiff = 0;
		
		for (int i = 0; i < cells.size(); ++i)
		{
			double CellsArea = 0;
			for (int j = 0; j < cells[i].size(); ++j)
			{
				double side[3];//存储三条边的长度;
				auto a = _RVD.vertex_(std::get<0>(cells[i][j]));
				auto b = _RVD.vertex_(std::get<1>(cells[i][j]));
				auto c = _RVD.vertex_(std::get<2>(cells[i][j]));
				side[0] = sqrt(pow(a.x() - b.x(), 2) + pow(a.y() - b.y(), 2) + pow(a.z() - b.z(), 2));
				side[1] = sqrt(pow(a.x() - c.x(), 2) + pow(a.y() - c.y(), 2) + pow(a.z() - c.z(), 2));
				side[2] = sqrt(pow(c.x() - b.x(), 2) + pow(c.y() - b.y(), 2) + pow(c.z() - b.z(), 2));
				double p = (side[0] + side[1] + side[2]) / 2;
				double area = sqrt(p * (p - side[0]) * (p - side[1]) * (p - side[2]));
				CellsArea += area;
			}
			AreaDiff += (CellsArea - TotArea / num) * (CellsArea - TotArea / num);
		}
		cout<<"CapVT AreaDiff: "<<AreaDiff<<endl;
		
		
	}

	void _CVT3D::calculate_CapCVTByGD(std::vector<BGAL::_Point3>& sites)
	{
		int num = sites.size();
		std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fgg
			= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
		{
			//for (int i = 0; i < num; ++i)
			//{
			//	BGAL::_Point3 p(X(i * 3), X(i * 3 + 1), X(i * 3 + 2));
			//	// 这里是不是需要将点投影到mesh表面有待考虑
			//	_sites[i] = p;
			//}
			_RVD.calculate_(_sites);
			const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RVD.get_cells_();
			double energy = 0;
			g.setZero();
			for (int i = 0; i < num; ++i)
			{
				for (int j = 0; j < cells[i].size(); ++j)
				{
					Eigen::VectorXd inte = BGAL::_Integral::integral_triangle3D(
						[&](BGAL::_Point3 p)
						{
							Eigen::VectorXd r(5);
							r(0) = _rho(p);
							r(1) = _rho(p) * ((_sites[i] - p).sqlength_());
							r(2) = 2 * _rho(p) * (_sites[i].x() - p.x());
							r(3) = 2 * _rho(p) * (_sites[i].y() - p.y());
							r(4) = 2 * _rho(p) * (_sites[i].z() - p.z());
							return r;
						}, _RVD.vertex_(std::get<0>(cells[i][j])), _RVD.vertex_(std::get<1>(cells[i][j])), _RVD.vertex_(std::get<2>(cells[i][j]))
							);
					energy += inte(1);
					g(i * 3) += inte(2);
					g(i * 3 + 1) += inte(3);
					g(i * 3 + 2) += inte(4);
				}
			}
			return energy;
		};
		
		//GD
		double learning_rate = 100;
		int maxiter = 30;
		int lastEng = 99999.0;

		auto BestSites = _sites;
		double MinLoss = 999999999.0;
		int NoChange = 0;
		while (1)
		{
			_para.epsilon = 1e-5;
			calculate_CapVT(_sites);
			
			Eigen::VectorXd iterX(num * 3);
			for (int i = 0; i < num; ++i)
			{
				iterX(i * 3) = _sites[i].x();
				iterX(i * 3 + 1) = _sites[i].y();
				iterX(i * 3 + 2) = _sites[i].z();
			}
			Eigen::VectorXd g(num * 3);
			double energy = fgg(iterX, g);
			cout << "  Energy before: " << energy ;
			if (energy < MinLoss)
			{
				BestSites = _sites;
				NoChange = 0;
			}
			else
			{
				NoChange++;
				if (NoChange > 5)
				{
					break;
				}
			}
			
			
			
			//if( energy < lastEng )
			//{
			//	//learning_rate = learning_rate * 1.1;
			//}
			//else
			//{
			//	learning_rate = learning_rate * 0.8;
			//}

			lastEng = energy;
			Eigen::VectorXd dX = -learning_rate * g;
			for (int i = 0; i < num; ++i)
			{
				iterX(i * 3) += dX(i * 3);
				iterX(i * 3 + 1) += dX(i * 3 + 1);
				iterX(i * 3 + 2) += dX(i * 3 + 2);
				_sites[i] = BGAL::_Point3(iterX(i * 3), iterX(i * 3 + 1), iterX(i * 3 + 2));
			}
			cout << "  Energy after: " << fgg(iterX, g) << endl;
			
			if (learning_rate < 1e-3)
			{
				break;
			}

			learning_rate *= 0.90;
			if (maxiter-- == 0)
			{
				break;
			}
			
		}

		_sites = BestSites;
		_para.epsilon = 1e-5;
		//calculate_CapVT(_sites);
		Eigen::VectorXd iterX(num * 3);
		for (int i = 0; i < num; ++i)
		{
			iterX(i * 3) = _sites[i].x();
			iterX(i * 3 + 1) = _sites[i].y();
			iterX(i * 3 + 2) = _sites[i].z();
		}
		Eigen::VectorXd g(num * 3);
		double energy = fgg(iterX, g);
		cout << "\nCVT  Energy Final: " << energy << endl;
	}


} // namespace BGAL
