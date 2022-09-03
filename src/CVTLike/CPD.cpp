#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "BGAL/CVTLike/CPD.h"
#include "BGAL/Algorithm/BOC/BOC.h"
#include "BGAL/Integral/Integral.h"
#include "BGAL/Optimization/LinearSystem/LinearSystem.h"

namespace BGAL
{
	_CPD3D::_CPD3D(const _ManifoldModel& model) : _model(model), _RPD(model), _para()
	{
		_rho = [](BGAL::_Point3& p)
		{
			return 1;
		};
		_para.is_show = true;
		_para.epsilon = 5e-5;
		_max_count = 50;
		_omt_eps = 1e-4;
		_pinvtoler = 1e-6;
		_hessian_eps = 1e-10;
	}
	_CPD3D::_CPD3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para) : _model(model), _RPD(model), _rho(rho), _para(para)
	{
		_max_count = 50;
		_omt_eps = 1e-4;
		_pinvtoler = 1e-6;
		_hessian_eps = 1e-10;
	}
	void _CPD3D::calculate_(const std::vector<double>& capacity, std::vector<_Point3> sites)
	{

		_max_count = 50;
		_omt_eps = 1e-4;
		_pinvtoler = 1e-4;
		_hessian_eps = 1e-5;
		_para.max_linearsearch = 20;

		_capacity = capacity;
		int num = _capacity.size();
		/*_sites.resize(num);
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
		}*/
		_sites = sites;
		_weights.resize(num, 0);
		_RPD.calculate_(_sites, _weights);
		std::function<bool(Eigen::SparseMatrix<double>& h)> cal_h
			= [&](Eigen::SparseMatrix<double>& h)
		{
			const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RPD.get_cells_();
			const std::vector<std::map<int, std::vector<std::pair<int, int>>>>& edges = _RPD.get_edges_();
			vector<Eigen::Triplet<double>> trilist;
			for (int i = 0; i < num; ++i)
			{
				double hii = 0;
				for (auto& kv : edges[i])
				{
					// 这一段临时写的，需要仔细打磨
					double len = 0;
					std::set<int> temp_p;
					for (auto& te : kv.second)
					{
						len += (_RPD.vertex_(te.first) - _RPD.vertex_(te.second)).length_();
						if (temp_p.find(te.first) == temp_p.end())
						{
							temp_p.insert(te.first);
						}
						else
						{
							temp_p.erase(te.first);
						}
						if (temp_p.find(te.second) == temp_p.end())
						{
							temp_p.insert(te.second);
						}
						else
						{
							temp_p.erase(te.second);
						}
					}
					std::vector<int> ps;
					for (auto& tp : temp_p)
					{
						ps.push_back(tp);
					}
					BGAL::_Point3 mid_p;
					// 对于退化情况，主要是薄板问题，特殊处理下，主要是代码能跑通，逻辑其实不对
					if (ps.size() < 2)
					{
						mid_p = (_sites[i] + _sites[kv.first]) * 0.5;
					}
					else
					{
						const BGAL::_Point3& sp = _RPD.vertex_(ps[0]);
						const BGAL::_Point3& tp = _RPD.vertex_(ps[1]);
						mid_p = (_RPD.vertex_(ps[0]) + _RPD.vertex_(ps[1])) * 0.5;
					}					
					double hij = -_rho(mid_p) * (len) * 0.5 / ((_sites[i] - _sites[kv.first]).length_());
					trilist.push_back(Eigen::Triplet<double>(i, kv.first, hij));
					hii -= hij;
				}
				trilist.push_back(Eigen::Triplet<double>(i, i, hii + _hessian_eps));
			}
			h.resize(num, num);
			h.setFromTriplets(trilist.begin(), trilist.end());
			return true;
		};
		std::function<bool(const Eigen::VectorXd& x, double& e, Eigen::VectorXd& g)> omt_fg
			= [&](const Eigen::VectorXd& x, double& e, Eigen::VectorXd& g)
		{
			e = 0;
			g.setZero();
			for (int i = 0; i < num; ++i)
			{
				_weights[i] = x(i);
			}
			_RPD.calculate_(_sites, _weights);
			if (_RPD.number_hidden_point_() > 0)
			{
				return false;
			}
			const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RPD.get_cells_();
			for (int i = 0; i < num; ++i)
			{
				for (int j = 0; j < cells[i].size(); ++j)
				{
					Eigen::VectorXd inte = BGAL::_Integral::integral_triangle3D(
						[&](BGAL::_Point3 p)
						{
							Eigen::VectorXd r(2);
							r(0) = _rho(p);
							r(1) = _rho(p) * ((_sites[i] - p).sqlength_());
							return r;
						}, _RPD.vertex_(std::get<0>(cells[i][j])), _RPD.vertex_(std::get<1>(cells[i][j])), _RPD.vertex_(std::get<2>(cells[i][j]))
							);
					e -= inte(1) - _weights[i] * inte(0);
					g(i) += inte(0);
				}
				g(i) -= _capacity[i];
				e -= _weights[i] * _capacity[i];
			}
			return true;
		};

		std::function<void()> update_w
			= [&]()
		{
			_weights.clear();
			_weights.resize(num, 0);
			Eigen::VectorXd iterW(num);
			iterW.setZero();
			for (int i = 0; i < num; ++i)
			{
				iterW(i) = _weights[i];
			}
			int count = 0;
			double e;
			Eigen::VectorXd g(num);
			omt_fg(iterW, e, g);
			while (1)
			{				
				
				if (g.norm() < _omt_eps)
					break;
				Eigen::SparseMatrix<double> hess(num, num);
				cal_h(hess);
				Eigen::VectorXd d = -BGAL::_LinearSystem::solve_ldlt(hess, g, _pinvtoler);
				double lambda = 1;
				double newe = e;
				while ((!omt_fg(iterW + lambda * d, newe, g)) || newe > e)
				{
					lambda *= 0.5;
				}
				e = newe;
				iterW = iterW + lambda * d;
				count++;
				if (count > _max_count)
					break;
			}
		};

		std::function<double(const Eigen::VectorXd& X, Eigen::VectorXd& g)> fg
			= [&](const Eigen::VectorXd& X, Eigen::VectorXd& g)
		{
			for (int i = 0; i < num; ++i)
			{
				BGAL::_Point3 p(X(i * 3), X(i * 3 + 1), X(i * 3 + 2));
				// 这里是不是需要将点投影到mesh表面有待考虑
				_sites[i] = p;
			}
			update_w();

			const std::vector<std::vector<std::tuple<int, int, int>>>& cells = _RPD.get_cells_();
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
						}, _RPD.vertex_(std::get<0>(cells[i][j])), _RPD.vertex_(std::get<1>(cells[i][j])), _RPD.vertex_(std::get<2>(cells[i][j]))
							);
					energy += inte(1) - _weights[i] * inte(0);
					g(i * 3) += inte(2);
					g(i * 3 + 1) += inte(3);
					g(i * 3 + 2) += inte(4);
				}
				energy += _weights[i] * _capacity[i];
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
		update_w();
		_RPD.calculate_(_sites, _weights);
	}
} // namespace BGAL
