#pragma once
#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
#include "BGAL/BaseShape/Triangle.h"
#include "BGAL/BaseShape/Line.h"
#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"
#include "BGAL/Tessellation3D/Tessellation3D.h"
#include "BGAL/Optimization/LBFGS/LBFGS.h"

namespace BGAL
{
	class _CPD3D
	{
	public:
		_CPD3D(const _ManifoldModel& model);
		_CPD3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para);
		void calculate_(const std::vector<double>& capacity, std::vector<_Point3> sites);
		const std::vector<_Point3>& get_sites() const
		{
			return _sites;
		}
		const std::vector<double>& get_weights() const
		{
			return _weights;
		}
		const _Restricted_Tessellation3D& get_RPD() const
		{
			return _RPD;
		}
	public:
		// 这些参数应该封装一下的，暂时没时间了，直接暴露出来
		const _ManifoldModel& _model;
		_Restricted_Tessellation3D _RPD;
		std::vector<double> _capacity{};
		std::vector<_Point3> _sites{};
		std::vector<double> _weights{};
		std::function<double(_Point3& p)> _rho;
		_LBFGS::_Parameter _para;
		int _max_count; // 牛顿法最大次数
		double _omt_eps; // 算OMT时停止条件
		double _pinvtoler; // 牛顿法求解线性方程组来求方向时的精度
		double _hessian_eps; // 海森矩阵是半正定的，保险起见加了epsilon
	};
} // namespace BGAL
