#pragma once
#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
#include "BGAL/BaseShape/Triangle.h"
#include "BGAL/BaseShape/Line.h"
#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"
#include "BGAL/Tessellation3D/Tessellation3D.h"
#include "BGAL/Optimization/LBFGS/LBFGS.h"


struct MyFaceCVT
{
	MyFaceCVT(Eigen::Vector3i a)
	{
		p = a;
	}
	MyFaceCVT(int a, int b, int c)
	{
		p.x() = a;
		p.y() = b;
		p.z() = c;
	}
	Eigen::Vector3i p;
	bool operator<(const MyFaceCVT& a) const
	{
		if (p.x() == a.p.x())
		{
			if (p.y() == a.p.y())
			{
				return p.z() > a.p.z();
			}
			return p.y() > a.p.y();
		}
		return p.x() > a.p.x();
	}
};
namespace BGAL
{
	class _CVT3D
	{
	public:
		_CVT3D(const _ManifoldModel& model);
		_CVT3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para);
		_CVT3D(const _ManifoldModel& model, std::function<double(_Point3& p)>& rho, _LBFGS::_Parameter para, std::string modelname);
		void calculate_(int site_num);
		void calculate_CapVT(std::vector<BGAL::_Point3>& sites);
		void _CVT3D::calculate_CapCVTByGD(std::vector<BGAL::_Point3>& sites);
		const std::vector<_Point3>& get_sites() const
		{
			return _sites;
		}
		const _Restricted_Tessellation3D& get_RVD() const
		{
			return _RVD;
		}
	public:
		const _ManifoldModel& _model;
		_Restricted_Tessellation3D _RVD;
		std::vector<_Point3> _sites{};
		std::function<double(_Point3& p)> _rho;
		_LBFGS::_Parameter _para;
		std::string _modelname;
		std::set<MyFaceCVT> rdtFaces;
	};
} // namespace BGAL
