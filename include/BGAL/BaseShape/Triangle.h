#pragma once
#include "Point.h"

namespace BGAL
{
	class _Triangle3
	{
	private:
		std::vector<_Point3> _points;
		bool _open_in;
		double _area;
	public:
		_Triangle3();
		_Triangle3(const std::vector<_Point3>& in_points);
		_Triangle3(const _Point3& p1, const _Point3& p2, const _Point3& p3);
		void start_();
		void insert_(const _Point3& in_p);
		void insert_(const double& in_x, const double& in_y, const double& in_z);
		void end_();
		const _Point3& operator[](const int& in_inx) const;
		const _Point3& point(const int& in_inx) const;
		double area_();
	};
}