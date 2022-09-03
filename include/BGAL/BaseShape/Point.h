#pragma once
#include <vector>
#include <iostream>
#include "BGAL/Algorithm/BOC/BOC.h"
#include <Eigen/Dense>
namespace BGAL 
{
	class _Point3;
	class _Point2;
	class _Point2 
	{
	private:
		std::vector<double> _value;
	public:
		_Point2();
		_Point2(const double& in_x, const double& in_y);
		_Point2(const _Point2& in_p);
		std::vector<double> get_value_() const;
		double x() const;
		double y() const;
		double sqlength_() const;
		double length_() const;
		_Point2 normalize_() const;
		_Point2& normalized_();
		friend std::ostream& operator<<(std::ostream& in_os, _Point2 in_p);
		double operator[](const int& in_inx) const;
		const _Point2& operator=(const _Point2& in_p);
		bool operator==(const _Point2& in_p) const;
		bool operator<(const _Point2& in_p) const;
		bool operator<=(const _Point2& in_p) const;
		bool operator>(const _Point2& in_p) const;
		bool operator>=(const _Point2& in_p) const;
		bool operator!=(const _Point2& in_p) const;
		_Point2 operator+(const _Point2& in_p) const;
		_Point2 operator-(const _Point2& in_p) const;
		//_Point2 operator-() const;
		_Point2& operator+=(const _Point2& in_p);
		_Point2& operator-=(const _Point2& in_p);
		_Point2 operator*(const _Point2& in_p) const;
		_Point2 operator*(const double& in_s) const;
		_Point2 operator/(const _Point2& in_p) const;
		_Point2 operator/(const double& in_s) const;
		_Point2& operator*=(const _Point2& in_p);
		_Point2& operator*=(const double& in_s);
		_Point2& operator/=(const _Point2& in_p);
		_Point2& operator/=(const double& in_s);
		double dot_(const _Point2& in_p) const;
		_Point3 cross_(const _Point2& in_p) const;
		static _Point2 intersection_two_line(const _Point2& v0, const double& d0, const _Point2& v1, const double& d1) 
		{
			double r0 = v0.x() * v1.y() - v0.y() * v1.x();
			double r1 = v0.y() * d1 - v1.y() * d0;
			double r2 = v1.x() * d0 - v0.x() * d1;
			_Point2 r(r1 / r0, r2 / r0);
			return r;
		}
	};

	class _Point3 
	{
	private:
		std::vector<double> _value;
	public:
		_Point3();
		_Point3(const double& in_x, const double& in_y, const double& in_z);
		_Point3(const _Point3& in_p);
		std::vector<double> get_value_() const;
		double x() const;
		double y() const;
		double z() const;
		double sqlength_() const;
		double length_() const;
		_Point3 normalize_() const;
		_Point3& normalized_();
		friend std::ostream& operator<<(std::ostream& in_os, const _Point3& in_p);
		double operator[](const int in_inx) const;
		double& operator[](const int in_inx);
		_Point3 rotate_(const Eigen::Matrix3d& RM) const;
		const _Point3& operator=(const _Point3& in_p);
		bool operator==(const _Point3& in_p) const;
		bool operator<(const _Point3& in_p) const;
		bool operator<=(const _Point3& in_p) const;
		bool operator>(const _Point3& in_p) const;
		bool operator>=(const _Point3& in_p) const;
		bool operator!=(const _Point3& in_p) const;
		_Point3 operator+(const _Point3& in_p) const;
		_Point3 operator-(const _Point3& in_p) const;
		_Point3 operator-() const;
		_Point3& operator+=(const _Point3& in_p);
		_Point3& operator-=(const _Point3& in_p);
		_Point3 operator*(const _Point3& in_p) const;
		_Point3 operator*(const double& in_s) const;
		friend _Point3 operator*(const double& in_s, const _Point3& in_p);
		_Point3 operator/(const _Point3& in_p) const;
		_Point3 operator/(const double& in_s) const;
		_Point3& operator*=(const _Point3& in_p);
		_Point3& operator*=(const double& in_s);
		_Point3& operator/=(const _Point3& in_p);
		_Point3& operator/=(const double& in_s);
		double dot_(const _Point3& in_p) const;
		_Point3 cross_(const _Point3& in_p) const;
		static _Point3 intersection_three_plane(const _Point3& v0, const double& d0, const _Point3& v1, const double& d1, const _Point3& v2, const double& d2) 
		{
			double det = v0.x() * (v1.y() * v2.z() - v1.z() * v2.y()) - v0.y() * (v1.x() * v2.z() - v1.z() * v2.x())
				+ v0.z() * (v1.x() * v2.y() - v1.y() * v2.x());
			double x = -d0 * (v1.y() * v2.z() - v1.z() * v2.y()) - v0.y() * (-d1 * v2.z() + v1.z() * d2)
				+ v0.z() * (-d1 * v2.y() + v1.y() * d2);
			x /= det;
			double y = v0.x() * (-d1 * v2.z() + v1.z() * d2) + d0 * (v1.x() * v2.z() - v1.z() * v2.x())
				+ v0.z() * (-d2 * v1.x() + d1 * v2.x());
			y /= det;
			double z = v0.x() * (-d2 * v1.y() + d1 * v2.y()) - v0.y() * (-d2 * v1.x() + d1 * v2.x())
				- d0 * (v1.x() * v2.y() - v1.y() * v2.x());
			z /= det;
			return _Point3(x, y, z);
		}
	};
}