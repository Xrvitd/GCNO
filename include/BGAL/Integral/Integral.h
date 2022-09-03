#pragma once
#include "Tetrahedron_arbq_rule.h"
#include "BGAL/BaseShape/Polygon.h"
#include <Eigen/Dense>
namespace BGAL 
{
	class _Integral 
	{
	public:
		template<class F>
		static Eigen::VectorXd integral_triangle(F f, const _Polygon& poly) 
		{
			if (poly.num_() != 3)
				throw "poly is not a triangle!";
			Eigen::VectorXd r1 = 1.0 / 30 * f(poly[1] * 0.5 + poly[2] * 0.5);
			Eigen::VectorXd r2 = 1.0 / 30 * f(poly[0] * 0.5 + poly[1] * 0.5);
			Eigen::VectorXd r3 = 1.0 / 30 * f(poly[0] * 0.5 + poly[2] * 0.5);
			Eigen::VectorXd r4 = 9.0 / 30 * f(poly[0] / 6.0 + poly[1] / 6.0 + poly[2] * 2.0 / 3.0);
			Eigen::VectorXd r5 = 9.0 / 30 * f(poly[2] / 6.0 + poly[0] / 6.0 + poly[1] * 2.0 / 3.0);
			Eigen::VectorXd r6 = 9.0 / 30 * f(poly[1] / 6.0 + poly[2] / 6.0 + poly[0] * 2.0 / 3.0);
			return poly.triangle_area_() * (r1 + r2 + r3 + r4 + r5 + r6);
		}
		template<class F>
		static Eigen::VectorXd integral_triangle(F f, const _Point2& p1, const _Point2& p2, const _Point2& p3) 
		{
			Eigen::VectorXd r1 = 1.0 / 30 * f(p2 * 0.5 + p3 * 0.5);
			Eigen::VectorXd r2 = 1.0 / 30 * f(p1 * 0.5 + p2 * 0.5);
			Eigen::VectorXd r3 = 1.0 / 30 * f(p1 * 0.5 + p3 * 0.5);
			Eigen::VectorXd r4 = 9.0 / 30 * f(p1 / 6.0 + p2 / 6.0 + p3 * 2.0 / 3.0);
			Eigen::VectorXd r5 = 9.0 / 30 * f(p3 / 6.0 + p1 / 6.0 + p2 * 2.0 / 3.0);
			Eigen::VectorXd r6 = 9.0 / 30 * f(p2 / 6.0 + p3 / 6.0 + p1 * 2.0 / 3.0);
			_Point2 v1 = p2 - p1;
			_Point2 v2 = p3 - p1;
			double area = v1.cross_(v2).length_() * 0.5;
			return area * (r1 + r2 + r3 + r4 + r5 + r6);
		}
		template<class F>
		static Eigen::VectorXd integral_triangle3D(F f, const _Point3& p1, const _Point3& p2, const _Point3& p3) 
		{
			Eigen::VectorXd r1 = 1.0 / 30 * f(p2 * 0.5 + p3 * 0.5);
			Eigen::VectorXd r2 = 1.0 / 30 * f(p1 * 0.5 + p2 * 0.5);
			Eigen::VectorXd r3 = 1.0 / 30 * f(p1 * 0.5 + p3 * 0.5);
			Eigen::VectorXd r4 = 9.0 / 30 * f(p1 / 6.0 + p2 / 6.0 + p3 * 2.0 / 3.0);
			Eigen::VectorXd r5 = 9.0 / 30 * f(p3 / 6.0 + p1 / 6.0 + p2 * 2.0 / 3.0);
			Eigen::VectorXd r6 = 9.0 / 30 * f(p2 / 6.0 + p3 / 6.0 + p1 * 2.0 / 3.0);
			_Point3 v1 = p2 - p1;
			_Point3 v2 = p3 - p1;
			double area = v1.cross_(v2).length_() * 0.5;
			return area * (r1 + r2 + r3 + r4 + r5 + r6);
		}
		template<class F>
		static Eigen::VectorXd integral_polygon(F f, const _Polygon& poly) 
		{
			std::vector<_Polygon> tris = poly.constrained_delaunay_triangulation_();
			Eigen::VectorXd r = integral_triangle(f, tris[0]);
			for (int i = 1; i < tris.size(); ++i) {
				r = r + integral_triangle(f, tris[i]);
			}
			return r;
		}
		template<class F>
		static Eigen::VectorXd integral_polygon_fast(F f, const _Polygon& poly) 
		{
			Eigen::VectorXd r = integral_triangle(f, poly[0], poly[1], poly[2]);
			for (int i = 2; i < poly.num_() - 1; ++i) {
				r = r + integral_triangle(f, poly[0], poly[i], poly[i + 1]);
			}
			return r;
		}
		template<class F>
		static Eigen::VectorXd integral_tetrahedron(F f, const _Point3& p1, const _Point3& p2, const _Point3& p3, const _Point3& p4) 
		{
			double node_xyz[3 * 4] = 
			{
				0.0, 0.0, 0.0,
				1.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
				0.0, 0.0, 1.0 
			};
			double node_xyz2[12];
			node_xyz2[0] = p1.x();
			node_xyz2[1] = p1.y();
			node_xyz2[2] = p1.z();
			node_xyz2[3] = p2.x();
			node_xyz2[4] = p2.y();
			node_xyz2[5] = p2.z();
			node_xyz2[6] = p3.x();
			node_xyz2[7] = p3.y();
			node_xyz2[8] = p3.z();
			node_xyz2[9] = p4.x();
			node_xyz2[10] = p4.y();
			node_xyz2[11] = p4.z();

			double* w;
			double* xyz;
			double* xyz2;
			int order_num = keast_order_num(4);
			xyz = new double[3 * order_num];
			xyz2 = new double[3 * order_num];
			w = new double[order_num];

			keast_rule(4, order_num, xyz, w);

			tetrahedron_reference_to_physical(node_xyz2, order_num, xyz, xyz2);

			double volume = tetrahedron_volume(node_xyz2);
			std::vector<_Point3> samples;
			for (int i = 0; i < order_num; ++i) 
			{
				samples.push_back(_Point3(xyz2[0 + i * 3], xyz2[1 + i * 3], xyz2[2 + i * 3]));
			}
			Eigen::VectorXd r = volume * w[0] * f(samples[0]);
			for (int i = 1; i < order_num; ++i) 
			{
				r = r + volume * w[i] * f(samples[i]);
			}

			delete[] w;
			delete[] xyz;
			delete[] xyz2;

			return r;
		}
	};
}