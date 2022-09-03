#pragma once
#include "BGAL/BaseShape/Point.h"

#include <fenv.h>
#if WIN32 || _WIN32
#include <corecrt_math.h>
#endif
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/determinant.h>
namespace BGAL
{
	class _Side2D
	{
	private:
		typedef CGAL::Simple_cartesian<double> SIDEK;
		typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> FK;
		typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;
		typedef CGAL::Cartesian_converter<SIDEK, EK> C2E;
		typedef CGAL::Cartesian_converter<SIDEK, FK> C2F;
		typedef SIDEK::Comparison_result Comparison_result;
		static inline _BOC::_Sign base_side1_(const double &p1x, const double &p1y, const double &w1,
											  const double &p2x, const double &p2y, const double &w2,
											  const double &qx, const double &qy);
		template <class RT>
		static inline RT side1_generic(const RT &p1x, const RT &p1y, const RT &w1,
									   const RT &p2x, const RT &p2y, const RT &w2,
									   const RT &qx, const RT &qy);
		static inline Comparison_result side1_filtered(const double &p1x, const double &p1y, const double &w1,
													   const double &p2x, const double &p2y, const double &w2,
													   const double &qx, const double &qy);

		static inline _BOC::_Sign base_side2_(const double &p1x, const double &p1y, const double &w1,
											  const double &p2x, const double &p2y, const double &w2,
											  const double &p3x, const double &p3y, const double &w3,
											  const double &q1x, const double &q1y,
											  const double &q2x, const double &q2y);
		template <class RT>
		static inline void side2_generic(RT &denom, RT &numer,
										 const RT &p1x, const RT &p1y, const RT &w1,
										 const RT &p2x, const RT &p2y, const RT &w2,
										 const RT &p3x, const RT &p3y, const RT &w3,
										 const RT &q1x, const RT &q1y,
										 const RT &q2x, const RT &q2y);
		static inline Comparison_result side2_filtered(const double &p1x, const double &p1y, const double &w1,
													   const double &p2x, const double &p2y, const double &w2,
													   const double &p3x, const double &p3y, const double &w3,
													   const double &q1x, const double &q1y,
													   const double &q2x, const double &q2y);
		static inline Comparison_result mult(Comparison_result s1, Comparison_result s2)
		{
			if (s1 == CGAL::EQUAL || s2 == CGAL::EQUAL)
			{
				return CGAL::EQUAL;
			}
			return (s1 == s2 ? CGAL::LARGER : CGAL::SMALLER);
		}

	public:
		static _BOC::_Sign side1_(const double &p1x, const double &p1y, const double &w1,
								  const double &p2x, const double &p2y, const double &w2,
								  const double &qx, const double &qy);
		static _BOC::_Sign side2_(const double &p1x, const double &p1y, const double &w1,
								  const double &p2x, const double &p2y, const double &w2,
								  const double &p3x, const double &p3y, const double &w3,
								  const double &q1x, const double &q1y,
								  const double &q2x, const double &q2y);
	};
} // namespace BGAL