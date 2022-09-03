#pragma once
#include "BGAL/BaseShape/Point.h"
#include "BGAL/Algorithm/BOC/BOC.h"
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
    class _Side3D
    {
    private:
        typedef CGAL::Simple_cartesian<double> K;
        typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> FK;
        typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;
        typedef CGAL::Cartesian_converter<K, EK> C2E;
        typedef CGAL::Cartesian_converter<K, FK> C2F;
        typedef K::Comparison_result Comparison_result;
        static inline K::Comparison_result mult(K::Comparison_result s1, K::Comparison_result s2);
        static inline _BOC::_Sign base_side1_(double p1x,
                                              double p1y,
                                              double p1z,
                                              double p2x,
                                              double p2y,
                                              double p2z,
                                              double w1,
                                              double w2,
                                              double qx,
                                              double qy,
                                              double qz);
        template <class RT>
        static RT side1_generic(
            const RT &p1x, const RT &p1y, const RT &p1z,
            const RT &p2x, const RT &p2y, const RT &p2z,
            const RT &w1, const RT &w2,
            const RT &qx, const RT &qy, const RT &qz);
        static Comparison_result side1_filtered(
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double w1, double w2,
            double qx, double qy, double qz);
        static _BOC::_Sign base_side2_(
            double p1x,
            double p1y,
            double p1z,
            double p2x,
            double p2y,
            double p2z,
            double p3x,
            double p3y,
            double p3z,
            double w1,
            double w2,
            double w3,
            double q1x,
            double q1y,
            double q1z,
            double q2x,
            double q2y,
            double q2z);
        template <class RT>
        static void side2_generic(
            RT &result,
            RT &denom,
            const RT &p1x, const RT &p1y, const RT &p1z,
            const RT &p2x, const RT &p2y, const RT &p2z,
            const RT &p3x, const RT &p3y, const RT &p3z,
            const RT &w1, const RT &w2, const RT &w3,
            const RT &q1x, const RT &q1y, const RT &q1z,
            const RT &q2x, const RT &q2y, const RT &q2z);
        static Comparison_result side2_filtered(
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double p3x, double p3y, double p3z,
            double w1, double w2, double w3,
            double q1x, double q1y, double q1z,
            double q2x, double q2y, double q2z);
        static _BOC::_Sign base_side3_(
            double p1x,
            double p1y,
            double p1z,
            double p2x,
            double p2y,
            double p2z,
            double p3x,
            double p3y,
            double p3z,
            double p4x,
            double p4y,
            double p4z,
            double w1,
            double w2,
            double w3,
            double w4,
            double q1x,
            double q1y,
            double q1z,
            double q2x,
            double q2y,
            double q2z,
            double q3x,
            double q3y,
            double q3z);
        template <class RT>
        static void side3_generic(
            RT &result,
            RT &denom,
            const RT &p1x, const RT &p1y, const RT &p1z,
            const RT &p2x, const RT &p2y, const RT &p2z,
            const RT &p3x, const RT &p3y, const RT &p3z,
            const RT &p4x, const RT &p4y, const RT &p4z,
            const RT &w1, const RT &w2, const RT &w3, const RT &w4,
            const RT &q1x, const RT &q1y, const RT &q1z,
            const RT &q2x, const RT &q2y, const RT &q2z,
            const RT &q3x, const RT &q3y, const RT &q3z);
        static Comparison_result side3_filtered(
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double p3x, double p3y, double p3z,
            double p4x, double p4y, double p4z,
            double w1, double w2, double w3, double w4,
            double q1x, double q1y, double q1z,
            double q2x, double q2y, double q2z,
            double q3x, double q3y, double q3z);

    public:
        static _BOC::_Sign side1_(
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double w1, double w2,
            double qx, double qy, double qz);
        static _BOC::_Sign side2_(
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double p3x, double p3y, double p3z,
            double w1, double w2, double w3,
            double q1x, double q1y, double q1z,
            double q2x, double q2y, double q2z);
        static _BOC::_Sign side3_(
            double p1x, double p1y, double p1z,
            double p2x, double p2y, double p2z,
            double p3x, double p3y, double p3z,
            double p4x, double p4y, double p4z,
            double w1, double w2, double w3, double w4,
            double q1x, double q1y, double q1z,
            double q2x, double q2y, double q2z,
            double q3x, double q3y, double q3z);
    };
} // namespace BGAL