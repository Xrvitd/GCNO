#include "BGAL/Tessellation2D/Side2D.h"

typedef CGAL::Simple_cartesian<double> SIDEK;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> FK;
typedef CGAL::Simple_cartesian<CGAL::MP_Float> EK;
typedef CGAL::Cartesian_converter<SIDEK, EK> C2E;
typedef CGAL::Cartesian_converter<SIDEK, FK> C2F;
typedef SIDEK::Comparison_result Comparison_result;

namespace BGAL
{
  inline _BOC::_Sign _Side2D::base_side1_(const double &p1x,
                                          const double &p1y,
                                          const double &w1,
                                          const double &p2x,
                                          const double &p2y,
                                          const double &w2,
                                          const double &qx,
                                          const double &qy)
  {
    double two = 2.00000000000000000000e+00;
    double p1x_bound = fabs(p1x);
    double p1y_bound = fabs(p1y);
    double w1_bound = fabs(w1);
    double p2x_bound = fabs(p2x);
    double p2y_bound = fabs(p2y);
    double w2_bound = fabs(w2);
    double qx_bound = fabs(qx);
    double qy_bound = fabs(qy);
    double t1 = p1x - p2x;
    double t1_bound = p1x_bound + p2x_bound;
    double t2 = p1y - p2y;
    double t2_bound = p1y_bound + p2y_bound;
    double t3 = two * qx - p1x - p2x;
    double t3_bound = two * qx_bound + p1x_bound + p2x_bound;
    double t4 = two * qy - p1y - p2y;
    double t4_bound = two * qy_bound + p1y_bound + p2y_bound;
    _BOC::_Sign int_tmp_result;
    if (fabs(t1 * t3 + t2 * t4 + w1 - w2) > (t1_bound * t3_bound + t2_bound * t4_bound + w1_bound + w2_bound) * 1.33226762955018784851e-15 && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0))
    {
      int_tmp_result = (((t1 * t3 + t2 * t4 + w1 - w2) < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE
                                                                                     : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    return int_tmp_result;
  }
  template <class RT>
  inline RT _Side2D::side1_generic(const RT &p1x, const RT &p1y, const RT &w1,
                                   const RT &p2x, const RT &p2y, const RT &w2,
                                   const RT &qx, const RT &qy)
  {
    static RT const two(2);
    return RT((p1x - p2x) * (two * qx - p1x - p2x) + (p1y - p2y) * (two * qy - p1y - p2y) + w1 - w2);
  }
  template <class RT>
  inline void _Side2D::side2_generic(RT &denom,
                                     RT &numer,
                                     const RT &p1x,
                                     const RT &p1y,
                                     const RT &w1,
                                     const RT &p2x,
                                     const RT &p2y,
                                     const RT &w2,
                                     const RT &p3x,
                                     const RT &p3y,
                                     const RT &w3,
                                     const RT &q1x,
                                     const RT &q1y,
                                     const RT &q2x,
                                     const RT &q2y)
  {
    RT two(2.00000000000000000000e+00);
    RT four(4.00000000000000000000e+00);
    RT t1 = p1x - p2x;
    RT t2 = p1y - p2y;
    RT t3 = p1x - p3x;
    RT t4 = p1y - p3y;
    RT t5 = q1x - q2x;
    RT t6 = q1y - q2y;
    denom = -two * t3 * t5 - two * t4 * t6;
    RT t7 = p1x + p2x;
    RT t8 = p1y + p2y;
    RT t9 = p1x + p3x;
    RT t10 = p1y + p3y;
    RT t11 = t3 * t9 + t4 * t10 + w3 - w1;
    RT t12 = two * t4 * (q1y * t5 - q1x * t6) - t5 * t11;
    RT t13 = -t6 * t11 - two * t3 * (q1y * t5 - q1x * t6);
    numer = two * (t1 * t12 + t2 * t13) - (t1 * t7 + t2 * t8 - w1 + w2) * denom;
  }
  inline _Side2D::Comparison_result _Side2D::side1_filtered(const double &p1x,
                                                            const double &p1y,
                                                            const double &w1,
                                                            const double &p2x,
                                                            const double &p2y,
                                                            const double &w2,
                                                            const double &qx,
                                                            const double &qy)
  {
    {
      C2F c2f;
      FK::Comparison_result res = sign(
          side1_generic(
              c2f(p1x), c2f(p1y), c2f(w1),
              c2f(p2x), c2f(p2y), c2f(w2),
              c2f(qx), c2f(qy)));
      if (is_certain(res))
      {
        return get_certain(res);
      }
    }
    C2E c2e;
    // Exact version (if filtered version failed)
    return sign(
        side1_generic(
            c2e(p1x), c2e(p1y), c2e(w1),
            c2e(p2x), c2e(p2y), c2e(w2),
            c2e(qx), c2e(qy)));
  }

  inline _BOC::_Sign _Side2D::base_side2_(const double &p1x,
                                          const double &p1y,
                                          const double &w1,
                                          const double &p2x,
                                          const double &p2y,
                                          const double &w2,
                                          const double &p3x,
                                          const double &p3y,
                                          const double &w3,
                                          const double &q1x,
                                          const double &q1y,
                                          const double &q2x,
                                          const double &q2y)
  {
    double two = 2.00000000000000000000e+00;
    double four = 4.00000000000000000000e+00;
    double p1x_bound = fabs(p1x);
    double p1y_bound = fabs(p1y);
    double w1_bound = fabs(w1);
    double p2x_bound = fabs(p2x);
    double p2y_bound = fabs(p2y);
    double w2_bound = fabs(w2);
    double p3x_bound = fabs(p3x);
    double p3y_bound = fabs(p3y);
    double w3_bound = fabs(w3);
    double q1x_bound = fabs(q1x);
    double q1y_bound = fabs(q1y);
    double q2x_bound = fabs(q2x);
    double q2y_bound = fabs(q2y);

    double t1 = p1x - p2x;
    double t1_bound = p1x_bound + p2x_bound;
    double t2 = p1y - p2y;
    double t2_bound = p1y_bound + p2y_bound;
    double t3 = p1x - p3x;
    double t3_bound = p1x_bound + p3x_bound;
    double t4 = p1y - p3y;
    double t4_bound = p1y_bound + p3y_bound;
    double t5 = q1x - q2x;
    double t5_bound = q1x_bound + q2x_bound;
    double t6 = q1y - q2y;
    double t6_bound = q1y_bound + q2y_bound;

    double denom = -two * t3 * t5 - two * t4 * t6;
    double denom_bound = two * t3_bound * t5_bound + two * t4_bound * t6_bound;
    _BOC::_Sign denom_sign;
    if (fabs(denom) > (denom_bound * 1.11022302462515654042e-15) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0))
    {
      denom_sign = ((denom < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    double t7 = p1x + p2x;
    double t7_bound = p1x_bound + p2x_bound;
    double t8 = p1y + p2y;
    double t8_bound = p1y_bound + p2y_bound;
    double t9 = p1x + p3x;
    double t9_bound = p1x_bound + p3x_bound;
    double t10 = p1y + p3y;
    double t10_bound = p1y_bound + p3y_bound;
    double t11 = t3 * t9 + t4 * t10 + w3 - w1;
    double t11_bound = t3_bound * t9_bound + t4_bound * t10_bound + w3_bound + w1_bound;

    double t12 = two * t4 * (q1y * t5 - q1x * t6) - t5 * t11;
    double t12_bound = two * t4_bound * (q1y_bound * t5_bound + q1x_bound * t6_bound) + t5_bound * t11_bound;
    double t13 = -t6 * t11 - two * t3 * (q1y * t5 - q1x * t6);
    double t13_bound = t6_bound * t11_bound + two * t3_bound * (q1y_bound * t5_bound + q1x_bound * t6_bound);

    double numer = two * (t1 * t12 + t2 * t13) - (t1 * t7 + t2 * t8 - w1 + w2) * denom;
    double numer_bound = two * (t1_bound * t12_bound + t2_bound * t13_bound) + (t1_bound * t7_bound + t2_bound * t8_bound + w1_bound + w2_bound) * denom_bound;

    _BOC::_Sign numer_sign;
    if (fabs(numer) > (numer_bound * 3.10862446895043831319e-15) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0))
    {
      numer_sign = ((numer < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    return ((denom_sign == numer_sign) ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE);
  }

  inline _Side2D::Comparison_result _Side2D::side2_filtered(const double &p1x,
                                                            const double &p1y,
                                                            const double &w1,
                                                            const double &p2x,
                                                            const double &p2y,
                                                            const double &w2,
                                                            const double &p3x,
                                                            const double &p3y,
                                                            const double &w3,
                                                            const double &q1x,
                                                            const double &q1y,
                                                            const double &q2x,
                                                            const double &q2y)
  {
    {
      C2F c2f;
      FK::FT denom;
      FK::FT numer;
      side2_generic(
          denom, numer,
          c2f(p1x), c2f(p1y), c2f(w1),
          c2f(p2x), c2f(p2y), c2f(w2),
          c2f(p3x), c2f(p3y), c2f(w3),
          c2f(q1x), c2f(q1y),
          c2f(q2x), c2f(q2y));
      FK::Comparison_result s1 = sign(numer);
      FK::Comparison_result s2 = sign(denom);

      if (is_certain(s1) && is_certain(s2))
      {
        return mult(get_certain(s1), get_certain(s2));
      }
    }
    C2E c2e;
    EK::FT denom;
    EK::FT numer;
    side2_generic(
        denom, numer,
        c2e(p1x), c2e(p1y), c2e(w1),
        c2e(p2x), c2e(p2y), c2e(w2),
        c2e(p3x), c2e(p3y), c2e(w3),
        c2e(q1x), c2e(q1y),
        c2e(q2x), c2e(q2y));
    return mult(sign(numer), sign(denom));
  }

  _BOC::_Sign _Side2D::side1_(const double &p1x,
                              const double &p1y,
                              const double &w1,
                              const double &p2x,
                              const double &p2y,
                              const double &w2,
                              const double &qx,
                              const double &qy)
  {
    _BOC::_Sign result = base_side1_(p1x, p1y, w1, p2x, p2y, w2, qx, qy);
    if (result == _BOC::_Sign::FaileD)
    {
      Comparison_result r = side1_filtered(p1x, p1y, w1, p2x, p2y, w2, qx, qy);
      result =
          (r == Comparison_result::ZERO ? _BOC::_Sign::ZerO : (r == Comparison_result::POSITIVE ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE));
    }
    return result;
  }

  _BOC::_Sign _Side2D::side2_(const double &p1x,
                              const double &p1y,
                              const double &w1,
                              const double &p2x,
                              const double &p2y,
                              const double &w2,
                              const double &p3x,
                              const double &p3y,
                              const double &w3,
                              const double &q1x,
                              const double &q1y,
                              const double &q2x,
                              const double &q2y)
  {
    _BOC::_Sign result = base_side2_(p1x, p1y, w1, p2x, p2y, w2, p3x, p3y, w3, q1x, q1y, q2x, q2y);
    if (result == _BOC::_Sign::FaileD)
    {
      Comparison_result r = side2_filtered(p1x, p1y, w1, p2x, p2y, w2, p3x, p3y, w3, q1x, q1y, q2x, q2y);
      result =
          (r == Comparison_result::ZERO ? _BOC::_Sign::ZerO : (r == Comparison_result::POSITIVE ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE));
    }
    return result;
  }
} // namespace BGAL
