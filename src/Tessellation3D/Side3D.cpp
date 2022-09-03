#include "BGAL/Tessellation3D/Side3D.h"

namespace BGAL
{

  _Side3D::K::Comparison_result _Side3D::mult(K::Comparison_result s1, K::Comparison_result s2)
  {
    if (s1 == CGAL::EQUAL || s2 == CGAL::EQUAL)
    {
      return CGAL::EQUAL;
    }
    return (s1 == s2 ? CGAL::LARGER : CGAL::SMALLER);
  }

  _BOC::_Sign _Side3D::base_side1_(double p1x,
                                   double p1y,
                                   double p1z,
                                   double p2x,
                                   double p2y,
                                   double p2z,
                                   double w1,
                                   double w2,
                                   double qx,
                                   double qy,
                                   double qz)
  {
    double p1x_bound = fabs(p1x);
    double p2x_bound = fabs(p2x);
    double v1x_bound;
    double v1x;
    v1x_bound = (p1x_bound + p2x_bound);
    v1x = (p1x - p2x);
    double p1y_bound = fabs(p1y);
    double p2y_bound = fabs(p2y);
    double v1y_bound;
    double v1y;
    v1y_bound = (p1y_bound + p2y_bound);
    v1y = (p1y - p2y);
    double p1z_bound = fabs(p1z);
    double p2z_bound = fabs(p2z);
    double v1z_bound;
    double v1z;
    v1z_bound = (p1z_bound + p2z_bound);
    v1z = (p1z - p2z);
    double qx_bound = fabs(qx);
    double v2x_bound;
    double v2x;
    v2x_bound = ((2.00000000000000000000e+00 * qx_bound) + (p1x_bound + p2x_bound));
    v2x = ((2.00000000000000000000e+00 * qx) - (p1x + p2x));
    double qy_bound = fabs(qy);
    double v2y_bound;
    double v2y;
    v2y_bound = ((2.00000000000000000000e+00 * qy_bound) + (p1y_bound + p2y_bound));
    v2y = ((2.00000000000000000000e+00 * qy) - (p1y + p2y));
    double qz_bound = fabs(qz);
    double v2z_bound;
    double v2z;
    v2z_bound = ((2.00000000000000000000e+00 * qz_bound) + (p1z_bound + p2z_bound));
    v2z = ((2.00000000000000000000e+00 * qz) - (p1z + p2z));
    _BOC::_Sign int_tmp_result;
    if (((fabs((((v1x * v2x) + (v1y * v2y)) + (v1z * v2z) + w1 - w2)) > ((((v1x_bound * v2x_bound) + (v1y_bound * v2y_bound)) + (v1z_bound * v2z_bound) + fabs(w1) + fabs(w2)) * 1.33226762955018784851e-15)) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0)))
    {
      int_tmp_result =
          (((((v1x * v2x) + (v1y * v2y)) + (v1z * v2z) + w1 - w2) < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE
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
  RT _Side3D::side1_generic(
      const RT &p1x, const RT &p1y, const RT &p1z,
      const RT &p2x, const RT &p2y, const RT &p2z,
      const RT &w1, const RT &w2,
      const RT &qx, const RT &qy, const RT &qz)
  {
    static RT const two(2);
    return (p1x - p2x) * (two * qx - p1x - p2x) + (p1y - p2y) * (two * qy - p1y - p2y) +
           (p1z - p2z) * (two * qz - p1z - p2z) + w1 - w2;
  }

  _Side3D::Comparison_result _Side3D::side1_filtered(double p1x,
                                                     double p1y,
                                                     double p1z,
                                                     double p2x,
                                                     double p2y,
                                                     double p2z,
                                                     double w1,
                                                     double w2,
                                                     double qx,
                                                     double qy,
                                                     double qz)
  {
    {
      C2F c2f;
      FK::Comparison_result res = sign(
          side1_generic(
              c2f(p1x), c2f(p1y), c2f(p1z),
              c2f(p2x), c2f(p2y), c2f(p2z),
              c2f(w1), c2f(w2),
              c2f(qx), c2f(qy), c2f(qz)));
      if (is_certain(res))
      {
        return get_certain(res);
      }
    }
    C2E c2e;
    // Exact version (if filtered version failed)
    return sign(
        side1_generic(
            c2e(p1x), c2e(p1y), c2e(p1z),
            c2e(p2x), c2e(p2y), c2e(p2z),
            c2e(w1), c2e(w2),
            c2e(qx), c2e(qy), c2e(qz)));
  }

  _BOC::_Sign _Side3D::base_side2_(double p1x,
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
                                   double q2z)
  {
    double q2z_bound = fabs(q2z);
    double q1z_bound = fabs(q1z);
    double t17_bound;
    double t17;
    t17_bound = (q2z_bound + q1z_bound);
    t17 = (q2z - q1z);
    double q2y_bound = fabs(q2y);
    double q1y_bound = fabs(q1y);
    double t18_bound;
    double t18;
    t18_bound = (q2y_bound + q1y_bound);
    t18 = (q2y - q1y);
    double q2x_bound = fabs(q2x);
    double q1x_bound = fabs(q1x);
    double t19_bound;
    double t19;
    t19_bound = (q2x_bound + q1x_bound);
    t19 = (q2x - q1x);
    double p1z_bound = fabs(p1z);
    double p3z_bound = fabs(p3z);
    double t20_bound;
    double t20;
    t20_bound = (p1z_bound + p3z_bound);
    t20 = (p1z - p3z);
    double p1y_bound = fabs(p1y);
    double p3y_bound = fabs(p3y);
    double t21_bound;
    double t21;
    t21_bound = (p1y_bound + p3y_bound);
    t21 = (p1y - p3y);
    double p1x_bound = fabs(p1x);
    double p3x_bound = fabs(p3x);
    double t22_bound;
    double t22;
    t22_bound = (p1x_bound + p3x_bound);
    t22 = (p1x - p3x);
    double denom_bound;
    double denom;
    denom_bound = (((t19_bound * t22_bound) + (t18_bound * t21_bound)) + (t17_bound * t20_bound));
    denom = (((t19 * t22) + (t18 * t21)) + (t17 * t20));
    _BOC::_Sign int_tmp_result;
    if (((fabs(denom) > (denom_bound * 1.11022302462515654042e-15)) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0)))
    {
      int_tmp_result = ((denom < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    double p2x_bound = fabs(p2x);
    double t28_bound;
    double t28;
    t28_bound = (p1x_bound + p2x_bound);
    t28 = (p1x - p2x);
    double t27_bound;
    double t27;
    t27_bound = (p1x_bound + (2 * q1x_bound));
    t27 = (p1x - (2 * q1x));
    double p2y_bound = fabs(p2y);
    double t26_bound;
    double t26;
    t26_bound = (p1y_bound + p2y_bound);
    t26 = (p1y - p2y);
    double t25_bound;
    double t25;
    t25_bound = (p1y_bound + (2 * q1y_bound));
    t25 = (p1y - (2 * q1y));
    double p2z_bound = fabs(p2z);
    double t24_bound;
    double t24;
    t24_bound = (p1z_bound + p2z_bound);
    t24 = (p1z - p2z);
    double t23_bound;
    double t23;
    t23_bound = (p1z_bound + (2 * q1z_bound));
    t23 = (p1z - (2 * q1z));
    double result_bound;
    double result;
    result_bound = (((((t28_bound * t19_bound) + (t26_bound * t18_bound)) + (t24_bound * t17_bound)) * ((((p3x_bound + t27_bound) * t22_bound) + ((p3y_bound + t25_bound) * t21_bound)) + ((p3z_bound + t23_bound) * t20_bound) + fabs(w1) + fabs(w3))) + ((((t28_bound * (p2x_bound + t27_bound)) + (t26_bound * (p2y_bound + t25_bound))) + (t24_bound * (p2z_bound + t23_bound)) + fabs(w1) + fabs(w2)) * denom_bound));
    result = (((((t28 * t19) + (t26 * t18)) + (t24 * t17)) * ((((p3x + t27) * t22) + ((p3y + t25) * t21)) + ((p3z + t23) * t20) - w1 + w3)) + ((((t28 * (-p2x - t27)) + (t26 * (-p2y - t25))) + (t24 * (-p2z - t23)) + w1 - w2) * denom));
    _BOC::_Sign int_tmp_result_FFWKCAA;
    if (((fabs(result) > (result_bound * 3.10862446895043831319e-15)) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0)))
    {
      int_tmp_result_FFWKCAA = ((result < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    return (int_tmp_result == int_tmp_result_FFWKCAA) ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE;
  }

  template <class RT>
  void _Side3D::side2_generic(RT &result,
                              RT &denom,
                              const RT &p1x,
                              const RT &p1y,
                              const RT &p1z,
                              const RT &p2x,
                              const RT &p2y,
                              const RT &p2z,
                              const RT &p3x,
                              const RT &p3y,
                              const RT &p3z,
                              const RT &w1,
                              const RT &w2,
                              const RT &w3,
                              const RT &q1x,
                              const RT &q1y,
                              const RT &q1z,
                              const RT &q2x,
                              const RT &q2y,
                              const RT &q2z)
  {
    static RT const two(2);
    const RT t17 = q2z - q1z;
    const RT t18 = q2y - q1y;
    const RT t19 = q2x - q1x;
    const RT t20 = p1z - p3z;
    const RT t21 = p1y - p3y;
    const RT t22 = p1x - p3x;
    denom = (t19 * t22 + t18 * t21 + t17 * t20);
    const RT t28 = p1x - p2x;
    const RT t27 = p1x - two * q1x;
    const RT t26 = p1y - p2y;
    const RT t25 = p1y - two * q1y;
    const RT t24 = p1z - p2z;
    const RT t23 = p1z - two * q1z;
    result =
        (t28 * t19 + t26 * t18 + t24 * t17) *
            ((p3x + t27) * t22 + (p3y + t25) * t21 + (p3z + t23) * t20 - w1 + w3) +
        (t28 * (-p2x - t27) + t26 * (-p2y - t25) + t24 * (-p2z - t23) + w1 - w2) * denom;
  }

  _Side3D::Comparison_result _Side3D::side2_filtered(double p1x,
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
                                                     double q2z)
  {
    {
      C2F c2f;
      FK::FT result;
      FK::FT denom;
      side2_generic(
          result, denom,
          c2f(p1x), c2f(p1y), c2f(p1z),
          c2f(p2x), c2f(p2y), c2f(p2z),
          c2f(p3x), c2f(p3y), c2f(p3z),
          c2f(w1), c2f(w2), c2f(w3),
          c2f(q1x), c2f(q1y), c2f(q1z),
          c2f(q2x), c2f(q2y), c2f(q2z));
      FK::Comparison_result s1 = sign(result);
      FK::Comparison_result s2 = sign(denom);

      if (is_certain(s1) && is_certain(s2))
      {
        return mult(get_certain(s1), get_certain(s2));
      }
    }
    // Exact version (if filtered version failed)
    C2E c2e;
    EK::FT result;
    EK::FT denom;
    side2_generic(
        result, denom,
        c2e(p1x), c2e(p1y), c2e(p1z),
        c2e(p2x), c2e(p2y), c2e(p2z),
        c2e(p3x), c2e(p3y), c2e(p3z),
        c2e(w1), c2e(w2), c2e(w3),
        c2e(q1x), c2e(q1y), c2e(q1z),
        c2e(q2x), c2e(q2y), c2e(q2z));
    return mult(sign(result), sign(denom));
  }

  _BOC::_Sign _Side3D::base_side3_(double p1x,
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
                                   double q3z)
  {
    double p4z_bound = fabs(p4z);
    double p1z_bound = fabs(p1z);
    double t46_bound;
    double t46;
    t46_bound = (p4z_bound + p1z_bound);
    t46 = (p4z - p1z);
    double p4y_bound = fabs(p4y);
    double p1y_bound = fabs(p1y);
    double t47_bound;
    double t47;
    t47_bound = (p4y_bound + p1y_bound);
    t47 = (p4y - p1y);
    double p4x_bound = fabs(p4x);
    double p1x_bound = fabs(p1x);
    double t48_bound;
    double t48;
    t48_bound = (p4x_bound + p1x_bound);
    t48 = (p4x - p1x);
    double t55_bound;
    double t55;
    t55_bound = (((t48_bound * (p1x_bound + p4x_bound)) + (t47_bound * (p1y_bound + p4y_bound))) + (t46_bound * (p1z_bound + p4z_bound))) + fabs(w1) + fabs(w4);
    t55 = (((t48 * (p1x + p4x)) + (t47 * (p1y + p4y))) + (t46 * (p1z + p4z))) + w1 - w4;
    double p3z_bound = fabs(p3z);
    double t49_bound;
    double t49;
    t49_bound = (p1z_bound + p3z_bound);
    t49 = (-p1z + p3z);
    double p3y_bound = fabs(p3y);
    double t50_bound;
    double t50;
    t50_bound = (p1y_bound + p3y_bound);
    t50 = (-p1y + p3y);
    double p3x_bound = fabs(p3x);
    double t51_bound;
    double t51;
    t51_bound = (p1x_bound + p3x_bound);
    t51 = (-p1x + p3x);
    double t54_bound;
    double t54;
    t54_bound = (((t51_bound * (p1x_bound + p3x_bound)) + (t50_bound * (p1y_bound + p3y_bound))) + (t49_bound * (p1z_bound + p3z_bound))) + fabs(w1) + fabs(w4);
    t54 = (((t51 * (p1x + p3x)) + (t50 * (p1y + p3y))) + (t49 * (p1z + p3z))) + w1 - w3;
    double q3z_bound = fabs(q3z);
    double q1z_bound = fabs(q1z);
    double t40_bound;
    double t40;
    t40_bound = (q3z_bound + q1z_bound);
    t40 = (q3z - q1z);
    double q3x_bound = fabs(q3x);
    double q1x_bound = fabs(q1x);
    double t42_bound;
    double t42;
    t42_bound = (q3x_bound + q1x_bound);
    t42 = (q3x - q1x);
    double q2z_bound = fabs(q2z);
    double t43_bound;
    double t43;
    t43_bound = (q2z_bound + q1z_bound);
    t43 = (q2z - q1z);
    double q2x_bound = fabs(q2x);
    double t45_bound;
    double t45;
    t45_bound = (q2x_bound + q1x_bound);
    t45 = (q2x - q1x);
    double t37_bound;
    double t37;
    t37_bound = ((t43_bound * t42_bound) + (t45_bound * t40_bound));
    t37 = ((t43 * t42) - (t45 * t40));
    double q3y_bound = fabs(q3y);
    double q1y_bound = fabs(q1y);
    double t41_bound;
    double t41;
    t41_bound = (q3y_bound + q1y_bound);
    t41 = (q3y - q1y);
    double q2y_bound = fabs(q2y);
    double t44_bound;
    double t44;
    t44_bound = (q2y_bound + q1y_bound);
    t44 = (q2y - q1y);
    double t38_bound;
    double t38;
    t38_bound = ((t44_bound * t40_bound) + (t43_bound * t41_bound));
    t38 = ((t44 * t40) - (t43 * t41));
    double t39_bound;
    double t39;
    t39_bound = ((t45_bound * t41_bound) + (t44_bound * t42_bound));
    t39 = ((t45 * t41) - (t44 * t42));
    double t53_bound;
    double t53;
    t53_bound = (2 * (((t38_bound * q1x_bound) + (t37_bound * q1y_bound)) + (t39_bound * q1z_bound)));
    t53 = (2 * (((t38 * q1x) + (t37 * q1y)) + (t39 * q1z)));
    double t34_bound;
    double t34;
    t34_bound = ((t48_bound * t37_bound) + (t47_bound * t38_bound));
    t34 = ((t48 * t37) - (t47 * t38));
    double t33_bound;
    double t33;
    t33_bound = ((t47_bound * t39_bound) + (t46_bound * t37_bound));
    t33 = ((t47 * t39) - (t46 * t37));
    double t32_bound;
    double t32;
    t32_bound = ((t46_bound * t38_bound) + (t48_bound * t39_bound));
    t32 = ((t46 * t38) - (t48 * t39));
    double denom_bound;
    double denom;
    denom_bound = (((t51_bound * t33_bound) + (t50_bound * t32_bound)) + (t49_bound * t34_bound));
    denom = (((t51 * t33) + (t50 * t32)) + (t49 * t34));
    _BOC::_Sign int_tmp_result;
    if (((fabs(denom) > (denom_bound * 2.44249065417534438893e-15)) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0)))
    {
      int_tmp_result = ((denom < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    double p2x_bound = fabs(p2x);
    double p2y_bound = fabs(p2y);
    double p2z_bound = fabs(p2z);
    double result_bound;
    double result;
    result_bound = ((((p1x_bound + p2x_bound) * ((((t33_bound * t54_bound) + (((t49_bound * t37_bound) + (t50_bound * t39_bound)) * t55_bound)) + (((t50_bound * t46_bound) + (t49_bound * t47_bound)) * t53_bound)) + (denom_bound * (p1x_bound + p2x_bound)))) + ((p1y_bound + p2y_bound) * ((((t32_bound * t54_bound) + (((t51_bound * t39_bound) + (t49_bound * t38_bound)) * t55_bound)) + (((t49_bound * t48_bound) + (t51_bound * t46_bound)) * t53_bound)) + (denom_bound * (p1y_bound + p2y_bound))))) + ((p1z_bound + p2z_bound) * ((((t34_bound * t54_bound) + (((t50_bound * t38_bound) + (t51_bound * t37_bound)) * t55_bound)) + (((t51_bound * t47_bound) + (t50_bound * t48_bound)) * t53_bound)) + (denom_bound * (p1z_bound + p2z_bound))))) + (fabs(w1) + fabs(w2)) * denom_bound;
    result = ((((p1x - p2x) * ((((t33 * t54) + (((t49 * t37) - (t50 * t39)) * t55)) + (((t50 * t46) - (t49 * t47)) * t53)) - (denom * (p1x + p2x)))) + ((p1y - p2y) * ((((t32 * t54) + (((t51 * t39) - (t49 * t38)) * t55)) + (((t49 * t48) - (t51 * t46)) * t53)) - (denom * (p1y + p2y))))) + ((p1z - p2z) * ((((t34 * t54) + (((t50 * t38) - (t51 * t37)) * t55)) + (((t51 * t47) - (t50 * t48)) * t53)) - (denom * (p1z + p2z))))) + (w1 - w2) * denom;
    _BOC::_Sign int_tmp_result_FFWKCAA;
    if (((fabs(result) > (result_bound * 4.44089209850062616169e-15)) && (fetestexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID) == 0)))
    {
      int_tmp_result_FFWKCAA = ((result < 0.00000000000000000000e+00) ? _BOC::_Sign::NegativE : _BOC::_Sign::PositivE);
    }
    else
    {
      feclearexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
      return _BOC::_Sign::FaileD;
    }
    return (int_tmp_result == int_tmp_result_FFWKCAA) ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE;
  }

  template <class RT>
  void _Side3D::side3_generic(RT &result,
                              RT &denom,
                              const RT &p1x,
                              const RT &p1y,
                              const RT &p1z,
                              const RT &p2x,
                              const RT &p2y,
                              const RT &p2z,
                              const RT &p3x,
                              const RT &p3y,
                              const RT &p3z,
                              const RT &p4x,
                              const RT &p4y,
                              const RT &p4z,
                              const RT &w1,
                              const RT &w2,
                              const RT &w3,
                              const RT &w4,
                              const RT &q1x,
                              const RT &q1y,
                              const RT &q1z,
                              const RT &q2x,
                              const RT &q2y,
                              const RT &q2z,
                              const RT &q3x,
                              const RT &q3y,
                              const RT &q3z)
  {
    static RT const two(2);

    const RT t46 = p4z - p1z;
    const RT t47 = p4y - p1y;
    const RT t48 = p4x - p1x;
    const RT t55 = t48 * (p1x + p4x) + t47 * (p1y + p4y) + t46 * (p1z + p4z) + w1 - w4;
    const RT t49 = -p1z + p3z;
    const RT t50 = -p1y + p3y;
    const RT t51 = -p1x + p3x;
    const RT t54 = t51 * (p1x + p3x) + t50 * (p1y + p3y) + t49 * (p1z + p3z) + w1 - w3;
    const RT t40 = q3z - q1z;
    const RT t42 = q3x - q1x;
    const RT t43 = q2z - q1z;
    const RT t45 = q2x - q1x;
    const RT t37 = t43 * t42 - t45 * t40;
    const RT t41 = q3y - q1y;
    const RT t44 = q2y - q1y;
    const RT t38 = t44 * t40 - t43 * t41;
    const RT t39 = t45 * t41 - t44 * t42;
    const RT t53 = two * (t38 * q1x + t37 * q1y + t39 * q1z);
    const RT t34 = t48 * t37 - t47 * t38;
    const RT t33 = t47 * t39 - t46 * t37;
    const RT t32 = t46 * t38 - t48 * t39;
    denom = t51 * t33 + t50 * t32 + t49 * t34;
    result =
        (p1x - p2x) * (t33 * t54 + (t49 * t37 - t50 * t39) * t55 + (t50 * t46 - t49 * t47) * t53 - denom * (p1x + p2x)) +
        (p1y - p2y) * (t32 * t54 + (t51 * t39 - t49 * t38) * t55 + (t49 * t48 - t51 * t46) * t53 - denom * (p1y + p2y)) +
        (p1z - p2z) * (t34 * t54 + (t50 * t38 - t51 * t37) * t55 + (t51 * t47 - t50 * t48) * t53 - denom * (p1z + p2z)) +
        (w1 - w2) * denom;
  }

  _Side3D::Comparison_result _Side3D::side3_filtered(double p1x,
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
                                                     double q3z)
  {
    {
      C2F c2f;
      FK::FT result;
      FK::FT denom;
      side3_generic(
          result, denom,
          c2f(p1x), c2f(p1y), c2f(p1z),
          c2f(p2x), c2f(p2y), c2f(p2z),
          c2f(p3x), c2f(p3y), c2f(p3z),
          c2f(p4x), c2f(p4y), c2f(p4z),
          c2f(w1), c2f(w2), c2f(w3), c2f(w4),
          c2f(q1x), c2f(q1y), c2f(q1z),
          c2f(q2x), c2f(q2y), c2f(q2z),
          c2f(q3x), c2f(q3y), c2f(q3z));
      FK::Comparison_result s1 = sign(result);
      FK::Comparison_result s2 = sign(denom);

      if (is_certain(s1) && is_certain(s2))
      {
        return mult(get_certain(s1), get_certain(s2));
      }
    }
    C2E c2e;
    // Exact version (if filtered version failed)
    EK::FT result;
    EK::FT denom;
    side3_generic(
        result, denom,
        c2e(p1x), c2e(p1y), c2e(p1z),
        c2e(p2x), c2e(p2y), c2e(p2z),
        c2e(p3x), c2e(p3y), c2e(p3z),
        c2e(p4x), c2e(p4y), c2e(p4z),
        c2e(w1), c2e(w2), c2e(w3), c2e(w4),
        c2e(q1x), c2e(q1y), c2e(q1z),
        c2e(q2x), c2e(q2y), c2e(q2z),
        c2e(q3x), c2e(q3y), c2e(q3z));
    return mult(sign(result), sign(denom));
  }

  _BOC::_Sign _Side3D::side1_(double p1x,
                              double p1y,
                              double p1z,
                              double p2x,
                              double p2y,
                              double p2z,
                              double w1,
                              double w2,
                              double qx,
                              double qy,
                              double qz)
  {
    _BOC::_Sign result = base_side1_(
        p1x, p1y, p1z,
        p2x, p2y, p2z,
        w1, w2,
        qx, qy, qz);

    if (result == _BOC::_Sign::FaileD)
    {
      Comparison_result r = side1_filtered(
          p1x, p1y, p1z,
          p2x, p2y, p2z,
          w1, w2,
          qx, qy, qz);
      result =
          (r == Comparison_result::ZERO ? _BOC::_Sign::ZerO : (r == Comparison_result::POSITIVE ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE));
    }
    return result;
  }

  _BOC::_Sign _Side3D::side2_(double p1x,
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
                              double q2z)
  {
    _BOC::_Sign result = base_side2_(
        p1x, p1y, p1z,
        p2x, p2y, p2z,
        p3x, p3y, p3z,
        w1, w2, w3,
        q1x, q1y, q1z,
        q2x, q2y, q2z);
    if (result == _BOC::_Sign::FaileD)
    {
      Comparison_result r = side2_filtered(
          p1x, p1y, p1z,
          p2x, p2y, p2z,
          p3x, p3y, p3z,
          w1, w2, w3,
          q1x, q1y, q1z,
          q2x, q2y, q2z);
      result =
          (r == Comparison_result::ZERO ? _BOC::_Sign::ZerO : (r == Comparison_result::POSITIVE ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE));
    }
    return result;
  }

  _BOC::_Sign _Side3D::side3_(double p1x,
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
                              double q3z)
  {
    _BOC::_Sign result = base_side3_(
        p1x, p1y, p1z,
        p2x, p2y, p2z,
        p3x, p3y, p3z,
        p4x, p4y, p4z,
        w1, w2, w3, w4,
        q1x, q1y, q1z,
        q2x, q2y, q2z,
        q3x, q3y, q3z);
    if (result == _BOC::_Sign::FaileD)
    {
      Comparison_result r = side3_filtered(
          p1x, p1y, p1z,
          p2x, p2y, p2z,
          p3x, p3y, p3z,
          p4x, p4y, p4z,
          w1, w2, w3, w4,
          q1x, q1y, q1z,
          q2x, q2y, q2z,
          q3x, q3y, q3z);
      result =
          (r == Comparison_result::ZERO ? _BOC::_Sign::ZerO : (r == Comparison_result::POSITIVE ? _BOC::_Sign::PositivE : _BOC::_Sign::NegativE));
    }
    return result;
  }

} // namespace BGAL
