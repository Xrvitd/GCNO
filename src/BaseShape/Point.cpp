#include "BGAL/BaseShape/Point.h"
//#include "../../include/BK_BaseShape/BKPoint.h"

namespace BGAL
{
  _Point2::_Point2()
  {
    _value.resize(2);
    _value[0] = _value[1] = 0;
  }
  _Point2::_Point2(const double &in_x, const double &in_y)
  {
    _value.resize(2);
    _value[0] = in_x;
    _value[1] = in_y;
  }
  _Point2::_Point2(const _Point2 &in_p)
  {
    _value.resize(2);
    _value[0] = in_p.x();
    _value[1] = in_p.y();
  }
  std::vector<double> _Point2::get_value_() const
  {
    return _value;
  }
  double _Point2::x() const
  {
    return _value[0];
  }
  double _Point2::y() const
  {
    return _value[1];
  }
  double _Point2::sqlength_() const
  {
    return _value[0] * _value[0] + _value[1] * _value[1];
  }
  double _Point2::length_() const
  {
    return sqrt(sqlength_());
  }
  _Point2 _Point2::normalize_() const
  {
    double __len = length_();
    _Point2 __r(_value[0] / __len, _value[1] / __len);
    return __r;
  }
  _Point2 &_Point2::normalized_()
  {
    double __len = length_();
    _value[0] /= __len;
    _value[1] /= __len;
    return *this;
  }
  double _Point2::operator[](const int &in_inx) const
  {
    if (in_inx > 1 || in_inx < 0)
      throw std::runtime_error("Beyond the index!");
    return _value[in_inx];
  }
  const _Point2 &_Point2::operator=(const _Point2 &in_p)
  {
    _value.resize(2);
    _value[0] = in_p.x();
    _value[1] = in_p.y();
    return *this;
  }
  bool _Point2::operator==(const _Point2 &in_p) const
  {
    return _BOC::sign_(_value[0] - in_p.x()) == _BOC::_Sign::ZerO && _BOC::sign_(_value[1] - in_p.y()) == _BOC::_Sign::ZerO;
  }
  bool _Point2::operator<(const _Point2 &in_p) const
  {
    if (_BOC::sign_(_value[0] - in_p.x()) == _BOC::_Sign::NegativE)
      return true;
    else if (_BOC::sign_(_value[0] - in_p.x()) == _BOC::_Sign::PositivE)
      return false;
    else if (_BOC::sign_(_value[1] - in_p.y()) == _BOC::_Sign::NegativE)
      return true;
    else
      return false;
  }
  bool _Point2::operator<=(const _Point2 &in_p) const
  {
    return this->operator<(in_p) || this->operator==(in_p);
  }
  bool _Point2::operator>(const _Point2 &in_p) const
  {
    return !(this->operator<=(in_p));
  }
  bool _Point2::operator>=(const _Point2 &in_p) const
  {
    return !(this->operator<(in_p));
  }
  bool _Point2::operator!=(const _Point2 &in_p) const
  {
    return !(this->operator==(in_p));
  }
  _Point2 _Point2::operator+(const _Point2 &in_p) const
  {
    _Point2 __r(_value[0] + in_p.x(), _value[1] + in_p.y());
    return __r;
  }
  _Point2 _Point2::operator-(const _Point2 &in_p) const
  {
    _Point2 __r(_value[0] - in_p.x(), _value[1] - in_p.y());
    return __r;
  }
  _Point2 &_Point2::operator+=(const _Point2 &in_p)
  {
    _value[0] += in_p.x();
    _value[1] += in_p.y();
    return *this;
  }
  _Point2 &_Point2::operator-=(const _Point2 &in_p)
  {
    _value[0] -= in_p.x();
    _value[1] -= in_p.y();
    return *this;
  }
  _Point2 _Point2::operator*(const _Point2 &in_p) const
  {
    _Point2 __r(_value[0] * in_p.x(), _value[1] * in_p.y());
    return __r;
  }
  _Point2 _Point2::operator*(const double &in_s) const
  {
    _Point2 __r(_value[0] * in_s, _value[1] * in_s);
    return __r;
  }
  _Point2 _Point2::operator/(const _Point2 &in_p) const
  {
    if (in_p.x() == 0 || in_p.y() == 0)
      throw std::runtime_error("Division by 0!");
    _Point2 __r(_value[0] / in_p.x(), _value[1] / in_p.y());
    return __r;
  }
  _Point2 _Point2::operator/(const double &in_s) const
  {
    if (in_s == 0)
      throw std::runtime_error("Division by 0!");
    _Point2 __r(_value[0] / in_s, _value[1] / in_s);
    return __r;
  }
  _Point2 &_Point2::operator*=(const _Point2 &in_p)
  {
    _value[0] *= in_p.x();
    _value[1] *= in_p.y();
    return *this;
  }
  _Point2 &_Point2::operator*=(const double &in_s)
  {
    _value[0] *= in_s;
    _value[1] *= in_s;
    return *this;
  }
  _Point2 &_Point2::operator/=(const _Point2 &in_p)
  {
    if (in_p.x() == 0 || in_p.y() == 0)
      throw std::runtime_error("Division by 0!");
    _value[0] /= in_p.x();
    _value[1] /= in_p.y();
    return *this;
  }
  _Point2 &_Point2::operator/=(const double &in_s)
  {
    if (in_s == 0)
      throw std::runtime_error("Division by 0!");
    _value[0] /= in_s;
    _value[1] /= in_s;
    return *this;
  }
  double _Point2::dot_(const _Point2 &in_p) const
  {
    return _value[0] * in_p.x() + _value[1] * in_p.y();
  }
  _Point3 _Point2::cross_(const _Point2 &in_p) const
  {
    _Point3 __r(0, 0, _value[0] * in_p.y() - _value[1] * in_p.x());
    return __r;
  }
  std::ostream &operator<<(std::ostream &in_os, _Point2 in_p)
  {
    in_os << in_p[0] << " " << in_p[1];
    return in_os;
  }
  std::ostream &operator<<(std::ostream &in_os, const _Point3 &in_p)
  {
    in_os << in_p.x() << " " << in_p.y() << " " << in_p.z();
    return in_os;
  }
  _Point3 operator*(const double &in_s, const _Point3 &in_p)
  {
    return in_p * in_s;
  }
  _Point3::_Point3()
  {
    _value.resize(3);
    _value[0] = _value[1] = _value[2] = 0;
  }
  _Point3::_Point3(const double &in_x, const double &in_y, const double &in_z)
  {
    _value.resize(3);
    _value[0] = in_x;
    _value[1] = in_y;
    _value[2] = in_z;
  }
  _Point3::_Point3(const _Point3 &in_p)
  {
    _value.resize(3);
    _value[0] = in_p.x();
    _value[1] = in_p.y();
    _value[2] = in_p.z();
  }
  std::vector<double> _Point3::get_value_() const
  {
    return _value;
  }
  double _Point3::x() const
  {
    return _value[0];
  }
  double _Point3::y() const
  {
    return _value[1];
  }
  double _Point3::z() const
  {
    return _value[2];
  }
  double _Point3::sqlength_() const
  {
    return _value[0] * _value[0] + _value[1] * _value[1] + _value[2] * _value[2];
  }
  double _Point3::length_() const
  {
    return sqrt(sqlength_());
  }
  _Point3 _Point3::normalize_() const
  {
    double __len = length_();
    _Point3 __r(_value[0] / __len, _value[1] / __len, _value[2] / __len);
    return __r;
  }
  _Point3 &_Point3::normalized_()
  {
    double __len = length_();
    _value[0] /= __len;
    _value[1] /= __len;
    _value[2] /= __len;
    return *this;
  }
  double _Point3::operator[](const int in_inx) const
  {
    if (in_inx > 2 || in_inx < 0)
      throw std::runtime_error("Beyond the index!");
    return _value[in_inx];
  }
  double &_Point3::operator[](const int in_inx)
  {
    if (in_inx > 2 || in_inx < 0)
      throw std::runtime_error("Beyond the index!");
    return _value[in_inx];
  }
  _Point3 _Point3::rotate_(const Eigen::Matrix3d &RM) const
  {
    return _Point3(RM(0, 0) * x() + RM(0, 1) * y() + RM(0, 2) * z(),
                   RM(1, 0) * x() + RM(1, 1) * y() + RM(1, 2) * z(),
                   RM(2, 0) * x() + RM(2, 1) * y() + RM(2, 2) * z());
  }
  const _Point3 &_Point3::operator=(const _Point3 &in_p)
  {
    _value.resize(3);
    _value[0] = in_p.x();
    _value[1] = in_p.y();
    _value[2] = in_p.z();
    return *this;
  }
  bool _Point3::operator==(const _Point3 &in_p) const
  {
    return _BOC::sign_(_value[0] - in_p.x()) == _BOC::_Sign::ZerO && _BOC::sign_(_value[1] - in_p.y()) == _BOC::_Sign::ZerO && _BOC::sign_(_value[2] - in_p.z()) == _BOC::_Sign::ZerO;
  }
  bool _Point3::operator<(const _Point3 &in_p) const
  {
    if (_BOC::sign_(_value[0] - in_p.x()) == _BOC::_Sign::NegativE)
      return true;
    else if (_BOC::sign_(_value[0] - in_p.x()) == _BOC::_Sign::PositivE)
      return false;
    else if (_BOC::sign_(_value[1] - in_p.y()) == _BOC::_Sign::NegativE)
      return true;
    else if (_BOC::sign_(_value[1] - in_p.y()) == _BOC::_Sign::PositivE)
      return false;
    else if (_BOC::sign_(_value[2] - in_p.z()) == _BOC::_Sign::NegativE)
      return true;
    else
      return false;
  }
  bool _Point3::operator<=(const _Point3 &in_p) const
  {
    return this->operator<(in_p) || this->operator==(in_p);
  }
  bool _Point3::operator>(const _Point3 &in_p) const
  {
    return !(this->operator<=(in_p));
  }
  bool _Point3::operator>=(const _Point3 &in_p) const
  {
    return !(this->operator<(in_p));
  }
  bool _Point3::operator!=(const _Point3 &in_p) const
  {
    return !(this->operator==(in_p));
  }
  _Point3 _Point3::operator+(const _Point3 &in_p) const
  {
    _Point3 __r(_value[0] + in_p.x(), _value[1] + in_p.y(), _value[2] + in_p.z());
    return __r;
  }
  _Point3 _Point3::operator-(const _Point3 &in_p) const
  {
    _Point3 __r(_value[0] - in_p.x(), _value[1] - in_p.y(), _value[2] - in_p.z());
    return __r;
  }
  _Point3 _Point3::operator-() const
  {
    _Point3 __r(-_value[0], _value[1], -_value[2]);
    return __r;
  }
  _Point3 &_Point3::operator+=(const _Point3 &in_p)
  {
    _value[0] += in_p.x();
    _value[1] += in_p.y();
    _value[2] += in_p.z();
    return *this;
  }
  _Point3 &_Point3::operator-=(const _Point3 &in_p)
  {
    _value[0] -= in_p.x();
    _value[1] -= in_p.y();
    _value[2] -= in_p.z();
    return *this;
  }
  _Point3 _Point3::operator*(const _Point3 &in_p) const
  {
    _Point3 __r(_value[0] * in_p.x(), _value[1] * in_p.y(), _value[2] * in_p.z());
    return __r;
  }
  _Point3 _Point3::operator*(const double &in_s) const
  {
    _Point3 __r(_value[0] * in_s, _value[1] * in_s, _value[2] * in_s);
    return __r;
  }
  _Point3 _Point3::operator/(const _Point3 &in_p) const
  {
    if (in_p.x() == 0 || in_p.y() == 0 || in_p.z() == 0)
      throw std::runtime_error("Division by 0!");
    _Point3 __r(_value[0] / in_p.x(), _value[1] / in_p.y(), _value[2] / in_p.z());
    return __r;
  }
  _Point3 _Point3::operator/(const double &in_s) const
  {
    if (in_s == 0)
      throw std::runtime_error("Division by 0!");
    _Point3 __r(_value[0] / in_s, _value[1] / in_s, _value[2] / in_s);
    return __r;
  }
  _Point3 &_Point3::operator*=(const _Point3 &in_p)
  {
    _value[0] *= in_p.x();
    _value[1] *= in_p.y();
    _value[2] *= in_p.z();
    return *this;
  }
  _Point3 &_Point3::operator*=(const double &in_s)
  {
    _value[0] *= in_s;
    _value[1] *= in_s;
    _value[2] *= in_s;
    return *this;
  }
  _Point3 &_Point3::operator/=(const _Point3 &in_p)
  {
    if (in_p.x() == 0 || in_p.y() == 0 || in_p.z() == 0)
      throw std::runtime_error("Division by 0!");
    _value[0] /= in_p.x();
    _value[1] /= in_p.y();
    _value[2] /= in_p.z();
    return *this;
  }
  _Point3 &_Point3::operator/=(const double &in_s)
  {
    if (in_s == 0)
      throw std::runtime_error("Division by 0!");
    _value[0] /= in_s;
    _value[1] /= in_s;
    _value[2] /= in_s;
    return *this;
  }
  double _Point3::dot_(const _Point3 &in_p) const
  {
    return _value[0] * in_p.x() + _value[1] * in_p.y() + _value[2] * in_p.z();
  }
  _Point3 _Point3::cross_(const _Point3 &in_p) const
  {
    _Point3 __r(_value[1] * in_p.z() - _value[2] * in_p.y(),
                _value[2] * in_p.x() - _value[0] * in_p.z(),
                _value[0] * in_p.y() - _value[1] * in_p.x());
    return __r;
  }
} // namespace BGAL