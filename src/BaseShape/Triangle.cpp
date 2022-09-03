#include "BGAL/BaseShape/Triangle.h"

namespace BGAL
{
  _Triangle3::_Triangle3()
  {
    _area = 0;
    _open_in = false;
  }
  _Triangle3::_Triangle3(const std::vector<_Point3> &in_points)
  {
    if (in_points.size() != 3)
      throw std::runtime_error("The number of points isn't right!");
    _points = in_points;
    _area = -1;
    _open_in = false;
  }
  _Triangle3::_Triangle3(const _Point3 &p1, const _Point3 &p2, const _Point3 &p3)
  {
    _points.resize(0);
    _points.push_back(p1);
    _points.push_back(p2);
    _points.push_back(p3);
    _area = -1;
    _open_in = false;
  }
  void _Triangle3::start_()
  {
    _points.clear();
    _open_in = true;
  }
  void _Triangle3::insert_(const _Point3 &in_p)
  {
    if (!_open_in)
      throw std::runtime_error("Cann't insert now!");
    if (_points.size() == 3)
      throw std::runtime_error("Too many points!");
    _points.push_back(in_p);
  }
  void _Triangle3::insert_(const double &in_x, const double &in_y, const double &in_z)
  {
    if (!_open_in)
      throw std::runtime_error("Cann't insert now!");
    if (_points.size() == 3)
      throw std::runtime_error("Too many points!");
    _Point3 __ip(in_x, in_y, in_z);
    _points.push_back(__ip);
  }
  void _Triangle3::end_()
  {
    _open_in = false;
    _area = -1;
  }
  const _Point3 &_Triangle3::operator[](const int &in_inx) const
  {
    if (in_inx < 0 || in_inx > _points.size() - 1)
      throw std::runtime_error("Beyond the index!");
    return _points[in_inx];
  }
  const _Point3 &_Triangle3::point(const int &in_inx) const
  {
    if (in_inx < 0 || in_inx > _points.size() - 1)
      throw std::runtime_error("Beyond the index!");
    return _points[in_inx];
  }
  double _Triangle3::area_()
  {
    if (_area >= 0)
      return _area;
    _Point3 v1 = _points[1] - _points[0];
    _Point3 v2 = _points[2] - _points[0];
    _area = abs(v1.cross_(v2).length_()) * 0.5;
    return _area;
  }
} // namespace BGAL