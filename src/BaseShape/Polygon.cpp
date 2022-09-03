//#include "BKPolygon.h"
//#include "../../include/BK_BaseShape/BKPolygon.h"
#include "BGAL/BaseShape/Polygon.h"

namespace BGAL
{
	_Polygon::_Polygon(const std::vector<_Point2> &in_points)
	{
		_points = in_points;
		if (_points.size() < 3)
			throw std::runtime_error("The number of points is too small!");
		double __maxx, __minx, __maxy, __miny;
		__maxx = __minx = _points[0][0];
		__maxy = __miny = _points[0][1];
		for (int i = 0; i < _points.size(); ++i)
		{
			if (_points[i][0] > __maxx)
				__maxx = _points[i][0];
			if (_points[i][0] < __minx)
				__minx = _points[i][0];
			if (_points[i][1] > __maxy)
				__maxy = _points[i][1];
			if (_points[i][1] < __miny)
				__miny = _points[i][1];
		}
		_bounding_box = std::make_pair(_Point2(__minx, __miny), _Point2(__maxx, __maxy));
	}
	void _Polygon::start_()
	{
		_points.clear();
		_open_in = true;
	}
	void _Polygon::insert_(const _Point2 &in_p)
	{
		if (!_open_in)
			throw std::runtime_error("Cann't insert now!");
		_points.push_back(in_p);
	}
	void _Polygon::insert_(const double &in_x, const double &in_y)
	{
		if (!_open_in)
			throw std::runtime_error("Cann't insert now!");
		_Point2 __ip(in_x, in_y);
		_points.push_back(__ip);
	}
	void _Polygon::end_()
	{
		_open_in = false;
		if (_points.size() < 3)
			throw std::runtime_error("The number of points is too small!");
		double __maxx, __minx, __maxy, __miny;
		__maxx = __minx = _points[0][0];
		__maxy = __miny = _points[0][1];
		for (int i = 0; i < _points.size(); ++i)
		{
			if (_points[i][0] > __maxx)
				__maxx = _points[i][0];
			if (_points[i][0] < __minx)
				__minx = _points[i][0];
			if (_points[i][1] > __maxy)
				__maxy = _points[i][1];
			if (_points[i][1] < __miny)
				__miny = _points[i][1];
		}
		_bounding_box = std::make_pair(_Point2(__minx, __miny), _Point2(__maxx, __maxy));
	}

	const _Point2 &_Polygon::operator[](const int &in_inx) const
	{
		if (in_inx < 0 || in_inx > _points.size() - 1)
			throw std::runtime_error("Beyond the index!");
		return _points[in_inx];
	}
	_Point2 &_Polygon::operator[](const int &in_inx)
	{
		if (in_inx < 0 || in_inx > _points.size() - 1)
			throw std::runtime_error("Beyond the index!");
		return _points[in_inx];
	}
	bool _Polygon::is_in_(const _Point2 &in_p) const
	{
		bool __res = false;
		for (int i = 0; i < _points.size(); ++i)
		{
			const _Point2 &__p1 = _points[i];
			const _Point2 &__p2 = _points[(i + _points.size() - 1) % _points.size()];
			if (_BOC::sign_((__p1 - in_p).cross_(__p2 - in_p).length_()) == _BOC::_Sign::ZerO && (_BOC::sign_((__p1 - in_p).dot_((__p2 - in_p))) == _BOC::_Sign::NegativE || _BOC::sign_((__p1 - in_p).dot_((__p2 - in_p))) == _BOC::_Sign::ZerO))
			{
				return false;
			}
			if (((_BOC::sign_(__p1.y() - in_p.y()) == _BOC::_Sign::PositivE) != (_BOC::sign_(__p2.y() - in_p.y()) == _BOC::_Sign::PositivE)) && _BOC::sign_(in_p.x() - (in_p.y() - __p1.y()) * (__p1.x() - __p2.x()) / (__p1.y() - __p2.y()) - __p1.x()) == _BOC::_Sign::NegativE)
			{
				__res = !__res;
			}
		}
		return __res;
	}
	_Point2 _Polygon::nearest_point_(const _Point2 &in_p)
	{
		double min_dis = (in_p - _points[0]).length_();
		_Point2 min_p = _points[0];
		for (int i = 0; i < num_(); ++i)
		{
			if (_BOC::sign_((_points[(i + 1) % num_()] - _points[i]).length_()) == _BOC::_Sign::ZerO)
			{
				double len = (_points[i] - in_p).length_();
				if (len < min_dis)
				{
					min_dis = len;
					min_p = _points[i];
				}
				continue;
			}
			double len = (_points[(i + 1) % num_()] - _points[i]).dot_(in_p - _points[i]);
			_Point2 lp;
			if (len < 0 || len > (_points[(i + 1) % num_()] - _points[i]).sqlength_())
			{
				double len1 = (in_p - _points[i]).length_();
				double len2 = (in_p - _points[(i + 1) % num_()]).length_();
				len = (len1 < len2 ? len1 : len2);
				lp = (len1 < len2 ? _points[i] : _points[(i + 1) % num_()]);
			}
			else
			{
				lp = _points[i] + (_points[(i + 1) % num_()] - _points[i]) * len / ((_points[(i + 1) % num_()] - _points[i]).sqlength_());
				len = (lp - in_p).length_();
			}
			if (len < min_dis)
			{
				min_dis = len;
				min_p = lp;
			}
		}
		return min_p;
	}
	double _Polygon::distance_to_boundary_(const _Point2 &in_p)
	{
		double min_dis = (in_p - _points[0]).length_();
		for (int i = 0; i < num_(); ++i)
		{
			if (_BOC::sign_((_points[(i + 1) % num_()] - _points[i]).length_()) == _BOC::_Sign::ZerO)
			{
				double len = (_points[i] - in_p).length_();
				if (len < min_dis)
				{
					min_dis = len;
				}
				continue;
			}
			double len = (_points[(i + 1) % num_()] - _points[i]).dot_(in_p - _points[i]);
			if (len < 0 || len > (_points[(i + 1) % num_()] - _points[i]).sqlength_())
			{
				double len1 = (in_p - _points[i]).length_();
				double len2 = (in_p - _points[(i + 1) % num_()]).length_();
				len = (len1 < len2 ? len1 : len2);
			}
			else
			{
				_Point2 lp = _points[i] + (_points[(i + 1) % num_()] - _points[i]) * len / ((_points[(i + 1) % num_()] - _points[i]).sqlength_());
				len = (lp - in_p).length_();
			}
			if (len < min_dis)
			{
				min_dis = len;
			}
		}
		return (is_in_(in_p) ? -min_dis : min_dis);
	}
	int _Polygon::intersection_with_linesegment_(const _Point2 &p1, const _Point2 &p2, std::vector<std::pair<int, _Point2>> &intersections) const
	{
		bool add_p1 = false;
		bool add_p2 = false;
		int p1_inx = -1;
		int p2_inx = -1;
		std::vector<bool> is_on_vertex;
		intersections.resize(0);
		for (int i = 0; i < _points.size(); ++i)
		{
			bool p1_on = false;
			bool p2_on = false;

			if (_BOC::sign_((_points[i] - p1).dot_(_points[(i + 1) % _points.size()] - p1) + ((_points[i] - p1).length_()) * ((_points[(i + 1) % _points.size()] - p1).length_())) == _BOC::_Sign::ZerO)
			{
				p1_on = true;
			}
			if (_BOC::sign_(((_points[i] - p2).dot_(_points[(i + 1) % _points.size()] - p2) + ((_points[i] - p2).length_()) * ((_points[(i + 1) % _points.size()] - p2).length_()))) == _BOC::_Sign::ZerO)
			{
				p2_on = true;
			}

			if (p1_on && p2_on)
			{
				intersections.push_back(std::make_pair(i, p1));
				intersections.push_back(std::make_pair(i, p2));
				return -2;
			}
			if (_BOC::sign_((p1 - p2).cross_(_points[i] - _points[(i + 1) % _points.size()]).length_()) == _BOC::_Sign::ZerO)
			{
				if (p1_on && (!p2_on))
				{
					intersections.resize(0);
					if (p1 == _points[(i + 1) % _points.size()])
					{
						intersections.push_back(std::make_pair((i + 1) % _points.size(), p1));
					}
					else
					{
						intersections.push_back(std::make_pair(i, p1));
					}
					return 1;
				}
				if (p2_on && (!p1_on))
				{
					intersections.resize(0);
					if (p2 == _points[i])
					{
						intersections.push_back(std::make_pair((i + _points.size() - 1) % _points.size(), p2));
					}
					else
					{
						intersections.push_back(std::make_pair(i, p2));
					}
					return 1;
				}
				continue;
			}

			if (p1_on)
			{
				if (intersections.size() == 1)
				{
					break;
				}
				if (add_p2)
				{
					intersections.resize(0);
					return 0;
				}
				add_p1 = true;
				p1_inx = i;
				if (p1 == _points[(i + 1) % _points.size()])
					p1_inx = (i + 1) % _points.size();
				continue;
			}
			if (p2_on)
			{
				if (intersections.size() == 1)
				{
					break;
				}
				if (add_p1)
				{
					intersections.resize(0);
					return 0;
				}
				add_p2 = true;
				p2_inx = i;
				if (p2 == _points[i])
					p2_inx = (i + _points.size() - 1) % _points.size();
				continue;
			}

			if (_BOC::sign_((p1 - _points[i]).dot_(p2 - _points[i]) + ((p1 - _points[i]).length_()) * ((p2 - _points[i]).length_())) == _BOC::_Sign::ZerO)
			{
				intersections.push_back(std::make_pair(i, _points[i]));
				is_on_vertex.push_back(true);
				continue;
			}
			if (_BOC::sign_((p1 - _points[(i + 1) % _points.size()]).dot_(p2 - _points[(i + 1) % _points.size()]) + ((p1 - _points[(i + 1) % _points.size()]).length_()) * ((p2 - _points[(i + 1) % _points.size()]).length_())) == _BOC::_Sign::ZerO)
			{
				continue;
			}

			double __c1 = (_points[i] - _points[(i + 1) % _points.size()]).cross_(_points[i] - p1).z();
			double __c2 = (_points[i] - _points[(i + 1) % _points.size()]).cross_(_points[i] - p2).z();
			double __c3 = (p2 - p1).cross_(p2 - _points[i]).z();
			double __c4 = (p2 - p1).cross_(p2 - _points[(i + 1) % _points.size()]).z();
			if (__c1 * __c2 <= 0 && __c3 * __c4 <= 0)
			{
				_Point2 base = _points[(i + 1) % _points.size()] - _points[i];
				double __d1 = base.cross_((_points[i] - p1)).length_();
				double __d2 = base.cross_((_points[i] - p2)).length_();
				_Point2 __inter = _Point2(p1 + (p2 - p1) * __d1 / (__d1 + __d2));
				intersections.push_back(std::make_pair(i, __inter));
				is_on_vertex.push_back(false);
			}
		}
		if (intersections.size() == 0)
		{
			if (add_p1)
			{
				intersections.push_back(std::make_pair(p1_inx, p1));
				return intersections.size();
			}
			if (add_p2)
			{
				intersections.push_back(std::make_pair(p2_inx, p2));
				return intersections.size();
			}
			return 0;
		}
		if (intersections.size() == 2)
		{
			if ((intersections[0].second - p1).length_() > (intersections[1].second - p1).length_())
			{
				std::pair<int, _Point2> tempp = intersections[0];
				bool tempb = is_on_vertex[0];
				intersections[0] = intersections[1];
				is_on_vertex[0] = is_on_vertex[1];
				intersections[1] = tempp;
				is_on_vertex[1] = tempb;
			}
			if (is_on_vertex[0])
			{
				intersections[0].first = (intersections[0].first + _points.size() - 1) % _points.size();
			}
			return 2;
		}
		if ((!is_in_(p1)) && (!is_in_(p2)))
		{
			intersections.resize(0);
			return 0;
		}
		if (is_in_(p2) && is_on_vertex[0])
		{
			intersections[0].first = (intersections[0].first + _points.size() - 1) % _points.size();
		}
		return intersections.size();
	}
	double _Polygon::area_() const
	{
		double sum = 0;
		std::vector<_Polygon> tris = constrained_delaunay_triangulation_();
		for (int i = 0; i < tris.size(); ++i)
		{
			sum += tris[i].triangle_area_();
		}
		return sum;
	}
	double _Polygon::triangle_area_() const
	{
		_Point2 v1 = _points[1] - _points[0];
		_Point2 v2 = _points[2] - _points[0];
		return v1.cross_(v2).length_() * 0.5;
	}
	double _Polygon::circumference_() const
	{
		double C = 0;
		for (int i = 0; i < num_(); ++i)
		{
			C += (_points[i] - _points[(i + 1) % num_()]).length_();
		}
		return C;
	}
} // namespace BGAL