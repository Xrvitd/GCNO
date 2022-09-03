#include "BGAL/Reconstruction/MarchingTetrahedra/MarchingTetrahedra.h"
#include <fstream>
namespace BGAL
{
  _Marching_Tetrahedra::_Marching_Tetrahedra()
      : _depth(1), _method(0)
  {
  }
  _Marching_Tetrahedra::_Marching_Tetrahedra(const std::pair<_Point3, _Point3> &in_boundingbox, const int &in_depth)
      : _boundingbox(in_boundingbox), _depth(in_depth), _method(0)
  {
    tiling_();
  }
  void _Marching_Tetrahedra::tiling_()
  {
    double box = 1.1 * std::max(_boundingbox.second.x() - _boundingbox.first.x(), _boundingbox.second.y() - _boundingbox.first.y());
    double l = 0.5 * box * (1 / sqrt(3.0) + 1.0);
    double e = l / pow(2, _depth - 1);
    double a = e * 0.25 * sqrt(2.0);
    _Point3 lu(-box * 0.25 * (1.0 / sqrt(3.0) + 1.0),
               box * 0.25 * (1.0 + sqrt(3.0)),
               _boundingbox.first.z() - 0.1 * (_boundingbox.second.z() - _boundingbox.first.z()) - 3 * a);
    int row = pow(2, _depth - 1);
    std::vector<std::vector<_Point3>> xy_vertices(row * 2 + 1);
    double delta_x = e * 0.5;
    double delta_y = e * 0.5 * sqrt(3.0);
    for (int i = 0; i < row; ++i)
    {
      for (int j = 0; j < row + 1 + i; ++j)
      {
        _Point3 p(lu.x() - delta_x * i + e * j, lu.y() - delta_y * i, lu.z());
        xy_vertices[i].push_back(p);
        _Point3 p2(p.x(), -p.y(), p.z());
        xy_vertices[row * 2 - i].push_back(p2);
      }
    }
    for (int i = 0; i < 2 * row + 1; ++i)
    {
      _Point3 p(-box * 0.5 * (1.0 / sqrt(3.0) + 1) + i * e, 0, lu.z());
      xy_vertices[row].push_back(p);
    }
    int height =
        (int)(((_boundingbox.second.z() + 0.1 * (_boundingbox.second.z() - _boundingbox.first.z())) - lu.z() - 2.0 * a) / (3.0 * a)) + 1;
    std::map<std::tuple<int, int, int>, int> from_sym_to_vertex;
    _tetra_vertices.clear();
    _tetras.clear();

    for (int i = 0; i < row; ++i)
    {
      for (int j = 0; j < row + i; ++j)
      {
        std::pair<int, int> p0, p1, p2;
        int l0 = (i + j) % 3;
        int l1 = (i + j + 1) % 3;
        int l2 = (i + j + 2) % 3;
        if (l0 == 0)
        {
          if (l1 == 1)
          {
            p0 = std::pair<int, int>(i, j);
            p1 = std::pair<int, int>(i, j + 1);
            p2 = std::pair<int, int>(i + 1, j + 1);
          }
          else
          {
            p0 = std::pair<int, int>(i, j);
            p2 = std::pair<int, int>(i, j + 1);
            p1 = std::pair<int, int>(i + 1, j + 1);
          }
        }
        else if (l1 == 0)
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(i, j);
            p0 = std::pair<int, int>(i, j + 1);
            p2 = std::pair<int, int>(i + 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(i, j);
            p0 = std::pair<int, int>(i, j + 1);
            p1 = std::pair<int, int>(i + 1, j + 1);
          }
        }
        else
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(i, j);
            p2 = std::pair<int, int>(i, j + 1);
            p0 = std::pair<int, int>(i + 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(i, j);
            p1 = std::pair<int, int>(i, j + 1);
            p0 = std::pair<int, int>(i + 1, j + 1);
          }
        }
        for (int h = 0; h < height; ++h)
        {
          std::tuple<int, int, int> v0 = std::make_tuple(p0.first, p0.second, h * 3);
          std::tuple<int, int, int> v1 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          std::tuple<int, int, int> v2 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          std::tuple<int, int, int> v3 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          v1 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          v2 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v3 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v1 = std::make_tuple(p2.first, p2.second, h * 3 + 5);
          v2 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          v3 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v3]));
        }
        if (l0 == 0)
        {
          if (l1 == 1)
          {
            p0 = std::pair<int, int>(2 * row - i, j);
            p1 = std::pair<int, int>(2 * row - i, j + 1);
            p2 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
          else
          {
            p0 = std::pair<int, int>(2 * row - i, j);
            p2 = std::pair<int, int>(2 * row - i, j + 1);
            p1 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
        }
        else if (l1 == 0)
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(2 * row - i, j);
            p0 = std::pair<int, int>(2 * row - i, j + 1);
            p2 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(2 * row - i, j);
            p0 = std::pair<int, int>(2 * row - i, j + 1);
            p1 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
        }
        else
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(2 * row - i, j);
            p2 = std::pair<int, int>(2 * row - i, j + 1);
            p0 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(2 * row - i, j);
            p1 = std::pair<int, int>(2 * row - i, j + 1);
            p0 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
        }
        for (int h = 0; h < height; ++h)
        {
          std::tuple<int, int, int> v0 = std::make_tuple(p0.first, p0.second, h * 3);
          std::tuple<int, int, int> v1 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          std::tuple<int, int, int> v2 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          std::tuple<int, int, int> v3 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          v1 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          v2 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v3 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v1 = std::make_tuple(p2.first, p2.second, h * 3 + 5);
          v2 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          v3 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v3]));
        }
      }
      for (int j = 0; j < row + i + 1; ++j)
      {
        std::pair<int, int> p0, p1, p2;
        int l0 = (i + j) % 3;
        int l1 = (i + j + 1) % 3;
        int l2 = (i + j + 2) % 3;
        if (l0 == 0)
        {
          if (l1 == 1)
          {
            p0 = std::pair<int, int>(i, j);
            p1 = std::pair<int, int>(i + 1, j);
            p2 = std::pair<int, int>(i + 1, j + 1);
          }
          else
          {
            p0 = std::pair<int, int>(i, j);
            p2 = std::pair<int, int>(i + 1, j);
            p1 = std::pair<int, int>(i + 1, j + 1);
          }
        }
        else if (l1 == 0)
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(i, j);
            p0 = std::pair<int, int>(i + 1, j);
            p2 = std::pair<int, int>(i + 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(i, j);
            p0 = std::pair<int, int>(i + 1, j);
            p1 = std::pair<int, int>(i + 1, j + 1);
          }
        }
        else
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(i, j);
            p2 = std::pair<int, int>(i + 1, j);
            p0 = std::pair<int, int>(i + 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(i, j);
            p1 = std::pair<int, int>(i + 1, j);
            p0 = std::pair<int, int>(i + 1, j + 1);
          }
        }
        for (int h = 0; h < height; ++h)
        {
          std::tuple<int, int, int> v0 = std::make_tuple(p0.first, p0.second, h * 3);
          std::tuple<int, int, int> v1 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          std::tuple<int, int, int> v2 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          std::tuple<int, int, int> v3 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          v1 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          v2 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v3 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v1 = std::make_tuple(p2.first, p2.second, h * 3 + 5);
          v2 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          v3 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v3]));
        }
        if (l0 == 0)
        {
          if (l1 == 1)
          {
            p0 = std::pair<int, int>(2 * row - i, j);
            p1 = std::pair<int, int>(2 * row - i - 1, j);
            p2 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
          else
          {
            p0 = std::pair<int, int>(2 * row - i, j);
            p2 = std::pair<int, int>(2 * row - i - 1, j);
            p1 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
        }
        else if (l1 == 0)
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(2 * row - i, j);
            p0 = std::pair<int, int>(2 * row - i - 1, j);
            p2 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(2 * row - i, j);
            p0 = std::pair<int, int>(2 * row - i - 1, j);
            p1 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
        }
        else
        {
          if (l0 == 1)
          {
            p1 = std::pair<int, int>(2 * row - i, j);
            p2 = std::pair<int, int>(2 * row - i - 1, j);
            p0 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
          else
          {
            p2 = std::pair<int, int>(2 * row - i, j);
            p1 = std::pair<int, int>(2 * row - i - 1, j);
            p0 = std::pair<int, int>(2 * row - i - 1, j + 1);
          }
        }
        for (int h = 0; h < height; ++h)
        {
          std::tuple<int, int, int> v0 = std::make_tuple(p0.first, p0.second, h * 3);
          std::tuple<int, int, int> v1 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          std::tuple<int, int, int> v2 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          std::tuple<int, int, int> v3 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p1.first, p1.second, h * 3 + 1);
          v1 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          v2 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v3 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v3]));
          v0 = std::make_tuple(p2.first, p2.second, h * 3 + 2);
          v1 = std::make_tuple(p2.first, p2.second, h * 3 + 5);
          v2 = std::make_tuple(p0.first, p0.second, h * 3 + 3);
          v3 = std::make_tuple(p1.first, p1.second, h * 3 + 4);
          if (from_sym_to_vertex.find(v0) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v0)][std::get<1>(v0)] + _Point3(0, 0, a) * std::get<2>(v0);
            from_sym_to_vertex[v0] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v1) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v1)][std::get<1>(v1)] + _Point3(0, 0, a) * std::get<2>(v1);
            from_sym_to_vertex[v1] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v2) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v2)][std::get<1>(v2)] + _Point3(0, 0, a) * std::get<2>(v2);
            from_sym_to_vertex[v2] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          if (from_sym_to_vertex.find(v3) == from_sym_to_vertex.end())
          {
            _Point3 p = xy_vertices[std::get<0>(v3)][std::get<1>(v3)] + _Point3(0, 0, a) * std::get<2>(v3);
            from_sym_to_vertex[v3] = _tetra_vertices.size();
            _tetra_vertices.push_back(p);
          }
          _tetras.push_back(std::make_tuple(from_sym_to_vertex[v0],
                                            from_sym_to_vertex[v2],
                                            from_sym_to_vertex[v1],
                                            from_sym_to_vertex[v3]));
        }
      }
    }
  }
} // namespace BGAL