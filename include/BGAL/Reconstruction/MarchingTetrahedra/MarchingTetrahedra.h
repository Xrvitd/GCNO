#pragma once
#include "BGAL/BaseShape/Point.h"
#include "BGAL/Model/ManifoldModel.h"
#include <tuple>
#include <functional>
#include <map>
#include <fstream>
namespace BGAL
{
  class _Marching_Tetrahedra
  {
  public:
    _Marching_Tetrahedra();
    _Marching_Tetrahedra(const std::pair<_Point3, _Point3> &in_boundingbox, const int &in_depth);
    template <class F>
    _ManifoldModel reconstruction_(const F &sign);
    void set_method_(const int &m)
    {
      _method = m;
    }

  private:
    void tiling_();

  private:
    std::pair<_Point3, _Point3> _boundingbox;
    int _depth;
    std::vector<_Point3> _tetra_vertices;
    std::vector<std::tuple<int, int, int, int>> _tetras;
    int _method;
  };
  template <class F>
  _ManifoldModel _Marching_Tetrahedra::reconstruction_(const F &sign)
  {
    std::vector<double> signs(_tetra_vertices.size());
    for (int i = 0; i < _tetra_vertices.size(); ++i)
    {
      signs[i] = sign(_tetra_vertices[i]).first;
    }
    std::vector<_Point3> _vertices;
    std::map<std::pair<int, int>, int> from_edge_to_vertex;
    std::vector<std::tuple<int, int, int>> tris;
    for (int i = 0; i < _tetras.size(); ++i)
    {
      std::vector<int> negetive;
      std::vector<int> positive;
      if (signs[std::get<0>(_tetras[i])] < 0)
        negetive.push_back(std::get<0>(_tetras[i]));
      else
        positive.push_back(std::get<0>(_tetras[i]));
      if (signs[std::get<1>(_tetras[i])] < 0)
        negetive.push_back(std::get<1>(_tetras[i]));
      else
        positive.push_back(std::get<1>(_tetras[i]));
      if (signs[std::get<2>(_tetras[i])] < 0)
        negetive.push_back(std::get<2>(_tetras[i]));
      else
        positive.push_back(std::get<2>(_tetras[i]));
      if (signs[std::get<3>(_tetras[i])] < 0)
        negetive.push_back(std::get<3>(_tetras[i]));
      else
        positive.push_back(std::get<3>(_tetras[i]));
      if (negetive.size() == 1)
      {
        std::vector<int> tri_v;
        int v0, v1;
        for (int j = 0; j < 3; ++j)
        {
          if (negetive[0] < positive[j])
          {
            v0 = negetive[0];
            v1 = positive[j];
          }
          else
          {
            v0 = positive[j];
            v1 = negetive[0];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
        }
        if (negetive[0] == std::get<3>(_tetras[i]) || negetive[0] == std::get<1>(_tetras[i]))
        {
          tris.push_back(std::make_tuple(tri_v[2], tri_v[1], tri_v[0]));
        }
        else
        {
          tris.push_back(std::make_tuple(tri_v[0], tri_v[1], tri_v[2]));
        }
      }
      else if (negetive.size() == 3)
      {
        std::vector<int> tri_v;
        int v0, v1;
        for (int j = 0; j < 3; ++j)
        {
          if (negetive[j] < positive[0])
          {
            v0 = negetive[j];
            v1 = positive[0];
          }
          else
          {
            v0 = positive[0];
            v1 = negetive[j];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
        }
        if (positive[0] == std::get<3>(_tetras[i]) || positive[0] == std::get<1>(_tetras[i]))
        {
          tris.push_back(std::make_tuple(tri_v[0], tri_v[1], tri_v[2]));
        }
        else
        {
          tris.push_back(std::make_tuple(tri_v[2], tri_v[1], tri_v[0]));
        }
      }
      else if (negetive.size() == 2)
      {
        std::vector<int> tri_v;
        if (negetive[1] == std::get<3>(_tetras[i]))
        {
          int v0, v1;
          if (negetive[1] < positive[0])
          {
            v0 = negetive[1];
            v1 = positive[0];
          }
          else
          {
            v0 = positive[0];
            v1 = negetive[1];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[0] < positive[0])
          {
            v0 = negetive[0];
            v1 = positive[0];
          }
          else
          {
            v0 = positive[0];
            v1 = negetive[0];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[0] < positive[1])
          {
            v0 = negetive[0];
            v1 = positive[1];
          }
          else
          {
            v0 = positive[1];
            v1 = negetive[0];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[1] < positive[1])
          {
            v0 = negetive[1];
            v1 = positive[1];
          }
          else
          {
            v0 = positive[1];
            v1 = negetive[1];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[0] == std::get<1>(_tetras[i]))
          {
            tris.push_back(std::make_tuple(tri_v[3], tri_v[2], tri_v[1]));
            tris.push_back(std::make_tuple(tri_v[3], tri_v[1], tri_v[0]));
          }
          else
          {
            tris.push_back(std::make_tuple(tri_v[0], tri_v[1], tri_v[2]));
            tris.push_back(std::make_tuple(tri_v[0], tri_v[2], tri_v[3]));
          }
        }
        else
        {
          int v0, v1;
          if (negetive[0] < positive[1])
          {
            v0 = negetive[0];
            v1 = positive[1];
          }
          else
          {
            v0 = positive[1];
            v1 = negetive[0];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[0] < positive[0])
          {
            v0 = negetive[0];
            v1 = positive[0];
          }
          else
          {
            v0 = positive[0];
            v1 = negetive[0];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[1] < positive[0])
          {
            v0 = negetive[1];
            v1 = positive[0];
          }
          else
          {
            v0 = positive[0];
            v1 = negetive[1];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (negetive[1] < positive[1])
          {
            v0 = negetive[1];
            v1 = positive[1];
          }
          else
          {
            v0 = positive[1];
            v1 = negetive[1];
          }
          if (from_edge_to_vertex.find(std::make_pair(v0, v1)) == from_edge_to_vertex.end())
          {
            from_edge_to_vertex[std::make_pair(v0, v1)] = _vertices.size();
            switch (_method)
            {
            case 0:
            {
              _vertices.push_back(sign((_tetra_vertices[v0] + _tetra_vertices[v1]) * 0.5).second);
              break;
            }
            case 1:
            {
              double len = fabs(signs[v0]) + fabs(signs[v1]);
              _vertices.push_back(
                  _tetra_vertices[v0] * fabs(signs[v1]) / len + _tetra_vertices[v1] * fabs(signs[v0] / len));
              break;
            }
            default:
              break;
            }
          }
          tri_v.push_back(from_edge_to_vertex[std::make_pair(v0, v1)]);
          if (positive[0] == std::get<1>(_tetras[i]))
          {
            tris.push_back(std::make_tuple(tri_v[0], tri_v[1], tri_v[2]));
            tris.push_back(std::make_tuple(tri_v[0], tri_v[2], tri_v[3]));
          }
          else
          {
            tris.push_back(std::make_tuple(tri_v[3], tri_v[2], tri_v[1]));
            tris.push_back(std::make_tuple(tri_v[3], tri_v[1], tri_v[0]));
          }
        }
      }
    }
    std::vector<_Model::_MFace> faces;
    for (int i = 0; i < tris.size(); ++i)
    {
      faces.push_back(_Model::_MFace(std::get<0>(tris[i]), std::get<1>(tris[i]), std::get<2>(tris[i])));
    }
    _ManifoldModel res_model(_vertices, faces);
    return res_model;
  }
} // namespace BGAL