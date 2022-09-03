#pragma once

#include <vector>
#include <queue>
#include <map>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
#include "BGAL/BaseShape/Triangle.h"
#include "BGAL/BaseShape/Line.h"
#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"

namespace BGAL
{
  class _Tessellation3D_Skeleton
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::FT Weight;
    typedef K::Point_3 Point;
    typedef K::Weighted_point_3 Weighted_point;
    typedef CGAL::Regular_triangulation_vertex_base_3<K> Vb0;
    typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
    typedef CGAL::Regular_triangulation_cell_base_3<K> Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
    typedef CGAL::Regular_triangulation_3<K, Tds> Rt;

  private:
    std::vector<std::set<int>> _neights;

  public:
    _Tessellation3D_Skeleton();
    _Tessellation3D_Skeleton(Rt &rt, const int &num_vertices);
    const std::set<int> &neight_(const int &in_i) const
    {
      return _neights[in_i];
    }
  };
  class _Restricted_Tessellation3D
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::FT Weight;
    typedef K::Point_3 Point;
    typedef K::Weighted_point_3 Weighted_point;
    typedef CGAL::Regular_triangulation_vertex_base_3<K> Vb0;
    typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
    typedef CGAL::Regular_triangulation_cell_base_3<K> Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
    typedef CGAL::Regular_triangulation_3<K, Tds> Rt;

  private:
    int _num_sites;
    std::vector<_Point3> _sites;
    std::vector<double> _weights;
    _ManifoldModel _model;
    std::vector<_Point3> _vertices;
    std::vector<std::vector<std::tuple<int, int, int>>> _cells;
    std::vector<std::map<int, std::vector<std::pair<int, int>>>> _edges;
    _Tessellation3D_Skeleton _skeleton;
    std::vector<bool> _is_hidden;

  private:
    struct _Symbolic_Point
    {
      std::set<int> _sym;
      int _site;
      std::vector<int> p;
      int flag;
      _Symbolic_Point() : flag(0)
      {
        _sym.clear();
        p.clear();
        _site = -1;
      }
      _Symbolic_Point(const std::set<int> &in_sym, const int &in_site)
          : _sym(in_sym), _site(in_site), flag(0)
      {
        p.clear();
      }
      _Symbolic_Point(const int &in1, const int &in2, const int &in3, const int &in_site)
          : _site(in_site), flag(0)
      {
        _sym.clear();
        _sym.insert(in1);
        _sym.insert(in2);
        _sym.insert(in3);
        p.clear();
      }
      void insert_(const int &in)
      {
        _sym.insert(in);
      }
      std::set<int> insec_(const _Symbolic_Point &in_sym) const
      {
        std::set<int> res;
        for (auto it = _sym.begin(); it != _sym.end(); ++it)
        {
          if (in_sym._sym.find(*it) != in_sym._sym.end())
          {
            res.insert(*it);
          }
        }
        return res;
      }
      int num_sites_() const
      {
        int num = 0;
        for (auto it = _sym.begin(); it != _sym.end(); ++it)
        {
          if (*it > 0)
            ++num;
        }
        return num;
      }
      void update_(const int &in_site, const _Symbolic_Point &adj_sym)
      {
        _sym = insec_(adj_sym);
        _sym.insert(in_site + 1);
      }
    };
    _BOC::_Sign side_(const int &ip1, const int &ip2, _Symbolic_Point &v);
    _Symbolic_Point insec_bisector_(const _Symbolic_Point &p1,
                                    const _Symbolic_Point &p2,
                                    const int &neigh,
                                    const int &center,
                                    const std::vector<int> &current_face_vid) const;
    void calculate_();

  public:
    _Restricted_Tessellation3D(const _ManifoldModel& in_model);
    _Restricted_Tessellation3D(const _ManifoldModel &in_model,
                               const std::vector<_Point3> &in_sites,
                               const std::vector<double> &in_weights);
    _Restricted_Tessellation3D(const _ManifoldModel &in_model, const std::vector<_Point3> &in_sites);
    void calculate_(const std::vector<_Point3> &in_sites);
    void calculate_(const std::vector<_Point3>& in_sites, const std::vector<double>& in_weights);
    int number_hidden_point_() const
    {
        int n = 0;
        for (int i = 0; i < _is_hidden.size(); ++i)
        {
            if (_is_hidden[i])
                n++;
        }
        return n;
    }
    int number_vertices_() const
    {
      return _vertices.size();
    }
    _Point3 vertex_(const int &id) const
    {
      if (id < 0 || id >= _vertices.size())
        throw std::runtime_error("Beyond the index!");
      return _vertices[id];
    }
    
    const  std::vector<_Point3>& get_sites_() const
    {
        return _sites;
    }
    const std::vector<std::vector<std::tuple<int, int, int>>> &get_cells_() const
    {
      return _cells;
    }
    const std::vector<std::map<int, std::vector<std::pair<int, int>>>> &get_edges_() const
    {
      return _edges;
    }
  };
} // namespace BGAL