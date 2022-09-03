#include "BGAL/Tessellation2D/Side2D.h"
#include "BGAL/Tessellation2D/Tessellation2D.h"

//std
#include <queue>

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_vertex_base_2<K> Vbase;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K, Vbase> Vb;
typedef CGAL::Regular_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Regular_triangulation_2<K, Tds> Regular;
typedef K::Point_2 Point;
typedef K::Weighted_point_2 Wpoint;
typedef Regular::Vertex_handle Vertex_handle;

namespace BGAL
{

  _Tessellation2D::_Tessellation2D(const _Polygon &in_boundary, const std::vector<_Point2> &in_sites)
      : _boundary(in_boundary), _sites(in_sites)
  {
    _num_sites = in_sites.size();
    _weights.resize(_num_sites, 0);
    calculate_();
  }

  _Tessellation2D::_Tessellation2D(const _Polygon &in_boundary,
                                   const std::vector<_Point2> &in_sites,
                                   const std::vector<double> &in_weights)
      : _boundary(in_boundary), _sites(in_sites), _weights(in_weights)
  {
    _num_sites = in_sites.size();
    calculate_();
  }

  void _Tessellation2D::calculate_(const _Polygon &in_boundary, const std::vector<_Point2> &in_sites)
  {
    _boundary = in_boundary;
    _sites = in_sites;
    _num_sites = in_sites.size();
    _weights.clear();
    _weights.resize(_num_sites, 0);
    calculate_();
  }

  void _Tessellation2D::calculate_(const _Polygon &in_boundary,
                                   const std::vector<_Point2> &in_sites,
                                   const std::vector<double> &in_weights)
  {
    _boundary = in_boundary;
    _sites = in_sites;
    _num_sites = in_sites.size();
    _weights = in_weights;
    calculate_();
  }

  //std::vector<BKHao::_Polygon> BKHao::_Tessellation2D::get_cell_polygons_() const
  //{
  //	return _cell_polygons;
  //}

  const std::vector<_Polygon> &_Tessellation2D::get_cell_polygons_() const
  {
    return _cell_polygons;
  }

  const std::vector<std::vector<std::pair<int, int>>> &_Tessellation2D::get_cells_() const
  {
    return _cells;
  }

  const _Point2 &_Tessellation2D::vertex_(const int &id) const
  {
    if (id < 0 || id >= _tessellation_vertices.size())
      throw std::runtime_error("Beyond the index!");
    return _tessellation_vertices[id];
  }

  int _Tessellation2D::num_hidden_() const
  {
    int count = 0;
    for (int i = 0; i < _num_sites; ++i)
    {
      if (_is_hidden[i])
        ++count;
    }
    return count;
  }
  
  int _Tessellation2D::num_vertex_() const
  {
      return _tessellation_vertices.size();
  }
  _BOC::_Sign _Tessellation2D::symbolic_point_site_(const int &site1,
                                                    const int &site2,
                                                    const std::pair<int, int> &symp,
                                                    const int &cur_edge)
  {
    _BOC::_Sign r;
    //std::cout << site1 << " " << site2 << " " << symp.first << "  " << symp.second << std::endl;
    if (symp.first < 0 && symp.second < 0)
    {
      if (symp.first - symp.second == -1)
        r = _Side2D::side1_(_sites[site1].x(),
                            _sites[site1].y(),
                            _weights[site1],
                            _sites[site2].x(),
                            _sites[site2].y(),
                            _weights[site2],
                            _boundary[-symp.first - 1].x(),
                            _boundary[-symp.first - 1].y());
      else if (symp.first - symp.second == 1)
        r = _Side2D::side1_(_sites[site1].x(),
                            _sites[site1].y(),
                            _weights[site1],
                            _sites[site2].x(),
                            _sites[site2].y(),
                            _weights[site2],
                            _boundary[-symp.second - 1].x(),
                            _boundary[-symp.second - 1].y());
      else
        r = _Side2D::side1_(_sites[site1].x(),
                            _sites[site1].y(),
                            _weights[site1],
                            _sites[site2].x(),
                            _sites[site2].y(),
                            _weights[site2],
                            _boundary[0].x(),
                            _boundary[0].y());
    }
    else
    {
      int site3 = symp.first > 0 ? symp.first - 1 : symp.second - 1;
      r = _Side2D::side2_(_sites[site1].x(),
                          _sites[site1].y(),
                          _weights[site1],
                          _sites[site2].x(),
                          _sites[site2].y(),
                          _weights[site2],
                          _sites[site3].x(),
                          _sites[site3].y(),
                          _weights[site3],
                          _boundary[cur_edge].x(),
                          _boundary[cur_edge].y(),
                          _boundary[(cur_edge + 1) % (_boundary.num_())].x(),
                          _boundary[(cur_edge + 1) % (_boundary.num_())].y());
    }
    return r;
  }

  _Point2 _Tessellation2D::convert_sym_to_point_(const std::set<int> &sym)
  {
    std::vector<int> s;
    std::vector<int> e;
    for (auto it = sym.begin(); it != sym.end(); ++it)
    {
      if ((*it) < 0)
        e.push_back(-(*it) - 1);
      else
        s.push_back((*it) - 1);
    }
    if (e.size() == 0)
    {
      _Point2 v0 = _sites[s[1]] - _sites[s[0]];
      _Point2 v1 = _sites[s[2]] - _sites[s[0]];
      double d0 = -0.5 * (v0.dot_((_sites[s[1]] + _sites[s[0]])) + (_weights[s[0]] - _weights[s[1]]));
      double d1 = -0.5 * (v1.dot_((_sites[s[2]] + _sites[s[0]])) + (_weights[s[0]] - _weights[s[2]]));
      return _Point2::intersection_two_line(v0, d0, v1, d1);
    }
    else if (e.size() == 1)
    {
      _Point2 v0 = _sites[s[1]] - _sites[s[0]];
      double d0 = -0.5 * (v0.dot_((_sites[s[1]] + _sites[s[0]])) + (_weights[s[0]] - _weights[s[1]]));
      _Point2 tv1 = _boundary[(e[0] + 1) % (_boundary.num_())] - _boundary[e[0]];
      _Point2 v1(-tv1.y(), tv1.x());
      double d1 = -v1.dot_(_boundary[e[0]]);
      return _Point2::intersection_two_line(v0, d0, v1, d1);
    }
    else
    {
      if (e[0] - e[1] == -1)
      {
        return _boundary[e[1]];
      }
      else if (e[0] - e[1] == 1)
      {
        return _boundary[e[0]];
      }
      else
      {
        return _boundary[0];
      }
    }
  }

  void _Tessellation2D::calculate_()
  {
    std::pair<_Point2, _Point2> bbox = _boundary.bounding_box_();
    _Point2 LBp = bbox.first + (bbox.first - bbox.second) * 2;
    _Point2 RTp = bbox.second + (bbox.second - bbox.first) * 2;
    double min_w = *(std::min_element(_weights.begin(), _weights.end()));
    std::vector<std::pair<Wpoint, int>> points;
    for (int i = 0; i < _num_sites; ++i)
    {
      points.push_back(std::make_pair(Wpoint(Point(_sites[i].x(), _sites[i].y()), _weights[i]), i));
    }
    points.push_back(std::make_pair(Wpoint(Point(LBp.x(), LBp.y()), min_w), -1));
    points.push_back(std::make_pair(Wpoint(Point(RTp.x(), LBp.y()), min_w), -2));
    points.push_back(std::make_pair(Wpoint(Point(RTp.x(), RTp.y()), min_w), -3));
    points.push_back(std::make_pair(Wpoint(Point(LBp.x(), RTp.y()), min_w), -4));
    Regular rt;
    rt.insert(points.begin(), points.end());
    _adjacent_edges.clear();
    _adjacent_edges.resize(_num_sites);
    for (auto it = rt.finite_faces_begin(); it != rt.finite_faces_end(); ++it)
    {
      int vi0 = it->vertex(0)->info();
      int vi1 = it->vertex(1)->info();
      int vi2 = it->vertex(2)->info();
      if (vi0 > -1)
        _adjacent_edges[vi0][vi1] = vi2;
      if (vi1 > -1)
        _adjacent_edges[vi1][vi2] = vi0;
      if (vi2 > -1)
        _adjacent_edges[vi2][vi0] = vi1;
    }
    std::vector<std::map<std::pair<int, int>, std::pair<int, int>>> boundary_edges(_num_sites);
    std::vector<std::map<std::pair<int, int>, std::pair<int, int>>> tessellation_edges(_num_sites);
    std::vector<std::vector<std::pair<std::pair<int, int>, int>>> boundary_cross(_num_sites);

    const _Point2 edge_center = (_boundary[0] + _boundary[1]) * 0.5;
    int nearest_site = rt.nearest_power_vertex(Point(edge_center.x(), edge_center.y()))->info();
    std::queue<std::pair<int, int>> Qse;
    Qse.push(std::make_pair(nearest_site, 0));
    std::set<std::pair<int, int>> site_edge_pair_is_visited;
    site_edge_pair_is_visited.insert(std::make_pair(nearest_site, 0));
    while (!Qse.empty())
    {
      const int cur_site = Qse.front().first;
      const int cur_edge = Qse.front().second;
      Qse.pop();
      std::pair<int, int>
          sym_ps = std::make_pair(-(cur_edge + _boundary.num_() - 1) % (_boundary.num_()) - 1, -cur_edge - 1);
      std::pair<int, int> sym_pt = std::make_pair(-cur_edge - 1, -(cur_edge + 1) % (_boundary.num_()) - 1);
      bool s_degraded = false;
      int s_degrad_site = -1;
      bool t_degraded = false;
      int t_degrad_site = -1;
      for (auto it = _adjacent_edges[cur_site].begin(); it != _adjacent_edges[cur_site].end(); ++it)
      {
        _BOC::_Sign side_s, side_t;
        if (it->first < 0)
        {
          side_s = _BOC::_Sign::PositivE;
          side_t = _BOC::_Sign::PositivE;
        }
        else
        {
          side_s = symbolic_point_site_(cur_site, it->first, sym_ps, cur_edge);
          side_t = symbolic_point_site_(cur_site, it->first, sym_pt, cur_edge);
        }
        if (side_s == _BOC::_Sign::NegativE && side_t == _BOC::_Sign::PositivE)
        {
          s_degraded = false;
          sym_ps = std::make_pair(-cur_edge - 1, it->first + 1);
        }
        else if (side_s == _BOC::_Sign::PositivE && side_t == _BOC::_Sign::NegativE)
        {
          t_degraded = false;
          sym_pt = std::make_pair(-cur_edge - 1, it->first + 1);
        }
        //else if ((side_s == _BOC::_Sign::NegativE && side_t == _BOC::_Sign::ZerO)
        //	|| (side_s == _BOC::_Sign::ZerO && side_t == _BOC::_Sign::NegativE))
        //{
        //	s_degraded = t_degraded = false;
        //	sym_ps = std::make_pair(0, 0);
        //	sym_pt = std::make_pair(-cur_edge - 1, it->first + 1);
        //	break;
        //}
        else if (side_s == _BOC::_Sign::NegativE && side_t == _BOC::_Sign::ZerO)
        {
          if (sym_pt.first < 0 && sym_pt.second < 0)
          {
            s_degraded = t_degraded = false;
            sym_ps = std::make_pair(0, 0);
            sym_pt = std::make_pair(-cur_edge - 1, it->first + 1);
            break;
          }
          else
          {
            s_degraded = false;
            sym_ps = std::make_pair(-cur_edge - 1, it->first + 1);
          }
        }
        else if (side_s == _BOC::_Sign::ZerO && side_t == _BOC::_Sign::NegativE)
        {
          if (sym_ps.first < 0 && sym_ps.second < 0)
          {
            s_degraded = t_degraded = false;
            sym_ps = std::make_pair(0, 0);
            sym_pt = std::make_pair(-cur_edge - 1, it->first + 1);
            break;
          }
          else
          {
            t_degraded = false;
            sym_pt = std::make_pair(-cur_edge - 1, it->first + 1);
          }
        }
        else if (side_s == _BOC::_Sign::ZerO && side_t == _BOC::_Sign::ZerO)
        {
          if ((_boundary[(cur_edge + 1) % (_boundary.num_())] - _boundary[cur_edge]).cross_(_sites[cur_site] - _sites[it->first]).z() < 0)
          {
            sym_ps = std::make_pair(0, 0);
            sym_pt = std::make_pair(0, 0);
            break;
          }
        }
        else if (side_s == _BOC::_Sign::NegativE && side_t == _BOC::_Sign::NegativE)
        {
          s_degraded = t_degraded = false;
          sym_ps = std::make_pair(0, 0);
          sym_pt = std::make_pair(0, 0);
          break;
        }
        else if (side_s == _BOC::_Sign::ZerO && side_t == _BOC::_Sign::PositivE)
        {
          s_degraded = true;
          s_degrad_site = it->first;
        }
        else if (side_s == _BOC::_Sign::PositivE && side_t == _BOC::_Sign::ZerO)
        {
          t_degraded = true;
          t_degrad_site = it->first;
        }
      }
      if (sym_ps.first != 0)
      {
        boundary_edges[cur_site][sym_ps] = sym_pt;
        if (sym_ps.first > 0 || sym_ps.second > 0)
        {
          boundary_cross[cur_site].push_back(std::make_pair(sym_ps, 1));
        }
        else if (s_degraded)
        {
          std::pair<int, int> new_sym;
          if (sym_ps.first - sym_ps.second == 1)
          {
            new_sym = std::make_pair(s_degrad_site + 1, sym_ps.second);
          }
          else if (sym_ps.first - sym_ps.second == -1)
          {
            new_sym = std::make_pair(s_degrad_site + 1, sym_ps.first);
          }
          else
          {
            new_sym = std::make_pair(s_degrad_site + 1, -1);
          }
          boundary_cross[cur_site].push_back(std::make_pair(new_sym, -1));
        }
        if (sym_pt.first > 0 || sym_pt.second > 0)
        {
          boundary_cross[cur_site].push_back(std::make_pair(sym_pt, 0));
        }
        else if (t_degraded)
        {
          std::pair<int, int> new_sym;
          if (sym_pt.first - sym_pt.second == 1)
          {
            new_sym = std::make_pair(t_degrad_site + 1, sym_pt.second);
          }
          else if (sym_pt.first - sym_pt.second == -1)
          {
            new_sym = std::make_pair(t_degrad_site + 1, sym_pt.first);
          }
          else
          {
            new_sym = std::make_pair(t_degrad_site + 1, -1);
          }
          boundary_cross[cur_site].push_back(std::make_pair(new_sym, -2));
        }
      }
      if (sym_ps.first == 0)
      {
        if (sym_pt.first != 0)
        {
          int site3 = sym_pt.first > 0 ? sym_pt.first - 1 : sym_pt.second - 1;
          if (site_edge_pair_is_visited.find(std::make_pair(site3, cur_edge)) == site_edge_pair_is_visited.end())
          {
            site_edge_pair_is_visited.insert(std::make_pair(site3, cur_edge));
            Qse.push(std::make_pair(site3, cur_edge));
          }
        }
      }
      else
      {
        if (sym_ps.first < 0 && sym_ps.second < 0)
        {
          if (site_edge_pair_is_visited.find(std::make_pair(cur_site,
                                                            (cur_edge + _boundary.num_() - 1) % (_boundary.num_()))) == site_edge_pair_is_visited.end())
          {
            site_edge_pair_is_visited.insert(std::make_pair(cur_site,
                                                            (cur_edge + _boundary.num_() - 1) % (_boundary.num_())));
            Qse.push(std::make_pair(cur_site, (cur_edge + _boundary.num_() - 1) % (_boundary.num_())));
          }
        }
        else
        {
          int site3 = sym_ps.first > 0 ? sym_ps.first - 1 : sym_ps.second - 1;
          if (site_edge_pair_is_visited.find(std::make_pair(site3, cur_edge)) == site_edge_pair_is_visited.end())
          {
            site_edge_pair_is_visited.insert(std::make_pair(site3, cur_edge));
            Qse.push(std::make_pair(site3, cur_edge));
          }
        }
        if (sym_pt.first < 0 && sym_pt.second < 0)
        {
          if (site_edge_pair_is_visited.find(std::make_pair(cur_site, (cur_edge + 1) % (_boundary.num_()))) == site_edge_pair_is_visited.end())
          {
            site_edge_pair_is_visited.insert(std::make_pair(cur_site, (cur_edge + 1) % (_boundary.num_())));
            Qse.push(std::make_pair(cur_site, (cur_edge + 1) % (_boundary.num_())));
          }
        }
        else
        {
          int site3 = sym_pt.first > 0 ? sym_pt.first - 1 : sym_pt.second - 1;
          if (site_edge_pair_is_visited.find(std::make_pair(site3, cur_edge)) == site_edge_pair_is_visited.end())
          {
            site_edge_pair_is_visited.insert(std::make_pair(site3, cur_edge));
            Qse.push(std::make_pair(site3, cur_edge));
          }
        }
      }
    }
    std::queue<std::pair<int, std::pair<int, int>>> Qss;
    std::vector<bool> site_visited(_num_sites, false);
    for (int i = 0; i < _num_sites; ++i)
    {
      if (boundary_cross[i].size() > 0)
      {
        for (int j = 0; j < boundary_cross[i].size(); ++j)
        {
          if (boundary_cross[i][j].second == 0 || boundary_cross[i][j].second == -2)
          {
            if (boundary_cross[i][j].second == -2)
            {
              int __bv = -boundary_cross[i][j].first.second - 1;
              std::pair<int, int>
                  __new_v = std::make_pair(-(__bv + _boundary.num_() - 1) % (_boundary.num_()) - 1, -__bv - 1);
              Qss.push(std::make_pair(i, __new_v));
            }
            else
            {
              Qss.push(std::make_pair(i, boundary_cross[i][j].first));
            }
            site_visited[i] = true;
            break;
          }
        }
      }
    }
    std::vector<std::map<std::pair<int, int>, int>> from_edge_start_to_adjacent(_num_sites);
    for (int i = 0; i < _num_sites; ++i)
    {
      tessellation_edges[i] = boundary_edges[i];
      for (auto it = boundary_edges[i].begin(); it != boundary_edges[i].end(); ++it)
      {
        from_edge_start_to_adjacent[i][it->first] = -1;
      }
    }
    while (!Qss.empty())
    {
      const int cur_site = Qss.front().first;
      const std::pair<int, int> start_vertex = Qss.front().second;
      Qss.pop();
      std::pair<int, int> next_vertex = start_vertex;
      int adj_site;
      do
      {
        if (next_vertex.first < 0 && next_vertex.second < 0)
        {
          int __vid = 0;
          if (next_vertex.first - next_vertex.second == 1)
          {
            __vid = next_vertex.second;
          }
          else if (next_vertex.first - next_vertex.second == -1)
          {
            __vid = next_vertex.first;
          }
          else
          {
            __vid = -1;
          }
          for (int i = 0; i < boundary_cross[cur_site].size(); ++i)
          {
            if ((boundary_cross[cur_site][i].second == -2) && (boundary_cross[cur_site][i].first.second == __vid))
            {
              adj_site = boundary_cross[cur_site][i].first.first;
              break;
            }
          }
        }
        else if (next_vertex.first < 0 || next_vertex.second < 0)
        {
          adj_site = next_vertex.first > 0 ? next_vertex.first : next_vertex.second;
        }
        else
        {
          adj_site = next_vertex.second;
        }
        bool is_done = false;
        for (int i = 0; i < boundary_cross[cur_site].size(); ++i)
        {
          if ((boundary_cross[cur_site][i].second == 1 || boundary_cross[cur_site][i].second == -1) && (boundary_cross[cur_site][i].first.first == adj_site || boundary_cross[cur_site][i].first.second == adj_site))
          {
            if (boundary_cross[cur_site][i].second == -1)
            {
              int __bv = -boundary_cross[cur_site][i].first.second - 1;
              std::pair<int, int>
                  __new_v = std::make_pair(-(__bv + _boundary.num_() - 1) % (_boundary.num_()) - 1, -__bv - 1);
              tessellation_edges[cur_site][next_vertex] = __new_v;
              from_edge_start_to_adjacent[cur_site][next_vertex] = (adj_site - 1);
              next_vertex = __new_v;
            }
            else
            {
              tessellation_edges[cur_site][next_vertex] = boundary_cross[cur_site][i].first;
              from_edge_start_to_adjacent[cur_site][next_vertex] = (adj_site - 1);
              next_vertex = boundary_cross[cur_site][i].first;
            }
            while (boundary_edges[cur_site].find(next_vertex) != boundary_edges[cur_site].end())
            {
              next_vertex = boundary_edges[cur_site][next_vertex];
            }
            is_done = true;
            break;
          }
        }
        if (!is_done)
        {
          int next_site = _adjacent_edges[cur_site][adj_site - 1];
          std::pair<int, int> new_vertex(adj_site, next_site + 1);
          tessellation_edges[cur_site][next_vertex] = new_vertex;
          from_edge_start_to_adjacent[cur_site][next_vertex] = adj_site - 1;
          next_vertex = new_vertex;
          if (!site_visited[next_site])
          {
            site_visited[next_site] = true;
            Qss.push(std::make_pair(next_site, std::make_pair(cur_site + 1, adj_site)));
          }
        }
      } while (next_vertex != start_vertex);
    }

    std::map<std::set<int>, int> from_sym_to_id;
    _tessellation_vertices.resize(0);
    _cells.resize(_num_sites);
    _is_hidden.clear();
    _is_hidden.resize(_num_sites, false);
    _cell_polygons.clear();
    _cell_polygons.resize(_num_sites);
    for (int i = 0; i < _num_sites; ++i)
    {
      _cells[i].resize(0);
      if (tessellation_edges[i].size() == 0)
      {
        _is_hidden[i] = true;
        continue;
      }
      _cell_polygons[i].start_();
      auto it = tessellation_edges[i].begin();
      std::pair<int, int> start = it->first;
      std::pair<int, int> next = start;
      do
      {
        std::pair<int, int> nextnext = tessellation_edges[i][next];
        std::set<int> _p;
        _p.insert(next.first);
        _p.insert(next.second);
        if (next.first < 0 && next.second < 0)
          _p.insert(1);
        else
          _p.insert(i + 1);
        int _adj = -1;
        _adj = from_edge_start_to_adjacent[i][next];

        if (from_sym_to_id.find(_p) != from_sym_to_id.end())
        {
          _cells[i].push_back(std::make_pair(from_sym_to_id[_p], _adj));
        }
        else
        {
          from_sym_to_id[_p] = _tessellation_vertices.size();
          _cells[i].push_back(std::make_pair(_tessellation_vertices.size(), _adj));
          _tessellation_vertices.push_back(convert_sym_to_point_(_p));
        }
        _cell_polygons[i].insert_(_tessellation_vertices[_cells[i].back().first]);
        next = tessellation_edges[i][next];
      } while (next != start);
      _cell_polygons[i].end_();
    }
  }

} // namespace BGAL
