#pragma once
#include "BGAL/BaseShape/Polygon.h"
#include <set>
#include <map>
namespace BGAL
{
  class _Tessellation2D
  {
  public:
    _Tessellation2D(const _Polygon &in_boundary, const std::vector<_Point2> &in_sites);
    _Tessellation2D(const _Polygon &in_boundary,
                    const std::vector<_Point2> &in_sites,
                    const std::vector<double> &in_weights);
    void calculate_(const _Polygon &in_boundary, const std::vector<_Point2> &in_sites);
    void calculate_(const _Polygon &in_boundary,
                    const std::vector<_Point2> &in_sites,
                    const std::vector<double> &in_weights);
    const std::vector<_Polygon> &get_cell_polygons_() const;
    const std::vector<std::vector<std::pair<int, int>>> &get_cells_() const;
    const _Point2 &vertex_(const int &id) const;
    int num_hidden_() const;
    int num_vertex_() const;
  private:
    void calculate_();
    _BOC::_Sign symbolic_point_site_(const int &site1,
                                     const int &site2,
                                     const std::pair<int, int> &symp,
                                     const int &cur_edge);
    _Point2 convert_sym_to_point_(const std::set<int> &sym);

  public:
  private:
    int _num_sites;
    _Polygon _boundary;
    std::vector<std::map<int, int>> _adjacent_edges;
    std::vector<_Point2> _sites;
    std::vector<double> _weights;
    std::vector<bool> _is_hidden;
    std::vector<_Point2> _tessellation_vertices;
    std::vector<std::vector<std::pair<int, int>>> _cells;
    std::vector<_Polygon> _cell_polygons;
  };
} // namespace BGAL