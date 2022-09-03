#include "BGAL/Geodesic/AbstractMethod.h"

namespace BGAL
{
  namespace Geodesic
  {
    _Abstract_Method::_Abstract_Method(const _ManifoldModel &in_model)
        : _model(in_model), _method(0)
    {
    }
    _Abstract_Method::_Abstract_Method(const _ManifoldModel &in_model, const std::map<int, double> &in_sources)
        : _model(in_model), _sources(in_sources), _method(0)
    {
    }
    void _Abstract_Method::execute_()
    {
      initialize_();
      implement_();
    }
    void _Abstract_Method::initialize_()
    {
      _distances.clear();
      _distances.resize(_model.number_vertices_(), std::numeric_limits<double>::max());
      _max_queue_length = 0;
      _max_result_depth = 0;
    }
    std::vector<double> _Abstract_Method::get_distances_() const
    {
      return _distances;
    }
  } // namespace Geodesic
} // namespace BGAL