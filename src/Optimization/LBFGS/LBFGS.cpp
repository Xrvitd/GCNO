#include "BGAL/Optimization/LBFGS/LBFGS.h"
namespace BGAL
{
  _LBFGS::_Parameter::_Parameter()
  {
    m = 6;
    max_iteration = 1000;
    max_linearsearch = 100;
    min_step = 1e-20;
    min_xtol = 1e-20;
    min_ftol = 1e-20;
    epsilon = 1e-6;
    is_show = false;
    wolfe = 0.9;
    max_time = -1.0;
  }
  _LBFGS::_Parameter::_Parameter(const _Parameter &in_parameter)
  {
    m = in_parameter.m;
    max_iteration = in_parameter.max_iteration;
    max_linearsearch = in_parameter.max_linearsearch;
    min_step = in_parameter.min_step;
    min_xtol = in_parameter.min_xtol;
    min_ftol = in_parameter.min_ftol;
    epsilon = in_parameter.epsilon;
    is_show = in_parameter.is_show;
    wolfe = in_parameter.wolfe;
    max_time = in_parameter.max_time;
  }
  _LBFGS::_LBFGS()
  {
    //_parameter = _Parameter();
  }
  _LBFGS::_LBFGS(const _Parameter &in_parameter)
      : _parameter(in_parameter)
  {
  }
} // namespace BGAL