#include "BGAL/Optimization/GradientDescent/GradientDescent.h"

namespace BGAL
{
  _GradientDescent::_Parameter::_Parameter()
  {
    max_iteration = 1000;
    max_linearsearch = 100;
    min_step = 1e-20;
    min_xtol = 1e-20;
    min_ftol = 1e-20;
    epsilon = 1e-6;
    is_show = false;
    wolfe = 0.9;
    init_step = 1.0;
  }
  _GradientDescent::_GradientDescent()
      : _parameter()
  {
  }
  _GradientDescent::_GradientDescent(const _Parameter &in_parameter)
      : _parameter(in_parameter)
  {
  }
} // namespace BGAL