#pragma once
#include <Eigen/Dense>
#include <time.h>
#include <vector>
#include <iostream>
namespace BGAL
{
  class _LBFGS
  {
  public:
    class _Parameter
    {
    public:
      int m;
      int max_iteration;
      int max_linearsearch;
      double min_step;
      double min_xtol;
      double min_ftol;
      double epsilon;
      bool is_show;
      double wolfe;
      double max_time;
      _Parameter();
      _Parameter(const _Parameter &in_parameter);
    };
    _Parameter _parameter;
    _LBFGS();
    _LBFGS(const _Parameter &in_parameter);
    template <class fun>
    int minimize(fun &f, Eigen::VectorXd &iterX);

  private:
    template <class fun>
    int linear_search_(fun &f,
                       double &fval,
                       Eigen::VectorXd &iterX,
                       Eigen::VectorXd &gradient,
                       double &step,
                       const Eigen::VectorXd &direction);
  };
  template <class fun>
  int _LBFGS::minimize(fun &f, Eigen::VectorXd &iterX)
  {
    const clock_t start_t = clock();
    const int n = iterX.size();
    std::vector<Eigen::VectorXd> s(_parameter.m);
    std::vector<Eigen::VectorXd> y(_parameter.m);
    Eigen::VectorXd alpha(_parameter.m);
    Eigen::VectorXd ys(_parameter.m);
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd oldX, oldGradient;
    double fval = f(iterX, gradient);
    if (_parameter.is_show)
    {
      std::cout << 0 << "\t" << 0 << "\t" << (clock() - start_t) * 1.0 / CLOCKS_PER_SEC << "\t" << gradient.norm() << "\t"
                << fval << std::endl;
    }
    Eigen::VectorXd direction = -gradient;
    int k = 0;
    int l = 0;
    int cursor = 0;
    double step = 1.0 / direction.norm();
    while (1)
    {
      if (gradient.norm() < _parameter.epsilon)
      {
        if (_parameter.is_show)
        {
          std::cout << "reach the gradient tolerance" << std::endl;
        }
        return k;
      }
      if (_parameter.max_time > 0 && _parameter.max_time < (clock() - start_t) * 1000.0 / CLOCKS_PER_SEC)
      {
        if (_parameter.is_show)
        {
          std::cout << "reach the max time" << std::endl;
        }
        return k;
      }
      if (_parameter.max_iteration != 0 && _parameter.max_iteration == k)
      {
        if (_parameter.is_show)
        {
          std::cout << "reach the max itertion time" << std::endl;
        }
        return k;
      }
      oldX = iterX;
      oldGradient = gradient;
      int num_linear = linear_search_(f, fval, iterX, gradient, step, direction);
      if (num_linear == _parameter.max_linearsearch)
      {
        if (_parameter.is_show)
        {
          std::cout << "reach the max linearsearch time" << std::endl;
        }
        return k;
      }
      if (num_linear < 0)
      {
        if (_parameter.is_show)
        {
          std::cout << "can't fine a right step" << std::endl;
        }
        return k;
      }
      l += num_linear;
      k++;
      if (_parameter.is_show)
      {
        std::cout << k << "\t" << l << "\t" << (clock() - start_t) * 1.0 / CLOCKS_PER_SEC << "\t" << gradient.norm()
                  << "\t" << fval << std::endl;
      }
      s[cursor] = iterX - oldX;
      y[cursor] = gradient - oldGradient;
      double _ys = y[cursor].dot(s[cursor]);
      double _yy = y[cursor].dot(y[cursor]);
      ys[cursor] = _ys;
      cursor = (cursor + 1) % _parameter.m;
      int incr, bound;
      if (k < _parameter.m)
      {
        incr = 0;
        bound = k;
      }
      else
      {
        incr = cursor;
        bound = _parameter.m;
      }
      direction = -gradient;
      int j = cursor;

      for (int i = 0; i < bound; i++)
      {
        j = (j + _parameter.m - 1) % _parameter.m;
        alpha(j) = s[j].dot(direction) / ys(j);
        direction = direction - alpha[j] * y[j];
      }
      direction = (_ys / _yy) * direction;
      for (int i = 0; i < bound; i++)
      {
        double beta = y[j].dot(direction) / ys(j);
        direction = direction + (alpha(j) - beta) * s[j];
        j = (j + 1) % _parameter.m;
      }
      step = 1.0;
    }
  }
  template <class fun>
  int _LBFGS::linear_search_(fun &f,
                             double &fval,
                             Eigen::VectorXd &iterX,
                             Eigen::VectorXd &gradient,
                             double &step,
                             const Eigen::VectorXd &direction)
  {
    if (step < 0)
    {
      std::cout << "error! step<0" << std::endl;
      throw std::runtime_error("error! step<0");
    }
    const double ifval = fval;
    const Eigen::VectorXd iX = iterX;
    const double idg = gradient.dot(direction);
    if (idg > 0)
    {
      std::cout << "error! direction is not a decline direction!" << std::endl;
      throw std::runtime_error("error! direction is not a decline direction!");
    }
    double step_l = 0;
    double step_u = std::numeric_limits<double>::infinity();
    int k = 1;
    while (1)
    {
      iterX = iX + step * direction;
      fval = f(iterX, gradient);
      if (fval > ifval)
      {
        step_u = step;
      }
      else
      {
        double dg = gradient.dot(direction);
        if (dg < _parameter.wolfe * idg)
        {
          step_l = step;
        }
        else
        {
          if (dg > -_parameter.wolfe * idg)
          {
            step_u = step;
          }
          else
          {
            break;
          }
        }
      }
      k++;
      if (k == _parameter.max_linearsearch)
      {
        iterX = iX;
        break;
      }
      if (step < _parameter.min_step)
      {
        k = -1;
        iterX = iX;
        break;
      }
      step = isinf(step_u) ? 2 * step : 0.5 * (step_l + step_u);
      //step /= 10;
    }
    return k;
  }

} // namespace BGAL