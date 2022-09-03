#include "BGAL/Draw/DrawPS.h"

namespace BGAL
{
  _PS::_PS(std::ostream &in_os)
      : os(in_os)
  {
  }
  void _PS::set_bbox_(const std::pair<_Point2, _Point2> &bbox)
  {
    double minX = bbox.first.x();
    double minY = bbox.first.y();
    double maxX = bbox.second.x();
    double maxY = bbox.second.y();
    set_bbox_(minX, minY, maxX, maxY);
  }
  void _PS::set_bbox_(const double &minX, const double &minY, const double &maxX, const double &maxY)
  {
    os.precision(18);
    os << "%!PS-Adobe-3.0 EPSF-3.0\n";
    os << "%%BoundingBox: " << minX * 28.34645 - 3 << " " << minY * 28.34645 - 3 << " " << maxX * 28.34645 + 3 << " "
       << maxY * 28.34645 + 3 << "\n";
    os << "%%BeginProlog\n";
    os << "/s 28.34645 def\n";
    os << "%%EndProlog\n";
    os << "%%Page : 1 1\n";
    os << "s s scale\n";
  }
  void _PS::end_()
  {
    os << "showpage\n";
  }
  void _PS::draw_point_(const _Point2 &p, double linewidth, double r, double g, double b)
  {
    os << linewidth << " setlinewidth\n";
    os << r << " " << g << " " << b << " setrgbcolor\n";
    os << "newpath\n";
    os << p.x() << " " << p.y() << " " << linewidth
       << " 0 360 arc fill\n";
  }
  void _PS::draw_text_(const _Point2 &pos, const std::string &info, double scale, double r, double g, double b)
  {
    os << "/Times-Roman findfont\n";
    os << scale << " scalefont\n";
    os << r << " " << g << " " << b << " setrgbcolor\n";
    os << "setfont\n";
    os << pos.x() << " " << pos.y() << " moveto\n";
    os << "(" << info << ") show\n";
  }
  void _PS::draw_polygon_(const _Polygon &P, double linewidth, double r, double g, double b, bool filled)
  {
    size_t n = P.num_();
    if (n == 0)
      return;

    os << 0.5 * linewidth << " setlinewidth\n";
    os << r << " " << g << " " << b << " setrgbcolor\n";

    os << P[n - 1].x() << " "
       << P[n - 1].y() << " newpath moveto\n";
    for (int i = 0; i < n; ++i)
    {
      os << P[i].x() << " " << P[i].y() << " lineto\n";
    }
    os << "stroke\n";

    if (filled)
    {
      os << P[n - 1].x() << " "
         << P[n - 1].y() << " newpath moveto\n";
      for (int i = 0; i < n; ++i)
      {
        os << P[i].x() << " " << P[i].y() << " lineto\n";
      }
      os << "fill\n\n";
    }
  }
  void _PS::draw_line_segment_(const _Point2 &start, const _Point2 &end, double linewidth, double r, double g, double b)
  {
    os << linewidth << " setlinewidth\n";
    os << r << " " << g << " " << b << " setrgbcolor\n";

    os << start.x() << " "
       << start.y() << " newpath moveto\n";
    os << end.x() << " "
       << end.y() << " lineto\n";
    os << "stroke\n\n";
  }
} // namespace BGAL