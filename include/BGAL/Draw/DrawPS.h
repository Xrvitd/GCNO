#pragma once
#include <iostream>
#include <string>
#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
namespace BGAL 
{
	class _PS 
	{
	private:
		std::ostream& os;
	public:
		_PS(std::ostream& in_os);
		void set_bbox_(const std::pair<_Point2, _Point2>& bbox);
		void set_bbox_(const double& minX, const double& minY, const double& maxX, const double& maxY);
		void end_();
		void draw_point_(const _Point2& p, double linewidth = 0.1, double r = 0, double g = 0, double b = 0);
		void draw_text_(const _Point2& pos, const std::string& info, double scale, double r = 0, double g = 0, double b = 0);
		void draw_polygon_(const _Polygon& P, double linewidth = 0.1, double r = 0, double g = 0, double b = 0, bool filled = false);
		void draw_line_segment_(const _Point2& start, const _Point2& end, double linewidth = 0.1, double r = 0, double g = 0, double b = 0);
	};
}