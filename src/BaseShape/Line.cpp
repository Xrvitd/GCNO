//#include "BKLine.h"
#include "BGAL/BaseShape/Line.h"
namespace BGAL
{
    _Segment3::_Segment3()
        : _s(_Point3(0, 0, 0)), _t(_Point3(0, 0, 0))
    {
    }
    _Segment3::_Segment3(const _Point3 &in_s, const _Point3 &in_t)
        : _s(in_s), _t(in_t)
    {
    }
} // namespace BGAL