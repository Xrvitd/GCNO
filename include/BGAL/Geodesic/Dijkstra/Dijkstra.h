#pragma once
#include "BGAL/Geodesic/AbstractMethod.h"
namespace BGAL 
{
	namespace Geodesic 
	{
		class _Dijkstra : public _Abstract_Method 
		{
		protected:
			virtual void initialize_();
			virtual void implement_();
		public:
			_Dijkstra(const _ManifoldModel& in_model, const std::map<int, double>& in_sources);
		protected:
			std::vector<std::tuple<int, int, int, double>> _result;
		};
	}
}