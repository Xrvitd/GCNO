#pragma once
#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"

namespace BGAL 
{
	namespace Geodesic 
	{
		class _Abstract_Method 
		{
		public:
			_Abstract_Method(const _ManifoldModel& in_model);
			_Abstract_Method(const _ManifoldModel& in_model, const std::map<int, double>& in_sources);
			virtual void execute_();
			virtual void initialize_();
			std::vector<double> get_distances_() const;
		protected:
			virtual void implement_() = 0;
		protected:
			const _ManifoldModel& _model;
			std::map<int, double> _sources;
			int _method;
			std::vector<double> _distances;
			int _max_queue_length;
			int _max_result_depth;
		};
	}
}