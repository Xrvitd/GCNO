#include "BGAL/Geodesic/Dijkstra/Dijkstra.h"
#include <queue>
namespace BGAL
{
	namespace Geodesic
	{
		void _Dijkstra::initialize_()
		{
			_Abstract_Method::initialize_();
		}
		void _Dijkstra::implement_()
		{
			struct _Event : std::tuple<double, int, int, int, int>
			{
				_Event(double dis, int self, int parent, int root, int level)
					: std::tuple<double, int, int, int, int>(dis, self, parent, root, level)
				{
				}
				bool operator<(const _Event &right) const
				{
					return std::get<0>(*this) < std::get<0>(right);
				}
				bool operator>(const _Event &right) const
				{
					return std::get<0>(*this) > std::get<0>(right);
				}
				double get_distance_() const
				{
					return std::get<0>(*this);
				}
				int get_id_() const
				{
					return std::get<1>(*this);
				}
				int get_parent_() const
				{
					return std::get<2>(*this);
				}
				int get_root_() const
				{
					return std::get<3>(*this);
				}
				int get_level_() const
				{
					return std::get<4>(*this);
				}
			};
			std::priority_queue<_Event, std::vector<_Event>, std::greater<_Event>> evtQue;
			for (auto it = _sources.begin(); it != _sources.end(); ++it)
			{
				_Event evt(it->second, it->first, -1, it->first, 0);
				evtQue.push(evt);
			}
			std::vector<bool> visited(_model.number_vertices_(), false);
			_result.clear();
			_result.resize(_model.number_vertices_());
			while (!evtQue.empty())
			{
				if (evtQue.size() > _max_queue_length)
					_max_queue_length = evtQue.size();
				_Event topEvt = evtQue.top();
				evtQue.pop();
				if (topEvt.get_level_() > _max_result_depth)
					_max_result_depth = topEvt.get_level_();
				if (visited[topEvt.get_id_()])
					continue;
				visited[topEvt.get_id_()] = true;
				_result[topEvt.get_id_()] = std::make_tuple(topEvt.get_parent_(), topEvt.get_root_(), topEvt.get_level_(), topEvt.get_distance_());
				for (auto veit = _model.ve_begin(topEvt.get_id_()); veit != _model.ve_end(topEvt.get_id_()); ++veit)
				{
					int v = (*veit)._id_right_vertex;
					if (visited[v])
						continue;
					_Event evt(topEvt.get_distance_() + (*veit).length_(), v, topEvt.get_id_(), topEvt.get_root_(), topEvt.get_level_() + 1);
					evtQue.push(evt);
				}
			}
			for (int i = 0; i < _model.number_vertices_(); ++i)
			{
				_distances[i] = std::get<3>(_result[i]);
			}
		}
		_Dijkstra::_Dijkstra(const _ManifoldModel &in_model, const std::map<int, double> &in_sources)
			: _Abstract_Method(in_model, in_sources)
		{
			_method = 1;
		}

	} // namespace Geodesic
} // namespace BGAL
