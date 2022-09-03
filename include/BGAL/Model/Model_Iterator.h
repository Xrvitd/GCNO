#pragma once
#include "Model.h"
#include "ManifoldModel.h"
namespace BGAL 
{
	class _Face_Iterator 
	{
	public:
		_Face_Iterator(const _Model* in_model, const int& in_cursor);
		_Face_Iterator& operator=(const _Face_Iterator& fit);
		const _Model::_MFace& operator*();
		_Face_Iterator& operator++();
		inline bool operator<(const _Face_Iterator& fit) const 
		{
			return _cursor < fit._cursor;
		}
		inline bool operator==(const _Face_Iterator& fit) const 
		{
			return _cursor == fit._cursor;
		}
		inline bool operator>(const _Face_Iterator& fit) const 
		{
			return !(*this < fit || *this == fit);
		}
		inline bool operator!=(const _Face_Iterator& fit) const 
		{
			return !(*this == fit);
		}
		inline bool operator!=(const int& fit) const 
		{
			return !(_cursor == fit);
		}
		inline bool operator<=(const _Face_Iterator& fit) const 
		{
			return (*this < fit || *this == fit);
		}
		inline bool operator>=(const _Face_Iterator& fit) const 
		{
			return !(*this < fit);
		}
		inline int id() const 
		{
			return _cursor;
		}
		~_Face_Iterator();
	protected:
		const _Model* _model;
		int _cursor;
	};
	class _Vertex_Iterator 
	{
	public:
		_Vertex_Iterator(const _Model* in_model, const int& in_cursor);
		_Vertex_Iterator& operator=(const _Vertex_Iterator& vit);
		const _Point3& operator*();
		_Vertex_Iterator& operator++();
		inline bool operator<(const _Vertex_Iterator& vit) const 
		{
			return _cursor < vit._cursor;
		}
		inline bool operator==(const _Vertex_Iterator& vit) const 
		{
			return _cursor == vit._cursor;
		}
		inline bool operator>(const _Vertex_Iterator& vit) const 
		{
			return !(*this < vit || *this == vit);
		}
		inline bool operator!=(const _Vertex_Iterator& vit) const 
		{
			return !(*this == vit);
		}
		inline bool operator!=(const int& vit) const 
		{
			return !(_cursor == vit);
		}
		inline bool operator<=(const _Vertex_Iterator& vit) const 
		{
			return (*this < vit || *this == vit);
		}
		inline bool operator>=(const _Vertex_Iterator& vit) const 
		{
			return !(*this < vit);
		}
		inline int id() const 
		{
			return _cursor;
		}
		~_Vertex_Iterator();
	protected:
		const _Model* _model;
		int _cursor;
	};
	class _Edge_Iterator 
	{
	public:
		_Edge_Iterator(const _ManifoldModel* in_model, const int& in_cursor);
		_Edge_Iterator& operator=(const _Edge_Iterator& eit);
		const _ManifoldModel::_MMEdge& operator*();
		_Edge_Iterator& operator++();
		inline bool operator<(const _Edge_Iterator& eit) const 
		{
			return _cursor < eit._cursor;
		}
		inline bool operator==(const _Edge_Iterator& eit) const 
		{
			return _cursor == eit._cursor;
		}
		inline bool operator>(const _Edge_Iterator& eit) const 
		{
			return !(*this < eit || *this == eit);
		}
		inline bool operator!=(const _Edge_Iterator& eit) const 
		{
			return !(*this == eit);
		}
		inline bool operator!=(const int& eit) const 
		{
			return !(_cursor == eit);
		}
		inline bool operator<=(const _Edge_Iterator& eit) const 
		{
			return (*this < eit || *this == eit);
		}
		inline bool operator>=(const _Edge_Iterator& eit) const 
		{
			return !(*this < eit);
		}
		inline int id() const 
		{
			return _cursor;
		}
		~_Edge_Iterator();
	protected:
		const _ManifoldModel* _model;
		int _cursor;
	};
	class _FV_Iterator 
	{
	public:
		_FV_Iterator(const _Model* in_model, const int& in_fid, const int& in_cursor);
		_FV_Iterator& operator=(const _FV_Iterator& fvit);
		const _Point3& operator*();
		_FV_Iterator& operator++();
		inline bool operator<(const _FV_Iterator& fvit) const 
		{
			return _cursor < fvit._cursor;
		}
		inline bool operator==(const _FV_Iterator& fvit) const 
		{
			return _cursor == fvit._cursor;
		}
		inline bool operator>(const _FV_Iterator& fvit) const 
		{
			return !(*this < fvit || *this == fvit);
		}
		inline bool operator!=(const _FV_Iterator& fvit) const 
		{
			return !(*this == fvit);
		}
		inline bool operator!=(const int& fvit) const 
		{
			return !(_cursor == fvit);
		}
		inline bool operator<=(const _FV_Iterator& fvit) const 
		{
			return (*this < fvit || *this == fvit);
		}
		inline bool operator>=(const _FV_Iterator& fvit) const 
		{
			return !(*this < fvit);
		}
		inline int id() const 
		{
			return _model->face_(_fid)[_cursor];
		}
		inline int fid() const 
		{
			return _fid;
		}
		~_FV_Iterator();
	protected:
		const _Model* _model;
		int _fid;
		int _cursor;
	};
	class _FE_Iterator 
	{
	public:
		_FE_Iterator(const _ManifoldModel* in_model, const int& in_fid, const int& in_cursor);
		_FE_Iterator& operator=(const _FE_Iterator& feit);
		const _ManifoldModel::_MMEdge& operator*();
		_FE_Iterator& operator++();
		inline bool operator<(const _FE_Iterator& feit) const 
		{
			return _cursor < feit._cursor;
		}
		inline bool operator==(const _FE_Iterator& feit) const 
		{
			return _cursor == feit._cursor;
		}
		inline bool operator>(const _FE_Iterator& feit) const 
		{
			return !(*this < feit || *this == feit);
		}
		inline bool operator!=(const _FE_Iterator& feit) const 
		{
			return !(*this == feit);
		}
		inline bool operator!=(const int& feit) const 
		{
			return !(_cursor == feit);
		}
		inline bool operator<=(const _FE_Iterator& feit) const 
		{
			return (*this < feit || *this == feit);
		}
		inline bool operator>=(const _FE_Iterator& feit) const 
		{
			return !(*this < feit);
		}
		inline int id() const 
		{
			return _eid;
		}
		inline int fid() const 
		{
			return _fid;
		}
		~_FE_Iterator();
	protected:
		const _ManifoldModel* _model;
		int _fid;
		int _cursor;
		int _eid;
	};
	class _FF_Iterator 
	{
	public:
		_FF_Iterator(const _ManifoldModel* in_model, const int& in_fid, const int& in_cursor);
		_FF_Iterator& operator=(const _FF_Iterator& ffit);
		const _Model::_MFace& operator*();
		_FF_Iterator& operator++();
		inline bool operator<(const _FF_Iterator& ffit) const 
		{
			return _cursor < ffit._cursor;
		}
		inline bool operator==(const _FF_Iterator& ffit) const 
		{
			return _cursor == ffit._cursor;
		}
		inline bool operator>(const _FF_Iterator& ffit) const 
		{
			return !(*this < ffit || *this == ffit);
		}
		inline bool operator!=(const _FF_Iterator& ffit) const 
		{
			return !(*this == ffit);
		}
		inline bool operator!=(const int& ffit) const 
		{
			return !(_cursor == ffit);
		}
		inline bool operator<=(const _FF_Iterator& ffit) const 
		{
			return (*this < ffit || *this == ffit);
		}
		inline bool operator>=(const _FF_Iterator& ffit) const 
		{
			return !(*this < ffit);
		}
		inline int id() const 
		{
			return _model->edge_(_eid)._id_face;
		}
		inline int fid() const 
		{
			return _fid;
		}
		~_FF_Iterator();
	protected:
		const _ManifoldModel* _model;
		int _fid;
		int _cursor;
		int _eid;
	};
	class _VV_Iterator 
	{
	public:
		_VV_Iterator(const _ManifoldModel* in_model, const int& in_vid, const int& in_cursor);
		_VV_Iterator& operator=(const _VV_Iterator& vvit);
		const _Point3& operator*();
		_VV_Iterator& operator++();
		inline bool operator<(const _VV_Iterator& vvit) const 
		{
			return _cursor < vvit._cursor;
		}
		inline bool operator==(const _VV_Iterator& vvit) const 
		{
			return _cursor == vvit._cursor;
		}
		inline bool operator>(const _VV_Iterator& vvit) const 
		{
			return !(*this < vvit || *this == vvit);
		}
		inline bool operator!=(const _VV_Iterator& vvit) const 
		{
			return !(*this == vvit);
		}
		inline bool operator!=(const int& vvit) const 
		{
			return !(_cursor == vvit);
		}
		inline bool operator<=(const _VV_Iterator& vvit) const 
		{
			return (*this < vvit || *this == vvit);
		}
		inline bool operator>=(const _VV_Iterator& vvit) const 
		{
			return !(*this < vvit);
		}
		inline int id() const 
		{
			return _model->edge_(_eid)._id_right_vertex;
		}
		inline int vid() const 
		{
			return _vid;
		}
		~_VV_Iterator();
	protected:
		const _ManifoldModel* _model;
		int _vid;
		int _cursor;
		int _eid;
	};
	class _VE_Iterator 
	{
	public:
		_VE_Iterator(const _ManifoldModel* in_model, const int& in_vid, const int& in_cursor);
		_VE_Iterator& operator=(const _VE_Iterator& veit);
		const _ManifoldModel::_MMEdge& operator*();
		_VE_Iterator& operator++();
		inline bool operator<(const _VE_Iterator& veit) const 
		{
			return _cursor < veit._cursor;
		}
		inline bool operator==(const _VE_Iterator& veit) const 
		{
			return _cursor == veit._cursor;
		}
		inline bool operator>(const _VE_Iterator& veit) const 
		{
			return !(*this < veit || *this == veit);
		}
		inline bool operator!=(const _VE_Iterator& veit) const 
		{
			return !(*this == veit);
		}
		inline bool operator!=(const int& veit) const 
		{
			return !(_cursor == veit);
		}
		inline bool operator<=(const _VE_Iterator& veit) const 
		{
			return (*this < veit || *this == veit);
		}
		inline bool operator>=(const _VE_Iterator& veit) const 
		{
			return !(*this < veit);
		}
		inline int id() const 
		{
			return _eid;
		}
		inline int vid() const 
		{
			return _vid;
		}
		~_VE_Iterator();
	protected:
		const _ManifoldModel* _model;
		int _vid;
		int _cursor;
		int _eid;
	};
	class _VF_Iterator 
	{
	public:
		_VF_Iterator(const _ManifoldModel* in_model, const int& in_vid, const int& in_cursor);
		_VF_Iterator& operator=(const _VF_Iterator& vfit);
		const _Model::_MFace& operator*();
		_VF_Iterator& operator++();
		inline bool operator<(const _VF_Iterator& vfit) const 
		{
			return _cursor < vfit._cursor;
		}
		inline bool operator==(const _VF_Iterator& vfit) const 
		{
			return _cursor == vfit._cursor;
		}
		inline bool operator>(const _VF_Iterator& vfit) const 
		{
			return !(*this < vfit || *this == vfit);
		}
		inline bool operator!=(const _VF_Iterator& vfit) const 
		{
			return !(*this == vfit);
		}
		inline bool operator!=(const int& vfit) const 
		{
			return !(_cursor == vfit);
		}
		inline bool operator<=(const _VF_Iterator& vfit) const 
		{
			return (*this < vfit || *this == vfit);
		}
		inline bool operator>=(const _VF_Iterator& vfit) const 
		{
			return !(*this < vfit);
		}
		inline int id() const 
		{
			return _model->edge_(_eid)._id_face;
		}
		inline int vid() const 
		{
			return _vid;
		}
		~_VF_Iterator();
	protected:
		const _ManifoldModel* _model;
		int _vid;
		int _cursor;
		int _eid;
	};
}