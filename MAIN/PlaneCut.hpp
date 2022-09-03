#pragma once
#include"MyHalfEdgeModel.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
struct FaceInfo2
{
	FaceInfo2() {}
	int nesting_level;
	bool in_domain() {
		return nesting_level % 2 == 1;
	}
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;
typedef CDT::Face_handle                                          Face_handle;
typedef CDT::Vertex_handle                                          Vertex_handle;
typedef K::Plane_3                                     Plane;
void
mark_domains(CDT& ct,
	Face_handle start,
	int index,
	std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1) {
		return;
	}
	std::list<Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()) {
		Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1) {
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++) {
				CDT::Edge e(fh, i);
				Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1) {
					if (ct.is_constrained(e)) border.push_back(e);
					else queue.push_back(n);
				}
			}
		}
	}
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void
mark_domains(CDT& cdt)
{
	for (CDT::Face_handle f : cdt.all_face_handles()) {
		f->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()) {
		CDT::Edge e = border.front();
		border.pop_front();
		Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1) {
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}

using namespace std;

#define My_MAX 0.1

class PlaneCut
{
public:
	MyBaseModel* CuttedMesh;
	double mx = My_MAX, mi = -My_MAX;
	vector<Plane> Planes;
public:
	PlaneCut();
	PlaneCut(Eigen::Vector3d center);
	PlaneCut(MyBaseModel* Mesh) :CuttedMesh(Mesh)
	{}
	//从一个特定的mesh开始
	PlaneCut(const char* filename);
	//从一个文件开始
	bool CutByPlane(Plane plane);
};

PlaneCut::PlaneCut()
{
	vector<Eigen::Vector3d> Points;
	vector<Eigen::Vector3i> Faces;
	Points.push_back(Eigen::Vector3d(mi, mi, mi));
	Points.push_back(Eigen::Vector3d(mi, mx, mi));
	Points.push_back(Eigen::Vector3d(mx, mx, mi));
	Points.push_back(Eigen::Vector3d(mx, mi, mi));

	Points.push_back(Eigen::Vector3d(mi, mi, mx));
	Points.push_back(Eigen::Vector3d(mi, mx, mx));
	Points.push_back(Eigen::Vector3d(mx, mx, mx));
	Points.push_back(Eigen::Vector3d(mx, mi, mx));

	Faces.push_back(Eigen::Vector3i(1, 3, 2));
	Faces.push_back(Eigen::Vector3i(1, 4, 3));

	Faces.push_back(Eigen::Vector3i(5, 6, 7));
	Faces.push_back(Eigen::Vector3i(5, 7, 8));

	Faces.push_back(Eigen::Vector3i(1, 2, 6));
	Faces.push_back(Eigen::Vector3i(1, 6, 5));

	Faces.push_back(Eigen::Vector3i(2, 3, 7));
	Faces.push_back(Eigen::Vector3i(2, 7, 6));

	Faces.push_back(Eigen::Vector3i(3, 4, 8));
	Faces.push_back(Eigen::Vector3i(3, 8, 7));

	Faces.push_back(Eigen::Vector3i(4, 1, 5));
	Faces.push_back(Eigen::Vector3i(4, 5, 8));

	for (int i = 0; i < Faces.size(); i++)
	{
		int t = Faces[i].y();
		Faces[i].y() = Faces[i].z();
		Faces[i].z() = t;
		Faces[i] -= Eigen::Vector3i(1, 1, 1).reverse();
	}

	CuttedMesh = new MyBaseModel(Points, Faces);
}//默认创建一个巨大的cube


PlaneCut::PlaneCut(Eigen::Vector3d center)
{
	vector<Eigen::Vector3d> Points;
	vector<Eigen::Vector3i> Faces;
	Points.push_back(Eigen::Vector3d(mi, mi, mi) + center);
	Points.push_back(Eigen::Vector3d(mi, mx, mi) + center);
	Points.push_back(Eigen::Vector3d(mx, mx, mi) + center);
	Points.push_back(Eigen::Vector3d(mx, mi, mi) + center);

	Points.push_back(Eigen::Vector3d(mi, mi, mx) + center);
	Points.push_back(Eigen::Vector3d(mi, mx, mx) + center);
	Points.push_back(Eigen::Vector3d(mx, mx, mx) + center);
	Points.push_back(Eigen::Vector3d(mx, mi, mx) + center);

	Faces.push_back(Eigen::Vector3i(1, 3, 2));
	Faces.push_back(Eigen::Vector3i(1, 4, 3));

	Faces.push_back(Eigen::Vector3i(5, 6, 7));
	Faces.push_back(Eigen::Vector3i(5, 7, 8));

	Faces.push_back(Eigen::Vector3i(1, 2, 6));
	Faces.push_back(Eigen::Vector3i(1, 6, 5));

	Faces.push_back(Eigen::Vector3i(2, 3, 7));
	Faces.push_back(Eigen::Vector3i(2, 7, 6));

	Faces.push_back(Eigen::Vector3i(3, 4, 8));
	Faces.push_back(Eigen::Vector3i(3, 8, 7));

	Faces.push_back(Eigen::Vector3i(4, 1, 5));
	Faces.push_back(Eigen::Vector3i(4, 5, 8));

	for (int i = 0; i < Faces.size(); i++)
	{
		int t = Faces[i].y();
		Faces[i].y() = Faces[i].z();
		Faces[i].z() = t;
		Faces[i] -= Eigen::Vector3i(1, 1, 1).reverse();
	}

	CuttedMesh = new MyBaseModel(Points, Faces);
}


PlaneCut::PlaneCut(const char* filename)
{
	MyBaseModel newModel;
	newModel.ReadObjFile(filename);
	*CuttedMesh = newModel;
}

bool PlaneCut::CutByPlane(Plane plane)
{
	double a, b, c, d;
	a = plane.a();
	b = plane.b();
	c = plane.c();
	d = plane.d();

	Planes.push_back(plane);

	return 1;


	vector<double> scalarField;
	vector<Eigen::Vector3d> vertics = CuttedMesh->GetVertices();
	map<Eigen::Vector3d, int> Vmap;


	Eigen::Vector3d P1 = vertics[0];
	K::Point_3 P1_Cgal(P1.x(), P1.y(), P1.z());
	bool side = plane.has_on_positive_side(P1_Cgal);
	bool conti = 0;
	for (int i = 1; i < vertics.size(); i++)
	{
		K::Point_3 Pnow(vertics[i].x(), vertics[i].y(), vertics[i].z());
		if (plane.has_on_positive_side(Pnow) != side)
		{
			conti = 1;
			break;
		}
	}
	if (!conti)
	{
		return 0;
	}




	//int id = 0;
	bool ifzheng = 0;
	for (auto p : vertics)
	{
		double value = a * p.x() + b * p.y() + c * p.z() + d;
		if (value > 0)
		{
			ifzheng = 1;
		}
		scalarField.push_back(value);
	}
	//计算标量场

	/*if (ifzheng)
	{
		return 1;
	}
	else
	{
		return 0;
	}*/




	auto isoline = CuttedMesh->ExtractIsoline(scalarField, 0);

	if (isoline.size() == 0)
	{
		return 0;
	}
	*CuttedMesh = CuttedMesh->SplitModelByIsoline(scalarField, 0).first;




	//切割模型
	auto vec = plane.orthogonal_vector();

	Eigen::Vector3d dir1(isoline[0][1].x() - isoline[0][0].x(), isoline[0][1].y() - isoline[0][0].y(), isoline[0][1].z() - isoline[0][0].z()), dir2(vec.x(), vec.y(), vec.z());
	dir1.normalize();
	dir2.normalize();
	Eigen::Vector3d dir3 = dir1.cross(dir2);
	dir3.normalize();
	//计算新的基向量

	struct isoPoints
	{
		Eigen::Vector3d locaInWorld;
		Eigen::Vector2d locaInPlane;
		int pId;
	};

	vector<isoPoints> isoPoint;
	Polygon_2 polygon;
	for (auto loop : isoline)
	{
		for (auto p : loop)
		{
			isoPoints np;
			np.locaInWorld = p;
			np.locaInPlane.x() = dir3.dot(p - isoline[0][0]);
			np.locaInPlane.y() = dir1.dot(p - isoline[0][0]);
			polygon.push_back(Point(np.locaInPlane.x(), np.locaInPlane.y()));
			isoPoint.push_back(np);
		}
	}
	//计算二维投影坐标

	CDT cdt;
	cdt.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);

	mark_domains(cdt);

	vector<Eigen::Vector2d> m_verts;
	vector<Eigen::Vector3i> m_faces;

	vector<Eigen::Vector3d> Last_verts = CuttedMesh->GetVertices();
	vector<Eigen::Vector3i> Last_faces = CuttedMesh->GetFaces();

	int Pnum = Last_verts.size();

	for (Face_handle f : cdt.finite_face_handles())
	{
		if (f->info().in_domain())
		{
			Eigen::Vector3i F_verts;
			map<int, bool> mp2d, mp3d;

			for (int i = 0; i < 3; i++)
			{
				Vertex_handle p = f->vertex(i);
				auto PPP = p->point();

				Eigen::Vector2d pt(PPP.x(), PPP.y());
				//F_verts.push_back(pt);
				double minDis = DBL_MAX;
				int Pid = -1;

				for (int j = 0; j < isoPoint.size(); j++)
				{
					double dis = sqrt((pt.x() - isoPoint[j].locaInPlane.x()) * (pt.x() - isoPoint[j].locaInPlane.x()) + (pt.y() - isoPoint[j].locaInPlane.y()) * (pt.y() - isoPoint[j].locaInPlane.y()));
					if (dis < minDis && mp2d.find(j) == mp2d.end())
					{
						minDis = dis;
						Pid = j;
					}
				}//二维转三维
				mp2d[Pid] = true;
				auto ThrDpoint = isoPoint[Pid].locaInWorld; // 3d  points

				minDis = DBL_MAX;
				Pid = -1;
				for (int j = 0; j < Last_verts.size(); j++)
				{
					auto v = Last_verts[j];
					double dis = sqrt((v.x() - ThrDpoint.x()) * (v.x() - ThrDpoint.x()) + (v.y() - ThrDpoint.y()) * (v.y() - ThrDpoint.y()) + (v.z() - ThrDpoint.z()) * (v.z() - ThrDpoint.z()));

					if (dis < minDis && mp3d.find(j) == mp3d.end())
					{
						minDis = dis;
						Pid = j;
					}
				}//三维查重
				mp3d[Pid] = true;
				if (i == 0)
					F_verts.x() = Pid;
				if (i == 1)
					F_verts.y() = Pid;
				if (i == 2)
					F_verts.z() = Pid;
				if (i == 2 && F_verts.z() == F_verts.y())
				{
					cout << f->vertex(1)->point() << endl;
					cout << f->vertex(2)->point() << endl;
				}

			}
			Last_faces.push_back(F_verts);
		}
	}
	//计算交界处 并合并重复的点

	*CuttedMesh = MyBaseModel(Last_verts, Last_faces);
	return 1;
}