#pragma once

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <Eigen\dense>

using namespace std;

class MyBaseModel
{
protected:
	vector<Eigen::Vector3d> m_verts;
	vector<Eigen::Vector3i> m_faces;

public:
	MyBaseModel() {}
	MyBaseModel(vector<Eigen::Vector3d> verts, vector<Eigen::Vector3i> faces)
		: m_verts(verts), m_faces(faces)
	{
	}
	vector<Eigen::Vector3d> GetVertices() const;
	vector<Eigen::Vector3i> GetFaces() const;
	void ReadObjFile(const char* filename);
	void WriteObjFile(const char* filename) const;
	vector<vector<Eigen::Vector3d>> ExtractIsoline(const vector<double>& scalarField, double val) const;
	void SaveIsoline(const char* filename, const vector<vector<Eigen::Vector3d>>&) const;
	pair<MyBaseModel, MyBaseModel> SplitModelByIsoline(const vector<double>& scalarField, double val) const;
	void WriteTexturedObjFile(const char* filename, const vector<pair<double, double>>& uvs) const;
	void WriteTexturedObjFile(const char* filename, const vector<double>& uvs) const;

	vector<int> VertsId2Component; 
};



void MyBaseModel::ReadObjFile(const char* filename)
{
	ifstream in(filename);

	char buf[1024];
	while (in.getline(buf, sizeof buf))
	{
		stringstream line(buf);
		string word;
		line >> word;
		if (word == "v")
		{
			double x, y, z;
			line >> x >> y >> z;
			m_verts.push_back(Eigen::Vector3d(x, y, z));
		}
		else if (word == "f")
		{
			vector<int> indices;
			int index;
			while (line >> index)
			{
				indices.push_back(index);
			}
			for (int i = 1; i <= indices.size() - 2; ++i)
			{
				m_faces.push_back(Eigen::Vector3i(indices[0]-1, indices[i]-1, indices[i+1]-1));
			}
		}
	}

	in.close();
}

void MyBaseModel::WriteObjFile(const char* filename) const
{
	ofstream out(filename);
	out << "# " << m_verts.size() << " vertices " << endl;
	out << "# " << m_faces.size() << " faces " << endl;

	for (auto v : m_verts)
	{
		out << "v " << setiosflags(ios::fixed) << setprecision(9) << v.transpose() << endl;
	}

	for (auto f : m_faces)
	{
		out << "f " << (f+Eigen::Vector3i(1,1,1)).transpose() << endl;
	}

	out.close();
}

vector<vector<Eigen::Vector3d>> MyBaseModel::ExtractIsoline(const vector<double>& scalarField, double val) const
{
	map<pair<int, int>, Eigen::Vector3d> pointsOnEdges;
	map<pair<int, int>, set<pair<int, int>> > fromOneEdgePoint2Neighbors;
	for (auto t : m_faces)
	{
		//get the ID of the triagle`s vertex
		int id1 = t.x();
		int id2 = t.y();
		int id3 = t.z();
		//get the coordinator of the Vertex
		auto v1 = m_verts[t.x()];
		auto v2 = m_verts[t.y()];
		auto v3 = m_verts[t.z()];
		vector<pair<int, int>> resPoints;//save the new points`s coordinate
		if ((scalarField[id1] > val) != (scalarField[id2] > val))
		{
			auto pos = pointsOnEdges.find(make_pair(min(id1, id2), max(id1, id2)));
			if (pos == pointsOnEdges.end())
			{
				double proportion = (scalarField[id2] - val) / (scalarField[id2] - scalarField[id1]);
				Eigen::Vector3d intersection = proportion * v1 + (1 - proportion) * v2;
				resPoints.push_back(make_pair(min(id1, id2), max(id1, id2)));
				pointsOnEdges[make_pair(min(id1, id2), max(id1, id2))] = intersection;
			}
			else
			{
				resPoints.push_back(pos->first);
			}
		}
		if ((scalarField[id2] > val) != (scalarField[id3] > val))
		{
			auto pos = pointsOnEdges.find(make_pair(min(id2, id3), max(id2, id3)));
			if (pos == pointsOnEdges.end())
			{
				double proportion = (scalarField[id3] - val) / (scalarField[id3] - scalarField[id2]);
				Eigen::Vector3d intersection = proportion * v2 + (1 - proportion) * v3;
				resPoints.push_back(make_pair(min(id3, id2), max(id3, id2)));
				pointsOnEdges[make_pair(min(id2, id3), max(id2, id3))] = intersection;
			}
			else
			{
				resPoints.push_back(pos->first);
			}
		}
		if ((scalarField[id3] > val) != (scalarField[id1] > val))
		{
			auto pos = pointsOnEdges.find(make_pair(min(id1, id3), max(id1, id3)));
			if (pos == pointsOnEdges.end())
			{
				double proportion = (scalarField[id1] - val) / (scalarField[id1] - scalarField[id3]);
				Eigen::Vector3d intersection = proportion * v3 + (1 - proportion) * v1;
				resPoints.push_back(make_pair(min(id1, id3), max(id1, id3)));
				pointsOnEdges[make_pair(min(id1, id3), max(id1, id3))] = intersection;
			}
			else
			{
				resPoints.push_back(pos->first);
			}
		}
		//Only store the situation that produces two intersections
		if (resPoints.size() == 2)
		{
			fromOneEdgePoint2Neighbors[resPoints[0]].insert(resPoints[1]);
			fromOneEdgePoint2Neighbors[resPoints[1]].insert(resPoints[0]);
		}
	}

	vector<vector<Eigen::Vector3d>> loops;

	while (!fromOneEdgePoint2Neighbors.empty())
	{
		vector<Eigen::Vector3d> loop;
		auto firstPoint = fromOneEdgePoint2Neighbors.begin()->first;
		auto prePoint = firstPoint;
		auto nxtPoint = *fromOneEdgePoint2Neighbors[prePoint].begin();

		fromOneEdgePoint2Neighbors[prePoint].erase(nxtPoint);
		if (fromOneEdgePoint2Neighbors[prePoint].empty())
			fromOneEdgePoint2Neighbors.erase(prePoint);

		fromOneEdgePoint2Neighbors[nxtPoint].erase(prePoint);

		if (fromOneEdgePoint2Neighbors[nxtPoint].empty())
			fromOneEdgePoint2Neighbors.erase(nxtPoint);
		loop.push_back(pointsOnEdges[prePoint]);

		while (nxtPoint != firstPoint)
		{
			prePoint = nxtPoint;
			nxtPoint = *fromOneEdgePoint2Neighbors[prePoint].begin();
			fromOneEdgePoint2Neighbors[prePoint].erase(nxtPoint);
			if (fromOneEdgePoint2Neighbors[prePoint].empty())
				fromOneEdgePoint2Neighbors.erase(prePoint);
			fromOneEdgePoint2Neighbors[nxtPoint].erase(prePoint);
			if (fromOneEdgePoint2Neighbors[nxtPoint].empty())
				fromOneEdgePoint2Neighbors.erase(nxtPoint);
			loop.push_back(pointsOnEdges[prePoint]);
		}
		loops.push_back(loop);
	}

	return loops;
}

void MyBaseModel::SaveIsoline(const char* filename, const vector<vector<Eigen::Vector3d>>& isoline) const
{
	ofstream out(filename);
	out << "g 3d_lines" << endl;
	int cnt(0);
	for (auto loop : isoline)
	{
		for (auto v : loop)
		{
			out << "v " << v.transpose() << endl;
		}		
		out << "l ";
		for (int i = 1; i <= loop.size(); ++i)
			out << i + cnt << " ";
		out << 1 + cnt << endl;

		cnt += loop.size();
	}
	out.close();
}

vector<Eigen::Vector3d> MyBaseModel::GetVertices() const
{
	return m_verts;
}

vector<Eigen::Vector3i> MyBaseModel::GetFaces() const
{
	return m_faces;
}

pair<MyBaseModel, MyBaseModel> MyBaseModel::SplitModelByIsoline(const vector<double>& scalarField, double val) const
{
	auto verts_copy = m_verts;
	vector<Eigen::Vector3i> facesInLarger;
	vector<Eigen::Vector3i> facesInLess;
	map<pair<int, int>, int> fromTwoIDs2IDofNewVert;

	set<int> uselessFaces;

	for (int index = 0; index < m_faces.size(); ++index)
	{
		int id1 = m_faces[index].x();
		int id2 = m_faces[index].y();
		int id3 = m_faces[index].z();
		auto v1 = m_verts[id1];
		auto v2 = m_verts[id2];
		auto v3 = m_verts[id3];
		bool flags[3] = { false };
		int cnt(0);
		if ((scalarField[id1] > val) != (scalarField[id2] > val))
		{
			flags[0] = true;
			++cnt;
			double proportion = (scalarField[id2] - val) / (scalarField[id2] - scalarField[id1]);
			Eigen::Vector3d intersection = proportion * v1 + (1 - proportion) * v2;
			if (fromTwoIDs2IDofNewVert.find(make_pair(min(id1, id2), max(id1, id2))) == fromTwoIDs2IDofNewVert.end())
			{
				verts_copy.push_back(intersection); // verts_copy.size() - 1;
				fromTwoIDs2IDofNewVert[make_pair(min(id1, id2), max(id1, id2))] = verts_copy.size() - 1;
			}
			else
			{

			}
			//resPoints.push_back(intersection);
		}
		if ((scalarField[id2] > val) != (scalarField[id3] > val))
		{
			flags[1] = true;
			++cnt;
			double proportion = (scalarField[id3] - val) / (scalarField[id3] - scalarField[id2]);
			Eigen::Vector3d intersection = proportion * v2 + (1 - proportion) * v3;
			if (fromTwoIDs2IDofNewVert.find(make_pair(min(id2, id3), max(id2, id3))) == fromTwoIDs2IDofNewVert.end())
			{
				verts_copy.push_back(intersection); // verts_copy.size() - 1;
				fromTwoIDs2IDofNewVert[make_pair(min(id2, id3), max(id2, id3))] = verts_copy.size() - 1;
			}
			else
			{

			}
		}
		if ((scalarField[id3] > val) != (scalarField[id1] > val))
		{
			flags[2] = true;
			++cnt;
			double proportion = (scalarField[id1] - val) / (scalarField[id1] - scalarField[id3]);
			Eigen::Vector3d intersection = proportion * v3 + (1 - proportion) * v1;
			if (fromTwoIDs2IDofNewVert.find(make_pair(min(id3, id1), max(id3, id1))) == fromTwoIDs2IDofNewVert.end())
			{
				verts_copy.push_back(intersection); // verts_copy.size() - 1;
				fromTwoIDs2IDofNewVert[make_pair(min(id3, id1), max(id3, id1))] = verts_copy.size() - 1;
			}
			else
			{

			}
		}
		if (cnt != 2)
			continue;
		uselessFaces.insert(index);
		if (flags[0] && flags[1])
		{
			int pre = fromTwoIDs2IDofNewVert[make_pair(min(id1, id2), max(id1, id2))];
			int middle = id2;
			int nxt = fromTwoIDs2IDofNewVert[make_pair(min(id2, id3), max(id2, id3))];
			if (scalarField[id2] > val)
			{
				facesInLarger.push_back(Eigen::Vector3i(pre, middle, nxt));
				facesInLess.push_back(Eigen::Vector3i(nxt, id3, id1));
				facesInLess.push_back(Eigen::Vector3i(pre, nxt, id1));
			}
			else
			{
				facesInLess.push_back(Eigen::Vector3i(pre, middle, nxt));
				facesInLarger.push_back(Eigen::Vector3i(nxt, id3, id1));
				facesInLarger.push_back(Eigen::Vector3i(pre, nxt, id1));
			}
		}
		if (flags[1] && flags[2])
		{
			int pre = fromTwoIDs2IDofNewVert[make_pair(min(id2, id3), max(id2, id3))];
			int middle = id3;
			int nxt = fromTwoIDs2IDofNewVert[make_pair(min(id1, id3), max(id1, id3))];
			if (scalarField[id3] > val)
			{
				facesInLarger.push_back(Eigen::Vector3i(pre, middle, nxt));
				facesInLess.push_back(Eigen::Vector3i(nxt, id1, id2));
				facesInLess.push_back(Eigen::Vector3i(pre, nxt, id2));
			}
			else
			{
				facesInLess.push_back(Eigen::Vector3i(pre, middle, nxt));
				facesInLarger.push_back(Eigen::Vector3i(nxt, id1, id2));
				facesInLarger.push_back(Eigen::Vector3i(pre, nxt, id2));
			}
		}
		if (flags[2] && flags[0])
		{
			int pre = fromTwoIDs2IDofNewVert[make_pair(min(id3, id1), max(id3, id1))];
			int middle = id1;
			int nxt = fromTwoIDs2IDofNewVert[make_pair(min(id1, id2), max(id1, id2))];
			if (scalarField[id1] > val)
			{
				facesInLarger.push_back(Eigen::Vector3i(pre, middle, nxt));
				facesInLess.push_back(Eigen::Vector3i(nxt, id2, id3));
				facesInLess.push_back(Eigen::Vector3i(pre, nxt, id3));
			}
			else
			{
				facesInLess.push_back(Eigen::Vector3i(pre, middle, nxt));
				facesInLarger.push_back(Eigen::Vector3i(nxt, id2, id3));
				facesInLarger.push_back(Eigen::Vector3i(pre, nxt, id3));
			}
		}
	}

	for (int index = 0; index < m_faces.size(); ++index)
	{
		if (uselessFaces.find(index) != uselessFaces.end())
			continue;
		double average = 1.0 / 3.0 * (scalarField[m_faces[index].x()] + scalarField[m_faces[index].y()] + scalarField[m_faces[index].z()]);
		if (average > val)
		{
			facesInLarger.push_back(m_faces[index]);
		}
		else
		{
			facesInLess.push_back(m_faces[index]);
		}
	}

	map<int, int> fromOldid2Newid;
	vector<Eigen::Vector3d> m_verts_less;

	for (int index = 0; index < facesInLess.size(); ++index)
	{
		auto t = facesInLess[index];
		int ids[] = { t.x(), t.y(), t.z() };
		for (int i = 0; i < 3; ++i)
		{
			int id = ids[i];
			if (fromOldid2Newid.find(id) == fromOldid2Newid.end())
			{
				int sz = fromOldid2Newid.size();
				fromOldid2Newid[id] = sz;
				m_verts_less.push_back(verts_copy[id]);
			}
		}
		facesInLess[index] = Eigen::Vector3i(fromOldid2Newid[ids[0]],
			fromOldid2Newid[ids[1]],
			fromOldid2Newid[ids[2]]);
	}

	MyBaseModel model_less(m_verts_less, facesInLess);

	vector<Eigen::Vector3d> m_verts_large;
	fromOldid2Newid.clear();

	for (int index = 0; index < facesInLarger.size(); ++index) 
	{
		auto t = facesInLarger[index];
		int ids[] = { t.x(),t.y(),t.z() };
		for (int i = 0; i < 3; i++) {
			int id = ids[i];
			if (fromOldid2Newid.find(id) == fromOldid2Newid.end()) {
				int sz = fromOldid2Newid.size();
				fromOldid2Newid[id] = sz;
				m_verts_large.push_back(verts_copy[id]);
			}
		}
		facesInLarger[index] = Eigen::Vector3i(fromOldid2Newid[ids[0]],
			fromOldid2Newid[ids[1]],
			fromOldid2Newid[ids[2]]);
	}
	MyBaseModel model_lager(m_verts_large, facesInLarger);
	
	return make_pair(model_less, model_lager);
}

void MyBaseModel::WriteTexturedObjFile(const char* filename, const vector<pair<double, double>>& uvs) const
{
	ofstream out(filename);
	out << "# " << m_verts.size() << " vertices " << endl;
	out << "# " << m_faces.size() << " faces " << endl;
	out << "mtllib defaultmaterial.mtl" << endl << "usemtl mydefault" << endl;
	for (auto v : m_verts)
	{
		out << "v " << setiosflags(ios::fixed) << setprecision(9) << v.transpose() << endl;
	}

	for (auto uv : uvs)
	{
		out << "vt " << setiosflags(ios::fixed) << setprecision(9) << uv.first << " " << uv.second << endl;
	}

	for (auto f : m_faces)
	{
		auto ids = (f + Eigen::Vector3i(1, 1, 1)).transpose();
		out << "f " << ids.x() << "/" << ids.x() << " " << ids.y() << "/" << ids.y() << " " << ids.z() << "/" << ids.z() << endl;
	}

	out.close();
}

void  MyBaseModel::WriteTexturedObjFile(const char* filename, const vector<double>& uvs) const
{
	ofstream out(filename);
	out << "# " << m_verts.size() << " vertices " << endl;
	out << "# " << m_faces.size() << " faces " << endl;
	out << "mtllib defaultmaterial.mtl" << endl << "usemtl mydefault" << endl;
	for (auto v : m_verts)
	{
		out << "v " << setiosflags(ios::fixed) << setprecision(9) << v.transpose() << endl;
	}

	for (auto uv : uvs)
	{
		out << "vt " << setiosflags(ios::fixed) << setprecision(9) << uv << " " << 0 << endl;
	}

	for (auto f : m_faces)
	{
		auto ids = (f + Eigen::Vector3i(1, 1, 1)).transpose();
		//cerr << ids << endl;
		out << "f " << ids.x() << "/" << ids.x() << " " << ids.y() << "/" << ids.y() << " " << ids.z() << "/" << ids.z() << endl;
	}

	out.close();
}