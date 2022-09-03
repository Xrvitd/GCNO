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
class MyPointCloudModel
{
protected:
	vector<Eigen::Vector3d> m_verts;
	vector<Eigen::Vector3d> m_normals;

public:
	MyPointCloudModel() {}
	MyPointCloudModel(vector<Eigen::Vector3d> verts, vector<Eigen::Vector3d> normals)
		: m_verts(verts), m_normals(normals)
	{
	}
	vector<Eigen::Vector3d> GetVertices() const;
	vector<Eigen::Vector3d> GetNormals() const;
	void ReadXYZFile(const char* filename,bool WithNor);
	void WriteXYZFile(const char* filename, bool WithNor) const;
	


};

void MyPointCloudModel::ReadXYZFile(const char* filename, bool WithNor)
{
	ifstream in(filename);
	if (WithNor)
	{
		Eigen::Vector3d Point, Normal;
		while (in >> Point.x() >> Point.y() >> Point.z() >> Normal.x() >> Normal.y() >> Normal.z())
		{
			m_verts.push_back(Point);
			m_normals.push_back(Normal);
		}

	}
	else
	{
		Eigen::Vector3d Point;
		while (in >> Point.x() >> Point.y() >> Point.z())
		{
			m_verts.push_back(Point);
		}
	}
	in.close();

}
void MyPointCloudModel::WriteXYZFile(const char* filename, bool WithNor) const
{
	ofstream out(filename);
	if(WithNor)
	{
		int n = m_verts.size();
		for (int i = 0; i < n; i++)
		{
			out << setiosflags(ios::fixed) << setprecision(9) << m_verts[i].transpose() <<" "<< m_normals[i].transpose() << endl;
		}
	}
	else
	{
		int n = m_verts.size();
		for (int i = 0; i < n; i++)
		{
			out << setiosflags(ios::fixed) << setprecision(9) << m_verts[i].transpose() << endl;
		}
	}
	out.close();
}




vector<Eigen::Vector3d> MyPointCloudModel::GetVertices() const
{
	return m_verts;
}

vector<Eigen::Vector3d> MyPointCloudModel::GetNormals() const
{
	return m_normals;
}