#pragma once
#include<vector>
#include<iostream>
#include<fstream>
using namespace std;
class Knn
{

public:
	int VerNum, k;
	vector<vector<int>> neighboor;
	Knn(string FILE, int n, int knum, bool IfSelf)
	{
		//k = knum;
		VerNum = n;
		for (int i = 0; i < VerNum; i++)
		{
			vector<int> tmp;
			neighboor.push_back(tmp);
		}

		ifstream in(FILE);
		//int iter=0;
		for (int i = 0; i < VerNum; i++)
		{
			int thisk = 0;
			in >> thisk;

			int q;
			vector<int> tmp;
			if (IfSelf)
			{
				for (int j = 0; j < thisk; j++)
				{
					in >> q;
					if (q != i)
						tmp.push_back(q);
				}
				//in >> q;
			}
			else
			{
				for (int j = 0; j < thisk; j++)
				{
					in >> q;
					tmp.push_back(q);
				}

			}

			neighboor[i] = tmp;
		}

	}

};

