#pragma once

#include <fstream>
#include <vector>
#include <string>
using namespace std;

class mydefoucs
{
public:
	struct Coordinate 
	{
		float x, y, z;
	};
	void Calculateparticledefocus(vector<Coordinate> Coordinate,vector<float> angles, vector<vector<float>> defocus);
	void WriteResultsToFile(const std::string& filePath, const std::vector<vector<vector<float>>>& particledefocus);
	//�������꣬�Ƕȣ���б����ÿһ��ͼƬ��defocus 
	vector<Coordinate> Readcoordinates(string fileName);
	vector<float> Readangle(string fileName);
	const vector<vector<vector<float>>>& getResults() 
		const {
		return particledefocus;
	}

	double PixelSize;
	private:
		vector<vector<vector<float>>> particledefocus;
};