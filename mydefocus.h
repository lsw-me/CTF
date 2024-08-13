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
	//传入坐标，角度，倾斜序列每一张图片的defocus 
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