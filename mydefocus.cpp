#include "mydefocus.h"
#include "exception.h"
#include "common.h"
#include <iomanip>
using namespace std;


vector<mydefoucs::Coordinate> mydefoucs::Readcoordinates(string fileName)  //读取坐标
{
    ifstream infile;
    infile.open(fileName.c_str());

    if (!infile.good())
    {
        throw ExceptionFileOpen(fileName.c_str());
    }
    vector<Coordinate> coordinates;
    Coordinate coord;

    while (infile >> coord.x >> coord.y >> coord.z) {
        coordinates.push_back(coord);
    }

    // 如果文件没有正确关闭（即最后一行没有三个浮点数），则抛出异常
    if (!infile.eof()) {
        throw std::runtime_error("Invalid format in the file.");
    }

    infile.close();

    return coordinates;
}

vector<float> mydefoucs::Readangle(string fileName)
{
    ifstream infile;
    infile.open(fileName.c_str());
    if (!infile.good())
    {
        throw ExceptionFileOpen(fileName.c_str());
    }
    vector<float> angles;
    
    float angle;
    while (infile >> angle) {
        angles.push_back(angle);
    }

    // 检查文件是否正常结束，即最后一行是否是一个有效的角度值
    if (!infile.eof()) {
        throw std::runtime_error("Invalid format in the file. Expected only one angle per line.");
    }
    infile.close();
    return angles;
}

void mydefoucs::Calculateparticledefocus(vector<Coordinate> Coordinate, vector<float> angles,vector<vector<float>> defocus)
{
    PixelSize=2.56;//单位是埃 
    particledefocus.resize(Coordinate.size());
    for (size_t i = 0; i < Coordinate.size(); ++i) {
        particledefocus[i].resize(angles.size(), std::vector<float>(4));
        for (size_t j = 0; j < angles.size(); ++j) {
            float angle = angles[j] * M_PI / 180.0f;;
            // 假设计算公式包含四个分量，这里仅作为示例，请替换为实际计算公式
            float result1 = ((Coordinate[i].x - 2048)* sin(angle) + Coordinate[i].z * cos(angle))+ defocus[j][0];  //计算后的defocus1 （埃）
            float result2 = ((Coordinate[i].x - 2048)* sin(angle) + Coordinate[i].z * cos(angle))+ defocus[j][1]; // 计算后的defocus2(埃)
            float result3 = defocus[j][2]; // azimuth of astigmatism in degrees  这个应该不会随位置改变吧（存疑）
            float result4 = defocus[j][3]; // phase shift in radians（大部分为0）

            particledefocus[i][j][0] = result1;
            particledefocus[i][j][1] = result2;
            particledefocus[i][j][2] = result3;
            particledefocus[i][j][3] = result4;
        }//这里有个问题，这个Coordinate[i].x * sin(angle) + Coordinate[i].z * cos(angle)直接变成埃
    }

}

void mydefoucs::WriteResultsToFile(const std::string& baseFilePath, const std::vector<vector<vector<float>>>& particledefocus)
{
    bool success = true;
    for (size_t i = 0; i < particledefocus.size(); ++i) {
        std::stringstream ss;
        ss << baseFilePath << "defocus_" << i << ".defcous";

        std::ofstream outputFile(ss.str());
        if (!outputFile.is_open()) {
            throw std::runtime_error("Failed to open output file.");
        }

        // 写入坐标（这里的坐标应该是输入中的坐标，不是计算结果的一部分）
      
        // 写入每个角度对应的 defocus 值
        for (size_t j = 0; j < particledefocus[i].size(); ++j) {
            outputFile << std::fixed << std::setprecision(6);
            outputFile << particledefocus[i][j][0] << " "  // result1
                << particledefocus[i][j][1] << " "  // result2
                << particledefocus[i][j][2] << " "  // result3
                << particledefocus[i][j][3] << "\n"; // result4
        }

      
        outputFile.close();
    }
   
    if (success) 
    {
       cout << "All files written successfully." << endl;
    }
    else {
       cout << "Some files failed to be written." << endl;
    }
}
