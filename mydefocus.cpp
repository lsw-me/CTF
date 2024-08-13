#include "mydefocus.h"
#include "exception.h"
#include "common.h"
#include <iomanip>
using namespace std;


vector<mydefoucs::Coordinate> mydefoucs::Readcoordinates(string fileName)  //��ȡ����
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

    // ����ļ�û����ȷ�رգ������һ��û�������������������׳��쳣
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

    // ����ļ��Ƿ����������������һ���Ƿ���һ����Ч�ĽǶ�ֵ
    if (!infile.eof()) {
        throw std::runtime_error("Invalid format in the file. Expected only one angle per line.");
    }
    infile.close();
    return angles;
}

void mydefoucs::Calculateparticledefocus(vector<Coordinate> Coordinate, vector<float> angles,vector<vector<float>> defocus)
{
    PixelSize=2.56;//��λ�ǰ� 
    particledefocus.resize(Coordinate.size());
    for (size_t i = 0; i < Coordinate.size(); ++i) {
        particledefocus[i].resize(angles.size(), std::vector<float>(4));
        for (size_t j = 0; j < angles.size(); ++j) {
            float angle = angles[j] * M_PI / 180.0f;;
            // ������㹫ʽ�����ĸ��������������Ϊʾ�������滻Ϊʵ�ʼ��㹫ʽ
            float result1 = ((Coordinate[i].x - 2048)* sin(angle) + Coordinate[i].z * cos(angle))+ defocus[j][0];  //������defocus1 ������
            float result2 = ((Coordinate[i].x - 2048)* sin(angle) + Coordinate[i].z * cos(angle))+ defocus[j][1]; // ������defocus2(��)
            float result3 = defocus[j][2]; // azimuth of astigmatism in degrees  ���Ӧ�ò�����λ�øı�ɣ����ɣ�
            float result4 = defocus[j][3]; // phase shift in radians���󲿷�Ϊ0��

            particledefocus[i][j][0] = result1;
            particledefocus[i][j][1] = result2;
            particledefocus[i][j][2] = result3;
            particledefocus[i][j][3] = result4;
        }//�����и����⣬���Coordinate[i].x * sin(angle) + Coordinate[i].z * cos(angle)ֱ�ӱ�ɰ�
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

        // д�����꣨���������Ӧ���������е����꣬���Ǽ�������һ���֣�
      
        // д��ÿ���Ƕȶ�Ӧ�� defocus ֵ
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
