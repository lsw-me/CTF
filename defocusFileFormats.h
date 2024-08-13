#pragma once
#include<iostream>
#include<vector>
#include <fstream>
#include"common.h"
#include "projectionSet.h"

using namespace std;
using namespace novaCTF;



//************     �˴������numofp��ʱȷ��Ϊnumber of picture  ������б������ͼƬ����



class DefocusFileFormat
{
public:

    static DefocusFileFormat* createFileFormat(string fileFormatName);
   virtual void read(string fileName, int numop) = 0;
   virtual void getValues(std::vector<std::vector<float>>& values, int numop) = 0;
   // virtual void writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units) = 0;
};

class CTFFind4 : public DefocusFileFormat
{
public:
    void read(string fileName,int numop);
   
    std::vector<std::vector<float>>& getOriginalValues();
    void getValues(std::vector<std::vector<float>>& values, int numop);
    void printValues(const std::vector<std::vector<float>>& values);

  

   // void writeWithShiftedDefocii(std::vector<std::vector<float>>& newValues, string fileName, ProjectionSet& projSet, std::string units);
private:
    void skipComments(ifstream& infile);

    std::vector<std::vector<float>> originalValues;
};
