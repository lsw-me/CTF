#pragma once

#include"mrcStack.h"
#include"mrcStack.h"
#include "parameterSetup.h"

class myCorrection
{
public:

    myCorrection(string inputname);
    ~myCorrection();

    void run(string inputdefocusname,string outputmrcname,string mask);
    void initVariables(string inputdefocusname);
    void computeFrequencyArray();
    void correctCTF(string outmrcname, string outmaskname);
    void initMatrixWithValue(std::vector<std::vector<float> >& matrix, float value);
    void computeCTF(std::vector<float>& ctfAmp, std::vector<float>& ctfFilter, float defocus1, float defocus2, float astigmatism, float phaseShift);
    void checkAstigmatism();
    void padProjection(std::vector<float>& projection, size_t dimX, size_t dimY);
    void cropProjection(std::vector<float>& projection, size_t dimX, size_t dimY);
    void get_particle_defocusFileValues(string fileName, int numop);//这个用来得到每个粒子的散焦文件，这个文件是我自己生成的 记住两个版本read的区别
private:

    
    string inputname;
    string defocusFilename;
    string outmask;
    string outmrc;
    ParameterSetup params;
    MRCStack* inputStack;
    //ProjectionSet* projSet;

    float pixelSize;
    float amplitude;
    float cs;
    float evk;

    size_t arraySizeX;
    size_t arraySizeY;

    string ctfCorrectionType;
    string defocusFileFormat;

    bool correctAstigmatism;

    int touying;

    std::vector<std::vector<float>> particledefocusFileValues;
    std::vector<std::vector<float> > frequencyArray;

};