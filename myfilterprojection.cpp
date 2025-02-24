#include "myfilterprojection.h"
#include "common.h"
#include "volumeIO.h"
#include <math.h>
#include <algorithm>
#include <float.h>
#include <sstream>
#include "exception.h"
#include <fstream>
#include "fftRoutines.h"

using namespace novaCTF;
using namespace std;



FilterProjections::FilterProjections(string inputname)//这是要求指定输入的mrc文件
{
    //params = aParams;
    inputStack = new MRCStack(inputname, true, false, params.KeepFilesOpen());//params.KeepFilesOpen()这个要不要改存疑
}

FilterProjections::~FilterProjections()
{
    delete inputStack;
}

void FilterProjections::run(string inputname)
{
    AngleFilename = inputname;
    initParameters();
    filter();
    string outputname = "E:\\VS studio\\proj\\project2\\test_filter_output\\filteroutput.mrc";
    VolumeIO::writeMRCStack(inputStack->getStackHeader(), filteredProjections, outputname, inputStack->getExtraData());
    cout << "done" << endl;
    //outputfilename要改一改
}

void FilterProjections::loadTiltAngles(string inputname)//这个要改 这里指定角度文件
{
    ifstream infile;
    infile.open(inputname.c_str());

    if (!infile.good())
    {
        throw ExceptionFileOpen(inputname);
    }

    if (params.StackOrientation() == novaCTF::VolumeRotation::ALONG_XY)//这个不知到有没有用
        numberOfAngles = inputStack->getResolution().z;
    else
        numberOfAngles = inputStack->getResolution().y;

    angles.resize(numberOfAngles);

    float degreesToRadiansFactor = M_PI / 180.0;

    for (size_t i = 0; i < numberOfAngles; i++)
    {
        infile >> angles[i];
        cout << angles[i]<<endl;
        angles[i] = -degreesToRadiansFactor * (angles[i]);
    }

    infile.close();
}

void FilterProjections::initParameters()//这个好像不用改  读取也没问题
{
    projRes.x = (int)inputStack->getResolution().x;
    projRes.y = (int)inputStack->getResolution().y;
    projRes.z = (int)inputStack->getResolution().z;

    // Set up padding: 10% of X size or minimum of 16, max of 50  设置填充：X 尺寸的 10% 或最小 16，最多 50
    npad = min(50, 2 * max(8, projRes.x / 20));
    paddedResX = projRes.x + npad + 2;				//2 is because of FFT routines here 2是因为这里的FFT例程

    loadTiltAngles(AngleFilename);

    // Set up radial weighting
    int irmax;
    int ifall;

    if (1)//判断是否使用用户设置的参数 这里改为 0 if (params.UseRadialFilter())
    {
        irmax = 0.3;
        ifall = 0.05;
        if (irmax == 0)  //this is little bit weird - this can happen if cutoff is smaller than 0.5 due to int cast
            irmax = projRes.x * 0.03;

        if (ifall == 0)
            ifall = projRes.x * 0.05;
    }
    else    // Default radial weighting parameters - no smooth cut-off at the edges  默认径向加权参数 - 边缘无平滑截止
    {
        irmax = inputStack->getResolution().x / 2 + 1;
        ifall = 0.0;
    }

    radialWeighting(irmax, ifall);

}

// Main loop over slices perpendicular to tilt axis
void FilterProjections::filter()
{
    filteredProjections.resize(projRes.z * projRes.y * projRes.x);

    inputStack->readAllProjections(&filteredProjections[0]);

    for (int sliceNumber = 0; sliceNumber < projRes.z; sliceNumber++)
    {
        currentRows.clear();

        for (int y = 0; y < projRes.y; y++)
        {
            for (int x = 0; x < projRes.x; x++)
            {
                currentRows.push_back(filteredProjections[x + y * projRes.x + sliceNumber * projRes.x * projRes.y]);
            }
            taperEndToStart(y);
        }

        transform();

        for (int y = 0; y < projRes.y; y++)
        {
            for (int x = 0; x < projRes.x; x++)
                filteredProjections[x + y * projRes.x + sliceNumber * projRes.x * projRes.y] = currentRows[x + y * paddedResX];
        }
    }
}

void FilterProjections::transform()
{
    // Apply forward Fourier transform
    FFTRoutines::many1DTransform(&currentRows[0], (int)projRes.x + (int)npad, (int)numberOfAngles, 0);

    // Apply Radial weighting
    int index = 0;
    for (unsigned int nv = 0; nv < numberOfAngles; nv++)
    {
        for (int i = 0; i < paddedResX; i++)
        {
            currentRows[index] = currentRows[index] * radialFilter[index];
            index++;
        }
    }

    // Apply inverse transform
    FFTRoutines::many1DTransform(&currentRows[0], projRes.x + npad, numberOfAngles, 1);
}


//当前行的投影数据进行频域滤波处理，首先进行正向傅里叶变换，然后应用径向权重滤波，最后进行逆向傅里叶变换得到滤波后的空间域数据
void FilterProjections::taperEndToStart(unsigned int projNumber)
{
    int nsum = 0;
    float xsum = 0.0f;
    size_t projOffset = projNumber * paddedResX;

    for (int ix = 0; ix <= min(2, (int)projRes.x - 1); ix++)
    {
        nsum++;
        xsum += currentRows[projOffset + ix];
    }
    float stmean = xsum / (float)nsum;

    nsum = 0;
    xsum = 0.0f;

    for (int ix = max(0, (int)(projRes.x - 3)); ix <= (int)(projRes.x - 1); ix++)
    {
        nsum++;
        xsum += currentRows[projOffset + ix];
    }

    float endmean = xsum / (float)nsum;

    for (int ipad = 0; ipad < npad; ipad++)
    {
        float value = (ipad + 1) / (npad + 1.0f);
        currentRows.push_back(value * stmean + (1.0f - value) * endmean);
    }

    currentRows.push_back(0.0f);
    currentRows.push_back(0.0f);

}

// Fake SIRT-like filter
// All the values here correspond to the values from tilt routine
// in IMOD package
// 假的类似 SIRT 的过滤器此处的所有值都对应于 IMOD 封装中倾斜例程中的值
//它的主要功能是对当前投影切片（projection）的两端进行渐变填充（tapering），以平滑投影图像的边缘
unsigned int FilterProjections::computeSirtIterations()
{
    int finalFakeSirtIterations = 0;

    if (0) //这里设置为0 if (params.UseFakeSirtIterations())
    {
        if (params.FakeSirtIterations() <= 0)
        {
            finalFakeSirtIterations = 0;
        }
        else if (params.FakeSirtIterations() > 30)
        {
            finalFakeSirtIterations = 27 + 0.6 * (params.FakeSirtIterations() - 30);
        }
        else if (params.FakeSirtIterations() > 15)
        {
            finalFakeSirtIterations = 15 + 0.8 * (params.FakeSirtIterations() - 15);
        }
        else
        {
            finalFakeSirtIterations = params.FakeSirtIterations();
        }
    }

    return (unsigned int)finalFakeSirtIterations;
}

// Set Radial Transform weighting
// Linear ramp plus Gaussian fall off
void FilterProjections::radialWeighting(float cutOff, float fallOff)
{
    std::vector<float> wincr;
    wincr.push_back(2.0f);
    wincr.push_back(1.0f / 1.5f);

    std::vector<float> wgtAngles;
    wgtAngles.resize(numberOfAngles);

    unsigned int fakeSirtIterations = computeSirtIterations();

    //constant values for fake SIRT filter - ugly but corresponds to the tilt routine from IMOD
    float fakeMatchAdd = 0.3f;
    float fakeAlpha = 0.00195f;

    for (unsigned int i = 0; i < numberOfAngles; i++)
    {
        wgtAngles[i] = angles[i];
    }

    std::vector<float> wgtAtten;
    wgtAtten.resize(numberOfAngles);

    int iend = paddedResX / 2;
    float stretch = float(projRes.x + npad) / projRes.x;

    int irmax = round(cutOff * stretch);
    int ifall = round(fallOff * stretch);
    float avgint = 1.0;
    float attensum = 0.0f;

    if (numberOfAngles > 1)
    {
        avgint = (wgtAngles[numberOfAngles - 1] - wgtAngles[0]) / (numberOfAngles - 1);
    }

    float atten;

    // Set up the attenuations for the weighting angles
    for (int iv = 0; iv < (int)numberOfAngles; iv++)
    {
        atten = 1.0;
        if (numberOfAngles > 1)
        {
            float sumint = 0;
            float wsum = 0.0f;
            for (int iw = 0; iw < 2; iw++)
            {
                if ((iv - iw) > 0)
                {
                    wsum = wsum + wincr[iw];
                    sumint = sumint + wincr[iw] * (wgtAngles[iv - iw] - wgtAngles[iv - iw - 1]);
                }

                if ((unsigned int)(iv + iw + 1) < numberOfAngles)
                {
                    wsum = wsum + wincr[iw];
                    sumint = sumint + wincr[iw] * (wgtAngles[iv + iw + 1] - wgtAngles[iv + iw]);
                }
            }
            atten = atten * (sumint / wsum) / avgint;
        }
        wgtAtten[iv] = atten;
    }

    // Set up linear ramp
    for (unsigned int iv = 0; iv < numberOfAngles; iv++)
    {
        // Get weighting from nearest weighting angle
        atten = 1.0f;
        if (numberOfAngles > 1)
        {
            float diffmin = FLT_MAX;
            for (int iw = 0; iw < (int)numberOfAngles; iw++)
            {
                float diff = fabs(angles[iv] - wgtAngles[iw]);

                if (diff < diffmin)
                {
                    diffmin = diff;
                    atten = wgtAtten[iw];
                }
            }
            //cout<< "angle: " << angles[iv] << ", atten: " << atten << endl;

        }
        attensum = attensum + atten;


        for (int i = 0; i < min(irmax, iend); i++)
        {
            float frequency = (float)(i) / (float)(paddedResX);

            // Value 0.2 based on Kak & Slaney book
            if (i == 0)
            {
                radialFilter.push_back(atten * 0.2f);
                radialFilter.push_back(atten * 0.2f);
            }
            else if (params.UseFakeSirtIterations() && frequency >= fakeAlpha)
            {
                float fakeSirtValue = atten * i * (1.0f - pow((1.0f - fakeAlpha / frequency), (fakeSirtIterations + fakeMatchAdd)));
                radialFilter.push_back(fakeSirtValue);
                radialFilter.push_back(fakeSirtValue);
            }
            else if (params.UseSteepRampFilter())
            {
                radialFilter.push_back(atten * (2 * i - 1));
                radialFilter.push_back(atten * 2 * i);
            }
            else
            {
                radialFilter.push_back(atten * i);
                radialFilter.push_back(atten * i);
            }

        }

        for (int i = irmax; i < iend; i++)
        {
            radialFilter.push_back(0.0f);		//filling to allocate for the final loop
            radialFilter.push_back(0.0f);
        }
    }

    attensum = attensum / numberOfAngles;
    cout << "Attenuation sum:" << attensum << endl;

    // Set up Gaussian
    for (int i = irmax; i < iend; i++)
    {
        float arg = (float)(i + 1 - irmax) / (float)(ifall);
        atten = exp(-arg * arg);
        int ibase = 0;
        for (unsigned int iv = 0; iv < numberOfAngles; iv++)
        {
            int radInd = ibase + 2 * irmax - 1;
            float radValue = radialFilter[radInd];
            float z = atten * radValue;
            radialFilter[ibase + 2 * i] = z;
            radialFilter[ibase + 2 * i + 1] = z;
            ibase = ibase + paddedResX;
        }
    }

}
