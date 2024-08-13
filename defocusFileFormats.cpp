#include"defocusFileFormats.h"
#include <sstream>
#include <iomanip>
#include "exception.h"



/* Assuming format from CTFFind4 (version 4.0.16):
 * First few lines are comments and are skipped
 * 1 column: tilt number
 * 2 column: defocus #1 in Angstroms
 * 3 column: defocus #2 in Angstroms
 * 4 column: azimuth of astigmatism in degrees
 * 5 column: phase shift in radians
 * 6 column: cross correlation
 * 7 column: spacing (in Angstroms) up to which CTF rings were fit successfully
 */


//第一版本


void CTFFind4::read(string fileName, int numop)    //更改成功，改为读取CTFFIND4格式的文件，并且跟numop（也就是投影数量读取内容）这些代表了 o点的散焦。
{
    ifstream infile;
    infile.open(fileName.c_str());

    if (!infile.good())
    {
        throw ExceptionFileOpen(fileName.c_str());
    }

    skipComments(infile);

    originalValues.resize(numop);

    for (int j=0; j <numop ;j++)
    {
        unsigned int projIndex = j;

        for (unsigned int i = 0; i < 7; i++)
        {
            float value;
            infile >> value;
            originalValues[projIndex].push_back(value);
        }
    }

    infile.close();
}  //这个版本的read是读取CTFFIN4C产生的文件


void CTFFind4::skipComments(ifstream& infile)
{
    std::string line;
    unsigned int linesToSkip = 0;

    while (getline(infile, line))
    {
        if (line[0] == '#')
            linesToSkip++;
        else
            break;
    }

    infile.seekg(0, infile.beg);

    for (unsigned int i = 0; i < linesToSkip; i++)
    {
        string line;
        getline(infile, line);
    }
}
std::vector<std::vector<float>>& CTFFind4::getOriginalValues() 
{
    return originalValues;
}

void CTFFind4::getValues(std::vector<std::vector<float>>& values, int numop) //getValues 结束 这里我直接去掉判断单位，因为CTFfind4的文件单位是埃
{

    // 假设 numop 代表投影的数量，且 originalValues 的大小与 numop 相同
    for (int projIndex = 0; projIndex < numop; projIndex++)
    {
        // 确保 originalValues[projIndex] 存在，否则需要添加适当的检查
        values[projIndex][0] = originalValues[projIndex][1] ; // defocus #1
        values[projIndex][1] = originalValues[projIndex][2] ; // defocus #2
        values[projIndex][2] = originalValues[projIndex][3];  // azimuth of astigmatism in degrees
        values[projIndex][3] = originalValues[projIndex][4];  // phase shift in radians
    }
}
//检查输出
void CTFFind4::printValues(const std::vector<std::vector<float>>& values) 
{
    for (unsigned int i = 0; i < values.size(); ++i)
    {
        for (unsigned int j = 0; j < values[i].size(); ++j)
        {
            std::cout << values[i][j];

            // 在每个元素之间输出一个空格，除了最后一个元素
            if (j != values[i].size() - 1)
            {
                std::cout << ' ';
            }
        }

        // 在每行结束时输出换行符
        std::cout << std::endl;
    }

}

DefocusFileFormat* DefocusFileFormat::createFileFormat(string fileFormatName)
{

    //if (fileFormatName == "imod")  //我把这些格式的都去掉了
    //   return new ImodCTFPlotter;
    if (fileFormatName == "ctffind4")
        return new CTFFind4;
   // if (fileFormatName == "gctf")
   //     return new GCTF;
    else
    {
        cout << "The following defocus file format is unknown: " << fileFormatName << std::endl;
        exit(EXIT_FAILURE);
    }
}


