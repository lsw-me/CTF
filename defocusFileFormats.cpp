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


//��һ�汾


void CTFFind4::read(string fileName, int numop)    //���ĳɹ�����Ϊ��ȡCTFFIND4��ʽ���ļ������Ҹ�numop��Ҳ����ͶӰ������ȡ���ݣ���Щ������ o���ɢ����
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
}  //����汾��read�Ƕ�ȡCTFFIN4C�������ļ�


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

void CTFFind4::getValues(std::vector<std::vector<float>>& values, int numop) //getValues ���� ������ֱ��ȥ���жϵ�λ����ΪCTFfind4���ļ���λ�ǰ�
{

    // ���� numop ����ͶӰ���������� originalValues �Ĵ�С�� numop ��ͬ
    for (int projIndex = 0; projIndex < numop; projIndex++)
    {
        // ȷ�� originalValues[projIndex] ���ڣ�������Ҫ����ʵ��ļ��
        values[projIndex][0] = originalValues[projIndex][1] ; // defocus #1
        values[projIndex][1] = originalValues[projIndex][2] ; // defocus #2
        values[projIndex][2] = originalValues[projIndex][3];  // azimuth of astigmatism in degrees
        values[projIndex][3] = originalValues[projIndex][4];  // phase shift in radians
    }
}
//������
void CTFFind4::printValues(const std::vector<std::vector<float>>& values) 
{
    for (unsigned int i = 0; i < values.size(); ++i)
    {
        for (unsigned int j = 0; j < values[i].size(); ++j)
        {
            std::cout << values[i][j];

            // ��ÿ��Ԫ��֮�����һ���ո񣬳������һ��Ԫ��
            if (j != values[i].size() - 1)
            {
                std::cout << ' ';
            }
        }

        // ��ÿ�н���ʱ������з�
        std::cout << std::endl;
    }

}

DefocusFileFormat* DefocusFileFormat::createFileFormat(string fileFormatName)
{

    //if (fileFormatName == "imod")  //�Ұ���Щ��ʽ�Ķ�ȥ����
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


