#include "volumeIO.h"
#include "parameterSetup.h"
#include "mrcHeader.h"
#include "geomHeader.h"
#include <fstream>
#include <iostream>
#include "exception.h"
#include <stdio.h>
#include <cstring>
#include <cstdio>
#include <cassert>

using namespace std;


void VolumeIO::writeMRCStack(MRCHeader header, std::vector<std::vector<float>>& data, string outputName, char* extraData)
{
    std::vector<float> linearData;//����MRCHeader �ṹ�� header����ά float ���� vector ���� data������ļ����ַ��� outputName���Լ���������ָ�� extraData��
    linearData.resize(data.size() * data[0].size());
    size_t k = 0;

    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[0].size(); j++)
        {
            linearData[k] = data[i][j];
            k++;
        }
    }

    writeMRCStack(header, linearData, outputName, extraData);
    //����ά vector ת��Ϊһά vector linearData��Ȼ�������һ�� writeMRCStack ��������һά���ݡ�
}

void VolumeIO::writeMRCStack(MRCHeader header, std::vector<float>& data, string outputName, char* extraData)
{
    //���������MRCHeader �ṹ�� header��һά float ���� vector ���� data������ļ����ַ��� outputName���Լ���������ָ�� extraData��
    std::ofstream outfile; 

    try
    {
        outfile.open(outputName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        //��һ������ļ��� outfile�������Զ����ƺͽض�ģʽ��ָ��������ļ� outputName��

        if (!outfile.good())
        {
            throw new ExceptionFileOpen(outputName);//���ļ���ʧ�ܣ��׳� ExceptionFileOpen �쳣����ӡ������Ϣ��
        }

        header.mode = 2;//���� MRCHeader �е� mode �ֶ�Ϊ 2����ʾ�������� 32 λ��������ʽ�洢��
        computeMinMaxMean(data, data.size(), header.dMin, header.dMax, header.dMean);//ʹ�� computeMinMaxMean ���������������ݵ���Сֵ�����ֵ��ƽ��ֵ�������µ� header �ṹ���С�

        outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));//�� header �ṹ��д�뵽�Ѵ򿪵�����ļ��С�

        //if (header.extra != 0)//�� header �еĶ������ݳ��ȣ�extra �ֶΣ���Ϊ�㣬�� extraData ָ�������д�뵽�ļ��С�
        //{
        //    outfile.write(reinterpret_cast<char*>(extraData), header.extra);
       // }

        writeProjections(outfile, data, header.nx * header.ny, header.nz);//���� writeProjections ���������ݰ��� header �е� nx��ny �� nz �����ά��˳��д�뵽�ļ��С�

        outfile.close();//�ر�����ļ�����

    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
    catch (ExceptionFileFormat& e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        cout << "Error while writing volume to disk!" << endl;
    }
    //�������̰����쳣�����������ļ���ʧ�ܻ������ļ���ʽ�쳣ʱ���Ჶ���Ӧ���쳣����ӡ������Ϣ��
}

void VolumeIO::writeProjections(std::ofstream& outfile, std::vector<float>& data, size_t projectionSize, unsigned int numberOfProjections)
{
    float* currentProjection = new float[projectionSize];

    for (unsigned int z = 0; z < numberOfProjections; z++)
    {
        for (size_t i = 0; i < projectionSize; i++)
        {
            currentProjection[i] = data[i + z * projectionSize];
        }

        outfile.write((const char*)currentProjection, projectionSize * sizeof(float));
    }

    delete[] currentProjection;

}


void VolumeIO::write(std::vector<float>& data, novaCTF::Vec3ui resolution, std::string outputVolumeFileName, int mode)
{
    novaCTF::VolumeRotation rotation = novaCTF::VolumeRotation::ALONG_XZ;
    write(data, resolution, outputVolumeFileName, mode, rotation);
}


//��ȡ����ջ�ļ���MRCͷ��Ϣ��������תģʽ���µ�����ߴ������е�������󽫵��������ͷ��Ϣд�뵽һ���µ���ʱ�������ļ��С�
void VolumeIO::writeHeader(std::string inputStackFileName, std::string outputVolumeFileName, novaCTF::Vec3ui volumeDimensions, int mode, novaCTF::VolumeRotation rotation)
{
    //����ջ�ļ���������Ϊ�ַ����� outputVolumeFileName: �������ļ���ģ�壬����Ϊ�ַ��� volumeDimensions: ����������ߴ磬����ΪnovaCTF::Vec3ui����ʾ(x, y, z)ά�ȡ�
    //mode: ������2 
    // rotation: ��תģʽ������novaCTF::VolumeRotationö�����ͣ�����ָ�����������������XZ���XY����ת��
    ifstream infile;
    infile.open(inputStackFileName.c_str(), std::ios::binary);

    if (!infile.good())
    {
        throw ExceptionFileOpen(inputStackFileName);
    }

    MRCHeader header;

    infile.read((char*)&header, sizeof(MRCHeader));

    if (rotation == novaCTF::VolumeRotation::ALONG_XZ)
    {
        float scalingFactor = header.cellDimX / (float)header.mx;

        header.nx = volumeDimensions.x;
        header.ny = volumeDimensions.z;
        header.nz = volumeDimensions.y;

        header.mx = volumeDimensions.x;
        header.my = volumeDimensions.z;
        header.mz = volumeDimensions.y;

        //we go from stack to volume here - if the stack was binned we assume the same bin factor for z-dimension as well
        header.cellDimY = (float)volumeDimensions.z * scalingFactor;
    }
    else //if (rotation==ALONG_XY)
    {
        float scalingFactor = header.cellDimX / (float)header.mx;

        header.nx = volumeDimensions.x;
        header.ny = volumeDimensions.y;
        header.nz = volumeDimensions.z;

        header.mx = volumeDimensions.x;
        header.my = volumeDimensions.y;
        header.mz = volumeDimensions.z;

        //we go from stack to volume here - if the stack was binned we assume the same bin factor for z-dimension as well
        header.cellDimZ = (float)volumeDimensions.z * scalingFactor;
    }

    std::ofstream outfile;
    outfile.open(generateTempName(outputVolumeFileName, mode).c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

    if (!outfile.good())
    {
        throw new ExceptionFileOpen(generateTempName(outputVolumeFileName, mode));
    }

    outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));       //�������ǸĲ��԰�������
    outfile.close();
    infile.close();
}

std::string VolumeIO::generateTempName(std::string filename, int mode)
{
    std::stringstream newFilename;

    newFilename << filename;

    if (mode != 2)
        newFilename << "_temp";

    return newFilename.str();

}

//������ݴ���Ĳ�����������ݽ���ֵ��ת��������ת��������ݱ��浽ָ��������ļ��С�������Ҫ��Ϊ���¼������裺
//����������Ǵ����ֵ��Ϊ0ǰ�������������
void VolumeIO::convertVolumeValues(std::string outputVolumeFileName, float min, float max, float mean, int mode)
{
    float newMin;
    float newMax;
    //modeΪ2ʱ������ԭ��Сֵ�����ֵ���䣻 
    switch (mode)
    {
    case 0: newMin = 0.f;
        newMax = 255.f;
        break;
    case 1: newMin = -32768.f;
        newMax = 32767.f;
        break;
    case 2: newMin = min;
        newMax = max;
        break;
    case 6: newMin = 0.f;
        newMax = 65535.f;
        break;
    default:	throw new ExceptionFileFormat();
        break;
    }

    float newMean = novaCTF::convertRange(mean, min, max, newMin, newMax);

    MRCHeader header;
    correctHeader(outputVolumeFileName, header, newMin, newMax, newMean, mode);
    //��ȡ����ļ���MRCͷ��Ϣ��MRCHeader�ṹ�壩��������correctHeader������MRCͷ����������ʹ�����µ�ֵ��Χƥ�䡣
    size_t sliceSize = header.nx * header.ny;
    //�����������ÿһ��Ƭ�Ĵ�С��sliceSize����
    //����modeѡ����Ӧ��ģ���������convertSliceData������
    // �ú������������Ƭ�ض�ȡԭʼ���ݲ�����ת�����µ�ֵ��Χ�ڣ�Ȼ��д�뵽����ļ��С����ڲ�ͬmodeֵ��ʹ�õ��������Ͳ�ͬ���ֱ�Ϊ
    switch (mode)
    {
    case 0: convertSliceData<unsigned char>(outputVolumeFileName, sliceSize, header.nz, min, max, newMin, newMax);
        break;
    case 1: convertSliceData<signed short>(outputVolumeFileName, sliceSize, header.nz, min, max, newMin, newMax);
        break;
    case 2: break;
    case 6: convertSliceData<unsigned short>(outputVolumeFileName, sliceSize, header.nz, min, max, newMin, newMax);
        break;
    default:	throw new ExceptionFileFormat();
        break;
    }
}

template <typename voxelType>
void VolumeIO::convertSliceData(std::string outputVolumeFileName, size_t sliceSize, unsigned int sliceNumber, float min, float max, float newMin, float newMax)
{
    ifstream originalFile;
    ofstream newFile;

    originalFile.open(generateTempName(outputVolumeFileName, 1).c_str(), std::ios_base::in | std::ios_base::binary);
    newFile.open(outputVolumeFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);

    if (!newFile.good())
    {
        throw ExceptionFileOpen(outputVolumeFileName);
    }

    if (!originalFile.good())
    {
        throw ExceptionFileOpen(generateTempName(outputVolumeFileName, 1).c_str());
    }

    originalFile.seekg(sizeof(MRCHeader));
    newFile.seekp(sizeof(MRCHeader));

    voxelType* convertedData = new voxelType[sliceSize];
    float* floatData = new float[sliceSize];

    for (unsigned int slice = 0; slice < sliceNumber; slice++)
    {
        originalFile.read((char*)(floatData), sliceSize * sizeof(float));

        for (size_t i = 0; i < sliceSize; i++)
        {
            convertedData[i] = static_cast<voxelType>(novaCTF::convertRange(floatData[i], min, max, newMin, newMax));
        }

        newFile.write((const char*)convertedData, sliceSize * sizeof(voxelType));
    }

    delete[] floatData;
    delete[] convertedData;

    newFile.close();
    originalFile.close();

    std::remove(generateTempName(outputVolumeFileName, 1).c_str());
}

void VolumeIO::correctHeader(std::string outputVolumeFileName, MRCHeader& header, float min, float max, float mean, int mode)
{
    fstream file;

    file.open(generateTempName(outputVolumeFileName, mode).c_str(), std::ios_base::in | std::ios_base::out | std::ios::binary);

    if (!file.good())
    {
        throw ExceptionFileOpen(generateTempName(outputVolumeFileName, mode));
    }

    file.seekg(0);
    file.read((char*)&header, sizeof(MRCHeader));

    header.dMin = min;
    header.dMax = max;
    header.dMean = mean;

    header.mode = mode;
    header.extra = 0;

    if (mode == 2)
    {
        file.seekp(0);
        file.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));
    }
    else
    {
        ofstream newFile;
        newFile.open(outputVolumeFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

        if (!newFile.good())
        {
            throw ExceptionFileOpen(outputVolumeFileName);
        }

        newFile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));
        newFile.close();
    }

    file.close();
}

//��һ�����������͵���������д�뵽ָ���ļ�����Ϊ������ݵ�һ����Ƭ
void VolumeIO::writeVolumeSliceInFloat(std::string outputVolumeFileName, std::vector<float> data, size_t sliceSize, int mode)
{
    //����Ҫд���ļ��ĸ��������ݵ���������Щ���ݴ���������ݵ�һ����Ƭ�� 
    std::ofstream outfile;
    outfile.open(generateTempName(outputVolumeFileName.c_str(), mode), fstream::out | fstream::binary | fstream::app);

    if (!outfile.good())
    {
        throw new ExceptionFileOpen(outputVolumeFileName);
    }

    outfile.write((const char*)&data[0], sliceSize * sizeof(float));
    //outfile.write()������data�����е��������ֽ�������ʽд�뵽�ļ��С�
    // ����ͨ��ȡdata�������׵�ַ��ת��Ϊconst char*���ͣ�Ȼ�����ÿ����������ռ���ֽ�����sizeof(float)���������ܵ���д����ֽ�����
    outfile.close();
}

void VolumeIO::writeVolumeSliceInFloat(std::string outputVolumeFileName, float& data, size_t sliceSize, int mode)
{
    std::ofstream outfile;
    outfile.open(generateTempName(outputVolumeFileName.c_str(), mode), fstream::out | fstream::binary | fstream::app);

    if (!outfile.good())
    {
        throw new ExceptionFileOpen(outputVolumeFileName);
    }

    outfile.write((const char*)(&data), sliceSize * sizeof(float));
    outfile.close();
}


void VolumeIO::computeMinMaxMean(std::vector<float>& data, size_t dataSize, float& min, float& max, float& mean)
{
    min = FLT_MAX;
    max = -FLT_MAX;
    mean = 0.0f;

    for (size_t i = 0; i < dataSize; i++)
    {
        min = std::min(min, data[i]);
        max = std::max(max, data[i]);
        mean += data[i];
    }

    mean /= dataSize;
}

void VolumeIO::write(std::vector<float>& data, novaCTF::Vec3ui resolution, std::string outputVolumeFileName, int mode, novaCTF::VolumeRotation rotation)
{
    std::ofstream outfile;

    try
    {
        outfile.open(outputVolumeFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

        if (!outfile.good())
        {
            throw new ExceptionFileOpen(outputVolumeFileName);
        }

        MRCHeader header;
        memset(&header, '\0', sizeof(MRCHeader));

        header.mode = mode;

        if (rotation == novaCTF::VolumeRotation::ALONG_XZ)
        {
            header.nx = resolution.x;
            header.ny = resolution.z;
            header.nz = resolution.y;

            header.mx = resolution.x;
            header.my = resolution.z;
            header.mz = resolution.y;

            header.cellDimX = (float)resolution.x;
            header.cellDimY = (float)resolution.z;
            header.cellDimZ = (float)resolution.y;
        }
        else //if (rotation==ALONG_XY)
        {
            header.nx = resolution.x;
            header.ny = resolution.y;
            header.nz = resolution.z;

            header.mx = resolution.x;
            header.my = resolution.y;
            header.mz = resolution.z;

            header.cellDimX = (float)resolution.x;
            header.cellDimY = (float)resolution.y;
            header.cellDimZ = (float)resolution.z;
        }

        header.map[0] = 'M';
        header.map[1] = 'A';
        header.map[2] = 'P';
        header.map[3] = ' ';

        header.mapC = 1;
        header.mapR = 2;
        header.mapS = 3;
        header.extra = 0;
        header.nint = 0;
        header.nreal = 0;

        float min, max, mean;

        computeMinMaxMean(data, (size_t)(resolution.x * resolution.y * resolution.z), min, max, mean);

        switch (header.mode)
        {
        case 0: header.dMin = 0.f;
            header.dMax = 255.f;
            header.dMean = novaCTF::convertRange(mean, min, max, header.dMin, header.dMax);
            break;
        case 1: header.dMin = -32768.f;
            header.dMax = 32767.f;
            header.dMean = novaCTF::convertRange(mean, min, max, header.dMin, header.dMax);
            break;
        case 2: header.dMin = min;
            header.dMax = max;
            header.dMean = mean;
            break;
        case 6: header.dMin = 0.f;
            header.dMax = 65535.f;
            header.dMean = novaCTF::convertRange(mean, min, max, header.dMin, header.dMax);
            break;
        default:	throw new ExceptionFileFormat();
            break;
        }

        outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));

        switch (mode)
        {
        case 0: writeData<unsigned char>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        case 1: writeData<signed short>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        case 2: writeData<float>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        case 6: writeData<unsigned short>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        default:
            throw new ExceptionFileFormat();
            break;
        }

        outfile.close();

    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
    catch (ExceptionFileFormat& e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        cout << "Error while writing volume to disk!" << endl;
    }
}


template <typename voxelType>
void VolumeIO::writeData(std::vector<float>& data, novaCTF::Vec3ui resolution, std::ofstream& outfile, novaCTF::VolumeRotation rotation, float oldMin, float oldMax, float newMin, float newMax)
{
    cout << "Writing the data..." << endl;
    if (rotation == novaCTF::VolumeRotation::ALONG_XZ)
        writeDataAlongXZPlane<voxelType>(data, resolution, outfile, oldMin, oldMax, newMin, newMax);
    else
        writeDataAlongXYPlane<voxelType>(data, resolution, outfile, oldMin, oldMax, newMin, newMax);
}

// slices along xz plane - classic eTomo reconstruction output
template <typename voxelType>
void VolumeIO::writeDataAlongXZPlane(std::vector<float>& data, novaCTF::Vec3ui resolution, std::ofstream& outfile, float oldMin, float oldMax, float newMin, float newMax)
{

    voxelType* dataToWrite = new voxelType[resolution.x * resolution.z];

    for (unsigned int y = 0; y < resolution.y; y++)
    {
        size_t bufferIndex = 0;

        for (size_t z = 0; z < resolution.z; z++)
        {
            for (size_t x = 0; x < resolution.x; x++)
            {
                dataToWrite[bufferIndex] = static_cast<voxelType>(novaCTF::convertRange(data[x + y * resolution.x + z * resolution.x * resolution.y], oldMin, oldMax, newMin, newMax));
                ++bufferIndex;
            }
        }

        outfile.write((const char*)dataToWrite, (size_t)resolution.x * (size_t)resolution.z * sizeof(voxelType));
    }

    delete[] dataToWrite;
}

template <typename voxelType>
void VolumeIO::writeDataAlongXYPlane(std::vector<float>& data, novaCTF::Vec3ui resolution, std::ofstream& outfile, float oldMin, float oldMax, float newMin, float newMax)
{

    voxelType* dataToWrite = new voxelType[resolution.x * resolution.y];

    for (size_t z = 0; z < resolution.z; z++)
    {
        size_t bufferIndex = 0;

        for (size_t y = 0; y < resolution.y; y++)
        {
            for (size_t x = 0; x < resolution.x; x++)
            {
                dataToWrite[bufferIndex] = static_cast<voxelType>(novaCTF::convertRange(data[x + y * resolution.x + z * resolution.x * resolution.y], oldMin, oldMax, newMin, newMax));
                ++bufferIndex;
            }
        }

        outfile.write((const char*)dataToWrite, (size_t)resolution.x * (size_t)resolution.y * sizeof(voxelType));
    }

    delete[] dataToWrite;
}

std::ostream& operator<<(std::ostream& oss, int mode)
{
    switch (mode)
    {
    case 2:
        oss << "gray scale (32 bit float)";
        break;
    case 0:
        oss << "gray scale (8 bit unsigned)";
        break;
    case 6:
        oss << "gray scale (16 bit unsigned)";
        break;
    case 1:
        oss << "gray scale (16 bit signed)";
        break;
    }
    return oss;
}



void VolumeIO::mywriteMRC( MRCHeader& header,  std::vector<float>& volume_data,  std::string& output_name,  char* extra_data = nullptr)
{
    // ȷ�����ݳߴ���� MRCHeader �е� nx, ny, nz
    assert(volume_data.size() == header.nx * header.ny * header.nz);

    std::ofstream mrc_file(output_name, std::ios::binary | std::ios::out | std::ios::trunc);

    if (!mrc_file.is_open())
    {
        throw std::runtime_error("Failed to open file for writing: " + output_name);
    }

    // ���� header �ı�Ҫ��Ϣ������ mode ��Ϊ 2 ��ʾ 32 λ������
    header.mode = 2;
    // ���㲢�������ݵ���Сֵ�����ֵ��ƽ��ֵ
    // �ڴ˼����Ѿ����� computeMinMaxMean ����
    computeMinMaxMean(volume_data, volume_data.size(), header.dMin, header.dMax, header.dMean);;

    // д�� MRC ͷ��
    mrc_file.write(reinterpret_cast<const char*>(&header), sizeof(MRCHeader));

    // д��������ݣ�������ڣ�
    if (header.extra != 0 && extra_data != nullptr)
    {
        mrc_file.write(extra_data, header.extra);
    }

    // ���������д���ļ�
    mrc_file.write(reinterpret_cast<const char*>(volume_data.data()), sizeof(float) * volume_data.size());

    // ����Ƿ�д��ɹ�
    if (!mrc_file)
    {
        throw std::runtime_error("Error occurred while writing volume data to disk.");
    }

    mrc_file.close();
}