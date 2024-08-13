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
    std::vector<float> linearData;//输入MRCHeader 结构体 header，二维 float 类型 vector 数据 data，输出文件名字符串 outputName，以及额外数据指针 extraData。
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
    //将二维 vector 转换为一维 vector linearData，然后调用另一个 writeMRCStack 函数处理一维数据。
}

void VolumeIO::writeMRCStack(MRCHeader header, std::vector<float>& data, string outputName, char* extraData)
{
    //输入参数：MRCHeader 结构体 header，一维 float 类型 vector 数据 data，输出文件名字符串 outputName，以及额外数据指针 extraData。
    std::ofstream outfile; 

    try
    {
        outfile.open(outputName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
        //打开一个输出文件流 outfile，尝试以二进制和截断模式打开指定的输出文件 outputName。

        if (!outfile.good())
        {
            throw new ExceptionFileOpen(outputName);//果文件打开失败，抛出 ExceptionFileOpen 异常，打印错误信息。
        }

        header.mode = 2;//设置 MRCHeader 中的 mode 字段为 2，表示数据是以 32 位浮点数形式存储。
        computeMinMaxMean(data, data.size(), header.dMin, header.dMax, header.dMean);//使用 computeMinMaxMean 函数计算输入数据的最小值、最大值和平均值，并更新到 header 结构体中。

        outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));//将 header 结构体写入到已打开的输出文件中。

        //if (header.extra != 0)//若 header 中的额外数据长度（extra 字段）不为零，则将 extraData 指向的内容写入到文件中。
        //{
        //    outfile.write(reinterpret_cast<char*>(extraData), header.extra);
       // }

        writeProjections(outfile, data, header.nx * header.ny, header.nz);//调用 writeProjections 函数将数据按照 header 中的 nx、ny 和 nz 定义的维度顺序写入到文件中。

        outfile.close();//关闭输出文件流。

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
    //整个过程包含异常处理，当出现文件打开失败或其他文件格式异常时，会捕获对应的异常并打印错误信息。
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


//读取输入栈文件的MRC头信息，依据旋转模式和新的体积尺寸对其进行调整，最后将调整后的新头信息写入到一个新的临时输出体积文件中。
void VolumeIO::writeHeader(std::string inputStackFileName, std::string outputVolumeFileName, novaCTF::Vec3ui volumeDimensions, int mode, novaCTF::VolumeRotation rotation)
{
    //输入栈文件名，类型为字符串。 outputVolumeFileName: 输出体积文件名模板，类型为字符串 volumeDimensions: 期望的体积尺寸，类型为novaCTF::Vec3ui，表示(x, y, z)维度。
    //mode: 这里是2 
    // rotation: 旋转模式，属于novaCTF::VolumeRotation枚举类型，这里指定了两种情况：沿着XZ轴或XY轴旋转。
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

    outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));       //这里总是改不对啊啊啊啊
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

//负责根据传入的参数对体积数据进行值域转换，并将转换后的数据保存到指定的输出文件中。函数主要分为以下几个步骤：
//这里的问题是传入的值都为0前边输出就有问题
void VolumeIO::convertVolumeValues(std::string outputVolumeFileName, float min, float max, float mean, int mode)
{
    float newMin;
    float newMax;
    //mode为2时，保持原最小值和最大值不变； 
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
    //读取输出文件的MRC头信息（MRCHeader结构体），并调用correctHeader函数对MRC头进行修正，使其与新的值域范围匹配。
    size_t sliceSize = header.nx * header.ny;
    //计算体积数据每一切片的大小（sliceSize）。
    //根据mode选择相应的模板参数调用convertSliceData函数，
    // 该函数负责逐个切片地读取原始数据并将其转换至新的值域范围内，然后写入到输出文件中。对于不同mode值，使用的数据类型不同，分别为
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

//将一个浮点数类型的数据向量写入到指定文件中作为体积数据的一个切片
void VolumeIO::writeVolumeSliceInFloat(std::string outputVolumeFileName, std::vector<float> data, size_t sliceSize, int mode)
{
    //包含要写入文件的浮点数数据的向量，这些数据代表体积数据的一个切片。 
    std::ofstream outfile;
    outfile.open(generateTempName(outputVolumeFileName.c_str(), mode), fstream::out | fstream::binary | fstream::app);

    if (!outfile.good())
    {
        throw new ExceptionFileOpen(outputVolumeFileName);
    }

    outfile.write((const char*)&data[0], sliceSize * sizeof(float));
    //outfile.write()方法将data向量中的数据以字节流的形式写入到文件中。
    // 这里通过取data向量的首地址并转换为const char*类型，然后乘以每个浮点数所占的字节数（sizeof(float)）来计算总的需写入的字节数。
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
    // 确保数据尺寸符合 MRCHeader 中的 nx, ny, nz
    assert(volume_data.size() == header.nx * header.ny * header.nz);

    std::ofstream mrc_file(output_name, std::ios::binary | std::ios::out | std::ios::trunc);

    if (!mrc_file.is_open())
    {
        throw std::runtime_error("Failed to open file for writing: " + output_name);
    }

    // 更新 header 的必要信息，例如 mode 设为 2 表示 32 位浮点数
    header.mode = 2;
    // 计算并更新数据的最小值、最大值和平均值
    // 在此假设已经有了 computeMinMaxMean 函数
    computeMinMaxMean(volume_data, volume_data.size(), header.dMin, header.dMax, header.dMean);;

    // 写入 MRC 头部
    mrc_file.write(reinterpret_cast<const char*>(&header), sizeof(MRCHeader));

    // 写入额外数据（如果存在）
    if (header.extra != 0 && extra_data != nullptr)
    {
        mrc_file.write(extra_data, header.extra);
    }

    // 将体积数据写入文件
    mrc_file.write(reinterpret_cast<const char*>(volume_data.data()), sizeof(float) * volume_data.size());

    // 检查是否写入成功
    if (!mrc_file)
    {
        throw std::runtime_error("Error occurred while writing volume data to disk.");
    }

    mrc_file.close();
}