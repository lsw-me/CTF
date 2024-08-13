#include <stdio.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <math.h>
#include <cstring>
#include <algorithm>

#include "exception.h"
#include "mrcStack.h"

using namespace std;


MRCStack::MRCStack(string aFileName, bool getRange, bool swapYZ, bool keepOpen)
{
	mKeepOpen = keepOpen; //初始化类成员变量 mKeepOpen 为传入参数 keepOpen 的值。这将决定文件是否在对象生命周期内保持打开状态。

	mStackName = aFileName;         //将构造函数参数 aFileName 赋值给成员变量 mStackName，表示要处理的 MRC 格式文件名

	mStackFile.open(mStackName.c_str(), std::ios::binary);//
	//使用 std::ios::binary 模式打开指定的 MRC 文件，并将结果赋值给 mStackFile 成员变量。如果打开失败，抛出一个 ExceptionFileOpen 异常，其中包含文件名。
	if (!mStackFile.good())
	{
		throw ExceptionFileOpen(mStackName);
	}

	mStackFile.read((char*)&mHeaderMRC, sizeof(MRCHeader)); //从文件中读取 MRC 头部信息（类型为 MRCHeader），并将数据存放在 mHeaderMRC 成员变量中。
	mStackFile.ignore(mHeaderMRC.extra);//忽略 MRC 文件头中的额外字节数量（由 mHeaderMRC.extra 指定）。

	if (swapYZ) //如果为真对 MRC 头部的部分字段进行 YZ 维度交换操作，以适应不同的数据排列方式
	{
		novaCTF::swap<int>(mHeaderMRC.ny, mHeaderMRC.nz);
		novaCTF::swap<int>(mHeaderMRC.my, mHeaderMRC.mz);
		novaCTF::swap<int>(mHeaderMRC.nyStart, mHeaderMRC.nzStart);
		novaCTF::swap<float>(mHeaderMRC.cellDimY, mHeaderMRC.cellDimZ);
		novaCTF::swap<float>(mHeaderMRC.cellAngleY, mHeaderMRC.cellAngleZ);
		novaCTF::swap<int>(mHeaderMRC.mapR, mHeaderMRC.mapS);
	}

	mInputMeanValue = mHeaderMRC.dMean;//初始化一些与 MRC 数据相关的成员变量，如平均值、投影数量、投影大小等。
	mNumberOfProjections = mHeaderMRC.nz;
	mProjectionSize = (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.ny;


	//writeOutHeader();

	if (getRange)  //若 getRange 参数为真，则根据 MRC 文件模式获取数据范围。调用 getDataRange 函数并使用模板参数匹配相应的数据类型。
	{
		switch (mHeaderMRC.mode)
		{
		case 0: getDataRange<unsigned char>(mStackFile);
			break;
		case 1: getDataRange<signed short>(mStackFile);
			break;
		case 2: getDataRange<float>(mStackFile);
			break;
		case 6: getDataRange<unsigned short>(mStackFile);
			break;
		default: throw ExceptionFileFormat();
			break;
		}
	}

	//header initialization - holds only for parallel geometry !!!
	mHeader = new GeomHeader();  //创建一个新的 GeomHeader 类型的对象 mHeader，并初始化其成员变量，这些变量主要与 MRC 数据的几何属性相关。
	mHeader->mWidth = mHeaderMRC.nx;
	mHeader->mHeight = mHeaderMRC.ny;
	mHeader->mDepth = mHeaderMRC.nz;

	mHeader->mSourcePosition = novaCTF::makeVec3f(mHeaderMRC.nx, mHeaderMRC.ny, mHeaderMRC.nx);
	mHeader->mDetectorPosition = novaCTF::makeVec3f(mHeaderMRC.nx, mHeaderMRC.ny, mHeaderMRC.nx);
	mHeader->mHorizontalPitch = novaCTF::makeVec3f(1.f, 1.f, 1.f);
	mHeader->mVerticalPitch = novaCTF::makeVec3f(1.f, 1.f, 1.f);  //这里不是很清楚


	switch (mHeaderMRC.mode)   //根据 MRC 文件模式设置 mSizeOfVoxelType 成员变量，表示每个体素的数据类型大小。
	{
	case 0: mSizeOfVoxelType = sizeof(unsigned char);
		break;
	case 1: mSizeOfVoxelType = sizeof(signed short);
		break;
	case 2: mSizeOfVoxelType = sizeof(float);
		break;
	case 6: mSizeOfVoxelType = sizeof(unsigned short);
		break;
	default: throw ExceptionFileFormat();
		break;
	}

	if (mHeaderMRC.extra != 0)//若 MRC 文件头中有额外数据，则分配内存空间存放这部分数据，并从文件中读取到 extraData 变量中；否则分配一个占位符字符。
	{
		extraData = new char[mHeaderMRC.extra];
		mStackFile.seekg(sizeof(MRCHeader));
		mStackFile.read(extraData, mHeaderMRC.extra);
	}
	else
	{
		extraData = new char[1];	//dummy
	}

	if (mKeepOpen)// 最后，根据 mKeepOpen 的值来决定是否关闭文件或移动文件指针至数据区开始位置。
	{
		mStackFile.seekg(sizeof(MRCHeader) + mHeaderMRC.extra);
	}
	else
		mStackFile.close();
}

MRCStack::~MRCStack()
{
	if (mKeepOpen)
		mStackFile.close();

	delete[] extraData;
}

void MRCStack::writeOutHeader()
{
	cout << "Number of columns, rows, sections" << mHeaderMRC.nx << ", " << mHeaderMRC.ny << ", " << mHeaderMRC.nz << endl;
	cout << "Map mode" << mHeaderMRC.mode << endl;
	cout << "Start columns, rows, sections" << mHeaderMRC.nxStart << ", " << mHeaderMRC.nyStart << ", " << mHeaderMRC.nzStart << endl;
	cout << "Grid size in x, y, z" << mHeaderMRC.mx << ", " << mHeaderMRC.my << ", " << mHeaderMRC.mz << endl;
	cout << "Cell dimensions in x, y, z" << mHeaderMRC.cellDimX << ", " << mHeaderMRC.cellDimY << ", " << mHeaderMRC.cellDimZ << endl;
	cout << "Cell angles in x, y, z (degrees)" << mHeaderMRC.cellAngleX << ", " << mHeaderMRC.cellAngleY << ", " << mHeaderMRC.cellAngleZ << endl;
	cout << "Axis corresponding to columns, rows and sections (1,2,3 for X,Y,Z)" << mHeaderMRC.mapC << ", " << mHeaderMRC.mapR << ", " << mHeaderMRC.mapS << endl;
	cout << "Origin on x, y, z" << mHeaderMRC.originX << ", " << mHeaderMRC.originY << ", " << mHeaderMRC.originZ << endl;
	cout << "Minimum density" << mHeaderMRC.dMin << endl;
	cout << "Maximum density" << mHeaderMRC.dMax << endl;
	cout << "Mean density" << mHeaderMRC.dMean << endl;
	cout << "Tilt angles - original" << mHeaderMRC.tiltangles[0] << ", " << mHeaderMRC.tiltangles[1] << ", " << mHeaderMRC.tiltangles[2] << endl;
	cout << "Tilt angles - current" << mHeaderMRC.tiltangles[3] << ", " << mHeaderMRC.tiltangles[4] << ", " << mHeaderMRC.tiltangles[5] << endl;
	cout << "Space group" << mHeaderMRC.ISPG << endl;
	cout << "Number of bytes used for symmetry data" << mHeaderMRC.NSYMBT << endl;
	cout << "Number of bytes in extended header" << mHeaderMRC.extra << endl;
	cout << "Creator ID" << mHeaderMRC.creatorID << endl;
	cout << "ID type" << mHeaderMRC.idtype << endl;
	cout << "Lens" << mHeaderMRC.lens << endl;
	cout << mHeaderMRC.nLabel << "labels:" << endl;
	for (int i = 0; i < mHeaderMRC.nLabel; i++)
		cout << mHeaderMRC.labels[i] << endl;
}

void MRCStack::readProjections(float* aData, unsigned int numberOfProjections, unsigned int sliceNumber)
{
	prepareFilePosition(sliceNumber * mProjectionSize * mSizeOfVoxelType);//
	//prepareFilePosition(sliceNumber * mProjectionSize * mSizeOfVoxelType);
	//mProjectionSize = (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.ny; 也就是图片大小
	// mSizeOfVoxelType 数据大小
	//该函数根据给定的切片编号sliceNumber计算并设置文件读取的当前位置，以便读取相应的投影数据。
	switch (mHeaderMRC.mode)
	{
	case 0: readSlices<unsigned char>(aData, mProjectionSize * (size_t)numberOfProjections);
		break;
	case 1: readSlices<signed short>(aData, mProjectionSize * (size_t)numberOfProjections);
		break;
	case 2: readSlices<float>(aData, mProjectionSize * (size_t)numberOfProjections);
		break;
	case 6: readSlices<unsigned short>(aData, mProjectionSize * (size_t)numberOfProjections);
		break;
	default: throw ExceptionFileFormat();
		break;
	}


	if (!mKeepOpen)
		mStackFile.close();
}

void MRCStack::readAllProjections(float* aData)
{
	prepareFilePosition();

	switch (mHeaderMRC.mode)
	{
	case 0: readSlices<unsigned char>(aData, mProjectionSize * mNumberOfProjections);
		break;
	case 1: readSlices<signed short>(aData, mProjectionSize * mNumberOfProjections);
		break;
	case 2: readSlices<float>(aData, mProjectionSize * mNumberOfProjections);
		break;
	case 6: readSlices<unsigned short>(aData, mProjectionSize * mNumberOfProjections);
		break;
	default: throw ExceptionFileFormat();
		break;
	}

	if (!mKeepOpen)
		mStackFile.close();
}

void MRCStack::prepareFilePosition(size_t offset)
{
	if (mKeepOpen)
	{
		if (!mStackFile.good())
		{
			throw ExceptionFileOpen(mStackName);
		}
	}
	else
	{
		mStackFile.open(mStackName.c_str(), std::ios::binary);

		if (!mStackFile.good())
		{
			throw ExceptionFileOpen(mStackName);
		}

		mStackFile.seekg(sizeof(MRCHeader) + mHeaderMRC.extra + offset);
		//函数移动文件流的读取指针到适当的位置：
		// 首先跳过MRC文件头的大小（sizeof(MRCHeader)）和额外的数据（mHeaderMRC.extra），然后移动到给定的偏移量offset处。

	}
}

/*template <typename voxelType>
float MRCStack::convertValue(voxelType rawValue, float minValue)
{
	float rawFloat = (float)rawValue;
	if(mLogarithmizeData)
	{
		return logf(mDataMax / (rawFloat - minValue));
	}
	else
	{
		return rawFloat;
	}
}*/

template <typename voxelType>
void MRCStack::readSlices(float* aData, size_t elementsToRead)//void MRCStack::readSlices(float* aData, size_t elementsToRead)
{
	//一个新的voxelType类型的动态数组data，大小为elementsToRead，用于存储从文件中读取的原始体素数据
	voxelType* data = new voxelType[elementsToRead];

	mStackFile.read((char*)(data), elementsToRead * sizeof(voxelType));

	//mStackFile.read方法从MRC文件中读取elementsToRead个体素数据到刚分配的data数组中
	

	//41次
	//cout << "进入循环" << endl;
	for (size_t i = 0; i < (size_t)elementsToRead; i++)
	{
		aData[i] = (float)data[i];
		//cout << i <<endl;
		//cout << aData[i] << endl;
	}
	delete[] data;
}

template <typename voxelType>
void MRCStack::getDataRange(ifstream& infile)
{
	mDataMax = -FLT_MAX;
	mDataMin = FLT_MAX;
	mDataMean = 0.0;

	for (size_t i = 0; i < (size_t)mHeaderMRC.ny * (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.nz; i++)
	{
		voxelType value;
		infile.read(reinterpret_cast<char*>(&value), sizeof(voxelType));
		mDataMax = std::max((float)value, mDataMax);
		mDataMin = std::min((float)value, mDataMin);
		mDataMean += (float)value;
	}

	mDataMean = mDataMean / ((float)mHeaderMRC.ny * (float)mHeaderMRC.nx * (float)mHeaderMRC.nz);
}

novaCTF::Vec3ui MRCStack::getResolution()
{
	return novaCTF::Vec3ui(mHeaderMRC.nx, mHeaderMRC.ny, mHeaderMRC.nz);
}

GeomHeader* MRCStack::getHeader()
{
	return mHeader;
}

float MRCStack::getStackMin()
{
	return mDataMin;
}

float MRCStack::getStackMax()
{
	return mDataMax;
}

float MRCStack::getInputMeanValue()
{
	return mInputMeanValue;
}

float MRCStack::getStackMean()
{
	return mDataMean;
}

unsigned int MRCStack::getNumberOfProjections()
{
	return mHeaderMRC.nz;
}

MRCHeader MRCStack::getStackHeader()
{
	return mHeaderMRC;
}

char* MRCStack::getExtraData()
{
	return extraData;
}

size_t MRCStack::getProjectionSize()
{
	return mProjectionSize;
}
