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
	mKeepOpen = keepOpen; //��ʼ�����Ա���� mKeepOpen Ϊ������� keepOpen ��ֵ���⽫�����ļ��Ƿ��ڶ������������ڱ��ִ�״̬��

	mStackName = aFileName;         //�����캯������ aFileName ��ֵ����Ա���� mStackName����ʾҪ����� MRC ��ʽ�ļ���

	mStackFile.open(mStackName.c_str(), std::ios::binary);//
	//ʹ�� std::ios::binary ģʽ��ָ���� MRC �ļ������������ֵ�� mStackFile ��Ա�����������ʧ�ܣ��׳�һ�� ExceptionFileOpen �쳣�����а����ļ�����
	if (!mStackFile.good())
	{
		throw ExceptionFileOpen(mStackName);
	}

	mStackFile.read((char*)&mHeaderMRC, sizeof(MRCHeader)); //���ļ��ж�ȡ MRC ͷ����Ϣ������Ϊ MRCHeader�����������ݴ���� mHeaderMRC ��Ա�����С�
	mStackFile.ignore(mHeaderMRC.extra);//���� MRC �ļ�ͷ�еĶ����ֽ��������� mHeaderMRC.extra ָ������

	if (swapYZ) //���Ϊ��� MRC ͷ���Ĳ����ֶν��� YZ ά�Ƚ�������������Ӧ��ͬ���������з�ʽ
	{
		novaCTF::swap<int>(mHeaderMRC.ny, mHeaderMRC.nz);
		novaCTF::swap<int>(mHeaderMRC.my, mHeaderMRC.mz);
		novaCTF::swap<int>(mHeaderMRC.nyStart, mHeaderMRC.nzStart);
		novaCTF::swap<float>(mHeaderMRC.cellDimY, mHeaderMRC.cellDimZ);
		novaCTF::swap<float>(mHeaderMRC.cellAngleY, mHeaderMRC.cellAngleZ);
		novaCTF::swap<int>(mHeaderMRC.mapR, mHeaderMRC.mapS);
	}

	mInputMeanValue = mHeaderMRC.dMean;//��ʼ��һЩ�� MRC ������صĳ�Ա��������ƽ��ֵ��ͶӰ������ͶӰ��С�ȡ�
	mNumberOfProjections = mHeaderMRC.nz;
	mProjectionSize = (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.ny;


	//writeOutHeader();

	if (getRange)  //�� getRange ����Ϊ�棬����� MRC �ļ�ģʽ��ȡ���ݷ�Χ������ getDataRange ������ʹ��ģ�����ƥ����Ӧ���������͡�
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
	mHeader = new GeomHeader();  //����һ���µ� GeomHeader ���͵Ķ��� mHeader������ʼ�����Ա��������Щ������Ҫ�� MRC ���ݵļ���������ء�
	mHeader->mWidth = mHeaderMRC.nx;
	mHeader->mHeight = mHeaderMRC.ny;
	mHeader->mDepth = mHeaderMRC.nz;

	mHeader->mSourcePosition = novaCTF::makeVec3f(mHeaderMRC.nx, mHeaderMRC.ny, mHeaderMRC.nx);
	mHeader->mDetectorPosition = novaCTF::makeVec3f(mHeaderMRC.nx, mHeaderMRC.ny, mHeaderMRC.nx);
	mHeader->mHorizontalPitch = novaCTF::makeVec3f(1.f, 1.f, 1.f);
	mHeader->mVerticalPitch = novaCTF::makeVec3f(1.f, 1.f, 1.f);  //���ﲻ�Ǻ����


	switch (mHeaderMRC.mode)   //���� MRC �ļ�ģʽ���� mSizeOfVoxelType ��Ա��������ʾÿ�����ص��������ʹ�С��
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

	if (mHeaderMRC.extra != 0)//�� MRC �ļ�ͷ���ж������ݣ�������ڴ�ռ����ⲿ�����ݣ������ļ��ж�ȡ�� extraData �����У��������һ��ռλ���ַ���
	{
		extraData = new char[mHeaderMRC.extra];
		mStackFile.seekg(sizeof(MRCHeader));
		mStackFile.read(extraData, mHeaderMRC.extra);
	}
	else
	{
		extraData = new char[1];	//dummy
	}

	if (mKeepOpen)// ��󣬸��� mKeepOpen ��ֵ�������Ƿ�ر��ļ����ƶ��ļ�ָ������������ʼλ�á�
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
	//mProjectionSize = (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.ny; Ҳ����ͼƬ��С
	// mSizeOfVoxelType ���ݴ�С
	//�ú������ݸ�������Ƭ���sliceNumber���㲢�����ļ���ȡ�ĵ�ǰλ�ã��Ա��ȡ��Ӧ��ͶӰ���ݡ�
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
		//�����ƶ��ļ����Ķ�ȡָ�뵽�ʵ���λ�ã�
		// ��������MRC�ļ�ͷ�Ĵ�С��sizeof(MRCHeader)���Ͷ�������ݣ�mHeaderMRC.extra����Ȼ���ƶ���������ƫ����offset����

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
	//һ���µ�voxelType���͵Ķ�̬����data����СΪelementsToRead�����ڴ洢���ļ��ж�ȡ��ԭʼ��������
	voxelType* data = new voxelType[elementsToRead];

	mStackFile.read((char*)(data), elementsToRead * sizeof(voxelType));

	//mStackFile.read������MRC�ļ��ж�ȡelementsToRead���������ݵ��շ����data������
	

	//41��
	//cout << "����ѭ��" << endl;
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
