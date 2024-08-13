#pragma once

#include <string>
#include "common.h"
#include "mrcHeader.h"
#include "geomHeader.h"
#include <fstream>
#include <vector>
#include <typeinfo>

using namespace std;

class MRCStack
{
public:

	MRCStack(string aFileName, bool getRange, bool swapYZ, bool keepOpen);

	~MRCStack();

	void readProjections(float* aData, unsigned int numberOfProjections);  //从堆栈中读取指定数量的投影数据到提供的浮点型缓冲区。
	void readProjections(float* aData, unsigned int numberOfProjections, unsigned int aIndex);//从指定索引开始读取投影数据。
	void readAllProjections(float* aData);//读取堆栈中的所有投影数据。

	void writeOutHeader();//输出堆栈头信息
	float getStackMin();
	float getStackMax();
	float getStackMean();//获取堆栈中数据的最小值、最大值和平均值
	float getInputMeanValue();//获取输入数据的平均值。
	unsigned int getNumberOfProjections();//获取堆栈中的投影数量。
	size_t getProjectionSize();//获取单个投影的数据大小（以字节为单位）。
	novaCTF::Vec3ui getResolution();//返回分辨率信息，类型为Vec3ui。

	GeomHeader* getHeader();//：获取几何头信息指针。
	MRCHeader getStackHeader();//获取MRC头信息
	char* getExtraData();

	ifstream mStackFile;

private:

	MRCHeader mHeaderMRC;
	char* extraData;

	size_t mProjectionSize;
	string mStackName;

	template <typename voxelType>
	void readSlices(float* aData, size_t elementsToRead);

	void prepareFilePosition(size_t offset = 0);

	template <typename voxelType>
	void getDataRange(ifstream& infile);

	GeomHeader* mHeader;					//MRC header

	float mDataMin;				//min density value in the stack
	float mDataMax;				//max denisty value in the stack
	float mDataMean;
	float mInputMeanValue;
	unsigned int mNumberOfProjections;
	size_t mSizeOfVoxelType;

	bool mKeepOpen;


};
