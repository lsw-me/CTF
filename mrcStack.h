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

	void readProjections(float* aData, unsigned int numberOfProjections);  //�Ӷ�ջ�ж�ȡָ��������ͶӰ���ݵ��ṩ�ĸ����ͻ�������
	void readProjections(float* aData, unsigned int numberOfProjections, unsigned int aIndex);//��ָ��������ʼ��ȡͶӰ���ݡ�
	void readAllProjections(float* aData);//��ȡ��ջ�е�����ͶӰ���ݡ�

	void writeOutHeader();//�����ջͷ��Ϣ
	float getStackMin();
	float getStackMax();
	float getStackMean();//��ȡ��ջ�����ݵ���Сֵ�����ֵ��ƽ��ֵ
	float getInputMeanValue();//��ȡ�������ݵ�ƽ��ֵ��
	unsigned int getNumberOfProjections();//��ȡ��ջ�е�ͶӰ������
	size_t getProjectionSize();//��ȡ����ͶӰ�����ݴ�С�����ֽ�Ϊ��λ����
	novaCTF::Vec3ui getResolution();//���طֱ�����Ϣ������ΪVec3ui��

	GeomHeader* getHeader();//����ȡ����ͷ��Ϣָ�롣
	MRCHeader getStackHeader();//��ȡMRCͷ��Ϣ
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
