#pragma once
#pragma once
#include "parameterSetup.h"
#include "mrcStack.h"
#include "projectionSet.h"
#include "microscopeGeometry.h"
#include "common.h"

using namespace novaCTF;

class CTF3d
{
public:

	struct voxelProjection {
		int startingIndex;
		int endingIndex;
		double position;
		int projectionImpact;
		bool fillWithProjectionValue;
	};

	CTF3d(string inputstackname, novaCTF::Vec3ui volumeDimensions);
	~CTF3d();

	void run();
private:

	void writeVolumeSlice();
	void backproject();
	void initVariables();
	void initAngles();

	// void setNeededSlices();
	void computeOneRow(int& volumeSliceIndex, int projectionIndex, int numberOfVoxels, double xPosition, float cosAngle, bool computeSliceStatistics, unsigned int defocusProjIndex, size_t volumeSliceOffset);
	void projectSlice();
	void computeProjectionBoundaries();

	void prepareCTF();
	void loadInputStacks();

	void setDataSize();
	void writeOutDefocusSlices(vector<vector<unsigned int> >& sliceDefocusSplits);

	string generateFilename(string originalFilename, unsigned int number);





	MRCStack* initialStack;

	std::vector<MRCStack*> inputStacks;
	ParameterSetup params;
	ProjectionSet* projSet;
	Geometry* microscopeGeometry;
	novaCTF::Vec3ui setvolumeDimensions;

	vector<vector<unsigned int> > sliceDefocusSplits;
	unsigned int numberOfStacks;
	string firstFileName;

	float globalMin;
	float globalMax;
	float globalMean;

	size_t volumeSliceSize;

	std::vector<std::vector<voxelProjection>> projectionBoundaries;

	unsigned int  slicesToProcessAtOnce;	//number of slices to process at once
	size_t processedDataSize;
	Vec3i projRes;
	unsigned int nviews;    	// number of input views, or number used  ������ͼ����ʹ�õ�����
	unsigned int sliceX;     	// width of output slice //����Ŀ��
	unsigned int sliceZ;    	// thickness of output slice made by PROJECT PROJECT�����������Ƭ��� ��ϣ����200������ֵ

	vector<float> angles;		// tilt angles
	float edgeFill;				// filling value for FS
	vector<float> cbet;			// cos angles in radians;
	vector<float> sbet;			// sin angles in radians;
	float scale;
	float addition;

	vector<vector<float> > currentRows;		// �洢�뵱ǰ��Ƭ��Ӧ��������ͼ�е�һ�л����
	vector<float> currentVolumeSlice;		// ���� xz ��Ƭ������Ƭ��;
	vector<int>	  listOfIndices;			// �����Ԥ��ָ���б�

	unsigned int reportIndex;



};
