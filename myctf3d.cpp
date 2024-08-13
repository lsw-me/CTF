#include "myctf3d.h"
#include "common.h"
#include "volumeIO.h"
#include <math.h>
#include <algorithm>
#include <float.h>
#include <sstream>
#include <fftw3.h>
#include "geometrySetup.h"
#include "geomHeader.h"
#include "myctf3d.h"

using namespace novaCTF;
using namespace std;



CTF3d::CTF3d(string inputname, novaCTF::Vec3ui volumeDimensions)
{
	//params = aParams;
	projSet = new ProjectionSetIdentity();
	setvolumeDimensions = volumeDimensions;
	firstFileName = inputname;
	initialStack = new MRCStack(firstFileName, true, true, false);
	string TiltAnglesFilename = "E:\\VS studio\\proj\\project2\\IS002_291013_005.tlt";//�Ƕ��ļ�
	projSet->init(initialStack->getNumberOfProjections());//get�������� mHeaderMRC.nz��
	microscopeGeometry = new Geometry(*initialStack, setvolumeDimensions, TiltAnglesFilename, 0.0f, 0,0.0f);//���ﴫ���volumeDimensions 
	cout << "done" << endl;
}

CTF3d::~CTF3d()
{
	for (unsigned int inputStackIndex = 0; inputStackIndex < numberOfStacks; inputStackIndex++)
		delete inputStacks[inputStackIndex];


	delete microscopeGeometry;
	delete initialStack;
	delete projSet;
}

void CTF3d::run()
{
	string outputname = "E:\\VS studio\\proj\\project2\\test_rec_output\\output.rec";
	VolumeIO::writeHeader(firstFileName, outputname, setvolumeDimensions, 2, novaCTF::VolumeRotation::ALONG_XZ);  // ��ת��xz
	cout << "writeHeader done" << endl;
	initVariables();
	cout << "init done " << endl;
	prepareCTF();
	cout << "prepareCTF done" << endl;
	loadInputStacks();
	cout << "loatStack done" << endl;
	computeProjectionBoundaries();
	cout << "cpmpute done" << endl;
	setDataSize();
	cout << "setdatasize done" << endl;
	backproject();                      
	cout << "wbp done" << endl;

	VolumeIO::convertVolumeValues(outputname, globalMin, globalMax, globalMean / (float)(200 * 200 * 200), 2);//
	cout << "volume done" << endl;
	//������ݴ���Ĳ�����������ݽ���ֵ��ת��������ת��������ݱ��浽ָ��������ļ��С�
}

void CTF3d::loadInputStacks()  //������Ϊnovactf���������ջ�ļ���������Щ�ļ�һ�����У����������ֻ��һ��У�����mrc�ļ���������б�����ǣ�������Ҫȥ�������Ż�
{
	//���������ص���MRC�ļ�����firstFileNameָ�����ļ�������������Ϊһ����ջʵ��ѹ��inputStacks�����С�ͬʱ�������ڴ˴�����ִ�в����ء�
	//if (!params.Use3DCTF())
	//{
		inputStacks.push_back(new MRCStack(firstFileName, false, false, true));
		return;
	//}

	//for (unsigned int inputStackIndex = 0; inputStackIndex < numberOfStacks; inputStackIndex++)
	//{
	//	inputStacks.push_back(new MRCStack(generateFilename(params.InputStackName(), inputStackIndex), false, false, params.KeepFilesOpen()));
	//}
}

void CTF3d::prepareCTF()
{

	sliceDefocusSplits.resize(initialStack->getNumberOfProjections());//����HeaderMRC.nzֵ������תͼ����nz��ֵΪ200

	size_t xzSliceSize = setvolumeDimensions.x*setvolumeDimensions.z;//���������xz����Ҳ��������y��


	//
	//��û�� CTF ������£�ֻ��ʹ��һ������ 0 ����Ƭ��ָ��Ψһʹ�õ������ջ
	//if (!params.Use3DCTF()) //�ж����Ҳ�����ӵ�
	// 
	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++)
	{
		unsigned int projIndex = it.second();
		sliceDefocusSplits[projIndex].resize(xzSliceSize);
		std::fill(sliceDefocusSplits[projIndex].begin(), sliceDefocusSplits[projIndex].end(), 0);
		//ѭ�������� sliceDefocusSplits��ά�����СΪ 200*40000
	}
		numberOfStacks = 1;
		currentRows.resize(1);
		return;

	//����ע�͵���ԭ�������������ջ����3dCTF Ȼ��������Ķ�ջ�Ǹ�����ά��Ϣ����У����Ľ��������ҿ���ֻ����һ��mrc���������������ֵ��������ά��ϢУ������
	//���Բ�ʹ������3dCTF
   /*
	if (params.DefocusStep() != 0.0f)
	{
		numberOfStacks = microscopeGeometry->computeNumberOfParts(params.DefocusStep(), params.PixelSize(), params.VolumeThicknessType());
		cout << "Number of input stacks based on DefocusStep, PixelSize and volume thickness (THICKNESS) is: " << numberOfStacks << endl;
	}
	else if (params.NumberOfInputStacks() != 0)
	{
		numberOfStacks = params.NumberOfInputStacks();
		cout << "Number of input stacks based on NumberOfInputStacks is: " << numberOfStacks << endl;
	}
	else
	{
		cout << "Either defocus step size (DefocusStep in nm) or number of input stacks (NumberOfInputStacks) has to be specified!" << endl;
		return;
	}

	currentRows.resize(numberOfStacks);

	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++)
	{
		unsigned int projIndex = it.second();
		microscopeGeometry->setProjectionGeometry(projIndex);
		sliceDefocusSplits[projIndex].reserve(xzSliceSize);
		microscopeGeometry->generateFocusGrid(sliceDefocusSplits[projIndex], numberOfStacks, params.VolumeThicknessType());
	}

	writeOutDefocusSlices(sliceDefocusSplits);
	*/

}


void CTF3d::initVariables()
{
	projRes.x = (int)initialStack->getResolution().x; 
	projRes.y = (int)initialStack->getResolution().y;
	projRes.z = (int)initialStack->getResolution().z;  // ��    ��ȡ��mrc�ļ���x y z�Ĵ�С

	nviews = projSet->getSize();  //�����ȡmsize  ͶӰ�����б��е�Ԫ���� ;


	sliceZ = setvolumeDimensions.z;  //sliceZ = params.VolumeDimensions().z�Ӳ��������л�ȡ�ع���������ĳߴ磨VolumeDimensions��������ֵ������sliceZ��sliceX�����������ģ����������Ҫ�ֶ�����sliceX��sliceZ��ֵ��
	sliceX = setvolumeDimensions.x;//����ά������ݵĳߴ磬Ҳ������XYZ����ά���ϵ������������������� �����ع������200*200*200

	volumeSliceSize = sliceX * sliceZ;//���㵥��������Ƭ�Ĵ�С�����ֽڼƣ���ͨ��sliceX��sliceZ��˲���ֵ��volumeSliceSize��

	addition = 0 * nviews;// ȫ������ƫ�Ʋ�������ͼ��������addition��������ȫ�����Ų���������ͼ��������scale  addition = params.ScalingOffset() * nviews
	scale = 0 / nviews;//����������ݵ�ƫ�ƺͳ߶ȣ�ȷ������������̵�׼ȷ�Ժ���Ч�ԡ� ����Ҫ��Ҫ����ô�Ļ�������

	initAngles();

	if (1)//�Ƿ�ʹ�ñ߽����  �ݶ��Ȳ�ʹ�� if (params.UseInputEdgeFill())
	{
		edgeFill = params.EdgeFill();
		cout << "Using input mean value: " << edgeFill << endl;
	}
	else
	{
		edgeFill = initialStack->getStackMean();
		cout << "Using computed mean value: " << edgeFill << endl;//������û����
	}

	globalMin = FLT_MAX;
	globalMax = -FLT_MAX;
	globalMean = 0.0;

}  //��

void CTF3d::initAngles()
{
	// Set up trigonometric tables, then convert angles to radians  �������Ǳ�Ȼ�󽫽Ƕ�ת��Ϊ����
	cbet.resize(nviews);  //�ֱ����ڴ洢ÿ��ͶӰ�ӽǶ�Ӧ������ֵ������ֵ��
	sbet.resize(nviews);
	angles.resize(nviews); //�Թ�

	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++)
		listOfIndices.push_back(it.second());  //����Ĳ��Ļ�û��� 


	float degreesToRadiansFactor = M_PI / 180.0; //����ת����ֵ����

	for (unsigned int i = 0; i < nviews; i++)
	{
		angles[i] = microscopeGeometry->getAngleInDegrees(listOfIndices[i]);
		float thetanv = angles[i] + 0.0f;//params.Offset().x=0.0f
		if (thetanv > 180.0)
			thetanv = thetanv - 360.0;

		if (thetanv <= -180.0)
			thetanv = thetanv + 360.0;

		cbet[i] = cos(thetanv * degreesToRadiansFactor);

		//Keep cosine from going to zero so it can be divided by �������Ҳ�Ϊ�㣬�������Ϳ��Գ���
		if (fabs(cbet[i]) < 1.e-6)
			cbet[i] = sign(cbet[i]) * 1.e-6;

		sbet[i] = -sin(thetanv * degreesToRadiansFactor);
		angles[i] = -degreesToRadiansFactor * (angles[i] + 0.0f);
		//�Ҽǵú����Ĳ�������������ʵ��ֵ��
	}


}

/* Precomputes projection boundaries for each row in xz slice
 * Assuming zero x-tilt, these values are same for all slices
 * Ԥ�ȼ���xz��Ƭ��ÿ�е�ͶӰ�߽�
 *����x��бΪ�㣬��������Ƭ����Щֵ����ͬ
 */
void CTF3d::computeProjectionBoundaries()
{
	projectionBoundaries.resize(nviews);//��ʼ��һ����ΪprojectionBoundaries�Ķ�ά���飬��СΪͶӰ��ͼ����nviews���� �����Ƭ�߶ȣ�sliceZ��

	float delxx = 0.0f;		// ��б��������ͼ�����ĵ�ƫ�ơ�Ĭ��ֵΪ��ƫ�ƻ���ת
	float yoffset = 0.0f; 	// �����Ƭ�����ݵĴ�ֱƫ�ơ���ʵ��zƫ��
	novaCTF::Vec2i subsetStart = novaCTF::Vec2i(0, 0);  //�����и���

	float xcenin = projRes.x / 2.0 + 0.5 - subsetStart.x; 	//������Ƭ����������

	float xoffAdj = 0 - (projRes.x / 2 + subsetStart.x - projRes.x / 2);//float xoffAdj = params.ZShift().x - (projRes.x / 2 + subsetStart.x - projRes.x / 2)
	//������� Z �᷽���ϵ�ʵ�ʺ���ƫ���������ͶӰͼ�����ĺ��Ӽ���ʼ������������ƫ����֮��Ĳ�ֵ�������ֵ��������У׼�����ͶӰ����ά�ع������е�λ�á�
	float xcen = sliceX / 2 + 0.5 + delxx + xoffAdj;	// �����Ƭ������x����
	float ycen = sliceZ / 2 + 0.5 + yoffset;			// �����Ƭ������y����
	
	for (unsigned int iv = 0; iv < nviews; iv++) //��������ͶӰ��ͼ
	{
		projectionBoundaries[iv].resize(sliceZ);//Ϊ�䴴��һ������sliceZ��Ԫ�ص�һά���飬���ڴ洢ͶӰ�߽���Ϣ��

		// Set view angle �����ӽ�
		float cosBeta = cbet[iv];
		float sinBeta = sbet[iv];  //����ÿ��ͶӰ�ӽǵĦ½ǵ����Һ�����ֵ

		for (unsigned int i = 0; i < sliceZ; i++)
		{
			float zz = (i + 1 - ycen);
			float zpart = zz * sinBeta + xcenin + delxx;

			float x = cosBeta;
			if (fabs(cosBeta) < 0.001)
				x = 0.001 * sign(cosBeta);

			float projStart = (1.0 - zpart) / x + xcen;
			float projEnd = (projRes.x - zpart) / x + xcen;

			if (projEnd < projStart)
				novaCTF::swap(projStart, projEnd);

			int projStartPixelIndex = projStart;

			if (projStartPixelIndex < projStart)
				projStartPixelIndex = projStartPixelIndex + 1;

			projStartPixelIndex = max(projStartPixelIndex - 1, 0);

			int projEndPixelIndex = projEnd;

			if (projEndPixelIndex == projEnd)
				projEndPixelIndex = projEndPixelIndex - 1;

			projEndPixelIndex = min(projEndPixelIndex - 1, (int)sliceX - 1);

			if (projStartPixelIndex <= projEndPixelIndex)
			{
				x = projStartPixelIndex - xcen + 1;
				projectionBoundaries[iv][i].position = zpart + x * cosBeta - 1;
				projectionBoundaries[iv][i].projectionImpact = projEndPixelIndex - projStartPixelIndex + 1;
				projectionBoundaries[iv][i].fillWithProjectionValue = true;

			}
			else
			{
				projEndPixelIndex = 0;
				projectionBoundaries[iv][i].fillWithProjectionValue = false;
			}

			projectionBoundaries[iv][i].startingIndex = projStartPixelIndex;
			projectionBoundaries[iv][i].endingIndex = projEndPixelIndex;

		}
	}
}

/* Actual projection of one or more slices
 * Goes over all slices in the batch and over all projections.
 * For the currently processed slice it takes every row and based on precomputed projection boundaries
 * adds to the voxels values from input projections
 * 
 * һ��������Ƭ��ʵ��ͶӰ���������е�������Ƭ������ͶӰ�����ڵ�ǰ�������Ƭ����ȡÿһ�У�������Ԥ�ȼ����ͶӰ�߽罫����ͶӰ�е�����ֵ��ӵ�������
 */
void CTF3d::projectSlice()
{
	int projectionIndex = 0; //projectionIndex�����ڼ�¼��ǰͶӰ��λ��������
	for (unsigned int sliceID = 0; sliceID < slicesToProcessAtOnce; sliceID++) 
		//sliceID��������ǰ����Ҫ�������Ƭ��slicesToProcessAtOnce����  ����slicesToProcessAtOnce=1  �������ѭ��Ҫѭ��һ��
	{
		size_t volumeSliceOffset = sliceID * volumeSliceSize;//���������Ƭ�����ݻ�������ƫ����volumeSliceOffset ѭ��һ�εĻ�������0

		for (unsigned int iv = 0; iv < nviews; iv++) //�������е�ͶӰ��ͼ��nviews����
		{
			int volumeIndex = volumeSliceSize * sliceID; //��ʼ��volumeIndex����ʾ��ǰ�������������ά��������е�������
			bool computeSliceStatistics = false; //����computeSliceStatistics��־�������ǰ��ͼ�����һ����ͼ��������Ϊtrue���Ա�ͳ����Ƭ���ݵ����ֵ����Сֵ��ƽ��ֵ��ͳ����Ϣ��
			if (iv == nviews - 1)
				computeSliceStatistics = true;

			for (unsigned int i = 0; i < sliceZ; i++)//��ÿ����Ƭ��Z���ϵ�ÿ��λ�ã�i��  sliceZ=200
			{
				if (projectionBoundaries[iv][i].fillWithProjectionValue) //�����λ����Ҫ���ͶӰֵ��projectionBoundaries[iv][i].fillWithProjectionValueΪtrue��
				{
					//projectionBoundaries ��һ����ά���飬��СΪͶӰ��ͼ����nviews���� �����Ƭ�߶ȣ�sliceZ��

					int startIndex = volumeIndex; 
					int endIndex = volumeIndex + projectionBoundaries[iv][i].startingIndex; //���ݱ߽���Ϣ�����һ����Ҫ���ͶӰֵ����������startIndex�ͽ�������endIndex��

					for (int ind = startIndex; ind < endIndex; ind++)
					{
						currentVolumeSlice[ind] = currentVolumeSlice[ind] + edgeFill;
						volumeIndex++;
						if (computeSliceStatistics)
						{
							currentVolumeSlice[ind] = currentVolumeSlice[ind] * scale + addition;
							cout << currentVolumeSlice[ind] << endl;
							globalMax = max(globalMax, currentVolumeSlice[ind]);
							globalMin = min(globalMin, currentVolumeSlice[ind]);
							globalMean += currentVolumeSlice[ind];
						}
					}
					//�����������ڵ�����������������currentVolumeSlice��������ۼӲ�����������computeSliceStatistics��־�������ֵ����Сֵ��ƽ��ֵ��
					

					//����computeOneRow������һ������ǰ�����ݣ��ú�������ͶӰ�߽硢�Ƕȡ�����ֵ����Ϣ���о�ȷ�ķ�ͶӰ���㡣
					
					computeOneRow(volumeIndex, projectionIndex, projectionBoundaries[iv][i].projectionImpact, projectionBoundaries[iv][i].position, cbet[iv], computeSliceStatistics, iv, volumeSliceOffset);
				
				}


				int startIndex2 = volumeIndex;
				int endIndex2 = volumeIndex + sliceX - projectionBoundaries[iv][i].endingIndex - 1;
				for (int ind = startIndex2; ind < endIndex2; ind++)
				{
					currentVolumeSlice[ind] = currentVolumeSlice[ind] + edgeFill;
					volumeIndex++;

					if (computeSliceStatistics)
					{
						currentVolumeSlice[ind] = currentVolumeSlice[ind] * scale + addition;
						globalMax = max(globalMax, currentVolumeSlice[ind]);
						globalMin = min(globalMin, currentVolumeSlice[ind]);
						globalMean += currentVolumeSlice[ind];
					}
				}

			}
			projectionIndex = projectionIndex + projRes.x;
		}
	}
	// ����������ɺ�currentVolumeSlice���齫������ǰ���δ������Ƭ����ͶӰ��������ά�������ݣ�
	// ͬʱ������˸���Ƭ��ȫ�����ֵ����Сֵ��ƽ��ֵ����Щ��Ϣ�����ں��������ݹ�һ������ʾ�Ȳ�����
}

/* Computes values for all voxels affected by the currently processed projection
 * It uses simple linear interpolation between two adjacent projection pixel values
 * �����ܵ�ǰ����ͶӰӰ����������ص�ֵ����ʹ����������ͶӰ����ֵ֮��ļ����Բ�ֵ
 * 
 * ��δ���ĺ��������д���ͶӰ���ݣ���ÿһ������Ӧ�����Բ�ֵ����ԭ������ά�ռ��е���ʵǿ�ȣ�����ѡ�ؼ��㷴ͶӰ�����в�����ͳ����Ϣ��
 */
void CTF3d::computeOneRow(int& volumeSliceIndex, int projectionIndex, int numberOfVoxels, double xPosition, float cosAngle, bool computeSliceStatistics, unsigned int defocusProjIndex, size_t volumeSliceOffset)
{
	for (int j = 0; j < numberOfVoxels; j++)
	{
		int pixelIndex = xPosition;   //xPosition����ͶӰͼ���ϵ���������pixelIndex��λ��������������֮��ı�������fraction
		float fraction = xPosition - pixelIndex;
		unsigned int defocusIndex = sliceDefocusSplits[defocusProjIndex][volumeSliceIndex - volumeSliceOffset];


		//��sliceDefocusSplits�����л�ȡ��ǰ����������Ƭ��Ӧ�Ľ�������defocusIndex           ����û�л�ý���ֵ
		currentVolumeSlice[volumeSliceIndex] = currentVolumeSlice[volumeSliceIndex] + (1.0f - fraction) * currentRows[defocusIndex][projectionIndex + pixelIndex] + fraction * currentRows[defocusIndex][projectionIndex + pixelIndex + 1];
		//ʹ�����Բ�ֵ���㵱ǰ���ص�ֵ����ϵ�ǰ������ͶӰͼ���϶�Ӧ���ص��������ֵ���ұ�����ֵ���Լ����ǵı�������fraction����������صĳ�����ͶӰֵ�����ۼӵ�currentVolumeSlice[volumeSliceIndex]
		if (computeSliceStatistics)
		{
			//���computeSliceStatisticsΪ�棬��ʾ��Ҫ���㵱ǰ��Ƭ��һЩͳ����Ϣ��
			//����ǰ���ص�ֵ����һ������ϵ��scale���ټ���һ��ƫ��ֵaddition��
			//����ȫ�����ֵglobalMax��ȫ����СֵglobalMin��
			//����ȫ��ƽ��ֵglobalMean������ǰ���ص�ֵ�ۼӽ�ȥ��
			currentVolumeSlice[volumeSliceIndex] = currentVolumeSlice[volumeSliceIndex] * scale + addition;
			globalMax = max(globalMax, currentVolumeSlice[volumeSliceIndex]);
			globalMin = min(globalMin, currentVolumeSlice[volumeSliceIndex]);
			globalMean += currentVolumeSlice[volumeSliceIndex];
		}
		//�ƶ�����һ�����أ���volumeSliceIndex����1��������xPositionΪ��һ�����ص�λ�ã���λ���ǵ�ǰ����λ�ü���һ����ͶӰ�Ƕ��йص�����cosAngle��
		volumeSliceIndex = volumeSliceIndex + 1;
		xPosition = xPosition + cosAngle;
	}
}

void CTF3d::setDataSize()  //�й��ڴ�Ĳ���
{

	float safeBuffer = 50.0f; //��MBΪ��λ-Ϊ�˰�ȫ�����ʵ�ʳ���ֻ��Ҫ��Լ20MB����ȣ�

	size_t parameterSetupSize = sizeof(ParameterSetup);//�洢ParameterSetup�ṹ��������ڴ�
	size_t mrcStackSize = sizeof(MRCStack) * numberOfStacks; //mrcStackSize���洢MRCStack��������������ڴ棬��������numberOfStacks  ��������*1
	size_t ctfClassSize = sizeof(CTF3d); //ctfClassSize���洢CTF3d���������������ڴ档
	size_t defocusGridSize = sizeof(unsigned int) * volumeSliceSize * nviews;//defocusGridSize���洢��ɢ��դ��defocus grid��������ڴ棬���С����volumeSliceSize��nviews
	size_t boundariesSize = sizeof(projectionBoundaries) + sizeof(voxelProjection) * projectionBoundaries.capacity();
	//�洢projectionBoundaries���������voxelProjection����������ڴ档

	float baseSum = parameterSetupSize + mrcStackSize + ctfClassSize + defocusGridSize + boundariesSize;
	baseSum = baseSum / 1000000.f + safeBuffer;
	//������ӣ������ϰ�ȫ���������õ������ڴ�������������MBΪ��λ����

	float memoryLimit = params.MemoryLimit() - baseSum;
	//�����û����õ��ڴ����ƣ�params.MemoryLimit()����ȥ�����ڴ������������õ�ʣ������ڴ������ݵ��ڴ棨memoryLimit����
	if (memoryLimit <= 0)
	{
		//���memoryLimitС�ڵ���0��˵���ڴ����ƹ��ͣ��޷����������Ĵ�����ʱ����slicesToProcessAtOnce����Ϊ1��reportIndex��Ϊ200�������������Ϣ��
		memoryLimit = 0;
		reportIndex = 200;
		slicesToProcessAtOnce = 1;
		processedDataSize = volumeSliceSize;

		if (params.MemoryLimit() != 0)
			cout << "Memory limit was too low!!!" << std::endl;

		cout << "The number of slices to be processed was set to 1." << std::endl << std::endl;

		return;  //����ֱ�ӷ��� ��������±ߵĲ���
	}

	size_t doubleVolumeSliceSize = 2 * sizeof(float) * volumeSliceSize;	//���ǵ�д��ʱ������Ҫ�������ʱ�ռ䣬������������Ƭ��˫����С��
	size_t loadingCurrentRowsSize = projRes.x * nviews * sizeof(float);//���ص�ǰ�����������ڴ�
	size_t currentRowsSize = projRes.x * nviews * numberOfStacks * sizeof(float);//�洢��ǰ�����������ڴ���ܺ�

	float sumForOneSlice = (doubleVolumeSliceSize + loadingCurrentRowsSize + currentRowsSize) / 1000000.f;
	//����ʣ���ڴ�memoryLimit�͵�����Ƭ������ڴ棬����һ���Կ��Դ������Ƭ������slicesToProcessAtOnce��  slicesToProcessAtOnce ����������1
	slicesToProcessAtOnce = floor(memoryLimit / sumForOneSlice);

	//ȷ��projRes.y����ʾY���ϵ���Ƭ�������ܱ�slicesToProcessAtOnce������������������������slicesToProcessAtOnce��
	while ((projRes.y % slicesToProcessAtOnce) != 0) 
	{
		slicesToProcessAtOnce--;
	}

	reportIndex = slicesToProcessAtOnce;
	while (reportIndex < 200) //�趨reportIndexΪslicesToProcessAtOnce��2���ݴη�������ȡ������ݣ���ȷ�������ڵ���200�����ڶ��ڱ��洦�����
	{
		reportIndex = 2 * reportIndex;
	}
	//��󣬸��ݴ������Ƭ��������processedDataSize
	processedDataSize = volumeSliceSize * slicesToProcessAtOnce;
	cout << "Number of slices processed at once was set to:" << slicesToProcessAtOnce << std::endl;

}

void CTF3d::backproject()
{
	currentVolumeSlice.resize(processedDataSize); //���ڴ洢��ǰ�����������Ƭ���ݣ���СΪprocessedDataSize  ������40000
	currentRows.resize(numberOfStacks);//���ڴ洢��ͶӰ���ݼ��ж�ȡ�ĵ�ǰͶӰ�����ݣ���numberOfStacks��Ԫ�أ�ÿ��Ԫ�ش�СΪͶӰ��ȳ���ͶӰ��������һ�δ������Ƭ����

	for (unsigned int inputStackIndex = 0; inputStackIndex < numberOfStacks; inputStackIndex++)  //����û��ʹ������3dCTF���� numberOfStacks�Ѿ������ó�1
		currentRows[inputStackIndex].resize(projRes.x * nviews * slicesToProcessAtOnce);  //currentRows ��ά�����С  1*40000

	//currentRows Ҳ��ǰ������ΪcurrentRows.resize��1��   

	for (int sliceNumber = 0; sliceNumber < projRes.y; sliceNumber += slicesToProcessAtOnce) //  projRes.y����mrc�ļ���ת�� projRes.y��41
	{
		cout << (size_t)sliceNumber << endl;
		//ѭ������sliceNumber��0��ʼ������ÿ������slicesToProcessAtOnce��ֱ���ﵽͶӰ���ݵĸ߶ȣ�projRes.y��Ϊֹ��
		if (sliceNumber % reportIndex == 0)
		{
			cout << "Started processing slice number " << (size_t)sliceNumber << "." << endl;
		}

		for (unsigned int inputStackIndex = 0; inputStackIndex < numberOfStacks; inputStackIndex++)
		{
			//����ÿ�������ջ����inputStackIndex = 0��numberOfStacks - 1��������inputStacks[inputStackIndex]->readProjections������
			// �������ջ�ж�ȡָ��������slicesToProcessAtOnce����ͶӰ��Ƭ���ݣ��������ݴ洢��currentRows[inputStackIndex]�����ָ��λ�ã�
			// ��ʼλ��ΪcurrentRows[inputStackIndex][0]

			inputStacks[inputStackIndex]->readProjections(&currentRows[inputStackIndex][0], slicesToProcessAtOnce, sliceNumber);

			//��δ�����������ڣ����ڵ� inputStackIndex ������ջ��
			// �Ӹ�ջ�ж�ȡ�� sliceNumber ��ʼ������ slicesToProcessAtOnce ����ƬͶӰ���ݣ�������Щ���ݴ��뵽 currentRows ��Ӧ�������С�
		}

		// clear out the slice �ڽ��з�ͶӰ����ǰ��ʹ��std::fill������currentVolumeSlice����ȫ�����㣬��׼����ű��μ�������
		fill(currentVolumeSlice.begin(), currentVolumeSlice.end(), 0.0f);
		//����projectSlice���������ݶ�ȡ��ͶӰ���ݽ��з�ͶӰ���㣬����������currentVolumeSlice�����С�

		projectSlice();

		//ʹ��VolumeIO��ľ�̬����writeVolumeSliceInFloat������ǰ����õ���������Ƭ���ݣ�currentVolumeSlice���Ը�������ʽд�뵽ָ��������ļ��У�
		// �ļ�����Դ��params.OutputFilename()��ͬʱָ�����ģʽΪparams.OutputMode()��
		string outputname = "E:\\VS studio\\proj\\project2\\test_rec_output\\output_4.18.6.rec";
		VolumeIO::writeVolumeSliceInFloat(outputname, currentVolumeSlice, processedDataSize, 2);//û����
	}
}

// Main loop over slices perpendicular to tilt axis  ��ֱ����б�����Ƭ�ϵ���ѭ��


string CTF3d::generateFilename(string originalFilename, unsigned int number)
{
	stringstream newFilename;
	newFilename << originalFilename << "_" << number;
	string ret = newFilename.str();
	return ret;
}

void CTF3d::writeOutDefocusSlices(vector<vector<unsigned int> >& sliceDefocusSplits) //�����������û���õ�����ΪprepareCTF()��������ԭ��
{
	if (!params.WriteOutDefocusSlices())
	{
		cout << "��������ж�" << endl;//����ж������ʱû�㶮
		return;
	}

	cout << "û�н�������ж�" << endl;
	vector<float> sliceData;  //��ʼ��һ������������sliceData�����ڴ洢���н�ɢ��Ƭ���ݡ�

	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++) //��������ͶӰ��ͶӰ��projSet�ĵ�������������ÿ��ͶӰ��
	{
		unsigned int projIndex = it.second();// ��ȡͶӰ����projIndex��
		for (size_t i = 0; i < sliceDefocusSplits[projIndex].size(); i++)
			sliceData.push_back(sliceDefocusSplits[projIndex][i]);
		//��ǰͶӰ��Ӧ�Ľ�ɢ��Ƭ���ݣ�sliceDefocusSplits[projIndex]�����Ԫ����ӵ�sliceData�����С�
	}

	VolumeIO::write(sliceData, Vec3ui(params.VolumeDimensions().x, params.VolumeDimensions().z, initialStack->getNumberOfProjections()), params.OutputFilename() + "_defocusSlices.mrc", 0, novaCTF::VolumeRotation::ALONG_XY);

}