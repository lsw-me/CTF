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
	string TiltAnglesFilename = "E:\\VS studio\\proj\\project2\\IS002_291013_005.tlt";//角度文件
	projSet->init(initialStack->getNumberOfProjections());//get函数返回 mHeaderMRC.nz，
	microscopeGeometry = new Geometry(*initialStack, setvolumeDimensions, TiltAnglesFilename, 0.0f, 0,0.0f);//这里传入的volumeDimensions 
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
	VolumeIO::writeHeader(firstFileName, outputname, setvolumeDimensions, 2, novaCTF::VolumeRotation::ALONG_XZ);  // 旋转从xz
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
	//负责根据传入的参数对体积数据进行值域转换，并将转换后的数据保存到指定的输出文件中。
}

void CTF3d::loadInputStacks()  //这里因为novactf会产生许多堆栈文件，利用这些文件一起进行校正，我这里只有一个校正后的mrc文件包含着倾斜序列们，这里需要去掉或者优化
{
	//函数将加载单个MRC文件（即firstFileName指定的文件），并将它作为一个堆栈实例压入inputStacks容器中。同时，函数在此处结束执行并返回。
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

	sliceDefocusSplits.resize(initialStack->getNumberOfProjections());//返回HeaderMRC.nz值，若旋转图像，则nz的值为200

	size_t xzSliceSize = setvolumeDimensions.x*setvolumeDimensions.z;//这个切面是xz方向也就是沿着y轴


	//
	//在没有 CTF 的情况下，只需使用一个带有 0 的切片来指向唯一使用的输入堆栈
	//if (!params.Use3DCTF()) //判断语句也可以扔掉
	// 
	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++)
	{
		unsigned int projIndex = it.second();
		sliceDefocusSplits[projIndex].resize(xzSliceSize);
		std::fill(sliceDefocusSplits[projIndex].begin(), sliceDefocusSplits[projIndex].end(), 0);
		//循环结束后 sliceDefocusSplits二维数组大小为 200*40000
	}
		numberOfStacks = 1;
		currentRows.resize(1);
		return;

	//这里注释掉的原因，他这里输入堆栈叫做3dCTF 然而我输入的堆栈是根据三维信息进行校正后的结果，因此我可以只输入一个mrc，但是里面的像素值是利用三维信息校正过的
	//所以不使用它的3dCTF
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
	projRes.z = (int)initialStack->getResolution().z;  // √    读取的mrc文件中x y z的大小

	nviews = projSet->getSize();  //这里读取msize  投影索引列表中的元素数 ;


	sliceZ = setvolumeDimensions.z;  //sliceZ = params.VolumeDimensions().z从参数设置中获取重构体素网格的尺寸（VolumeDimensions），并赋值给变量sliceZ和sliceX。根据上下文，这里可能需要手动调整sliceX和sliceZ的值。
	sliceX = setvolumeDimensions.x;//是三维体积数据的尺寸，也就是在XYZ三个维度上的像素数量或体素数量 我想重构结果是200*200*200

	volumeSliceSize = sliceX * sliceZ;//计算单个体素切片的大小（以字节计），通过sliceX和sliceZ相乘并赋值给volumeSliceSize。

	addition = 0 * nviews;// 全局缩放偏移参数和视图数量计算addition，并根据全局缩放参数除以视图数量计算scale  addition = params.ScalingOffset() * nviews
	scale = 0 / nviews;//负责调整数据的偏移和尺度，确保后续计算过程的准确性和有效性。 这里要不要改怎么改还有问题

	initAngles();

	if (1)//是否使用边界填充  暂定先不使用 if (params.UseInputEdgeFill())
	{
		edgeFill = params.EdgeFill();
		cout << "Using input mean value: " << edgeFill << endl;
	}
	else
	{
		edgeFill = initialStack->getStackMean();
		cout << "Using computed mean value: " << edgeFill << endl;//到这里没问题
	}

	globalMin = FLT_MAX;
	globalMax = -FLT_MAX;
	globalMean = 0.0;

}  //√

void CTF3d::initAngles()
{
	// Set up trigonometric tables, then convert angles to radians  设置三角表，然后将角度转换为弧度
	cbet.resize(nviews);  //分别用于存储每个投影视角对应的余弦值和正弦值。
	sbet.resize(nviews);
	angles.resize(nviews); //对勾

	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++)
		listOfIndices.push_back(it.second());  //这里改不改还没清楚 


	float degreesToRadiansFactor = M_PI / 180.0; //度数转弧度值因子

	for (unsigned int i = 0; i < nviews; i++)
	{
		angles[i] = microscopeGeometry->getAngleInDegrees(listOfIndices[i]);
		float thetanv = angles[i] + 0.0f;//params.Offset().x=0.0f
		if (thetanv > 180.0)
			thetanv = thetanv - 360.0;

		if (thetanv <= -180.0)
			thetanv = thetanv + 360.0;

		cbet[i] = cos(thetanv * degreesToRadiansFactor);

		//Keep cosine from going to zero so it can be divided by 保持余弦不为零，这样它就可以除以
		if (fabs(cbet[i]) < 1.e-6)
			cbet[i] = sign(cbet[i]) * 1.e-6;

		sbet[i] = -sin(thetanv * degreesToRadiansFactor);
		angles[i] = -degreesToRadiansFactor * (angles[i] + 0.0f);
		//我记得后续的部分中这里是有实际值的
	}


}

/* Precomputes projection boundaries for each row in xz slice
 * Assuming zero x-tilt, these values are same for all slices
 * 预先计算xz切片中每行的投影边界
 *假设x倾斜为零，则所有切片的这些值都相同
 */
void CTF3d::computeProjectionBoundaries()
{
	projectionBoundaries.resize(nviews);//初始化一个名为projectionBoundaries的二维数组，大小为投影视图数（nviews）× 输出切片高度（sliceZ）

	float delxx = 0.0f;		// 倾斜轴与输入图像中心的偏移。默认值为无偏移或旋转
	float yoffset = 0.0f; 	// 输出切片中数据的垂直偏移。事实上z偏移
	novaCTF::Vec2i subsetStart = novaCTF::Vec2i(0, 0);  //这里有更改

	float xcenin = projRes.x / 2.0 + 0.5 - subsetStart.x; 	//输入切片的中心坐标

	float xoffAdj = 0 - (projRes.x / 2 + subsetStart.x - projRes.x / 2);//float xoffAdj = params.ZShift().x - (projRes.x / 2 + subsetStart.x - projRes.x / 2)
	//计算的是 Z 轴方向上的实际横向偏移量与基于投影图像中心和子集起始点计算出的理论偏移量之间的差值。这个差值可能用于校准或调整投影在三维重构过程中的位置。
	float xcen = sliceX / 2 + 0.5 + delxx + xoffAdj;	// 输出切片的中心x坐标
	float ycen = sliceZ / 2 + 0.5 + yoffset;			// 输出切片的中心y坐标
	
	for (unsigned int iv = 0; iv < nviews; iv++) //遍历所有投影视图
	{
		projectionBoundaries[iv].resize(sliceZ);//为其创建一个包含sliceZ个元素的一维数组，用于存储投影边界信息。

		// Set view angle 设置视角
		float cosBeta = cbet[iv];
		float sinBeta = sbet[iv];  //计算每个投影视角的β角的余弦和正弦值

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
 * 一个或多个切片的实际投影覆盖批次中的所有切片和所有投影。对于当前处理的切片，它取每一行，并基于预先计算的投影边界将输入投影中的体素值添加到体素中
 */
void CTF3d::projectSlice()
{
	int projectionIndex = 0; //projectionIndex：用于记录当前投影的位置索引。
	for (unsigned int sliceID = 0; sliceID < slicesToProcessAtOnce; sliceID++) 
		//sliceID：遍历当前批次要处理的切片（slicesToProcessAtOnce）。  这里slicesToProcessAtOnce=1  所以这个循环要循环一次
	{
		size_t volumeSliceOffset = sliceID * volumeSliceSize;//计算体积切片在数据缓冲区的偏移量volumeSliceOffset 循环一次的话这里是0

		for (unsigned int iv = 0; iv < nviews; iv++) //遍历所有的投影视图（nviews个）
		{
			int volumeIndex = volumeSliceSize * sliceID; //初始化volumeIndex，表示当前处理的体素在三维体积数据中的索引。
			bool computeSliceStatistics = false; //设置computeSliceStatistics标志，如果当前视图是最后一个视图，则设置为true，以便统计切片数据的最大值、最小值、平均值等统计信息。
			if (iv == nviews - 1)
				computeSliceStatistics = true;

			for (unsigned int i = 0; i < sliceZ; i++)//对每个切片在Z轴上的每个位置（i）  sliceZ=200
			{
				if (projectionBoundaries[iv][i].fillWithProjectionValue) //如果该位置需要填充投影值（projectionBoundaries[iv][i].fillWithProjectionValue为true）
				{
					//projectionBoundaries 是一个二维数组，大小为投影视图数（nviews）× 输出切片高度（sliceZ）

					int startIndex = volumeIndex; 
					int endIndex = volumeIndex + projectionBoundaries[iv][i].startingIndex; //根据边界信息计算第一个需要填充投影值的体素索引startIndex和结束索引endIndex。

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
					//遍历该区间内的所有体素索引，对currentVolumeSlice数组进行累加操作，并根据computeSliceStatistics标志更新最大值、最小值和平均值。
					

					//调用computeOneRow函数进一步处理当前行数据，该函数根据投影边界、角度、余弦值等信息进行精确的反投影计算。
					
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
	// 整个过程完成后，currentVolumeSlice数组将包含当前批次处理的切片经反投影计算后的三维体素数据，
	// 同时计算出了该切片的全局最大值、最小值和平均值。这些信息可用于后续的数据归一化或显示等操作。
}

/* Computes values for all voxels affected by the currently processed projection
 * It uses simple linear interpolation between two adjacent projection pixel values
 * 计算受当前处理投影影响的所有体素的值，它使用两个相邻投影像素值之间的简单线性插值
 * 
 * 这段代码的核心是逐行处理投影数据，对每一个体素应用线性插值法还原其在三维空间中的真实强度，并可选地计算反投影过程中产生的统计信息。
 */
void CTF3d::computeOneRow(int& volumeSliceIndex, int projectionIndex, int numberOfVoxels, double xPosition, float cosAngle, bool computeSliceStatistics, unsigned int defocusProjIndex, size_t volumeSliceOffset)
{
	for (int j = 0; j < numberOfVoxels; j++)
	{
		int pixelIndex = xPosition;   //xPosition计算投影图像上的像素索引pixelIndex和位于两个相邻像素之间的比例因子fraction
		float fraction = xPosition - pixelIndex;
		unsigned int defocusIndex = sliceDefocusSplits[defocusProjIndex][volumeSliceIndex - volumeSliceOffset];


		//从sliceDefocusSplits数组中获取当前体素所在切片对应的焦深索引defocusIndex           这里没有获得焦深值
		currentVolumeSlice[volumeSliceIndex] = currentVolumeSlice[volumeSliceIndex] + (1.0f - fraction) * currentRows[defocusIndex][projectionIndex + pixelIndex] + fraction * currentRows[defocusIndex][projectionIndex + pixelIndex + 1];
		//使用线性插值计算当前体素的值：结合当前体素在投影图像上对应像素的左边像素值和右边像素值，以及它们的比例因子fraction，计算出体素的初步反投影值，并累加到currentVolumeSlice[volumeSliceIndex]
		if (computeSliceStatistics)
		{
			//如果computeSliceStatistics为真，表示需要计算当前切片的一些统计信息：
			//将当前体素的值乘以一个缩放系数scale，再加上一个偏置值addition。
			//更新全局最大值globalMax和全局最小值globalMin。
			//更新全局平均值globalMean，将当前体素的值累加进去。
			currentVolumeSlice[volumeSliceIndex] = currentVolumeSlice[volumeSliceIndex] * scale + addition;
			globalMax = max(globalMax, currentVolumeSlice[volumeSliceIndex]);
			globalMin = min(globalMin, currentVolumeSlice[volumeSliceIndex]);
			globalMean += currentVolumeSlice[volumeSliceIndex];
		}
		//移动到下一个体素：将volumeSliceIndex递增1，并更新xPosition为下一个体素的位置，该位置是当前体素位置加上一个与投影角度有关的增量cosAngle。
		volumeSliceIndex = volumeSliceIndex + 1;
		xPosition = xPosition + cosAngle;
	}
}

void CTF3d::setDataSize()  //有关内存的操作
{

	float safeBuffer = 50.0f; //以MB为单位-为了安全起见，实际程序只需要大约20MB（库等）

	size_t parameterSetupSize = sizeof(ParameterSetup);//存储ParameterSetup结构体所需的内存
	size_t mrcStackSize = sizeof(MRCStack) * numberOfStacks; //mrcStackSize：存储MRCStack对象数组所需的内存，数量等于numberOfStacks  单张输入*1
	size_t ctfClassSize = sizeof(CTF3d); //ctfClassSize：存储CTF3d类对象自身所需的内存。
	size_t defocusGridSize = sizeof(unsigned int) * volumeSliceSize * nviews;//defocusGridSize：存储焦散格栅（defocus grid）所需的内存，其大小基于volumeSliceSize和nviews
	size_t boundariesSize = sizeof(projectionBoundaries) + sizeof(voxelProjection) * projectionBoundaries.capacity();
	//存储projectionBoundaries和相关联的voxelProjection数组所需的内存。

	float baseSum = parameterSetupSize + mrcStackSize + ctfClassSize + defocusGridSize + boundariesSize;
	baseSum = baseSum / 1000000.f + safeBuffer;
	//需求相加，并加上安全缓冲区，得到基础内存消耗总量（以MB为单位）。

	float memoryLimit = params.MemoryLimit() - baseSum;
	//根据用户设置的内存限制（params.MemoryLimit()）减去基础内存消耗总量，得到剩余可用于处理数据的内存（memoryLimit）。
	if (memoryLimit <= 0)
	{
		//如果memoryLimit小于等于0，说明内存限制过低，无法进行正常的处理。这时，将slicesToProcessAtOnce设置为1，reportIndex设为200，并输出警告信息。
		memoryLimit = 0;
		reportIndex = 200;
		slicesToProcessAtOnce = 1;
		processedDataSize = volumeSliceSize;

		if (params.MemoryLimit() != 0)
			cout << "Memory limit was too low!!!" << std::endl;

		cout << "The number of slices to be processed was set to 1." << std::endl << std::endl;

		return;  //这里直接返回 不会进行下边的操作
	}

	size_t doubleVolumeSliceSize = 2 * sizeof(float) * volumeSliceSize;	//考虑到写入时可能需要额外的临时空间，这里估算体积切片的双倍大小。
	size_t loadingCurrentRowsSize = projRes.x * nviews * sizeof(float);//加载当前行数据所需内存
	size_t currentRowsSize = projRes.x * nviews * numberOfStacks * sizeof(float);//存储当前行数据所需内存的总和

	float sumForOneSlice = (doubleVolumeSliceSize + loadingCurrentRowsSize + currentRowsSize) / 1000000.f;
	//根据剩余内存memoryLimit和单个切片所需的内存，计算一次性可以处理的切片数量（slicesToProcessAtOnce）  slicesToProcessAtOnce 这是设置是1
	slicesToProcessAtOnce = floor(memoryLimit / sumForOneSlice);

	//确保projRes.y（表示Y轴上的切片总数）能被slicesToProcessAtOnce整除，如果不能整除，则减少slicesToProcessAtOnce。
	while ((projRes.y % slicesToProcessAtOnce) != 0) 
	{
		slicesToProcessAtOnce--;
	}

	reportIndex = slicesToProcessAtOnce;
	while (reportIndex < 200) //设定reportIndex为slicesToProcessAtOnce的2的幂次方（向上取最近的幂），确保它大于等于200，用于定期报告处理进度
	{
		reportIndex = 2 * reportIndex;
	}
	//最后，根据处理的切片数量更新processedDataSize
	processedDataSize = volumeSliceSize * slicesToProcessAtOnce;
	cout << "Number of slices processed at once was set to:" << slicesToProcessAtOnce << std::endl;

}

void CTF3d::backproject()
{
	currentVolumeSlice.resize(processedDataSize); //用于存储当前处理的体素切片数据，大小为processedDataSize  这里是40000
	currentRows.resize(numberOfStacks);//用于存储从投影数据集中读取的当前投影行数据，有numberOfStacks个元素，每个元素大小为投影宽度乘以投影数量乘以一次处理的切片数量

	for (unsigned int inputStackIndex = 0; inputStackIndex < numberOfStacks; inputStackIndex++)  //这里没有使用它的3dCTF所以 numberOfStacks已经被设置成1
		currentRows[inputStackIndex].resize(projRes.x * nviews * slicesToProcessAtOnce);  //currentRows 二维数组大小  1*40000

	//currentRows 也在前边设置为currentRows.resize（1）   

	for (int sliceNumber = 0; sliceNumber < projRes.y; sliceNumber += slicesToProcessAtOnce) //  projRes.y由于mrc文件旋转过 projRes.y是41
	{
		cout << (size_t)sliceNumber << endl;
		//循环变量sliceNumber从0开始递增，每次增加slicesToProcessAtOnce，直到达到投影数据的高度（projRes.y）为止。
		if (sliceNumber % reportIndex == 0)
		{
			cout << "Started processing slice number " << (size_t)sliceNumber << "." << endl;
		}

		for (unsigned int inputStackIndex = 0; inputStackIndex < numberOfStacks; inputStackIndex++)
		{
			//对于每个输入堆栈（从inputStackIndex = 0到numberOfStacks - 1），调用inputStacks[inputStackIndex]->readProjections方法，
			// 从输入堆栈中读取指定数量（slicesToProcessAtOnce）的投影切片数据，并将数据存储在currentRows[inputStackIndex]数组的指定位置，
			// 起始位置为currentRows[inputStackIndex][0]

			inputStacks[inputStackIndex]->readProjections(&currentRows[inputStackIndex][0], slicesToProcessAtOnce, sliceNumber);

			//这段代码的意义在于，对于第 inputStackIndex 个输入栈，
			// 从该栈中读取从 sliceNumber 开始的连续 slicesToProcessAtOnce 个切片投影数据，并将这些数据存入到 currentRows 对应的容器中。
		}

		// clear out the slice 在进行反投影计算前，使用std::fill函数将currentVolumeSlice数组全部清零，以准备存放本次计算结果。
		fill(currentVolumeSlice.begin(), currentVolumeSlice.end(), 0.0f);
		//调用projectSlice函数，根据读取的投影数据进行反投影计算，将结果存放在currentVolumeSlice数组中。

		projectSlice();

		//使用VolumeIO类的静态方法writeVolumeSliceInFloat，将当前计算得到的体素切片数据（currentVolumeSlice）以浮点数形式写入到指定的输出文件中，
		// 文件名来源于params.OutputFilename()，同时指定输出模式为params.OutputMode()。
		string outputname = "E:\\VS studio\\proj\\project2\\test_rec_output\\output_4.18.6.rec";
		VolumeIO::writeVolumeSliceInFloat(outputname, currentVolumeSlice, processedDataSize, 2);//没问题
	}
}

// Main loop over slices perpendicular to tilt axis  垂直于倾斜轴的切片上的主循环


string CTF3d::generateFilename(string originalFilename, unsigned int number)
{
	stringstream newFilename;
	newFilename << originalFilename << "_" << number;
	string ret = newFilename.str();
	return ret;
}

void CTF3d::writeOutDefocusSlices(vector<vector<unsigned int> >& sliceDefocusSplits) //这个函数这里没有用到，因为prepareCTF()这个步骤的原因
{
	if (!params.WriteOutDefocusSlices())
	{
		cout << "进入这个判断" << endl;//这个判断语句暂时没搞懂
		return;
	}

	cout << "没有进入这个判断" << endl;
	vector<float> sliceData;  //初始化一个浮点数向量sliceData，用于存储所有焦散切片数据。

	for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++) //遍历所有投影（投影集projSet的迭代器），对于每个投影：
	{
		unsigned int projIndex = it.second();// 获取投影索引projIndex。
		for (size_t i = 0; i < sliceDefocusSplits[projIndex].size(); i++)
			sliceData.push_back(sliceDefocusSplits[projIndex][i]);
		//当前投影对应的焦散切片数据（sliceDefocusSplits[projIndex]）逐个元素添加到sliceData向量中。
	}

	VolumeIO::write(sliceData, Vec3ui(params.VolumeDimensions().x, params.VolumeDimensions().z, initialStack->getNumberOfProjections()), params.OutputFilename() + "_defocusSlices.mrc", 0, novaCTF::VolumeRotation::ALONG_XY);

}