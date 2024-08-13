#include "ctfCorrection.h"
#include "exception.h"
#include <sstream>
#include "fftRoutines.h"
#include "volumeIO.h"
#include "defocusFileFormats.h"

CTFCorrection::CTFCorrection(ParameterSetup& aParams)
{
    params = aParams;
    projSet = new ProjectionSetIdentity();
    inputStack = new MRCStack(params.InputStackName(), false, false, true);
    projSet->init(inputStack->getNumberOfProjections());
}
CTFCorrection::~CTFCorrection()
{
    delete inputStack;
}
void CTFCorrection::run()
{
    initVariables();
    computeFrequencyArray();
    checkAstigmatism();
    correctCTF();
}
void CTFCorrection::initVariables()
{
    defocusFileFormat = params.DefocusFileFormat(); //ctffind4
    ctfCorrectionType = params.CtfCorrectionType();//

    pixelSize = params.PixelSize();
    amplitude = params.Amplitude();
    cs = params.Cs();
    evk = params.Evk();

    //Convert spherical aberration term from mm to Angstroms
    cs = cs * (1.0e7);

    //Convert from nm to Angstroms
    pixelSize = pixelSize * 10.0f;

    //readDefocusFile();
    defocusFileValues.resize(inputStack->getNumberOfProjections());
    DefocusFileFormat* defFile = DefocusFileFormat::createFileFormat(defocusFileFormat);
    defFile->read(params.DefocusFile(), *projSet);
    defFile->getValues(defocusFileValues, *projSet, "microns");

    correctAstigmatism = params.CorrectAstigmatism();
}

void CTFCorrection::computeFrequencyArray()//计算输入图像的频率数组   这个函数不用更改
{
    size_t resX = inputStack->getResolution().x;//读取输入图像的x y值
    size_t resY = inputStack->getResolution().y;//

    // For non-squared images take the larger dimension
    if (resX != resY)
    {
        resX = max(resX, resY);
        resY = max(resX, resY);
    }                                     //如果x y 不同则保留最大值为x y  

    std::vector<float> xArray;    
    std::vector<float> yArray;//两个浮点数向量：xArray和yArray，分别用于存储x轴和y轴方向上的频率值

    for (int i = -floor(resX / 2); i < floor(resX / 2) + (resX % 2); i++)
        xArray.push_back(i / (resX * pixelSize));              //使用两个循环生成x轴方向的频率值：  floor 用于对浮点数进行向下取整操作。
                                                             // 循环范围从 - floor(resX / 2) 到 floor(resX / 2) + (resX % 2)，确保覆盖整个图像的中心区域，包括奇数尺寸时的中心像素。
                                                             // 对于每个i值，将它除以(resX * pixelSize)得到频率值，并将结果添加到xArray中。
                                                             //样的方式生成y轴方向的频率值。
    for (int i = -floor(resY / 2); i < floor(resY / 2) + (resY % 2); i++)
        yArray.push_back(i / (resY * pixelSize));

    frequencyArray.resize(xArray.size());          //初始化一个二维动态数组frequencyArray，其大小为xArray的长度。

    for (size_t x = 0; x < xArray.size(); x++)
        for (size_t y = 0; y < yArray.size(); y++)
            frequencyArray[x].push_back(sqrt(xArray[x] * xArray[x] + yArray[y] * yArray[y])); 
    //再次使用嵌套循环遍历xArray和yArray，计算每个坐标点对应的频率值（根据欧几里得距离公式），并将结果添加到frequencyArray[x]中。
    arraySizeX = xArray.size();
    arraySizeY = yArray.size(); //设置arraySizeX和arraySizeY分别为xArray和yArray的长度。
}

/*If astigmatism correction is required it performs check whether the values are available.
* If they are not, than the warning is written out and the correction is performed without
* the astigmatism correction.
* 如果需要散光校正，它会检查这些值是否可用。如果不是，则写出警告并在没有的情况下执行更正。散光校正
*/
void CTFCorrection::checkAstigmatism() //这个也不用改
{
    if (!correctAstigmatism)
        return;

    float diff = 0.0f;

    for (unsigned int i = 0; i < defocusFileValues.size(); i++)
    {
        diff = diff + fabs(defocusFileValues[i][0] - defocusFileValues[i][1]);
    }

    if (diff < 0.0001)
    {
        cout << "The values necessary for astigmatism correction are not present in the defocus file!!!" << std::endl;
        cout << "The CTF correction will be performed without the astigmatism correction!!!" << std::endl << std::endl;
        correctAstigmatism = false;
    }
}

void CTFCorrection::correctCTF()
{
    std::vector<std::vector<float>> correctedProjections;//存储经过校正后的投影数据，类型为二维浮点数向量。
    std::vector<float> ctfAmp; //用于存储计算出的对比传递函数（Contrast Transfer Function, CTF）的振幅部分。
    std::vector<float> ctfFilter;//用于存储计算出的对比传递函数的滤波器，可能代表相位或振幅矫正系数。

    correctedProjections.resize(inputStack->getNumberOfProjections());//根据输入Stack中投影的数量调整correctedProjections的大小。

    for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++) //使用迭代器it遍历projSet中的所有投影，并获取每个投影在集合中的索引projIndex
    {
        unsigned int projIndex = it.second();
        correctedProjections[projIndex].resize(inputStack->getProjectionSize());
        inputStack->readProjections(&correctedProjections[projIndex][0], 1, projIndex);
        //从输入Stack中读取投影数据到correctedProjections[projIndex]中，按需分配内存空间。
        padProjection(correctedProjections[projIndex], inputStack->getResolution().x, inputStack->getResolution().y);
        //调用padProjection函数来扩展或填充投影图像以适应特定分辨率要求。
        computeCTF(ctfAmp, ctfFilter, defocusFileValues[projIndex][0], defocusFileValues[projIndex][1], defocusFileValues[projIndex][2], defocusFileValues[projIndex][3]);
        //使用提供的焦距参数（四个值分别对应不同的像差项），调用computeCTF函数计算当前投影对应的CTF信息。
        if (ctfCorrectionType == "phaseflip")
        {
            if (!correctAstigmatism)
                FFTRoutines::real2DTransform(arraySizeX, arraySizeY, correctedProjections[projIndex], ctfFilter);
            else
                FFTRoutines::complex2DTransform(arraySizeX, arraySizeY, correctedProjections[projIndex], ctfFilter);
        }
        else //multiplication
        {
            if (!correctAstigmatism)
                FFTRoutines::real2DTransform(arraySizeX, arraySizeY, correctedProjections[projIndex], ctfAmp);
            else
                FFTRoutines::complex2DTransform(arraySizeX, arraySizeY, correctedProjections[projIndex], ctfAmp);
        }

        cropProjection(correctedProjections[projIndex], inputStack->getResolution().x, inputStack->getResolution().y);
        //调用cropProjection函数去除之前填充的额外像素区域，恢复至原始分辨率大小。
    }

    VolumeIO::writeMRCStack(inputStack->getStackHeader(), correctedProjections, params.OutputFilename(), inputStack->getExtraData());

}

void CTFCorrection::computeCTF(std::vector<float>& ctfAmp, std::vector<float>& ctfFilter, float defocus1, float defocus2, float astigmatism, float phaseShift)
{
    // Convert defocii from microns to Angstroms 将输入的两个焦距值从微米（microns）转换为埃（Angstroms）。
    defocus1 = defocus1 * (1.0e4);
    defocus2 = defocus2 * (1.0e4);


    // Calculate electron wavelength 使用普朗克常数、光速、电子静止能量、加速电压以及电子电荷等物理常数来计算电子在指定加速电压下的波长。
    // Most of the equations here are from Mindell and Grigorieff (2003) 

    double h = 6.62606957e-34;
    double c = 299792458;
    float eRest = 511000;
    float v = evk * 1000;
    float eCharge = 1.602e-19;
    double lambda = (c * h) / sqrt(((2 * eRest * v) + (v * v)) * (eCharge * eCharge)) * (pow(10.0, 10));

    // Calculate astigmatic defocus array 初始化并计算非球面像差焦距数组：

    // Initialize defocus array
    std::vector<std::vector<float> > defocusArray; //二维向量 散焦数组

    // Calculate center of image  //计算中心
    float centerX = arraySizeX / 2.0f - 1.0f;
    float centerY = arraySizeY / 2.0f - 1.0f;
    
    /*arraySizeX 表示图像在 x 轴方向上的像素数量。将 arraySizeX 除以 2.0f 是为了得到图像宽度的一半，这通常会得到图像的中心点所在的列索引（假设图像索引从0开始）。
然而，在许多编程环境和图像处理库中，图像坐标的原点位于左上角，且坐标系是包含0在内的整数索引。所以，如果要找到图像中心的实际像素位置（即，不考虑浮点数部分，
只取最接近的整数像素），需要对结果进行调整。这里 - 1.0f 就是为了做这个调整，它确保了当 arraySizeX 是偶数时，
centerX 指向的是紧邻中心点左侧的像素；如果是奇数，则指向的就是确切的中间像素。这样做的原因是，对于一个大小为 n 的数组，
实际可访问的索引是从 0 到 n-1，因此数组中心点的实际索引应是 (n-1)/2（向下取整）。
总结起来，该行代码用于确定图像在 x 轴方向上的中心像素坐标，并假定图像索引是从0开始计数的。*/



    // For no astigmatism 如果没有 astigmatism
    //这段代码是用来计算图像在不同位置上的焦距值，以应对是否存在像散（astigmatism）的情况。
    // 如果不存在像散（!correctAstigmatism），则所有位置的焦距都使用平均焦距 defocus1；若存在像散，则需要为每个像素位置单独计算焦距
    if (!correctAstigmatism)
    {
        initMatrixWithValue(defocusArray, defocus1); //则所有位置使用平均焦距。
    }
    else	// For astigmatism
    {
        // Precalculate some numbers
        float defSum = defocus1 + defocus2;//两个焦距之和
        float defDiff = defocus1 - defocus2;//两个焦距之差

        // Calculate the astigmatic defocus and store in array
        defocusArray.resize(arraySizeX);    //根据图像大小调整 defocusArray 的尺寸，并逐个计算每个像素点的焦距。

        /*使用两个嵌套循环遍历图像的所有像素坐标(x, y)。
          计算当前像素相对于图像中心(centerX, centerY) 的角度 angle，这里使用了反三角函数 atan2 来确定角度，并转换为度数。
          根据这个角度以及预先计算好的 defSum 和 defDiff、像散度 astigmatism 计算该像素位置对应的焦距值，并存储到 defocusArray[x][y] 中。
        */

        for (size_t x = 0; x < arraySizeX; x++)
        {
            defocusArray[x].resize(arraySizeY);

            for (size_t y = 0; y < arraySizeY; y++)
            {
                float angle = (-atan2((centerY - (y + 1)), ((x + 1) - centerX))) * 180 / M_PI;
                defocusArray[x][y] = (defSum + defDiff * cosd(2.0f * (angle - astigmatism))) / 2.0f;
            }
        }
    }

    // CTF calculation

    // Calculate weighting factors
    float w = sqrt(1.0f - amplitude * amplitude);//振幅

    // Calculate phase array 

    ctfFilter.resize(arraySizeX * arraySizeY);// 用于初始化一个名为ctfFilter的一维浮点数向量。
    //大小调整为arraySizeX和arraySizeY（分别代表图像在x轴和y轴上的像素数量）的乘积，它通过线性索引表示了二维图像的所有像素
    std::fill(ctfFilter.begin(), ctfFilter.end(), 1.0f);

    /*这行代码使用C++标准库中的std::fill函数来填充ctfFilter向量的所有元素。ctfFilter.begin()返回向量的第一个元素的迭代器，
      ctfFilter.end()返回向量最后一个元素之后的位置的迭代器，因此这个区间包含了向量的所有元素。
      std::fill会将指定范围内的所有元素都赋值为给定的值，在这里就是1.0f（单精度浮点数）。
      这样就将ctfFilter中所有的元素初始值设为了1.0f。
      在后续的计算中，这些值会被用来根据对比传递函数（CTF）计算结果更新成相应的滤波器系数。
    */

    ctfAmp.resize(arraySizeX * arraySizeY);//ctfAmp 向量将具有与输入图像相同数量的元素，可以用来存储每个像素位置对应的CTF振幅值。


    //历输入图像在 x 和 y 轴上的所有像素位置。其主要目的是计算每个像素点的对比传递函数（CTF）的振幅值，并将其存储到 ctfAmp 数组中；
    //同时根据计算结果更新滤波器数组 ctfFilter

    for (size_t x = 0; x < arraySizeX; x++)
    {
        for (size_t y = 0; y < arraySizeY; y++)
        {
            float phase = ((M_PI * lambda * (frequencyArray[x][y] * frequencyArray[x][y])) * (defocusArray[x][y] - 0.5 * (lambda * lambda) * (frequencyArray[x][y] * frequencyArray[x][y]) * cs)) + phaseShift;
            // Calculate CTF
            // Has to be transposed here for the case of astigmatism!!! 对于散光的情况，这里必须换位!!
            ctfAmp[y + x * arraySizeY] = (w * sin(phase)) + (amplitude * cos(phase));//存储CTF函数计算后的结果 
            //将计算得出的结果存储到 ctfAmp 数组中。由于这里采用了一维数组表示原本的二维图像数据，因此需要通过转换公式 y + x * arraySizeY 计算一维数组中的线性索引。
            if (ctfAmp[y + x * arraySizeY] <= 0)
                ctfFilter[y + x * arraySizeY] = -1.0f;
            //如果计算得到的振幅值小于等于0，
            // 说明该像素点的 CTF 贡献可能较小或为负数，
            // 因此将对应位置的 ctfFilter 设置为 -1.0f，
            // 这可能是为了在后续处理中忽略这些像素或以某种方式过滤它们。
        }
    }
}


//initMatrixWithValue函数会创建一个与arraySizeX和arraySizeY大小相匹配的二维浮点数数组，
// 并将其中所有的元素均初始化为给定的value。
void CTFCorrection::initMatrixWithValue(std::vector<std::vector<float> >& matrix, float value)
{
    matrix.resize(arraySizeX);

    for (size_t i = 0; i < matrix.size(); i++)
    {
        matrix[i].resize(arraySizeY);
        std::fill(matrix[i].begin(), matrix[i].end(), value);
    }
}


//非正方形图像的投影数据进行填充，以确保图像变为正方形。具体做法是根据较大边的尺寸对较小边进行线性插值填充。
void CTFCorrection::padProjection(std::vector<float>& projection, size_t dimX, size_t dimY)
{
    // For squared images do nothing
    if (dimX == dimY) //如果是正方形不做处理
        return;

    // number of pixels used to compute start and end mean for row/column
    unsigned int meanSize = 10; //定义变量meanSize为用于计算行/列起始和结束平均值的像素数量

    // For non-squared images take the larger dimension 
    // 确定较大的维度padDim、较小维度的起始位置padStart以及分别存储较大和较小维度的变量largeDim和smallDim。
    size_t padDim = max(dimX, dimY);
    size_t padStart = min(dimX, dimY);
    size_t largeDim = max(dimX, dimY);
    size_t smallDim = min(dimX, dimY);

    bool padX; //根据哪个维度较大确定填充方向标志变量padX，如果dimX较大则向x轴填充，否则向y轴填充。

    if (largeDim == dimX)
        padX = false;
    else
        padX = true;

    std::vector<float> padProjection;     //新的浮点数向量padProjection，大小为padDim*padDim，用于存放填充后的投影数据。
    padProjection.resize(padDim * padDim);

    for (size_t i = 0; i < largeDim; i++)  //对于较大的维度largeDim，遍历每个“行”或“列”：
    {
        float startMean = 0;
        float endMean = 0;     //该行或列的起始部分和结束部分的平均值。原始投影数据复制到新数组中对应的位置。

        size_t x = i;

        for (size_t j = 0; j < smallDim; j++) //对较小维度 smallDim 进行遍历
        {
            size_t y = j;
            if (padX)
            {
                x = j;
                y = i;
            }
            //定义索引变量 x 为当前的行/列编号 i。根据填充方向标志 padX，调整 x 和另一个索引变量 y 的值。如果 padX 为真（表示需要在 x 轴方向填充），则交换 x 和 y 的值。

            if (j < meanSize)
                startMean += projection[x + y * dimX]; //计算并累加起始部分像素值到 startMean   
             
            if (j > (smallDim - meanSize))
                endMean += projection[x + y * dimX];

            padProjection[x + y * padDim] = projection[x + y * dimX];
        }

        startMean /= meanSize;                                                                        //这部分不太明白
        endMean /= meanSize;

        unsigned int dist_counter = 0;
        for (size_t j = padStart; j < padDim; j++)
        {
            size_t y = j;
            if (padX)
            {
                x = j;
                y = i;
            }

            float dist = dist_counter / (padDim - padStart - 1);
            padProjection[x + y * padDim] = dist * startMean + (1.0f - dist) * endMean;
            dist_counter++;
        }
    }

    projection.resize(padDim * padDim);

    std::copy(padProjection.begin(), padProjection.end(), projection.begin()); //将填充好的padProjection数据复制回原输入向量projection中，同时调整projection的大小为正方形尺寸。
}


void CTFCorrection::cropProjection(std::vector<float>& projection, size_t dimX, size_t dimY)
{
    if (dimX == dimY)
        return;   //首先检查输入的投影是否已经是正方形（即宽度dimX和高度dimY相等），如果是，则不需要做任何处理直接返回。

    // For non-squared images take the larger dimension
    size_t largeDim = max(dimX, dimY);
    size_t smallDim = min(dimX, dimY);  //计算较大的维度largeDim和较小维度smallDim，分别存储较大和较小的宽高值。

    std::vector<float> croppedProjection;
    croppedProjection.resize(dimX * dimY);  //新的浮点数向量croppedProjection，大小为dimX*dimY，用于存放裁剪后的投影数据。

    bool padX;

    if (largeDim == dimX)       //维度较大确定填充方向标志变量padX，如果dimX较大则不需调整坐标，否则需要交换x和y坐标。
        padX = false;
    else
        padX = true;

    for (size_t i = 0; i < largeDim; i++)  //对于较大的维度largeDim，遍历每个“行”或“列
    {
        size_t x = i;  //定义索引变量 x 为当前的行/列编号 i。

        for (size_t j = 0; j < smallDim; j++)//再次遍历较小维度smallDim
        {
            size_t y = j;// 定义索引变量 y 为当前的行 / 列编号 j
            if (padX)
            {                 //需要根据padX交换坐标，则更新 x 和 y 的值。
                x = j;
                y = i;
            }

            croppedProjection[x + y * dimX] = projection[x + y * largeDim]; //将原始投影数据中位于正确位置的像素复制到新数组croppedProjection对应的位置上。
        }

    }

    projection.resize(dimX * dimY);

    std::copy(croppedProjection.begin(), croppedProjection.end(), projection.begin()); 
    //最后，将裁剪好的croppedProjection数据复制回原输入向量projection中，并调整projection的大小为裁剪后的正方形尺寸。
}