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

void CTFCorrection::computeFrequencyArray()//��������ͼ���Ƶ������   ����������ø���
{
    size_t resX = inputStack->getResolution().x;//��ȡ����ͼ���x yֵ
    size_t resY = inputStack->getResolution().y;//

    // For non-squared images take the larger dimension
    if (resX != resY)
    {
        resX = max(resX, resY);
        resY = max(resX, resY);
    }                                     //���x y ��ͬ�������ֵΪx y  

    std::vector<float> xArray;    
    std::vector<float> yArray;//����������������xArray��yArray���ֱ����ڴ洢x���y�᷽���ϵ�Ƶ��ֵ

    for (int i = -floor(resX / 2); i < floor(resX / 2) + (resX % 2); i++)
        xArray.push_back(i / (resX * pixelSize));              //ʹ������ѭ������x�᷽���Ƶ��ֵ��  floor ���ڶԸ�������������ȡ��������
                                                             // ѭ����Χ�� - floor(resX / 2) �� floor(resX / 2) + (resX % 2)��ȷ����������ͼ����������򣬰��������ߴ�ʱ���������ء�
                                                             // ����ÿ��iֵ����������(resX * pixelSize)�õ�Ƶ��ֵ�����������ӵ�xArray�С�
                                                             //���ķ�ʽ����y�᷽���Ƶ��ֵ��
    for (int i = -floor(resY / 2); i < floor(resY / 2) + (resY % 2); i++)
        yArray.push_back(i / (resY * pixelSize));

    frequencyArray.resize(xArray.size());          //��ʼ��һ����ά��̬����frequencyArray�����СΪxArray�ĳ��ȡ�

    for (size_t x = 0; x < xArray.size(); x++)
        for (size_t y = 0; y < yArray.size(); y++)
            frequencyArray[x].push_back(sqrt(xArray[x] * xArray[x] + yArray[y] * yArray[y])); 
    //�ٴ�ʹ��Ƕ��ѭ������xArray��yArray������ÿ��������Ӧ��Ƶ��ֵ������ŷ����þ��빫ʽ�������������ӵ�frequencyArray[x]�С�
    arraySizeX = xArray.size();
    arraySizeY = yArray.size(); //����arraySizeX��arraySizeY�ֱ�ΪxArray��yArray�ĳ��ȡ�
}

/*If astigmatism correction is required it performs check whether the values are available.
* If they are not, than the warning is written out and the correction is performed without
* the astigmatism correction.
* �����Ҫɢ��У������������Щֵ�Ƿ���á�������ǣ���д�����沢��û�е������ִ�и�����ɢ��У��
*/
void CTFCorrection::checkAstigmatism() //���Ҳ���ø�
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
    std::vector<std::vector<float>> correctedProjections;//�洢����У�����ͶӰ���ݣ�����Ϊ��ά������������
    std::vector<float> ctfAmp; //���ڴ洢������ĶԱȴ��ݺ�����Contrast Transfer Function, CTF����������֡�
    std::vector<float> ctfFilter;//���ڴ洢������ĶԱȴ��ݺ������˲��������ܴ�����λ���������ϵ����

    correctedProjections.resize(inputStack->getNumberOfProjections());//��������Stack��ͶӰ����������correctedProjections�Ĵ�С��

    for (ProjectionSet::iterator it = projSet->begin(); it != projSet->end(); it++) //ʹ�õ�����it����projSet�е�����ͶӰ������ȡÿ��ͶӰ�ڼ����е�����projIndex
    {
        unsigned int projIndex = it.second();
        correctedProjections[projIndex].resize(inputStack->getProjectionSize());
        inputStack->readProjections(&correctedProjections[projIndex][0], 1, projIndex);
        //������Stack�ж�ȡͶӰ���ݵ�correctedProjections[projIndex]�У���������ڴ�ռ䡣
        padProjection(correctedProjections[projIndex], inputStack->getResolution().x, inputStack->getResolution().y);
        //����padProjection��������չ�����ͶӰͼ������Ӧ�ض��ֱ���Ҫ��
        computeCTF(ctfAmp, ctfFilter, defocusFileValues[projIndex][0], defocusFileValues[projIndex][1], defocusFileValues[projIndex][2], defocusFileValues[projIndex][3]);
        //ʹ���ṩ�Ľ���������ĸ�ֵ�ֱ��Ӧ��ͬ������������computeCTF�������㵱ǰͶӰ��Ӧ��CTF��Ϣ��
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
        //����cropProjection����ȥ��֮ǰ���Ķ����������򣬻ָ���ԭʼ�ֱ��ʴ�С��
    }

    VolumeIO::writeMRCStack(inputStack->getStackHeader(), correctedProjections, params.OutputFilename(), inputStack->getExtraData());

}

void CTFCorrection::computeCTF(std::vector<float>& ctfAmp, std::vector<float>& ctfFilter, float defocus1, float defocus2, float astigmatism, float phaseShift)
{
    // Convert defocii from microns to Angstroms ���������������ֵ��΢�ף�microns��ת��Ϊ����Angstroms����
    defocus1 = defocus1 * (1.0e4);
    defocus2 = defocus2 * (1.0e4);


    // Calculate electron wavelength ʹ�����ʿ˳��������١����Ӿ�ֹ���������ٵ�ѹ�Լ����ӵ�ɵ������������������ָ�����ٵ�ѹ�µĲ�����
    // Most of the equations here are from Mindell and Grigorieff (2003) 

    double h = 6.62606957e-34;
    double c = 299792458;
    float eRest = 511000;
    float v = evk * 1000;
    float eCharge = 1.602e-19;
    double lambda = (c * h) / sqrt(((2 * eRest * v) + (v * v)) * (eCharge * eCharge)) * (pow(10.0, 10));

    // Calculate astigmatic defocus array ��ʼ��������������������飺

    // Initialize defocus array
    std::vector<std::vector<float> > defocusArray; //��ά���� ɢ������

    // Calculate center of image  //��������
    float centerX = arraySizeX / 2.0f - 1.0f;
    float centerY = arraySizeY / 2.0f - 1.0f;
    
    /*arraySizeX ��ʾͼ���� x �᷽���ϵ������������� arraySizeX ���� 2.0f ��Ϊ�˵õ�ͼ���ȵ�һ�룬��ͨ����õ�ͼ������ĵ����ڵ�������������ͼ��������0��ʼ����
Ȼ����������̻�����ͼ������У�ͼ�������ԭ��λ�����Ͻǣ�������ϵ�ǰ���0���ڵ��������������ԣ����Ҫ�ҵ�ͼ�����ĵ�ʵ������λ�ã����������Ǹ��������֣�
ֻȡ��ӽ����������أ�����Ҫ�Խ�����е��������� - 1.0f ����Ϊ���������������ȷ���˵� arraySizeX ��ż��ʱ��
centerX ָ����ǽ������ĵ��������أ��������������ָ��ľ���ȷ�е��м����ء���������ԭ���ǣ�����һ����СΪ n �����飬
ʵ�ʿɷ��ʵ������Ǵ� 0 �� n-1������������ĵ��ʵ������Ӧ�� (n-1)/2������ȡ������
�ܽ����������д�������ȷ��ͼ���� x �᷽���ϵ������������꣬���ٶ�ͼ�������Ǵ�0��ʼ�����ġ�*/



    // For no astigmatism ���û�� astigmatism
    //��δ�������������ͼ���ڲ�ͬλ���ϵĽ���ֵ����Ӧ���Ƿ������ɢ��astigmatism���������
    // �����������ɢ��!correctAstigmatism����������λ�õĽ��඼ʹ��ƽ������ defocus1����������ɢ������ҪΪÿ������λ�õ������㽹��
    if (!correctAstigmatism)
    {
        initMatrixWithValue(defocusArray, defocus1); //������λ��ʹ��ƽ�����ࡣ
    }
    else	// For astigmatism
    {
        // Precalculate some numbers
        float defSum = defocus1 + defocus2;//��������֮��
        float defDiff = defocus1 - defocus2;//��������֮��

        // Calculate the astigmatic defocus and store in array
        defocusArray.resize(arraySizeX);    //����ͼ���С���� defocusArray �ĳߴ磬���������ÿ�����ص�Ľ��ࡣ

        /*ʹ������Ƕ��ѭ������ͼ���������������(x, y)��
          ���㵱ǰ���������ͼ������(centerX, centerY) �ĽǶ� angle������ʹ���˷����Ǻ��� atan2 ��ȷ���Ƕȣ���ת��Ϊ������
          ��������Ƕ��Լ�Ԥ�ȼ���õ� defSum �� defDiff����ɢ�� astigmatism ���������λ�ö�Ӧ�Ľ���ֵ�����洢�� defocusArray[x][y] �С�
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
    float w = sqrt(1.0f - amplitude * amplitude);//���

    // Calculate phase array 

    ctfFilter.resize(arraySizeX * arraySizeY);// ���ڳ�ʼ��һ����ΪctfFilter��һά������������
    //��С����ΪarraySizeX��arraySizeY���ֱ����ͼ����x���y���ϵ������������ĳ˻�����ͨ������������ʾ�˶�άͼ�����������
    std::fill(ctfFilter.begin(), ctfFilter.end(), 1.0f);

    /*���д���ʹ��C++��׼���е�std::fill���������ctfFilter����������Ԫ�ء�ctfFilter.begin()���������ĵ�һ��Ԫ�صĵ�������
      ctfFilter.end()�����������һ��Ԫ��֮���λ�õĵ��������������������������������Ԫ�ء�
      std::fill�Ὣָ����Χ�ڵ�����Ԫ�ض���ֵΪ������ֵ�����������1.0f�������ȸ���������
      �����ͽ�ctfFilter�����е�Ԫ�س�ʼֵ��Ϊ��1.0f��
      �ں����ļ����У���Щֵ�ᱻ�������ݶԱȴ��ݺ�����CTF�����������³���Ӧ���˲���ϵ����
    */

    ctfAmp.resize(arraySizeX * arraySizeY);//ctfAmp ����������������ͼ����ͬ������Ԫ�أ����������洢ÿ������λ�ö�Ӧ��CTF���ֵ��


    //������ͼ���� x �� y ���ϵ���������λ�á�����ҪĿ���Ǽ���ÿ�����ص�ĶԱȴ��ݺ�����CTF�������ֵ��������洢�� ctfAmp �����У�
    //ͬʱ���ݼ����������˲������� ctfFilter

    for (size_t x = 0; x < arraySizeX; x++)
    {
        for (size_t y = 0; y < arraySizeY; y++)
        {
            float phase = ((M_PI * lambda * (frequencyArray[x][y] * frequencyArray[x][y])) * (defocusArray[x][y] - 0.5 * (lambda * lambda) * (frequencyArray[x][y] * frequencyArray[x][y]) * cs)) + phaseShift;
            // Calculate CTF
            // Has to be transposed here for the case of astigmatism!!! ����ɢ��������������뻻λ!!
            ctfAmp[y + x * arraySizeY] = (w * sin(phase)) + (amplitude * cos(phase));//�洢CTF���������Ľ�� 
            //������ó��Ľ���洢�� ctfAmp �����С��������������һά�����ʾԭ���Ķ�άͼ�����ݣ������Ҫͨ��ת����ʽ y + x * arraySizeY ����һά�����е�����������
            if (ctfAmp[y + x * arraySizeY] <= 0)
                ctfFilter[y + x * arraySizeY] = -1.0f;
            //�������õ������ֵС�ڵ���0��
            // ˵�������ص�� CTF ���׿��ܽ�С��Ϊ������
            // ��˽���Ӧλ�õ� ctfFilter ����Ϊ -1.0f��
            // �������Ϊ���ں��������к�����Щ���ػ���ĳ�ַ�ʽ�������ǡ�
        }
    }
}


//initMatrixWithValue�����ᴴ��һ����arraySizeX��arraySizeY��С��ƥ��Ķ�ά���������飬
// �����������е�Ԫ�ؾ���ʼ��Ϊ������value��
void CTFCorrection::initMatrixWithValue(std::vector<std::vector<float> >& matrix, float value)
{
    matrix.resize(arraySizeX);

    for (size_t i = 0; i < matrix.size(); i++)
    {
        matrix[i].resize(arraySizeY);
        std::fill(matrix[i].begin(), matrix[i].end(), value);
    }
}


//��������ͼ���ͶӰ���ݽ�����䣬��ȷ��ͼ���Ϊ�����Ρ����������Ǹ��ݽϴ�ߵĳߴ�Խ�С�߽������Բ�ֵ��䡣
void CTFCorrection::padProjection(std::vector<float>& projection, size_t dimX, size_t dimY)
{
    // For squared images do nothing
    if (dimX == dimY) //����������β�������
        return;

    // number of pixels used to compute start and end mean for row/column
    unsigned int meanSize = 10; //�������meanSizeΪ���ڼ�����/����ʼ�ͽ���ƽ��ֵ����������

    // For non-squared images take the larger dimension 
    // ȷ���ϴ��ά��padDim����Сά�ȵ���ʼλ��padStart�Լ��ֱ�洢�ϴ�ͽ�Сά�ȵı���largeDim��smallDim��
    size_t padDim = max(dimX, dimY);
    size_t padStart = min(dimX, dimY);
    size_t largeDim = max(dimX, dimY);
    size_t smallDim = min(dimX, dimY);

    bool padX; //�����ĸ�ά�Ƚϴ�ȷ����䷽���־����padX�����dimX�ϴ�����x����䣬������y����䡣

    if (largeDim == dimX)
        padX = false;
    else
        padX = true;

    std::vector<float> padProjection;     //�µĸ���������padProjection����СΪpadDim*padDim�����ڴ�������ͶӰ���ݡ�
    padProjection.resize(padDim * padDim);

    for (size_t i = 0; i < largeDim; i++)  //���ڽϴ��ά��largeDim������ÿ�����С����С���
    {
        float startMean = 0;
        float endMean = 0;     //���л��е���ʼ���ֺͽ������ֵ�ƽ��ֵ��ԭʼͶӰ���ݸ��Ƶ��������ж�Ӧ��λ�á�

        size_t x = i;

        for (size_t j = 0; j < smallDim; j++) //�Խ�Сά�� smallDim ���б���
        {
            size_t y = j;
            if (padX)
            {
                x = j;
                y = i;
            }
            //������������ x Ϊ��ǰ����/�б�� i��������䷽���־ padX������ x ����һ���������� y ��ֵ����� padX Ϊ�棨��ʾ��Ҫ�� x �᷽����䣩���򽻻� x �� y ��ֵ��

            if (j < meanSize)
                startMean += projection[x + y * dimX]; //���㲢�ۼ���ʼ��������ֵ�� startMean   
             
            if (j > (smallDim - meanSize))
                endMean += projection[x + y * dimX];

            padProjection[x + y * padDim] = projection[x + y * dimX];
        }

        startMean /= meanSize;                                                                        //�ⲿ�ֲ�̫����
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

    std::copy(padProjection.begin(), padProjection.end(), projection.begin()); //�����õ�padProjection���ݸ��ƻ�ԭ��������projection�У�ͬʱ����projection�Ĵ�СΪ�����γߴ硣
}


void CTFCorrection::cropProjection(std::vector<float>& projection, size_t dimX, size_t dimY)
{
    if (dimX == dimY)
        return;   //���ȼ�������ͶӰ�Ƿ��Ѿ��������Σ������dimX�͸߶�dimY��ȣ�������ǣ�����Ҫ���κδ���ֱ�ӷ��ء�

    // For non-squared images take the larger dimension
    size_t largeDim = max(dimX, dimY);
    size_t smallDim = min(dimX, dimY);  //����ϴ��ά��largeDim�ͽ�Сά��smallDim���ֱ�洢�ϴ�ͽ�С�Ŀ��ֵ��

    std::vector<float> croppedProjection;
    croppedProjection.resize(dimX * dimY);  //�µĸ���������croppedProjection����СΪdimX*dimY�����ڴ�Ųü����ͶӰ���ݡ�

    bool padX;

    if (largeDim == dimX)       //ά�Ƚϴ�ȷ����䷽���־����padX�����dimX�ϴ�����������꣬������Ҫ����x��y���ꡣ
        padX = false;
    else
        padX = true;

    for (size_t i = 0; i < largeDim; i++)  //���ڽϴ��ά��largeDim������ÿ�����С�����
    {
        size_t x = i;  //������������ x Ϊ��ǰ����/�б�� i��

        for (size_t j = 0; j < smallDim; j++)//�ٴα�����Сά��smallDim
        {
            size_t y = j;// ������������ y Ϊ��ǰ���� / �б�� j
            if (padX)
            {                 //��Ҫ����padX�������꣬����� x �� y ��ֵ��
                x = j;
                y = i;
            }

            croppedProjection[x + y * dimX] = projection[x + y * largeDim]; //��ԭʼͶӰ������λ����ȷλ�õ����ظ��Ƶ�������croppedProjection��Ӧ��λ���ϡ�
        }

    }

    projection.resize(dimX * dimY);

    std::copy(croppedProjection.begin(), croppedProjection.end(), projection.begin()); 
    //��󣬽��ü��õ�croppedProjection���ݸ��ƻ�ԭ��������projection�У�������projection�Ĵ�СΪ�ü���������γߴ硣
}