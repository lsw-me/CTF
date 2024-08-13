#include "fftRoutines.h"
#include <fftw3.h>
#include "common.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "volumeIO.h"
#include <opencv2/opencv.hpp>

using namespace std;

void FFTRoutines::real1DTransform(size_t size, std::vector<float>& data, std::vector<float>& fft)
{

	double* fftIn;
	fftw_complex* fftOut;
	fftw_plan plan_forward;

	fftIn = (double*)fftw_malloc(sizeof(double) * size);

	for (unsigned int i = 0; i < size; i++)
	{
		fftIn[i] = (double)data[i];
	}

	size_t fftSize = (size / 2) + 1;

	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);

	plan_forward = fftw_plan_dft_r2c_1d(size, fftIn, fftOut, FFTW_ESTIMATE);

	fftw_execute(plan_forward);

	if (size % 2 == 0) /* N is even */
		fft.push_back(fftOut[size / 2][0]);

	for (int k = (size + 1) / 2 - 1; k >= 1; k--)  // (k < N/2 rounded up)
		fft.push_back(fftOut[k][0]);

	for (int k = 1; k < (int)(size + 1) / 2; k++)  // (k < N/2 rounded up)
		fft.push_back(fftOut[k][0]);
	if (size % 2 == 0) // N is even
		fft.push_back(fftOut[size / 2][0]);  // Nyquist freq.

	fftw_destroy_plan(plan_forward);

	fftw_free(fftOut);
	fftw_free(fftIn);

	return;

}

double FFTRoutines::logarithmizeValue(double value, bool logarithmizeData)
{
	if (logarithmizeData && value >= 1.0)
		return log(value);
	else
		return value;
}


/*经过傅里叶变换之后的频域数据进行逆fftshift操作，即将频谱从中心对齐的布局恢复到边缘对齐的布局。
这种操作常用于图像处理和信号处理领域，在进行傅里叶变换前后保持数据在空间域中的排列一致性
接受四个参数   输出数组，用于存储逆fftshift后的结果。  输入数组，存储了已经中心对齐的频域数据。 

图像的宽度（或数组在X轴方向的长度  。 图像的高度（或数组在Y轴方向的长度）。


fftshift操作保证了逆傅里叶变换后生成的空间域数据与原始空间域数据在排列上完全一致
*/
void inversefftshift(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)  //反向 FFT 偏移
{
	for (size_t i = 0; i < ydim / 2; i++)  //第一部分频谱左上四分之一部分（索引1->4区域），通过计算新的输出索引(outIndex1和outIndex2)和输入索引(inIndex1和inIndex2)，
		                                   //将输入数组in中的数据复制到输出数组out中，实现频谱的逆fftshift操作。
	{
		for (size_t j = 0; j < xdim / 2; j++)
		{
			size_t outIndex1 = (j + xdim / 2) + (i + ydim / 2) * xdim;
			size_t inIndex1 = j + i * xdim;
			size_t outIndex2 = j + i * xdim;
			size_t inIndex2 = (j + xdim / 2) + (i + ydim / 2) * xdim;
			out[inIndex1] = in[outIndex1];	//1->4
			out[inIndex2] = in[outIndex2];	//4->1
		}
	}

	for (size_t i = 0; i < ydim / 2; i++)  //对频谱右上四分之一部分（索引2->3区域），执行类似的逆fftshift操作。 
	{
		for (size_t j = xdim / 2; j < xdim; j++)
		{
			size_t outIndex1 = (j - xdim / 2) + (i + ydim / 2) * xdim;
			size_t inIndex1 = j + i * xdim;
			size_t outIndex2 = j + i * xdim;
			size_t inIndex2 = (j - xdim / 2) + (i + ydim / 2) * xdim;
			out[inIndex1] = in[outIndex1];	//2->3
			out[inIndex2] = in[outIndex2];	//3->2
		}
	}
}


void fftshift(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)
{
	for (size_t i = 0; i < ydim / 2; i++)
	{
		for (size_t j = 0; j < xdim / 2; j++)
		{
			size_t outIndex1 = (j + xdim / 2) + (i + ydim / 2) * xdim;
			size_t inIndex1 = j + i * xdim;
			size_t outIndex2 = j + i * xdim;
			size_t inIndex2 = (j + xdim / 2) + (i + ydim / 2) * xdim;
			out[outIndex1] = in[inIndex1];	//1->4
			out[outIndex2] = in[inIndex2];	//4->1
		}
	}

	for (size_t i = 0; i < ydim / 2; i++)
	{
		for (size_t j = xdim / 2; j < xdim; j++)
		{
			size_t outIndex1 = (j - xdim / 2) + (i + ydim / 2) * xdim;
			size_t inIndex1 = j + i * xdim;
			size_t outIndex2 = j + i * xdim;
			size_t inIndex2 = (j - xdim / 2) + (i + ydim / 2) * xdim;
			out[outIndex1] = in[inIndex1];	//2->3
			out[outIndex2] = in[inIndex2];	//3->2
		}
	}
}


//和写的那个displayFFTPowerSpectrum 功能类似，计算功率谱的方法也差不多
void FFTRoutines::computePowerSpectrum(std::vector<float>& powerSpectrum, fftw_complex* fftOut, size_t sizeX, size_t nyh, bool logarithmizeData)
{
	size_t k = 0;
	for (size_t j = 0; j < nyh; j++)
	{
		for (size_t i = 0; i < sizeX; i++)
		{
			powerSpectrum[k] = logarithmizeValue(fftOut[j + i * nyh][0] * fftOut[j + i * nyh][0] + fftOut[j + i * nyh][1] * fftOut[j + i * nyh][1], logarithmizeData);
			k++;
		}
	}

	for (int j = nyh - 2; j > 0; j--)
	{
		powerSpectrum[k] = fftOut[j][0] * fftOut[j][0] + fftOut[j][1] * fftOut[j][1];
		k++;
		for (int i = sizeX - 1; i > 0; i--)
		{
			powerSpectrum[k] = logarithmizeValue(fftOut[j + i * nyh][0] * fftOut[j + i * nyh][0] + fftOut[j + i * nyh][1] * fftOut[j + i * nyh][1], logarithmizeData);
			k++;
		}
	}
}  

void FFTRoutines::maskFFT(fftw_complex* fftOut, std::vector<double>& mask, size_t sizeX, size_t nyh)
{
	size_t k = 0;
	for (size_t j = 0; j < nyh; j++)
	{
		for (size_t i = 0; i < sizeX; i++)
		{
			if (mask[k] != 1.0f)
			{
				fftOut[j + i * nyh][0] = fftOut[j + i * nyh][0] * mask[k];
				fftOut[j + i * nyh][1] = fftOut[j + i * nyh][1] * mask[k];
			}
			k++;
		}
	}
}


//二维FFT变换结果进行滤波操作。函数接收四个参数： 
//fftw_complex* fftOut          一个指向fftw_complex类型数组的指针，它存储了经过二维快速傅里叶变换（FFT）后得到的数据，其中每个元素代表一个复数，通常包含实部和虚部。
//std::vector<float>& filter   这是一个引用类型的浮点数向量，表示滤波器系数序列，用于调整变换后频谱的各个频率分量。
//size_t sizeX                 输入信号在x方向上的维度大小，即FFT后数据矩阵的列数  
//size_t nyh                   输入信号在y方向上的一半加一，因为在二维FFT的结果中，由于对称性，下半部分可以通过上半部分来获得，所以通常只需要处理上半部分频谱。
void FFTRoutines::filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh)
{
	size_t k = 0;
	for (size_t j = 0; j < nyh; j++)                          //遍历整个二维频谱数据，如果filter 值不为1，系数分别乘以原频谱对应位置的实部和虚部 
	{
		for (size_t i = 0; i < sizeX; i++)
		{
			if (filter[k] != 1.0f)
			{
				fftOut[j + i * nyh][0] = fftOut[j + i * nyh][0] * filter[k];
				fftOut[j + i * nyh][1] = fftOut[j + i * nyh][1] * filter[k];
			}
			k++;
		}
	}
}

void FFTRoutines::normalizeValues(std::vector<float>& normalizedData, std::vector<double>& originalData, size_t dataSize, novaCTF::DataStats& dataStats)
{
	for (size_t i = 0; i < dataSize; i++)
	{
		normalizedData[i] = originalData[i] / (dataSize);

		dataStats.mean += normalizedData[i];
		dataStats.max = max(dataStats.max, normalizedData[i]);
		dataStats.min = min(dataStats.min, normalizedData[i]);
	}
}

void FFTRoutines::complex2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter)
{

	//对给定的一维实数数组（表示为二维图像数据）进行2D复数离散傅里叶变换(DFT)，应用一个过滤器（filter），然后再进行逆变换，将结果归一化后写回原始数据数组
	fftw_complex* fftOut;
	fftw_complex* fftIn;                //为输入和输出数据分配内存，将它们转换为fftw_complex类型的指针fftOut和fftIn。
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);
	fftIn = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);

	for (size_t i = 0; i < sizeX * sizeY; i++) //实数数组data的内容复制到fftIn数组的实部，虚部初始化为0.0，将实数数据转换为复数形式
	{
		fftIn[i][0] = (double)data[i];
		fftIn[i][1] = 0.0;
	}

	fftw_plan plan_forward = fftw_plan_dft_2d(sizeX, sizeY, fftIn, fftOut, -1, FFTW_ESTIMATE); //创建一个2D DFT计划（plan_forward），使用fftw_plan_dft_2d函数，执行方向为-1（表示正变换），并执行该计划
	fftw_execute(plan_forward);

	/*输出fftout的值
	for (size_t i = 0; i < sizeX; ++i)
	{
		for (size_t j = 0; j < sizeY; ++j)
		{
			size_t index = i * sizeX + j;
			cout<<"output["<<j<<","<< i<<"] =("<<fftOut[index][0]<<","<< fftOut[index][1] << ")" << std::endl;
		}
	}
	*/
	std::vector<float> isFilter;
	isFilter.resize(filter.size());

	inversefftshift(isFilter, filter, sizeX, sizeY);
	filterFFT(fftOut, isFilter, sizeX, sizeY);


	// displayFFTPowerSpectrum(fftOut, sizeX, sizeY);
	// 我想在这里查看fftout的值，最好是能够输出图片

	fftw_plan plan_backward = fftw_plan_dft_2d(sizeX, sizeY, fftOut, fftIn, 1, FFTW_ESTIMATE);//创建一个2D IDFT计划（plan_backward），执行方向为1（表示逆变换），并执行该计划
	fftw_execute(plan_backward);
	//这里看fftout的值
	fftw_destroy_plan(plan_backward);
	fftw_destroy_plan(plan_forward);
	fftw_free(fftOut);

	for (size_t i = 0; i < data.size(); i++)
	{
		data[i] = fftIn[i][0] / data.size();
	}

	fftw_free(fftIn);

}

void FFTRoutines::real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter)
{
	std::vector<double> fftIn;
	size_t nyh;
	fftw_complex* fftOut;
	fftw_plan plan_backward;
	fftw_plan plan_forward;

	size_t sliceSize = sizeX * sizeY;
	fftIn.resize(sizeX * sizeY);

	for (size_t i = 0; i < sliceSize; i++)
	{
		fftIn[i] = (double)data[i];
	}

	nyh = (sizeY / 2) + 1;
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * nyh);

	plan_forward = fftw_plan_dft_r2c_2d(sizeX, sizeY, &fftIn[0], fftOut, FFTW_ESTIMATE);
	fftw_execute(plan_forward);


	std::vector<float> isFilter;
	isFilter.resize(filter.size());


	inversefftshift(isFilter, filter, sizeX, sizeY);
	filterFFT(fftOut, isFilter, sizeX, nyh);

	plan_backward = fftw_plan_dft_c2r_2d(sizeX, sizeY, fftOut, &fftIn[0], FFTW_ESTIMATE);
	fftw_execute(plan_backward);

	fftw_destroy_plan(plan_backward);
	fftw_destroy_plan(plan_forward);
	fftw_free(fftOut);

	for (size_t i = 0; i < data.size(); i++)
	{
		data[i] = fftIn[i] / data.size();
	}
}

void FFTRoutines::real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<double>& mask, novaCTF::DataStats& dataStats)
{
	std::vector<double> fftIn;
	size_t nyh;
	fftw_complex* fftOut;
	fftw_plan plan_backward;
	fftw_plan plan_forward;

	size_t sliceSize = sizeX * sizeY;
	fftIn.resize(sizeX * sizeY);

	for (size_t i = 0; i < sliceSize; i++)
	{
		fftIn[i] = (double)data[i];
	}

	nyh = (sizeY / 2) + 1;
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * nyh);

	plan_forward = fftw_plan_dft_r2c_2d(sizeX, sizeY, &fftIn[0], fftOut, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	maskFFT(fftOut, mask, sizeX, nyh);

	plan_backward = fftw_plan_dft_c2r_2d(sizeX, sizeY, fftOut, &fftIn[0], FFTW_ESTIMATE);
	fftw_execute(plan_backward);

	normalizeValues(data, fftIn, sliceSize, dataStats);

	fftw_destroy_plan(plan_backward);
	fftw_destroy_plan(plan_forward);
	fftw_free(fftOut);

}

void FFTRoutines::many1DTransform(float* data, int sizeX, int sizeY, int direction)
{
	size_t nx = sizeX;
	size_t ny = sizeY;

	int nxpad = 2 * (nx / 2 + 1);
	float invScale = (float)(1.0f / sqrt((double)nx));

	fftwf_plan plan = nullptr;
	if (direction == 0)
	{
		plan = fftwf_plan_many_dft_r2c(1, &sizeX, ny, data, NULL, 1, nxpad, (fftwf_complex*)data, NULL, 1, nxpad / 2, FFTW_ESTIMATE);
	}
	else if (direction == 1)
	{
		plan = fftwf_plan_many_dft_c2r(1, &sizeX, ny, (fftwf_complex*)data, NULL, 1, nxpad / 2, data, NULL, 1, nxpad, FFTW_ESTIMATE);
	}

	fftwf_execute(plan);

	normalize(data, invScale, (size_t)nxpad * ny);

	fftwf_destroy_plan(plan);
}


/* Normalize the given number of real elements by the scaling factor */
void FFTRoutines::normalize(float* array, float scale, size_t dataSize)
{
	size_t i;
	for (i = 0; i < dataSize; i++)
		array[i] *= scale;
}

void FFTRoutines::displayFFTPowerSpectrum(const fftw_complex* fftOut, size_t sizeX, size_t sizeY)
{
	// 计算功率谱密度
	cv::Mat powerSpectrum(sizeY, sizeX, CV_32F);
	for (size_t i = 0; i < sizeX * sizeY; ++i) {
		float real = fftOut[i][0];
		float imag = fftOut[i][1];
		double magnitudeSquared = real * real + imag * imag; // 幅度的平方即为功率谱密度
		powerSpectrum.at<float>(i / sizeX, i % sizeX) = static_cast<float>(magnitudeSquared);
		// static_cast<float>(magnitudeSquared) 类型转换操作，使用 static_cast 关键字将 magnitudeSquared 从其原始类型（转换为 float 类型。这是因为 powerSpectrum 矩阵存储的是 float 类型的元素，所以需要将计算得到的 double 类型值转换为匹配的类型。
		//.at<float> 是cv::Mat 类提供的成员函数，用于访问矩阵中指定位置的元素  (i / sizeX, i % sizeX) 用来访问矩阵中的指定元素
		//powerSpectrum 是之前定义的 cv::Mat 类型变量，代表一个二维浮点数矩阵，用于存储计算出的功率谱密度。
	}

	// 对功率谱密度进行对数尺度处理，使其更适合人眼观察
	cv::log(powerSpectrum, powerSpectrum); // 取自然对数
	cv::normalize(powerSpectrum, powerSpectrum, 0, 255, cv::NORM_MINMAX); // 归一化到0-255范围内

	// 转换为8位灰度图像并进行fftshift操作
	cv::Mat logPowerSpectrum8U;
	powerSpectrum.convertTo(logPowerSpectrum8U, CV_8U);
	cv::Mat logPowerSpectrumShifted;
	cv::flip(logPowerSpectrum8U, logPowerSpectrumShifted, -1);

	//调整大小
	cv::Mat resizedImage;
	cv::resize(logPowerSpectrumShifted, resizedImage, cv::Size(400, 400));
	// 显示图像
	cv::imshow("2D FFT Power Spectrum Density", resizedImage);
	cv::waitKey(0);
}

