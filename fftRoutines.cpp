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


/*��������Ҷ�任֮���Ƶ�����ݽ�����fftshift����������Ƶ�״����Ķ���Ĳ��ָֻ�����Ե����Ĳ��֡�
���ֲ���������ͼ������źŴ��������ڽ��и���Ҷ�任ǰ�󱣳������ڿռ����е�����һ����
�����ĸ�����   ������飬���ڴ洢��fftshift��Ľ����  �������飬�洢���Ѿ����Ķ����Ƶ�����ݡ� 

ͼ��Ŀ�ȣ���������X�᷽��ĳ���  �� ͼ��ĸ߶ȣ���������Y�᷽��ĳ��ȣ���


fftshift������֤���渵��Ҷ�任�����ɵĿռ���������ԭʼ�ռ�����������������ȫһ��
*/
void inversefftshift(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)  //���� FFT ƫ��
{
	for (size_t i = 0; i < ydim / 2; i++)  //��һ����Ƶ�������ķ�֮һ���֣�����1->4���򣩣�ͨ�������µ��������(outIndex1��outIndex2)����������(inIndex1��inIndex2)��
		                                   //����������in�е����ݸ��Ƶ��������out�У�ʵ��Ƶ�׵���fftshift������
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

	for (size_t i = 0; i < ydim / 2; i++)  //��Ƶ�������ķ�֮һ���֣�����2->3���򣩣�ִ�����Ƶ���fftshift������ 
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


//��д���Ǹ�displayFFTPowerSpectrum �������ƣ����㹦���׵ķ���Ҳ���
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


//��άFFT�任��������˲����������������ĸ������� 
//fftw_complex* fftOut          һ��ָ��fftw_complex���������ָ�룬���洢�˾�����ά���ٸ���Ҷ�任��FFT����õ������ݣ�����ÿ��Ԫ�ش���һ��������ͨ������ʵ�����鲿��
//std::vector<float>& filter   ����һ���������͵ĸ�������������ʾ�˲���ϵ�����У����ڵ����任��Ƶ�׵ĸ���Ƶ�ʷ�����
//size_t sizeX                 �����ź���x�����ϵ�ά�ȴ�С����FFT�����ݾ��������  
//size_t nyh                   �����ź���y�����ϵ�һ���һ����Ϊ�ڶ�άFFT�Ľ���У����ڶԳ��ԣ��°벿�ֿ���ͨ���ϰ벿������ã�����ͨ��ֻ��Ҫ�����ϰ벿��Ƶ�ס�
void FFTRoutines::filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh)
{
	size_t k = 0;
	for (size_t j = 0; j < nyh; j++)                          //����������άƵ�����ݣ����filter ֵ��Ϊ1��ϵ���ֱ����ԭƵ�׶�Ӧλ�õ�ʵ�����鲿 
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

	//�Ը�����һάʵ�����飨��ʾΪ��άͼ�����ݣ�����2D������ɢ����Ҷ�任(DFT)��Ӧ��һ����������filter����Ȼ���ٽ�����任���������һ����д��ԭʼ��������
	fftw_complex* fftOut;
	fftw_complex* fftIn;                //Ϊ�����������ݷ����ڴ棬������ת��Ϊfftw_complex���͵�ָ��fftOut��fftIn��
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);
	fftIn = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY);

	for (size_t i = 0; i < sizeX * sizeY; i++) //ʵ������data�����ݸ��Ƶ�fftIn�����ʵ�����鲿��ʼ��Ϊ0.0����ʵ������ת��Ϊ������ʽ
	{
		fftIn[i][0] = (double)data[i];
		fftIn[i][1] = 0.0;
	}

	fftw_plan plan_forward = fftw_plan_dft_2d(sizeX, sizeY, fftIn, fftOut, -1, FFTW_ESTIMATE); //����һ��2D DFT�ƻ���plan_forward����ʹ��fftw_plan_dft_2d������ִ�з���Ϊ-1����ʾ���任������ִ�иüƻ�
	fftw_execute(plan_forward);

	/*���fftout��ֵ
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
	// ����������鿴fftout��ֵ��������ܹ����ͼƬ

	fftw_plan plan_backward = fftw_plan_dft_2d(sizeX, sizeY, fftOut, fftIn, 1, FFTW_ESTIMATE);//����һ��2D IDFT�ƻ���plan_backward����ִ�з���Ϊ1����ʾ��任������ִ�иüƻ�
	fftw_execute(plan_backward);
	//���￴fftout��ֵ
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
	// ���㹦�����ܶ�
	cv::Mat powerSpectrum(sizeY, sizeX, CV_32F);
	for (size_t i = 0; i < sizeX * sizeY; ++i) {
		float real = fftOut[i][0];
		float imag = fftOut[i][1];
		double magnitudeSquared = real * real + imag * imag; // ���ȵ�ƽ����Ϊ�������ܶ�
		powerSpectrum.at<float>(i / sizeX, i % sizeX) = static_cast<float>(magnitudeSquared);
		// static_cast<float>(magnitudeSquared) ����ת��������ʹ�� static_cast �ؼ��ֽ� magnitudeSquared ����ԭʼ���ͣ�ת��Ϊ float ���͡�������Ϊ powerSpectrum ����洢���� float ���͵�Ԫ�أ�������Ҫ������õ��� double ����ֵת��Ϊƥ������͡�
		//.at<float> ��cv::Mat ���ṩ�ĳ�Ա���������ڷ��ʾ�����ָ��λ�õ�Ԫ��  (i / sizeX, i % sizeX) �������ʾ����е�ָ��Ԫ��
		//powerSpectrum ��֮ǰ����� cv::Mat ���ͱ���������һ����ά�������������ڴ洢������Ĺ������ܶȡ�
	}

	// �Թ������ܶȽ��ж����߶ȴ���ʹ����ʺ����۹۲�
	cv::log(powerSpectrum, powerSpectrum); // ȡ��Ȼ����
	cv::normalize(powerSpectrum, powerSpectrum, 0, 255, cv::NORM_MINMAX); // ��һ����0-255��Χ��

	// ת��Ϊ8λ�Ҷ�ͼ�񲢽���fftshift����
	cv::Mat logPowerSpectrum8U;
	powerSpectrum.convertTo(logPowerSpectrum8U, CV_8U);
	cv::Mat logPowerSpectrumShifted;
	cv::flip(logPowerSpectrum8U, logPowerSpectrumShifted, -1);

	//������С
	cv::Mat resizedImage;
	cv::resize(logPowerSpectrumShifted, resizedImage, cv::Size(400, 400));
	// ��ʾͼ��
	cv::imshow("2D FFT Power Spectrum Density", resizedImage);
	cv::waitKey(0);
}

