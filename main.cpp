#include<iostream>
#include<vector>
#include <fstream>
#include <string>
#include"defocusFileFormats.h"
#include "mydefocus.h"
#include "mrcStack.h"
#include"myCorrection.h"
#include"fftRoutines.h"
#include"myfilterprojection.h"
#include"projectionset.h"
#include"myctf3d.h"
using namespace std;

#include <opencv2\opencv.hpp>

using namespace cv;

int test_opencv()
{
	Mat img;
	img = imread("E:\\VS studio\\proj\\opcv test\\RTTT.png");
	if (img.empty())
	{
		cout << "��ȷ��ͼ���ļ����Ƿ���ȷ" << endl;
		return 0;
	}
	imshow("test", img);
	waitKey(0);
	return 0;
}

int main()
{
	//����opencv
   // test_opencv();

	
	
	int numop = 41;
	CTFFind4 test;
	std::vector<std::vector<float>> values(numop, std::vector<float>(4));
	string inputCTFout = "D:\\�о���\\Data\\Fudan\\Position01\\Position_01_01_output.txt";
	test.read(inputCTFout,numop);
	vector<vector<float>> check=test.getOriginalValues();
	
	//��ȡ��ʼɢ���ļ� �Ĳ���
	/*
	for (int i = 0; i < check.size(); i++)
	{
		for (int j = 0; j < check[i].size(); j++)
		{
			cout << check[i][j]; // �Ƴ� endl��ֻ���ڲ�ѭ�����������
			if (j != check[i].size() - 1) 
			{ // ������ǵ�ǰ�е����һ��Ԫ�أ������һ���ո�
				cout << ' ';
			}
		}
		cout << endl;
	}
	*/

	//���Դӳ�ʼ�ļ���ȡǰ�ĸ�ֵ
    test.getValues(values, numop);
	//���values������
	test.printValues(values);

	//����ȡ����
	string inputCoords = "D:\\�о���\\Data\\Fudan\\01coord.txt";
	mydefoucs reader;
	auto coords=reader.Readcoordinates(inputCoords);
	cout << "Number of particles: " << coords.size() << '\n';

	// ���������������������
	/*for (const auto& coord : coords)
	{
		std::cout << "Coordinates: (" << coord.x << ", " << coord.y << ", " << coord.z << ")\n";
	}*/

	//����ȡ�Ƕ�
	string anglesFile = "D:\\�о���\\Data\\Fudan\\Position01\\Position_01_01.tlt";

	auto angles = reader.Readangle(anglesFile);
	cout << "Number of angles: " << angles.size() << '\n';

	/*
	for (size_t i = 0; i < angles.size(); ++i)
	{
		std::cout << "Angle " << i + 1 << ": " << angles[i] << '\n';
	}
	*/

	//���ÿ�����ӵļ�����   ��ȷ ����洢��ÿ��������ÿ���Ƕ��µ�4��ֵ��
	reader.Calculateparticledefocus(coords, angles, values);
	const auto& computedResults = reader.getResults();    
	//��ӡ���
	int Coordnum = computedResults.size();
	cout << Coordnum <<endl;

	/*
	for (size_t i = 0; i < computedResults.size(); ++i) {
		// ��ӡ������Ϣ
		cout << "Coordinates: (" << coords[i].x << ", " << coords[i].y << ", " << coords[i].z << ")" << endl;

		for (size_t j = 0; j < computedResults[i].size(); ++j) {
			cout << "Angle: " << angles[j] << " degrees" << endl;

			for (size_t k = 0; k < computedResults[i][j].size(); ++k) {
				cout << "Computed result " << k + 1 << ": " << computedResults[i][j][k] << endl;
			}
			cout << endl; // �ָ���ͬ�ǶȵĽ��
		}
		cout << "-----------------------------" << endl; // �ָ���ͬ����ļ�����
	}*/


	//���԰������ӷֱ�дɢ���ļ�
	string outpath = "D:\\�о���\\ʵ��\\FudanTest\\0101\\defcous\\";
	reader.WriteResultsToFile(outpath, computedResults);
// ���϶��������ж�ȡ��ʼɢ��������ÿ�������ڲ�ͬ�Ƕ��µ�ɢ�����洢���� ��




//�����ǽ���CTFУ���������У�����stack�ļ���
	const string baseDir = "D:\\�о���\\ʵ��\\FudanTest\\0101\\";
	const string inputDir = baseDir + "corp_mrc\\";
	const string defocusDir = baseDir + "defcous\\";
	const string outputMRCDir = baseDir + "ctf_mrc\\";
	const string outputMaskDir = baseDir + "mask\\";
	for (int i = 0; i < Coordnum; ++i)
	{
		string inputName = inputDir + "corp_"+ to_string(i) + ".mrcs";
		string inputDefocusName = defocusDir +"defocus_"+to_string(i) + ".defcous";
		string outputMRC = outputMRCDir + "output_" + to_string(i) + ".mrcs";
		string outputMask = outputMaskDir + "mask_" + to_string(i) + ".mrcs";

		myCorrection test2(inputName);
		test2.run(inputDefocusName, outputMRC, outputMask);
		cout << "done"<< i << endl;
		// ����ʵ����Ҫ��������Ҫ���ÿ�������Ƿ�ɹ�
	}



	/* ��Ϊimod�ж������������ӵ����
    string inputmrcname = "E:\\VS studio\\proj\\project2\\test_ctf_output\\output.mrc";   //δ���з�����ת��
    //string inputmrcname = "E:\\VS studio\\proj\\project2\\test_ctf_output\\output_flipped.mrc";
	novaCTF::Vec3ui volumeDimensions(200, 200, 200);
    CTF3d test4(inputmrcname,volumeDimensions);
	test4.run();
	*/



	return 0;
}
