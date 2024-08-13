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
		cout << "请确认图像文件名是否正确" << endl;
		return 0;
	}
	imshow("test", img);
	waitKey(0);
	return 0;
}

int main()
{
	//测试opencv
   // test_opencv();

	
	
	int numop = 41;
	CTFFind4 test;
	std::vector<std::vector<float>> values(numop, std::vector<float>(4));
	string inputCTFout = "D:\\研究生\\Data\\Fudan\\Position01\\Position_01_01_output.txt";
	test.read(inputCTFout,numop);
	vector<vector<float>> check=test.getOriginalValues();
	
	//读取初始散焦文件 的测试
	/*
	for (int i = 0; i < check.size(); i++)
	{
		for (int j = 0; j < check[i].size(); j++)
		{
			cout << check[i][j]; // 移除 endl，只在内部循环结束后添加
			if (j != check[i].size() - 1) 
			{ // 如果不是当前行的最后一个元素，则添加一个空格
				cout << ' ';
			}
		}
		cout << endl;
	}
	*/

	//测试从初始文件中取前四个值
    test.getValues(values, numop);
	//输出values的内容
	test.printValues(values);

	//检查读取坐标
	string inputCoords = "D:\\研究生\\Data\\Fudan\\01coord.txt";
	mydefoucs reader;
	auto coords=reader.Readcoordinates(inputCoords);
	cout << "Number of particles: " << coords.size() << '\n';

	// 遍历并输出所有坐标数据
	/*for (const auto& coord : coords)
	{
		std::cout << "Coordinates: (" << coord.x << ", " << coord.y << ", " << coord.z << ")\n";
	}*/

	//检查读取角度
	string anglesFile = "D:\\研究生\\Data\\Fudan\\Position01\\Position_01_01.tlt";

	auto angles = reader.Readangle(anglesFile);
	cout << "Number of angles: " << angles.size() << '\n';

	/*
	for (size_t i = 0; i < angles.size(); ++i)
	{
		std::cout << "Angle " << i + 1 << ": " << angles[i] << '\n';
	}
	*/

	//检查每个粒子的计算结果   正确 结果存储了每个坐标在每个角度下的4个值。
	reader.Calculateparticledefocus(coords, angles, values);
	const auto& computedResults = reader.getResults();    
	//打印输出
	int Coordnum = computedResults.size();
	cout << Coordnum <<endl;

	/*
	for (size_t i = 0; i < computedResults.size(); ++i) {
		// 打印坐标信息
		cout << "Coordinates: (" << coords[i].x << ", " << coords[i].y << ", " << coords[i].z << ")" << endl;

		for (size_t j = 0; j < computedResults[i].size(); ++j) {
			cout << "Angle: " << angles[j] << " degrees" << endl;

			for (size_t k = 0; k < computedResults[i][j].size(); ++k) {
				cout << "Computed result " << k + 1 << ": " << computedResults[i][j][k] << endl;
			}
			cout << endl; // 分隔不同角度的结果
		}
		cout << "-----------------------------" << endl; // 分隔不同坐标的计算结果
	}*/


	//测试按照粒子分别写散焦文件
	string outpath = "D:\\研究生\\实验\\FudanTest\\0101\\defcous\\";
	reader.WriteResultsToFile(outpath, computedResults);
// 以上都是来进行读取初始散焦，计算每个粒子在不同角度下的散焦并存储保存 。




//下面是进程CTF校正并且输出校正后的stack文件。
	const string baseDir = "D:\\研究生\\实验\\FudanTest\\0101\\";
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
		// 根据实际需要，可能需要检查每个操作是否成功
	}



	/* 因为imod有独立功能所以扔掉这块
    string inputmrcname = "E:\\VS studio\\proj\\project2\\test_ctf_output\\output.mrc";   //未进行方向旋转的
    //string inputmrcname = "E:\\VS studio\\proj\\project2\\test_ctf_output\\output_flipped.mrc";
	novaCTF::Vec3ui volumeDimensions(200, 200, 200);
    CTF3d test4(inputmrcname,volumeDimensions);
	test4.run();
	*/



	return 0;
}
