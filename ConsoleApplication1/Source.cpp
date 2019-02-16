#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<random>
#include<Eigen/LU>
#include<Eigen/Core>
#include<Eigen/Dense>
#include"function.h"
#include <other/narimhd.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "ItkImageIO.h"

typedef unsigned short	PixelType;
const unsigned int	Dimension = 3;
typedef itk::Image< PixelType, Dimension >	ImageType;
typedef itk::ImageFileWriter< ImageType > WriterType;

double function(double r){

	//正規分布
	double A = 667;
	double sigma = 0.906;
	double result = A * exp(- (r*r) / 2.0 / pow(sigma,2.0));
	
	return result;
}

int main(int argc, char * argv[]){


	//setting
	std::string output_path = /*argv[1]*/"E:/git/TFRecord_example/input/shift/train/"/*"E:/from_kubo/ConsoleApplication1/x64/Release/output/shift/"*/;

	int data_size = /*atoi(argv[2])*/10000;

	nari::mhd mhdI;

	const int patch_side_size = 15;
	const int patch_size = patch_side_size * patch_side_size * patch_side_size;
	int hani = (int)patch_side_size / 2;
	double range = M_PI / 2;

	//input_image
	Eigen::MatrixXd image_vec;
	image_vec.resize(patch_size, 1);

	//rnd_theta
	std::random_device seed_theta;
	std::mt19937_64 mt_theta(seed_theta());
	std::uniform_real_distribution<double> dist2(-M_PI / 2, M_PI/2);
	Eigen::VectorXd rnd_theta; rnd_theta.resize(data_size, 1);
	for (int i = 0; i < data_size; i++)
		rnd_theta(i) = dist2(mt_theta);

	//rnd_phi
	std::random_device seed_phi;
	std::mt19937_64 mt_phi(seed_phi());
	std::uniform_real_distribution<double> dist3(0, 1.0);
	Eigen::VectorXd rnd_phi; rnd_phi.resize(data_size, 1);
	for (int i = 0; i < data_size; i = i + 2){
		rnd_phi(i) = asin(dist3(mt_phi));
		rnd_phi(i+1) = -asin(dist3(mt_phi));
	}

	//rnd_delta
	std::random_device seed_delta;
	std::mt19937_64 mt_delta(seed_delta());
	std::uniform_real_distribution<double> dist4(-2.688, 2.688);
	Eigen::MatrixXd rnd_delta = Eigen::MatrixXd::Zero(3, data_size);
	for (int i = 0; i < data_size; i++) {
		double x = dist4(mt_delta); double y = dist4(mt_delta); double z = dist4(mt_delta);
		rnd_delta(0, i) = x;
		rnd_delta(1, i) = y;
		rnd_delta(2, i) = z;
		//std::cout << rnd_delta << std::endl;
		//system("pause");
	}


	//make_loop
	for (int loop = 0; loop < data_size; loop++){

		//rotation_ini
		Eigen::MatrixXd Ry;	Ry.resize(3, 3);
		Ry << cos(rnd_phi(loop)), 0, sin(rnd_phi(loop)),
			0, 1, 0,
			-sin(rnd_phi(loop)), 0, cos(rnd_phi(loop));

		Eigen::MatrixXd Rz; Rz.resize(3, 3);
		Rz << cos(rnd_theta(loop)), -sin(rnd_theta(loop)), 0,
			sin(rnd_theta(loop)), cos(rnd_theta(loop)), 0,
			0, 0, 1;

		Eigen::MatrixXd rotation_temp;
		rotation_temp = Rz * Ry;
		//rotation_temp = Rz ;

		//rotation
		Eigen::MatrixXd rotation = rotation_temp.leftCols(1);
		
		//rnd_noise
		int mean = 0;
		double var = 22.2;
		std::random_device seed_noise;
		std::mt19937_64 engine(seed_noise());
		std::normal_distribution<double> dist(mean, var);
		Eigen::VectorXd noise; noise.resize(patch_size, 1);
		for (int i = 0; i < patch_size; i++)
			noise(i) = dist(engine);


		// make_structure
		int i = 0;
		for (int z = -hani; z <= hani; z++){
			for (int y = -hani; y <= hani; y++){
				for (int x = -hani; x <= hani; x++){
					Eigen::MatrixXd X;
					X.resize(3, 1);
					X << x, y, z;

					Eigen::MatrixXd X_dash;
					X_dash.resize(3, 1);
					X_dash = X - rotation * (rotation.transpose() * X);
					double r = sqrt(X_dash.squaredNorm());
					//double r = sqrt(X.squaredNorm());
					//image_vec(i) = function(r);
					image_vec(i) = function(r) + noise(i);
					i++;
				}
			}
		}



		//input_matrix
		int ajd = (patch_side_size - 1) / 2;	//入力画像の中心座標算出
		Eigen::MatrixXd In; In.resize(Dimension, patch_size);	//出力画像のサイズ
		for (int k = 0; k < patch_side_size; k++) {
			for (int j = 0; j < patch_side_size; j++) {
				for (int i = 0; i < patch_side_size; i++) {
					In(0, i + patch_side_size * j + patch_side_size * patch_side_size * k) = i;
					In(1, i + patch_side_size * j + patch_side_size * patch_side_size * k) = j;
					In(2, i + patch_side_size * j + patch_side_size * patch_side_size * k) = k;
				}
			}
		}

		//Eigen::MatrixXd image_vec; image_vec.resize(Dimension, patch_size);
		//Eigen::Vector3d X;
		//for (int i = 0; i < patch_size; i++) {
		//	X = In.col(i);
		//	image_vec.col(i) = X - rotation * (rotation.transpose() * X);
		//}

		Eigen::MatrixXd Inout; Inout.resize(Dimension, patch_size);
		Eigen::Vector3d in;
		for (int i = 0; i < patch_size; i++) {
			in = In.col(i);
			Inout.col(i) = in - rnd_delta.col(loop);
		}

		//calcu_coordination
		Eigen::MatrixXi temp;
		temp.resize(Dimension, patch_size);
		for (int i = 0; i < Dimension; i++)
			for (int j = 0; j < patch_size; j++)
				temp(i, j) = (int)Inout(i, j);

		//vec_copy
		std::vector<double> v;
		v.resize(image_vec.size());
		Eigen::VectorXd::Map(&v[0], image_vec.size()) = image_vec;

		//interpolation		
		Eigen::MatrixXd output = Eigen::MatrixXd::Zero(patch_size, 1);
		for (unsigned int coordi = 0; coordi < patch_size; coordi++)
		{
			//縁処理
			int x0 = int_max(0, int_min(patch_side_size - 1, (int)temp(0, coordi)));
			int y0 = int_max(0, int_min(patch_side_size - 1, (int)temp(1, coordi)));
			int z0 = int_max(0, int_min(patch_side_size - 1, (int)temp(2, coordi)));
			int x1 = int_max(0, int_min(patch_side_size - 1, (int)(temp(0, coordi) + 1)));
			int y1 = int_max(0, int_min(patch_side_size - 1, (int)(temp(1, coordi) + 1)));
			int z1 = int_max(0, int_min(patch_side_size - 1, (int)(temp(2, coordi) + 1)));


			double xd = (Inout(0, coordi) - x0) / (x1 - x0 + 1e-5);
			double yd = (Inout(1, coordi) - y0) / (y1 - y0 + 1e-5);
			double zd = (Inout(2, coordi) - z0) / (z1 - z0 + 1e-5);


			double c00 = value(x0, y0, z0, v) * (1.0 - xd)
				+ value(x1, y0, z0, v) * xd;
			double c01 = value(x0, y0, z1, v) * (1.0 - xd)
				+ value(x1, y0, z1, v) * xd;
			double c10 = value(x0, y1, z0, v) * (1.0 - xd)
				+ value(x1, y1, z0, v) * xd;
			double c11 = value(x0, y1, z1, v) * (1.0 - xd)
				+ value(x1, y1, z1, v) * xd;

			double c0 = c00 * (1.0 - yd) + c10 * yd;
			double c1 = c01 * (1.0 - yd) + c11 * yd;
			double c = c0 * (1.0 - zd) + c1 * zd;

			output(coordi) = c;

			//std::cout << coordi << std::endl;
			//std::cout << "c00 " << c00 << std::endl;
			//std::cout << "c01 " << c01 << std::endl;
			//std::cout << "c10 " << c10 << std::endl;
			//std::cout << "c11 " << c11 << std::endl;
			//std::cout << "de " << value(x0, y0, z1, v) << std::endl;
			//std::cout << "de " << value(x1, y0, z1, v) << std::endl;
			//std::cout << "xd " << xd << std::endl;;
			//std::cout << "yd " << yd << std::endl;;
			//std::cout << "zd " << zd << std::endl;;
			//std::cout << "c " << c << std::endl;
			//
			//system("pause");
		}

		//int out_patch_side = 9;
		//Eigen::MatrixXd crop = Eigen::MatrixXd::Zero(out_patch_side, out_patch_side, out_patch_side);

		//for(int x = 3; x < out_patch_side; x++)
		//	for(int y = 3; y < out_patch_side; y++)
		//		for(int z = 3; z < out_patch_side; z++)
		//		{
		//			int position = x + out_patch_side * y + out_patch_side * out_patch_side * z;
		//			crop[x][y][z] = output[position];			
		//		}

		
		std::string theta = std::to_string(rnd_theta(loop)*180.0/M_PI);
		std::string phi = std::to_string(rnd_phi(loop)*180.0/M_PI/*0*/);

		std::vector<double> image_std(output.data(), output.data() + output.size());

		ImageIO<3> imageio;
		imageio.SetSize(0, patch_side_size);
		imageio.SetSize(1, patch_side_size);
		imageio.SetSize(2, patch_side_size);
		imageio.SetSpacing(0, 0.885);
		imageio.SetSpacing(1, 0.885);
		imageio.SetSpacing(2, 1);
		imageio.Write(image_std, output_path + "output_" + theta + "_" + phi + "_learn.mhd");

		write_raw_and_txt(output, output_path + "output_" + theta + "_" + phi + "_learn");
		std::ofstream outputfile(output_path + "filename.txt", std::ios::app);
		outputfile << output_path + "output_" + theta + "_" + phi + "_learn" + ".mhd" << std::endl;
		outputfile.close();

		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(output_path + "output_" + theta + "_" + phi + "_learn");
	

	}
	return	0;
}