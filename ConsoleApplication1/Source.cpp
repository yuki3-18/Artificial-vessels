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

	//ê≥ãKï™ïz
	double A = 667;
	double sigma = 0.906;
	double result = A * exp(- (r*r) / 2.0 / pow(sigma,2.0));
	
	return result;
}

int main(int argc, char * argv[]) {


	//setting
	std::string output_path = /*argv[1]*//*"E:/git/TFRecord_example/in/2axis/noise/train/"*/"E:/from_kubo/ConsoleApplication1/x64/Release/output/Glow/";

	int data_size = /*atoi(argv[2])*/10;

	nari::mhd mhdI;

	const int patch_side_size = 8;
	const int patch_size = patch_side_size * patch_side_size * patch_side_size;
	int hani = (int)patch_side_size / 2;
	double range = M_PI / 2;

	//input_image
	Eigen::MatrixXd image_vec;
	image_vec.resize(patch_size, 1);

	//rnd_theta
	std::random_device seed_theta;
	std::mt19937_64 mt_theta(seed_theta());
	std::uniform_real_distribution<double> dist2(-M_PI / 2, M_PI / 2);
	Eigen::VectorXd rnd_theta; rnd_theta.resize(data_size, 1);
	for (int i = 0; i < data_size; i++)
		rnd_theta(i) = dist2(mt_theta);

	//rnd_phi
	std::random_device seed_phi;
	std::mt19937_64 mt_phi(seed_phi());
	std::uniform_real_distribution<double> dist3(0, 1.0);
	Eigen::VectorXd rnd_phi; rnd_phi.resize(data_size, 1);
	for (int i = 0; i < data_size; i = i + 2) {
		rnd_phi(i) = asin(dist3(mt_phi));
		rnd_phi(i + 1) = -asin(dist3(mt_phi));
	}

	//rnd_delta
	std::random_device seed_delta;
	std::mt19937_64 mt_delta(seed_delta());
	std::uniform_real_distribution<double> dist4(-2.188, 2.188);
	//Eigen::MatrixXd delta; delta = Eigen::MatrixXd::Zero(3, data_size);
	//for (int i = 0; i < data_size; i++) {
	//	double dx = dist4(mt_delta);	double dy = dist4(mt_delta);	double dz = dist4(mt_delta);
	//	Eigen::Vector3d d; d << dx, dy, dz;
	//	//delta.leftCols(i + 1) = d;
	//	delta.col(i) = d;
	//	std::cout << "i=" << i << std::endl;
	//	std::cout << "d=" << d/*x << dy << dz */<< std::endl;
	//	std::cout << "delta=" << delta << std::endl;
	//}

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
		//rotation_temp = Rz * Ry;
		rotation_temp = Rz ;

		//rotation
		Eigen::MatrixXd rotation = rotation_temp.leftCols(1);
		//std::cout << "loop=" << loop << std::endl;

		//shift
		//Eigen::MatrixXd shift = delta.col(loop);
		int dx = dist4(mt_delta);	int dy = dist4(mt_delta);	int dz = dist4(mt_delta);
		//std::cout << "dx=" << dx <<"dy=" << dy << "dz=" << dz << std::endl;
		
		//rnd_noise
		int mean = 0;
		double var = 22.2;
		std::random_device seed_noise;
		std::mt19937_64 engine(seed_noise());
		std::normal_distribution<double> dist(mean, var);
		Eigen::VectorXd noise; noise.resize(patch_size, 1);
		for (int i = 0; i < patch_size; i++)
			noise(i) = dist(engine);

		//
		int i = 0;
		for (int z = -hani + dz; z <= hani + dz; z++) {
			for (int y = -hani + dy; y <= hani + dy; y++){
				for (int x = -hani + dx; x <= hani + dx; x++){
					Eigen::MatrixXd X;
					X.resize(3, 1);
					X << x , y , z;

					Eigen::MatrixXd X_dash;
					X_dash.resize(3, 1);
					X_dash = rotation * (rotation.transpose() * X) - X;
					double r = sqrt(X_dash.squaredNorm());
					image_vec(i) = function(r);
					//image_vec(i) = function(r) + noise(i);
					i++;
				}
			}
		}
		
		std::string theta = std::to_string(rnd_theta(loop)*180.0/M_PI);
		std::string phi = std::to_string(/*rnd_phi(loop)*180.0/M_PI*/0);

		std::vector<double> image_std(image_vec.data(), image_vec.data() + image_vec.size());

		ImageIO<3> imageio;
		imageio.SetSize(0, 8);
		imageio.SetSize(1, 8);
		imageio.SetSize(2, 8);
		imageio.SetSpacing(0, 0.885);
		imageio.SetSpacing(1, 0.885);
		imageio.SetSpacing(2, 1);
		imageio.Write(image_std, output_path + "output_" + theta + "_" + phi + "_learn.mhd");

		write_raw_and_txt(image_vec, output_path + "output_" + theta + "_" + phi + "_learn");
		std::ofstream outputfile(output_path + "filename.txt", std::ios::app);
		outputfile << output_path + "output_" + theta + "_" + phi + "_learn" + ".mhd" << std::endl;
		outputfile.close();

		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(output_path + "output_" + theta + "_" + phi + "_learn");

	

	}
	//system("pause");
	return	0;
}