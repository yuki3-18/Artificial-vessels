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

	//���K���z
	double A = 667;
	double sigma = 0.906;
	double result = A * exp(- (r*r) / 2.0 / pow(sigma,2.0));
	
	return result;
}

int main(int argc, char * argv[]){


	//setting
	std::string output_path = argv[1]/*"E:/git/TFRecord_example/in/2axis/noise/train/"*//*"E:/from_kubo/ConsoleApplication1/x64/Release/output/"*/;

	int data_size = atoi(argv[2])/*10*/;

	nari::mhd mhdI;

	const int patch_side_size = 9;
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

		//
		int i = 0;
		for (int z = -hani; z <= hani; z++){
			for (int y = -hani; y <= hani; y++){
				for (int x = -hani; x <= hani; x++){
					Eigen::MatrixXd X;
					X.resize(3, 1);
					X << x, y, z;

					Eigen::MatrixXd X_dash;
					X_dash.resize(3, 1);
					X_dash = X - rotation*(rotation.transpose() * X);
					double r = sqrt(X_dash.squaredNorm());
					//image_vec(i) = function(r);
					image_vec(i) = function(r) + noise(i);
					//image_vec(i) = function(r + delta);
					i++;
				}
			}
		}
		
		std::string theta = std::to_string(rnd_theta(loop)*180.0/M_PI);
		std::string phi = std::to_string(rnd_phi(loop)*180.0/M_PI/*0*/);

		std::vector<double> image_std(image_vec.data(), image_vec.data() + image_vec.size());

		ImageIO<3> imageio;
		imageio.SetSize(0, 9);
		imageio.SetSize(1, 9);
		imageio.SetSize(2, 9);
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
	return	0;
}