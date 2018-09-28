#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "Eigen/Core"
//#include "common_debug.h"



long long get_file_size(const std::string filename)
{
	/*
	�t�@�C���̃T�C�Y���擾����v���O����
	filename : �t�@�C����
	*/

	FILE *fp;
	struct _stat64 st;
	if (fopen_s(&fp, filename.c_str(), "rb") != 0){
		std::cerr << "Cannot open file: " << filename << std::endl;
		std::abort();
	}
	_fstat64(_fileno(fp), &st);
	fclose(fp);
	return st.st_size;
}

template< class T >
void read_vector(std::vector<T> &v, const std::string filename)
{
	std::cerr << "loading file: " << filename << std::endl;
	/*
	raw�摜��ǂݍ����vector�Ɋi�[
	v : �i�[����x�N�^�[
	filename : �t�@�C����
	*/
	auto num = get_file_size(filename) / sizeof(T);
	FILE *fp;
	if (fopen_s(&fp, filename.c_str(), "rb") != 0){
		std::cerr << "Cannot open file: " << filename << std::endl;
		std::abort();
	}
	v.resize(num);
	fread(v.data(), sizeof(T), num, fp);

	fclose(fp);
}

template< class T >
void write_vector(std::vector<T> &v, const std::string filename)
{
	std::cerr << "writing file: " << filename << std::endl;
	/*
	vector��raw�摜�ɏ�������
	v : �i�[����x�N�^�[
	filename : �ۑ��ꏊ��΃p�X
	*/
	FILE *fp;
	if (fopen_s(&fp, filename.c_str(), "wb") != 0){
		std::cerr << "Cannot open file: " << filename << std::endl;
		std::abort();
	}
	fwrite(v.data(), sizeof(T), v.size(), fp);
	fclose(fp);
}

template< class T >
void write_vector_if_not_exist(std::vector<T> &v, const std::string filename)
{
	std::cerr << "writing file: " << filename << std::endl;
	/*
	vector��raw�摜�ɏ�������
	v : �i�[����x�N�^�[
	filename : �ۑ��ꏊ��΃p�X
	*/
	FILE *fp;

	if (fopen_s(&fp, filename.c_str(), "wb") != 0){
		std::cerr << "Cannot open file: " << filename << std::endl;
		return;
	}
	fwrite(v.data(), sizeof(T), v.size(), fp);
	fclose(fp);
}

template<typename T>
void load_raw_to_eigen(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& eigen, std::string filename, int row_size)
{
	std::vector<T> v;
	auto num = get_file_size(filename) / sizeof(T);
	FILE *fp;
	if (fopen_s(&fp, filename.c_str(), "rb") != 0){
		std::cerr << "Cannot open file: " << filename << std::endl;
		std::abort();
	}
	v.resize(num);
	fread(v.data(), sizeof(T), num, fp);

	fclose(fp);

	size_t col_size = v.size() / row_size;

	eigen = eigen.setZero(row_size, col_size);


	for (long int y = 0; y < row_size; y++)
	{
		for (long int x = 0; x < col_size; x++)
		{
			eigen(y, x) = v[x * row_size + y];
		}
	}
}



template<typename T>
void write_raw_and_txt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& data, std::string filename)
{
	/*
	eigen��matrix��raw�摜�Ƃ��ĕۑ�
	�o�͂�raw�摜�ƃT�C�Y�ƌ^���L���ꂽtxt�t�@�C��
	data : �ۑ�����f�[�^
	filename : �g���q�̑O�܂ł̃p�X
	*/

	size_t rows = data.rows();
	size_t cols = data.cols();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Data;
	Data = data;


	/*std::ofstream fs1(filename + ".txt");
	fs1 << "rows = " << rows << std::endl;
	fs1 << "cols = " << cols << std::endl;
	fs1 << typeid(data).name() << std::endl;
	fs1.close();*/

	std::vector<T> save_data(rows * cols);
	Data.resize(data.rows()*data.cols(), 1);
	for (size_t i = 0; i < save_data.size(); i++)
		save_data[i] = Data(i, 0);
	write_vector(save_data, filename + ".raw");
	Data.resize(rows, cols);
}