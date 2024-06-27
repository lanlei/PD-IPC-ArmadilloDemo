#include "RondomValue.h"

double RondomValue::getRandomValue()
{
	//std::rand
	std::random_device rd;
	std::mt19937 mt(rd());

	return rd();
}

int RondomValue::getIntRandomValue(int max, int min)
{
	std::random_device rd;  // 将用于为随机数引擎获得种子
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> u(min, max);	
	return u(gen);
}

double RondomValue::getGaussianDistribution(double mean, double dev)
{
	std::random_device rd;
	std::mt19937 gen{ rd() };
	// 值最可能接近平均
	// 标准差影响生成的值距离平均数的分散
	std::normal_distribution<> d{ mean, dev};

	return d(gen);
}

void RondomValue::getGaussianDistribution(double mean, double dev, Eigen::VectorXd & value)
{
	std::random_device rd;
	std::mt19937 gen{ rd() };
	// 值最可能接近平均
	// 标准差影响生成的值距离平均数的分散
	std::normal_distribution<> d{ mean, dev };

	for(int i = 0; i < value.size(); i++)
		value.data()[i] = d(gen);
}
