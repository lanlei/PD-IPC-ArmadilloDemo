#pragma once
#ifndef  __RONDOMVALUE_H__
#define __RONDOMVALUE_H__
//===============================================================
//Summary:
//          This File...
//FileName:
//          RondomValue.h
//Remarks:
//          ...
//Date:
//          2018/4/20
//===============================================================
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <Eigen/Eigen>

namespace RondomValue
{
	double getRandomValue();
	int getIntRandomValue(int max = 0, int min = 0);
	double getGaussianDistribution(double mean, double dev);
	void getGaussianDistribution(double mean, double dev, Eigen::VectorXd& value);
}

#endif //__RONDOMVALUE_H__