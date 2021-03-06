#pragma once

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
using namespace std;


class MonteCarloProbability
{
public:
	MonteCarloProbability();
	~MonteCarloProbability();
	// return a uniformly distributed random number
	double uniformRandom()
	{
		return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
	}
	// return a normally distributed random number
	double normalRandom()
	{
		double u1 = uniformRandom();
		double u2 = uniformRandom();
		return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
	}

	//int main()
	//{
	//	for (int i = 0; i<100; i++)
	//		cout << " x_i = " << normalRandom() << endl;
	//	return 0;
	//}

};

