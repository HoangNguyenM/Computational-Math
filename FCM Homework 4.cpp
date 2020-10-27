#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>
#include <chrono>

using namespace std;

// variable m is the number of intervals in [a,b]
// accuracy is the desired accuracy, we use 0.00001 for this program
// at the beginning of the program, please change the function "func" and "integral_value" as desired
// also change a and b for the interval
// if the program has error, reduce accuracy for Gauss-Legendre method

double pi = 3.14159265358979323846;
int m;
double accuracy;



// return function value
double func(double x)
{
	//return exp(x);
	//return exp(sin(2 * x)) * cos(2 * x);
	//return tanh(x);
	//return x * cos(2 * pi * x);
	return (x + 1 / x);
}

// return the true integral value
double integral_value(void)
{
	//return exp(3) - 1;
	//return (-1 + exp(sqrt(3) / 2)) / 2;
	//return log(cosh(1) / cosh(2));
	//return -1 / (2 * pow(pi, 2));
	return (pow(2.5, 2) - pow(0.1, 2)) / 2 + log(2.5 / 0.1);
}


// this function return the endpoints of the m intervals from a to b
vector<double> create_interval(double a, double b, int m)
{
	vector<double> vect = {};
	for (int i = 0; i <= m; i++)
	{
		vect.push_back(a + i * (b - a) / m);
	}

	return vect;
}

// approximate the integral using the composite midpoint rule
double midpoint(double a, double b, int m)
{
	vector<double> vec = create_interval(a, b, m);

	double h = (b - a) / m;
	double sum = 0;
	for (int i = 0; i < vec.size() - 1; i++)
	{
		sum = sum + func(vec[i] + h / 2);
	}
	
	sum = sum * h;

	return sum;
}


// approximate the integral using the composite trapezoidal rule
double trap(double a, double b, int m)
{
	vector<double> vec = create_interval(a, b, m);

	double h = (b - a) / m;
	double sum = 0;
	for (int i = 0; i < vec.size() - 1; i++)
	{
		sum = sum + func(vec[i]) + func(vec[i + 1]);
	}

	sum = sum * h / 2;

	return sum;
}


// approximate the integral using the Gauss-Legendre method
double gl(double a, double b, int m)
{
	vector<double> vec = create_interval(a, b, m);

	double h = (b - a) / m;
	double sum = 0;
	for (int i = 0; i < vec.size() - 1; i++)
	{
		double z0 = (-(vec[i + 1] - vec[i]) / sqrt(3) + vec[i] + vec[i + 1]) / 2;
		double z1 = ((vec[i + 1] - vec[i]) / sqrt(3) + vec[i] + vec[i + 1]) / 2;
		sum = sum + func(z0) + func(z1);
	}

	sum = sum * h / 2;
	
	return sum;
}

int main()
{
	ofstream file1;
	ofstream file2;
	ofstream file3;

	double a = 0.1;
	double b = 2.5;

	file1.open("midpoint rule.csv", ios::out);
	file2.open("trap rule.csv", ios::out);
	file3.open("gauss-legendre.csv", ios::out);

	file1 << "m" << ", " << "1" << ", " << "3" << ", " << "9" << ", " << "27" << ", " << "81" << ", " << "243" << "\n" << "integral values" << ", ";
	file2 << "m" << ", " << "1" << ", " << "2" << ", " << "4" << ", " << "8" << ", " << "16" << ", " << "32" << "\n" << "integral values" << ", ";
	file3 << "m" << ", " << "1" << ", " << "2" << ", " << "4" << ", " << "8" << ", " << "16" << ", " << "32" << "\n" << "integral values" << ", ";

	// approximate the integral using midpoint rule
	m = 1;
	accuracy = 0.00001;
	double track = 0;
	double e = 1;
	vector<double> integral = {};
	vector<double> error = {};

		// keep measuring the integral until it reaches the desired accuracy
	while (abs(e) > accuracy)
	{
		integral.push_back(midpoint(a, b, m));
		e = integral_value() - integral[track];
		error.push_back(e);
		m = 3 * m;
		track = track + 1;
	}

		// print out the integral values
	for (int i = 0; i < 6; i++)
	{
		file1 << integral[i] << ", ";
	}

	file1 << "\n" << "error" << ", ";

		// print out the errors
	for (int i = 0; i < 6; i++)
	{
		file1 << error[i] << ", ";
	}

	file1 << "\n" << "desired m" << ", ";

		// print out the required m to obtain the desired accuracy
	file1 << m / 3;

	// approximate the integral using trapezoidal rule
	m = 1;
	accuracy = 0.00001;
	track = 0;
	e = 1;
	integral = {};
	error = {};

		// keep measuring the integral until it reaches the desired accuracy
	while (abs(e) > accuracy)
	{
		integral.push_back(trap(a, b, m));
		e = integral_value() - integral[track];
		error.push_back(e);
		m = 2 * m;
		track = track + 1;
	}

		// print out the integral values
	for (int i = 0; i < 6; i++)
	{
		file2 << integral[i] << ", ";
	}

	file2 << "\n" << "error" << ", ";

		// print out the errors
	for (int i = 0; i < 6; i++)
	{
		file2 << error[i] << ", ";
	}

	file2 << "\n" << "desired m" << ", ";

		// print out the required m to obtain the desired accuracy
	file2 << m / 2;

	// approximate the integral using Gauss-Legendre method
	m = 1;
	accuracy = 0.00001;
	track = 0;
	e = 1;
	integral = {};
	error = {};

		// keep measuring the integral until it reaches the desired accuracy
	while (abs(e) > accuracy)
	{
		integral.push_back(gl(a, b, m));
		e = integral_value() - integral[track];
		error.push_back(e);
		m = 2 * m;
		track = track + 1;
	}

		// print out the integral values
	for (int i = 0; i < 6; i++)
	{
		file3 << integral[i] << ", ";
	}

	file3 << "\n" << "error" << ", ";

		// print out the errors
	for (int i = 0; i < 6; i++)
	{
		file3 << error[i] << ", ";
	}

	file3 << "\n" << "desired m" << ", ";

		// print out the required m to obtain the desired accuracy
	file3 << m / 2;



	file1.close();
	file2.close();
	file3.close();


	return 0;
}
