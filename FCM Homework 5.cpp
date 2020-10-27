#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>
#include <chrono>

// The four methods are run in the main programs
// The function F(t) can be changed below in the function "func"
// Also the derivative F'(t) needs to be changed in the function "derivative"

using namespace std;

int lambda;
double ynot;


// return function F value
double func(double x)
{
	return sin(x);
	//return exp(x);
	//return 0;
}


// return the derivative of the function F
double derivative(double x)
{
	return cos(x);
	//return exp(x);
	//return 0;
}

// return the function f value
double smallf(double y, double t)
{
	return lambda * (y - func(t)) + derivative(t);
}


// return the true value of y at t
double true_value(double t)
{
	return (ynot - func(0)) * exp(lambda * t) + func(t);
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

// calculate the norm difference, here vec1 and vec2 have the same size
double norm(vector<double> vec1, vector<double> vec2)
{
	double max = 0;
	for (int i = 0; i < vec1.size(); i++)
	{
		double x = abs(vec1[i] - vec2[i]);
		if (x > max)
		{
			max = x;
		}
	}

	return max;
}


// Evaluate Adams Bashforth two-step method
void AB2(double a, double b, int m)
{
	vector<double> y = {};
	vector<double> truey = {};
	vector<double> f = {};
	vector<double> t = create_interval(a, b, m);
	double h = (b - a) / m;
	
	// process the initial conditions
	y.push_back(ynot);
	y.push_back(true_value(t[1]));
	truey.push_back(ynot);
	truey.push_back(true_value(t[1]));
	f.push_back(smallf(ynot, t[0]));
	f.push_back(smallf(y[1], t[1]));

	// evaluate y using the method
	for (int i = 2; i <= m; i++)
	{
		y.push_back(y[i - 1] + h * (3 * f[i - 1] - f[i - 2]) / 2);
		f.push_back(smallf(y[i], t[i]));
		truey.push_back(true_value(t[i]));
	}

	// print out the values
	ofstream file1;
	ofstream file2;
	ofstream file3;
	ofstream file4;


	file1.open("AB2 - y.csv", ios::out);
	file2.open("AB2 - true y.csv", ios::out);
	file3.open("AB2 - errors.csv", ios::out);
	file4.open("tn.csv", ios::out);

	for (int i = 0; i <= m; i++)
	{
		file1 << y[i] << "\n";
		file2 << truey[i] << "\n";
		file4 << t[i] << "\n";
	}

	file3 << "First step local error" << ", " << abs(y[2] - truey[2]) << "\n";
	file3 << "Final global error" << ", " << abs(y[m] - truey[m]) << "\n";
	file3 << "Maximum global error" << ", " << norm(y, truey);
	
	file1.close();
	file2.close();
	file3.close();
	file4.close();

}

// Evaluate Runge Kutta four-stage method
void RK4(double a, double b, int m)
{
	vector<double> y = {};
	vector<double> truey = {};
	vector<double> t = create_interval(a, b, m);
	double h = (b - a) / m;

	// process the initial conditions
	y.push_back(ynot);
	truey.push_back(ynot);

	// evaluate y using the method
	for (int i = 1; i <= m; i++)
	{
		double K1 = smallf(y[i - 1], t[i - 1]);
		double K2 = smallf(y[i - 1] + h * K1 / 2, t[i - 1] + h / 2);
		double K3 = smallf(y[i - 1] + h * K2 / 2, t[i - 1] + h / 2);
		double K4 = smallf(y[i - 1] + h * K3, t[i]);
		y.push_back(y[i - 1] + h * (K1 / 6 + K2 / 3 + K3 / 3 + K4 / 6));

		truey.push_back(true_value(t[i]));
	}

	// print out the values
	ofstream file1;
	ofstream file2;
	ofstream file3;


	file1.open("RK4 - y.csv", ios::out);
	file2.open("RK4 - true y.csv", ios::out);
	file3.open("RK4 - errors.csv", ios::out);

	for (int i = 0; i <= m; i++)
	{
		file1 << y[i] << "\n";
		file2 << truey[i] << "\n";
	}

	file3 << "First step local error" << ", " << abs(y[1] - truey[1]) << "\n";
	file3 << "Final global error" << ", " << abs(y[m] - truey[m]) << "\n";
	file3 << "Maximum global error" << ", " << norm(y, truey);

	file1.close();
	file2.close();
	file3.close();

}


// Evaluate Adams Moulton one-step method
void AM(double a, double b, int m)
{
	vector<double> y = {};
	vector<double> truey = {};
	vector<double> t = create_interval(a, b, m);
	double h = (b - a) / m;

	// process the initial conditions
	y.push_back(ynot);
	truey.push_back(ynot);

	// evaluate y using the method
	for (int i = 1; i <= m; i++)
	{
		double temp = smallf(y[i - 1], t[i - 1]) - lambda * func(t[i]) + derivative(t[i]);
		temp = temp * h / 2 + y[i - 1];
		temp = temp / (1 - h * lambda / 2);
		y.push_back(temp);

		truey.push_back(true_value(t[i]));
	}

	// print out the values
	ofstream file1;
	ofstream file2;
	ofstream file3;


	file1.open("AM - y.csv", ios::out);
	file2.open("AM - true y.csv", ios::out);
	file3.open("AM - errors.csv", ios::out);

	for (int i = 0; i <= m; i++)
	{
		file1 << y[i] << "\n";
		file2 << truey[i] << "\n";
	}

	file3 << "First step local error" << ", " << abs(y[1] - truey[1]) << "\n";
	file3 << "Final global error" << ", " << abs(y[m] - truey[m]) << "\n";
	file3 << "Maximum global error" << ", " << norm(y, truey);

	file1.close();
	file2.close();
	file3.close();

}


// Evaluate BDF two-step method
void BDF2(double a, double b, int m)
{
	vector<double> y = {};
	vector<double> truey = {};
	vector<double> t = create_interval(a, b, m);
	double h = (b - a) / m;

	// process the initial conditions
	y.push_back(ynot);
	y.push_back(true_value(t[1]));
	truey.push_back(ynot);
	truey.push_back(true_value(t[1]));

	// evaluate y using the method
	for (int i = 2; i <= m; i++)
	{
		double temp = derivative(t[i]) - lambda * func(t[i]);
		temp = temp * h * 2 / 3 + 4 * y[i - 1] / 3 - y[i - 2] / 3;
		temp = temp / (1 - 2 * h * lambda / 3);
		y.push_back(temp);

		truey.push_back(true_value(t[i]));
	}

	// print out the values
	ofstream file1;
	ofstream file2;
	ofstream file3;


	file1.open("BDF2 - y.csv", ios::out);
	file2.open("BDF2 - true y.csv", ios::out);
	file3.open("BDF2 - errors.csv", ios::out);

	for (int i = 0; i <= m; i++)
	{
		file1 << y[i] << "\n";
		file2 << truey[i] << "\n";
	}

	file3 << "First step local error" << ", " << abs(y[2] - truey[2]) << "\n";
	file3 << "Final global error" << ", " << abs(y[m] - truey[m]) << "\n";
	file3 << "Maximum global error" << ", " << norm(y, truey);

	file1.close();
	file2.close();
	file3.close();

}


int main()
{
	// The problem is evaluated in the interval [a,b] with stepsize h = (b-a)/m
	// m is the number of intervals
	// ynot is y0, which is given
	
	double a = 0;
	double b = 10;
	int m = 100;
	lambda = -10;
	ynot = 1;


	AB2(a, b, m);
	RK4(a, b, m);
	AM(a, b, m);
	BDF2(a, b, m);

	return 0;
}