#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <chrono>

// for this program, please change only the variables in the main program to perform the task
// first run the matrixcreation line to create the matrix, then solve the system in another language
// after that, disable the matrixcreation line and enable the calculation line to obtain the interpolated values
// the function can be changed in the function "function" below

using namespace std;

// n + 1 is the number of mesh points
// c + 1 is the number of points we want to measure the interpolation values
int n;
int c;

// define a function to generate a vector of size "count" following uniform distribution [d1, d2]
vector<double> create(int d1, int d2, int count)
{
	default_random_engine generator;
	uniform_real_distribution<double> distribution(d1, d2);
	vector<double> vec = {};
	for (int j = 0; j < count; j++)
	{
		generator.seed(chrono::system_clock::now().time_since_epoch().count());
		vec.push_back(distribution(generator));
	}
	return vec;
}

// The function can be changed
double function(double x)
{
	//return (1 / (1 + 25 * pow(x, 2)));
	return (x - 2) * (x + 2) * (x - 4) * (x + 4);
}

// order the mesh points in increasing order
vector<double> ordering(vector<double> vec)
{
	for (int i = 0; i < vec.size() - 1; i++)
	{
		for (int j = i + 1; j < vec.size(); j++)
		{
			if (vec[j] < vec[i])
			{
				float b = vec[j];
				vec[j] = vec[i];
				vec[i] = b;
			}
		}
	}

	return vec;
}

// output the system of equations, second derivative of s0 and sn are specified
void matrix1(vector<double> h, vector<double> func, double s0, double sn)
{
	ofstream file3;
	ofstream file4;
	file3.open("matrix.csv", ios::out);
	file4.open("d.csv", ios::out);

	vector<double> mu = {};
	vector<double> lambda = {};
	vector<double> d = {};

	// push back the first element of all the above vectors to match the enumeration of h
	mu.push_back(0);
	lambda.push_back(0);
	d.push_back(0);

	// calculate values of mu, lambda and d from i = 1 to i = n - 1
	for (int i = 1; i < (func.size() - 1); i++)
	{
		mu.push_back(h[i] / (h[i] + h[i + 1]));
		lambda.push_back(1 - mu[i]);
		double temp = 6 * ((func[i + 1] - func[i]) / h[i + 1] - (func[i] - func[i - 1]) / h[i]) / (h[i] + h[i + 1]);

		// check for two special case i = 1 and i = n - 1 (func.size is n + 1)
		if (i == 1)
		{
			temp = temp - mu[i] * s0;
		}

		if (i == func.size() - 2)
		{
			temp = temp - lambda[i] * sn;
		}

		d.push_back(temp);
		file4 << d[i] << "\n";
	}

	// print out the matrix (n - 1) x (n - 1)
	// track1 is the number of 0 at the beginning of a line
	// track2 is the number of 0 at the end of a line
	for (int i = 1; i < (func.size() - 1); i++)
	{
		int track1 = i - 2;
		int track2 = func.size() - i - 3;
		while (track1 > 0)
		{
			file3 << "0,";
			track1 = track1 - 1;
		}

		if (i == 1)
		{
			file3 << "2," << lambda[i] << ", ";
		}
		else
		{
			if (i == func.size() - 2)
			{
				file3 << mu[i] << ", " << "2,";
			}
			else
			{
				file3 << mu[i] << ", " << "2," << lambda[i] << ", ";
			}
		}

		while (track2 > 0)
		{
			file3 << "0,";
			track2 = track2 - 1;
		}

		file3 << "\n";
	}

	file3.close();
	file4.close();
}

// output the system of equations, first derivative of s0 and sn are specified
void matrix2(vector<double> h, vector<double> func, double s0, double sn)
{
	ofstream file3;
	ofstream file4;
	file3.open("matrix.csv", ios::out);
	file4.open("d.csv", ios::out);

	vector<double> mu = {};
	vector<double> lambda = {};
	vector<double> d = {};

	// push back the first element of all the above vectors to match the enumeration of h
	mu.push_back(0);
	lambda.push_back(0);

	// add the additional equation for first derivative bound
	d.push_back(s0 - (func[1] - func[0]) / h[1]);

	// calculate values of mu, lambda and d from i = 1 to i = n - 1
	for (int i = 1; i < (func.size() - 1); i++)
	{
		mu.push_back(h[i] / (h[i] + h[i + 1]));
		lambda.push_back(1 - mu[i]);
		double temp = 6 * ((func[i + 1] - func[i]) / h[i + 1] - (func[i] - func[i - 1]) / h[i]) / (h[i] + h[i + 1]);

		d.push_back(temp);
	}

	// add the additional equation for first derivative bound
	d.push_back(sn - (func[func.size() - 1] - func[func.size() - 2]) / h[func.size() - 1]);

	// print out d from i = 0 to i = n
	for (int i = 0; i < func.size(); i++)
	{
		file4 << d[i] << "\n";
	}

	// print out the matrix (n + 1) x (n + 1)
	// track1 is the number of 0 at the beginning of a line
	// track2 is the number of 0 at the end of a line
	for (int i = 0; i < func.size(); i++)
	{
		int track1 = i - 1;
		int track2 = func.size() - i - 2;
		while (track1 > 0)
		{
			file3 << "0,";
			track1 = track1 - 1;
		}

		if (i == 0)
		{
			file3 << -h[1] / 3 << ", " << -h[1] / 6 << ", ";
		}
		else
		{
			if (i == (func.size() - 1))
			{
				file3 << h[n] / 6 << ", " << h[n] / 3 << ", ";
			}
			else
			{
				file3 << mu[i] << ", " << "2," << lambda[i] << ", ";
			}
		}

		while (track2 > 0)
		{
			file3 << "0,";
			track2 = track2 - 1;
		}

		file3 << "\n";
	}

	file3.close();
	file4.close();
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

// this function divide the intervals of the mesh points
vector<double> divide_h(vector<double> vec)
{
	vector<double> vect = {};
	vect.push_back(vec[0]);
	for (int i = 1; i < vec.size(); i++)
	{
		vect.push_back((vec[i] + vec[i - 1]) / 2);
		vect.push_back(vec[i]);
	}

	return vect;
}

void matrixcreation(bool m, bool divide, double s0, double sn, bool task2)
{
	ofstream file1;
	ofstream file2;

	vector<double> vect = {};

	if (task2 == false)
	{
		// check to create matrix for h or h/2
		if (divide == false)
		{
			// create a vector of random mesh points
			vect = create(-5, 5, n + 1);
			vect = ordering(vect);
		}
		else
		{
			// perform the following code only for h/2
			ifstream file3("input1.csv");

			// read the already created mesh points

			vector<double> vectemp = {};
			for (int i = 0; i <= n; i++)
			{
				double temp = 0;
				file3 >> temp;
				vectemp.push_back(temp);
			}

			// create the midpoint in each interval
			vect = divide_h(vectemp);

			file3.close();
		}
	}
	else
	{
		// do this for task 2
		ifstream file4("input1.csv");

		// read the already created mesh points

		vector<double> vectemp = {};
		for (int i = 0; i <= n; i++)
		{
			double temp = 0;
			file4 >> temp;
			vectemp.push_back(temp);
		}

		vect = vectemp;
		file4.close();

	}


	file1.open("input.csv", ios::out);
	file2.open("functionvalue.csv", ios::out);
	vector<double> funcvalue = {};

	if (task2 == false)
	{
		// evaluate the function values f

		for (int i = 0; i < vect.size(); i++)
		{
			file1 << vect[i] << "\n";
			funcvalue.push_back(function(vect[i]));
			file2 << funcvalue[i] << "\n";
		}
	}
	else
	{
		// do this for task 2
		ifstream file5("functionvalue1.csv");

		// read the already created function values

		for (int i = 0; i <= n; i++)
		{
			double temp = 0;
			file5 >> temp;
			funcvalue.push_back(temp);
		}

		file5.close();
	}


	// make vector h i, the length of the intervals between the mesh points
	vector<double> h = {};
	h.push_back(0);

	for (int i = 1; i < vect.size(); i++)
	{
		h.push_back(vect[i] - vect[i - 1]);
	}


	// create the system of equations to solve
	if (m == true)
	{
		matrix1(h, funcvalue, s0, sn);
	}
	else
	{
		matrix2(h, funcvalue, s0, sn);
	}


	file1.close();
	file2.close();
}


void calculation(bool m, bool task2)
{

	// import the result from matlab
	ifstream file1("input.csv");
	ifstream file2("functionvalue.csv");
	ifstream file3("s.csv");

	// this file is to write the interpolated values
	ofstream file4;
	ofstream file5;
	ofstream file6;
	ofstream file7;
	ofstream file8;
	ofstream file9;
	ofstream file10;
	ofstream file11;

	file4.open("values.csv", ios::out);
	file5.open("interpolated values.csv", ios::out);
	file6.open("true values.csv", ios::out);
	file7.open("error.csv", ios::out);
	file8.open("first der.csv", ios::out);
	file9.open("true first der.csv", ios::out);
	file10.open("second der.csv", ios::out);
	file11.open("true second der.csv", ios::out);


	vector<double> vect = {};
	vector<double> func = {};
	vector<double> s = {};



	// read the mesh points and function values

	for (int i = 0; i <= n; i++)
	{
		double temp1 = 0;
		double temp2 = 0;
		file1 >> temp1;
		file2 >> temp2;
		vect.push_back(temp1);
		func.push_back(temp2);
	}


	// read solution from matlab
	// divide into two situations: if m = True, second derivatives are specified, if m = False, first derivatives are specified
	if (m == true)
	{
		s.push_back(0);
		for (int i = 1; i < (vect.size() - 1); i++)
		{
			double temp = 0;
			file3 >> temp;
			s.push_back(temp);
		}
		s.push_back(0);
	}
	else
	{
		for (int i = 0; i < (vect.size()); i++)
		{
			double temp = 0;
			file3 >> temp;
			s.push_back(temp);
		}
	}


	// make vector h i, the length of the intervals between the mesh points
	vector<double> h = {};
	h.push_back(0);

	for (int i = 1; i < vect.size(); i++)
	{
		h.push_back(vect[i] - vect[i - 1]);
	}

	// calculate vector gamma
	vector<double> gamma = {};

	for (int i = 1; i < vect.size(); i++)
	{
		gamma.push_back((func[i] - func[i - 1]) / h[i] - h[i] * (s[i] - s[i - 1]) / 6);
	}

	// calculate vector tiltgamma
	vector<double> tiltgamma = {};

	for (int i = 1; i < vect.size(); i++)
	{
		tiltgamma.push_back(func[i - 1] - s[i - 1] * pow(h[i], 2) / 6);
	}

	// create a vector of variables we want to measure the interpolation values
	vector<double> values = {};
	for (int i = 0; i <= c; i++)
	{
		values.push_back(vect[0] + i * (vect[vect.size() - 1] - vect[0]) / c);
	}

	// create a vector of interpolated values and true values
	vector<double> p = {};
	vector<double> truevalues = {};
	vector<double> firstder = {};
	vector<double> secondder = {};
	vector<double> truefirst = {};
	vector<double> truesecond = {};
	

	for (int i = 0; i <= c; i++)
	{
		// find the interval in which the value is
		int track = 1;
		for (int j = 1; j < vect.size(); j++)
		{
			if ((values[i] >= vect[j - 1]) & (values[i] <= vect[j]))
			{
				track = j;
			}
		}

		// calculate the interpolated values
		double temp = s[track - 1] * pow((vect[track] - values[i]), 3) / (6 * h[track]);
		temp = temp + s[track] * pow((values[i] - vect[track - 1]), 3) / (6 * h[track]);
		temp = temp + gamma[track - 1] * (values[i] - vect[track - 1]) + tiltgamma[track - 1];
		p.push_back(temp);
		truevalues.push_back(function(values[i]));

		// calculate the first order derivatives
		temp = - s[track - 1] * pow((vect[track] - values[i]), 2) / (2 * h[track]);
		temp = temp + s[track] * pow((values[i] - vect[track - 1]), 2) / (2 * h[track]);
		temp = temp + gamma[track - 1];
		firstder.push_back(temp);
		truefirst.push_back(4 * pow(values[i], 3) - 40 * values[i]);

		// calculate the second order derivatives
		temp = s[track - 1] * (vect[track] - values[i]) / h[track];
		temp = temp + s[track] * (values[i] - vect[track - 1]) / h[track];
		secondder.push_back(temp);
		truesecond.push_back(12 * pow(values[i], 2) - 40);

		file4 << values[i] << "\n";
		file5 << p[i] << "\n";
		file6 << truevalues[i] << "\n";
		file8 << firstder[i] << "\n";
		file9 << truefirst[i] << "\n";
		file10 << secondder[i] << "\n";
		file11 << truesecond[i] << "\n";
	}

	// calculate the norm f(x) - s(x)
	double error = norm(p, truevalues);
	file7 << error;

	// calculate the function values for task 2
	if (task2 == true)
	{
		ofstream file12;
		ofstream file13;

		file12.open("task2 - f(t).csv", ios::out);
		file13.open("task2 - G(t).csv", ios::out);

		for (int i = 0; i <= c; i++)
		{
			// calculate the values of f(t)
			file12 << p[i] + values[i] * firstder[i] << "\n";

			// calculate the values of G(t)
			file13 << exp(-values[i] * p[i]) << "\n";
		}


		file12.close();
		file13.close();
	}

	file1.close();
	file2.close();
	file3.close();
	file4.close();
	file5.close();
	file6.close();
	file7.close();
	file8.close();
	file9.close();
	file10.close();
	file11.close();
}


int main()
{
	// c could be changed to increase the number of values
	// n could be changed. The larger n is, the more h converge to 0
	// bound could be changed, if bound = true, second derivatives are specified, if bound = false, first derivatives are specified
	// divide could be changed, divide = false means calculation for h, divide = true means calculation for h/2
	// task2 = true to do task2, task2 = false otherwise

	c = 100;
	n = 20;
	bool bound = true;
	bool divide = false;
	bool task2 = false;
	double s0 = 0;
	double sn = 0;

	// perform the following code to create the mesh points and the system of equations
	matrixcreation(bound, divide, s0, sn, task2);

	if (divide == true)
	{
		n = 2 * n;
	}

	// perform the following code after getting the result from Matlab
	// run the below code or the matrixcreation code only one at a time
	//calculation(bound, task2);


	return 0;
}

