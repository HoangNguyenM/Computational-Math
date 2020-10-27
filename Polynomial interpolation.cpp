// This program includes some polynomial interpolation methods including Barycentric form 1, Barycentric form 2 and Newton form on different types of mesh points
// The program also includes some examples of certain given functions

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <chrono>

using namespace std;

float pi = 3.14159265358979323846;
int n;
int c;
int deg;


// define a function to generate a vector of size "count" following uniform distribution [d1, d2]
vector<float> create(int d1, int d2, int count)
{
	default_random_engine generator;
	uniform_real_distribution<float> distribution(d1, d2);
	vector<float> vec = {};
	for (int j = 0; j < count; j++)
	{
		generator.seed(chrono::system_clock::now().time_since_epoch().count());
		vec.push_back(distribution(generator));
	}
	return vec;
}

// creating Chebyshev mesh of kind 1 and 2, m is the size parameter of the vector
// variable "type" decides Chebyshev kind. type = true means kind 1, type = false means kind 2
vector<float> cheb(int m, bool type)
{
	vector<float> vec = {};

	if (type == true)
	{
		// type 1
		for (int i = 0; i < m; i++)
		{
			float x = cos((2 * i + 1) * pi / (2 * m + 2));
			vec.push_back(x);
		}
	}
	else
	{
		// type 2
		for (int i = 0; i < m; i++)
		{
			float x = cos(i * pi / m);
			vec.push_back(x);
		}
	}

	return vec;

}

// define the functions
double function1(float x)
{
	double y = pow(x - 2, 9);
	return y;
}

double function2(float x, int d)
{
	double y = 1;
	for (int i = 1; i <= d; i++)
	{
		y = y * (x - i);
	}
	return y;
}

double function3(vector<float> vec, float x)
{
	double product1 = 1;
	double product2 = 1;
	for (int i = 0; i < vec.size() - 1; i++)
	{
		product1 = product1 * (x - vec[i]);
		product2 = product2 * (vec[vec.size() - 1] - vec[i]);
	}

	return product1 / product2;
}

double function4(float x)
{
	return (1 / (1 + 25 * pow(x, 2)));
}

// function 5 is divided into two parts and are interpolated separately
double function5_x(float x)
{
	if ((-1 <= x) & (x <= -1 / 3))
	{
		return (- 3 * x / 2 - 1 / 2);
	}
	else
	{
		return 0;
	}
}

double function5_y(float x)
{
	if ((-1 / 3 <= x) & (x <= 1))
	{
		return (3 * x / 2 + 1 / 2);
	}
	else
	{
		return 0;
	}
}

// make a function to choose among the functions of interest
// x is the variable, f takes values 1, 2, 4, d is the degree parameter for certain functions
double choose(float x, int f, int d)
{
	switch (f)
	{
	case 1:
		return function1(x);
	case 2:
		return function2(x, d);
	case 4:
		return function4(x);
	default:
		return 0;
	}
}


// barycentric 1 routine, calculate gammas
vector<double> barycentric1(vector<float> vec)
{
	vector<double> gamma = {};
	for (int i = 0; i < vec.size(); i++)
	{
		// calculate gamma i and save as variable a
		double a = 1;
		for (int j = 0; j < vec.size(); j++)
		{
			if (i != j)
			{
				a = a * (vec[i] - vec[j]);
			}
		}

		a = 1 / a;
		gamma.push_back(a);
	}

	return gamma;
}

// this routine return the values of interpolated polynomial of Barycentric 1 form
// "values" is the vector of values we need to compute the polynomial at, "vec" is the vector of x i
vector<double> barycompute1(vector<float> values, vector<float> vec, vector<double> y)
{
	vector<double> p = {};

	// run through all values we need to evaluate
	for (int i = 0; i < values.size(); i++)
	{
		double omega = 1;
		int track = 0;
		// calculate omega n+1
		for (int j = 0; j < vec.size(); j++)
		{
			omega = omega * (values[i] - vec[j]);
			if (values[i] == vec[j])
			{
				track = j;
			}
		}

		// if the value taken is equal to one of xi, omega is 0 and the return value should be the value of the function
		if (omega == 0)
		{
			p.push_back(y[track]);
		}
		else
		{
			vector<double> gamma = barycentric1(vec);
			double x = 0;
			for (int j = 0; j < vec.size(); j++)
			{
				x = x + (gamma[j] * y[j] / (values[i] - vec[j]));
			}
			x = x * omega;

			p.push_back(x);
		}
	}

	return p;
}


// this routine evaluates beta i depending on the mesh
// variable "type" decides type of the mesh, type = 1 means Chebyshev type 1, type = 2 means Chebyshev type 2, type = 3 means uniform
vector<double> barycentric2(vector<float> vec, int type)
{
	vector<double> beta = {};
	if (type == 1)
	{
		// Chebyshev first kind formula for beta i
		for (int i = 0; i < vec.size(); i++)
		{
			beta.push_back(pow(-1, i) * sin((2 * i + 1) * pi / (2 * n + 2)));
		}
	}
	else
	{
		if (type == 2)
		{
			// Chebyshev second kind formula for beta i
			for (int i = 0; i < vec.size(); i++)
			{
				if ((i == 0) || (i == vec.size() - 1))
				{
					beta.push_back(pow(-1, i) / 2);
				}
				else
				{
					beta.push_back(pow(-1, i));
				}
			}
		}
		else
		{
			// uniform, beta is the same as gamma i from barycentric 1
			beta = barycentric1(vec);
		}
	}

	return beta;
}

// this routine return the values of interpolated polynomial of Barycentric 2 form
// "values" is the vector of values we need to compute the polynomial at, "vec" is the vector of x i
vector<double> barycompute2(vector<float> values, vector<float> vec, vector<double> y, int type)
{
	vector<double> p = {};
	vector<double> beta = barycentric2(vec, type);

	for (int i = 0; i < values.size(); i++)
	{
		// check if the value is equal to one of the x i, if it is, return the value of the function at x i
		bool check = true;
		int track = 0;
		for (int j = 0; j < vec.size(); j++)
		{
			if (values[i] == vec[j])
			{
				check = false;
				track = j;
			}
		}

		if (check == true)
		{
			// calculate the numerator and denominator of Barycentric form 2
			double num = 0;
			double den = 0;
			for (int j = 0; j < vec.size(); j++)
			{
				num = num + (y[j] * beta[j] / (values[i] - vec[j]));
				den = den + (beta[j] / (values[i] - vec[j]));
			}

			p.push_back(num / den);
		}
		else
		{
			p.push_back(y[track]);
		}
	}

	return p;
}

// this routine measures divided differences required for Newton form
vector<double> divided_diff(vector<float> vec, vector<double> y)
{
	vector<double> gamma = barycentric1(vec);
	vector<double> d;
	for (int i = 0; i < vec.size(); i++)
	{
		d.push_back(y[i] * gamma[i]);
	}

	vector<double> fdiff;

	for (int i = 0; i < vec.size(); i++)
	{
		double s = 0;
		for (int j = 0; j < i; j++)
		{
			s = s + d[j];
		}
		fdiff.push_back(s);
	}

	return fdiff;
}

// Adapting Horner's rule to evaluate the polynomial
vector<double> horner(vector<float> values, vector<float> vec, vector<double> y)
{
	vector<double> p = {};
	vector<double> diff = divided_diff(vec, y);
	for (int i = 0; i < values.size(); i++)
	{
		double product = 1;
		double x = 0;
		for (int j = 0; j < vec.size(); j++)
		{
			x = x + diff[j] * product;
			product = product * (values[i] - vec[j]);
		}

		p.push_back(x);
	}

	return p;
}


vector<float> ordering(vector<float> vec, int x)
{
	// increasing order
	if (x == 1)
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

	// decreasing order
	else if (x == 2)
	{
		for (int i = 0; i < vec.size() - 1; i++)
		{
			for (int j = i + 1; j < vec.size(); j++)
			{
				if (vec[j] > vec[i])
				{
					float b = vec[j];
					vec[j] = vec[i];
					vec[i] = b;
				}
			}
		}

		return vec;
	}

	// Leja order
	else
	{
		// find max for the first element
		for (int j = 1; j < vec.size(); j++)
		{
			if (abs(vec[j]) > abs(vec[0]))
			{
				float b = vec[j];
				vec[j] = vec[0];
				vec[0] = b;
			}
		}

		// sort the rest of the vector
		// make a vector p full of one's to keep track of the products
		vector<double> p = {};
		for (int i = 0; i < vec.size(); i++)
		{
			p.push_back(1);
		}

		for (int i = 1; i < vec.size() - 1; i++)
		{
			for (int j = i; j < vec.size(); j++)
			{
				p[j] = p[j] * abs(vec[j] - vec[i - 1]);
			}

			// find the element satisfies Leja order at the current position
			double pmax = p[i];
			int track = i;

			for (int j = i + 1; j < vec.size(); j++)
			{
				if (p[j] > pmax)
				{
					pmax = p[j];
					track = j;
				}
			}

			// swap the elements

			if (track != i)
			{
				// swap x
				float b = vec[i];
				vec[i] = vec[track];
				vec[track] = b;

				// swap p
				double c = p[i];
				p[i] = p[track];
				p[track] = c;
			}

		}

		return vec;
	}
}


// this routine evaluates the maximum of k(x,n,1)
double kxn1(vector<float> values, vector<float> vec)
{
	vector<double> contnumber = {};

	for (int i = 0; i < values.size(); i++)
	{
		// evaluate l i
		vector<double> l = {};
		bool check = false;

		// first check if the value evaluated is one of the x i, if it is, return 1 as the condition number
		for (int j = 0; j < vec.size(); j++)
		{
			if (values[i] == vec[j])
			{
				check = true;
			}
		}

		if (check == true)
		{
			contnumber.push_back(1);
		}
		else
		{
			vector<double> gamma = barycentric1(vec);
			double product = 1;

			// evaluate the product (x - x i)
			for (int j = 0; j < vec.size(); j++)
			{
				product = product * abs(values[i] - vec[j]);
			}

			// calculate l i
			for (int j = 0; j < vec.size(); j++)
			{
				l.push_back(gamma[j] * product / (values[i] - vec[j]));
			}

			double sum = 0;
			for (int j = 0; j < l.size(); j++)
			{
				sum = sum + abs(l[j]);
			}

			contnumber.push_back(sum);
		}
	}

	// find the maximum of kxn1
	double max = 0;
	for (int i = 0; i < contnumber.size(); i++)
	{
		if (contnumber[i] > max)
		{
			max = contnumber[i];
		}
	}
	
	return max;
}

// this routine evaluates k(x,n,y)
// p is the vector of the values of the interpolated polynomial with respect to the points of interest in vector "values"
double kxny(vector<float> values, vector<float> vec, vector<double> y, vector<double> p)
{
	vector<double> contnumber = {};

	for (int i = 0; i < values.size(); i++)
	{
		// evaluate l i
		vector<double> l = {};
		bool check = false;
		int track = 0;

		// first check if the value evaluated is one of the x i, if it is, return 1 as the condition number
		for (int j = 0; j < vec.size(); j++)
		{
			if (values[i] == vec[j])
			{
				check = true;
				track = j;
			}
		}

		if (check == true)
		{
			for (int j = 0; j < vec.size(); j++)
			{
				if (track == j)
				{
					l.push_back(1);
				}
				else
				{
					l.push_back(0);
				}
			}
		}
		else
		{
			vector<double> gamma = barycentric1(vec);
			double product = 1;

			// evaluate the product (x - x i)
			for (int j = 0; j < vec.size(); j++)
			{
				product = product * abs(values[i] - vec[j]);
			}

			// calculate l i
			for (int j = 0; j < vec.size(); j++)
			{
				l.push_back(gamma[j] * product / (values[i] - vec[j]));
			}
		}

		// after obtaining l i, we evaluate the conditioning number
		double sum = 0;
		for (int j = 0; j < vec.size(); j++)
		{
			sum = sum + abs(l[j] * y[j]);
		}

		sum = sum / abs(p[i]);
		contnumber.push_back(sum);
	}

	// find the maximum of kxny
	double max = 0;
	for (int i = 0; i < contnumber.size(); i++)
	{
		if (contnumber[i] > max)
		{
			max = contnumber[i];
		}
	}

	return max;
}

// create a perturbation of y
vector<double> perturb(vector<double> y)
{
	vector<double> pert = {};
	vector<float> p = create(-1, 1, y.size());
	for (int i = 0; i < y.size(); i++)
	{
		pert.push_back(y[i] + abs(y[i]) * p[i] * pow(10, -15));
	}

	return pert;
}

// calculate the max norm of r(x), here vec1 and vec2 have the same size
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


int main()
{

	ofstream file1;
	ofstream file2;
	ofstream file3;
	ofstream file4;
	ofstream file5;
	ofstream file6;
	ofstream file7;
	ofstream file8;
	ofstream file9;
	ofstream file10;

	// n is the amount of mesh points
	// deg is the degree parameter used in certain functions
	// c is the count of times we evaluate the interpolated functions
	n = 20;
	deg = 30;
	c = 30;


	file1.open("input.txt", ios::out);
	file2.open("functionvalue.txt", ios::out);
	file3.open("barycentric1.txt", ios::out);
	file4.open("barycentric2.txt", ios::out);
	file5.open("newton1.txt", ios::out);
	file6.open("newton2.txt", ios::out);
	file7.open("newton3.txt", ios::out);
	file8.open("conditioning.txt", ios::out);
	file9.open("error.txt", ios::out);
	file10.open("convergence.txt", ios::out);


	// generate a vector of points
	// "create" function generates uniform mesh, "cheb" function generates chebyshev points
	vector<float> vect = create(-1, 1, n);
	vector<double> f = {};
	vector<double> fpert = {};
	
	// evaluate the function values
	// change the value in the "choose" function below to fit the function (1, 2, or 4)
	// manually replace "choose" function by "function3(vect, vect[i])" for function 3
	// manually replace "choose" function by "function5_x(vect[i])" or "function5_y(vect[i])" for parametric function
	for (int i = 0; i < n; i++)
	{
		f.push_back(choose(vect[i], 1, deg));
	}

	// generate a vector of points to calculate the values of interpolated polynomial
	// use "value = cheb(c, true)" for function 5 to maintain the "value" vector for both x(z) and y(z)
	vector<float> value = create(-1, 1, c);
	vector<double> funcvalue = {};

	// print out the value of the function
	// change the value in the "choose" function below to fit the function (1, 2 or 4)
	// manually replace "choose" function by "function3(vect, value[i])" for function 3
	// manually replace "choose" function by "function5_x(value[i])" or "function5_y(value[i])" for parametric function
	for (int i = 0; i < c; i++)
	{
		file1 << value[i] << "\n";
		funcvalue.push_back(choose(value[i], 1, deg));
		file2 << funcvalue[i] << "\n";
	}

	// print out kxn1
	file8 << "kxn1 is " << kxn1(value, vect) << "\n";

	// measuring barycentric form 1
	vector<double> p = barycompute1(value, vect, f);
	vector<double> ppert = {};
	file10 << "barycentric 1 convergence norm is " << norm(p, funcvalue) << "\n";


	for (int i = 0; i < c; i++)
	{
		file3 << p[i] << "\n";
	}

	// print kxny for barycentric 1
	file8 << "barycentric 1 kxny is " << kxny(value, vect, f, p) << "\n";

	// create 100 perturbations of y to measure the norm of r(x)
	file9 << "					barycentric 1 r(x)\n";
	for (int i = 0; i < 100; i++)
	{
		fpert = perturb(f);
		ppert = barycompute1(value, vect, fpert);
		double r = norm(p, ppert);
		file9 << r << "\n";
	}

	// measuring barycentric form 2
	p = barycompute2(value, vect, f, 3);
	file10 << "barycentric 2 convergence norm is " << norm(p, funcvalue) << "\n";

	for (int i = 0; i < c; i++)
	{
		file4 << p[i] << "\n";
	}

	// print kxny for barycentric 2
	file8 << "barycentric 2 kxny is " << kxny(value, vect, f, p) << "\n";

	// create 100 perturbations of y to measure the norm of r(x)
	file9 << "					barycentric 2 r(x)\n";
	for (int i = 0; i < 100; i++)
	{
		fpert = perturb(f);
		ppert = barycompute2(value, vect, fpert, 3);
		double r = norm(p, ppert);
		file9 << r << "\n";
	}

	//  measuring newton form in increasing order mesh

	vect = ordering(vect, 1);

	p = horner(value, vect, f);
	file10 << "newton form - increasing - convergence norm is " << norm(p, funcvalue) << "\n";

	for (int i = 0; i < c; i++)
	{
		file5 << p[i] << "\n";
	}

	// print kxny for newton form - increasing order
	file8 << "newton - increasing - kxny is " << kxny(value, vect, f, p) << "\n";

	// create 100 perturbations of y to measure the norm of r(x)
	file9 << "					newton - increasing - r(x)\n";
	for (int i = 0; i < 100; i++)
	{
		fpert = perturb(f);
		ppert = horner(value, vect, fpert);
		double r = norm(p, ppert);
		file9 << r << "\n";
	}

	//  measuring newton form in decreasing order mesh

	vect = ordering(vect, 2);

	p = horner(value, vect, f);
	file10 << "newton form - decreasing - convergence norm is " << norm(p, funcvalue) << "\n";

	for (int i = 0; i < c; i++)
	{
		file6 << p[i] << "\n";
	}

	// print kxny for newton form - decreasing order
	file8 << "newton - decreasing - kxny is " << kxny(value, vect, f, p) << "\n";

	// create 100 perturbations of y to measure the norm of r(x)
	file9 << "					newton - decreasing - r(x)\n";
	for (int i = 0; i < 100; i++)
	{
		fpert = perturb(f);
		ppert = horner(value, vect, fpert);
		double r = norm(p, ppert);
		file9 << r << "\n";
	}

	//  measuring newton form in leja order mesh

	vect = ordering(vect, 3);

	p = horner(value, vect, f);
	file10 << "newton form - leja - convergence norm is " << norm(p, funcvalue) << "\n";

	for (int i = 0; i < c; i++)
	{
		file7 << p[i] << "\n";
	}

	// print kxny for newton form - leja order
	file8 << "newton - leja - kxny is " << kxny(value, vect, f, p) << "\n";

	// create 100 perturbations of y to measure the norm of r(x)
	file9 << "					newton - leja - r(x)\n";
	for (int i = 0; i < 100; i++)
	{
		fpert = perturb(f);
		ppert = horner(value, vect, fpert);
		double r = norm(p, ppert);
		file9 << r << "\n";
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

	return 0;
}
