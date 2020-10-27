// This program contains C++ code on floating point system representation
// including constructing mathematical operators on floating points such as translation from real numbers to floating points, rounding, addition, accumulation
// perform error tests on the above operations and error bound

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <random>
#include <fstream>
#include <chrono>

using namespace std;

// precision is the variable for precision, i is the variable for vector size
// we will also use countmax to represent the number of times a certain program runs
// you can increase i and countmax for a problem to obtain a larger amount of observations
// But please do not use too large number to avoid overflow
// Also, the program might take multiple times to achieve the final results
// The run time should not exceed 2 seconds if the program successfully runs.

int precision, i;

// check the sign of a certain number
int sign_check(double x)
{
	if (x >= 0)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

// define floating point number class
class floating
{
	public:
		int mantissa;
		int exp;
		int remainder;
};

// define a function to generate a vector of type "double" following uniform distribution [d1, d2]
vector<double> generate(int d1, int d2)
{
	default_random_engine generator;
	uniform_real_distribution<double> distribution(d1, d2);
	vector<double> vec;
	for (int j = 0; j < i; j++)
	{
		generator.seed(chrono::system_clock::now().time_since_epoch().count());
		vec.push_back(distribution(generator));
	}
	return vec;
}

// routine 1
floating round_to_near(int m, int e, int remainder)
{
	if (remainder >= 5) {
		m = m + 1;
		if (m == pow(10, precision)) {
			m = m / 10;
			e = e + 1;
		}
	}
	floating y;
	y.mantissa = m;
	y.exp = e;
	y.remainder = 0;
	return y;
}

// routine 2
floating translate(double x, bool round)
{
	int sign = sign_check(x);
	double a = abs(x);
	int r = 0;
	int e = 0;
	floating y;
	if (a >= pow(10, precision))
	{
		int m = floor(a);
		while (m >= pow(10, precision))
		{
			r = m % 10;
			m = m / 10;
			e = e + 1;
		}

		if (round == true)
		{
			y = round_to_near(m, e, r);
			y.mantissa = sign * y.mantissa;
			return y;
		}
		else
		{
			m = m * sign;
			y.mantissa = m;
			y.exp = e;
			y.remainder = r;
			return y;
		}
	}
	else
	{
		while (a < pow(10, precision - 1))
		{
			a = a * 10;
			e = e - 1;
		}
		r = floor(10 * (a - floor(a)));
		int m = floor(a);
		if (round == true)
		{
			y = round_to_near(m, e, r);
			y.mantissa = y.mantissa * sign;
			return y;
		}
		else
		{
			y.mantissa = sign * m;
			y.exp = e;
			y.remainder = r;
			return y;
		}
	}
}

// routine 3
floating addition(floating x, floating y, bool round)
{
	long int a, b, c;
	int m, e;
	int	r = 0;
	int shift = x.exp - y.exp;
	if (shift > 0)
	{
		a = x.mantissa * pow(10, shift);
		b = y.mantissa;
		e = y.exp;
	}
	else
	{
		a = y.mantissa * pow(10, -shift);
		b = x.mantissa;
		e = x.exp;
	}
	c = a + b;
	int sign = sign_check(c);
	c = abs(c);
	if (c >= pow(10, precision))
	{
		while (c >= pow(10, precision))
		{
			r = c % 10;
			c = c / 10;
			e = e + 1;
		}
	}
	else
	{
		while (c < pow(10, precision - 1))
		{
			c = c * 10;
			e = e - 1;
		}
	}
	m = static_cast<int>(c);
	floating z;

	if (round == true)
	{
		z = round_to_near(m, e, r);
		z.mantissa = z.mantissa * sign;
	}
	else
	{
		z.mantissa = sign * m;
		z.exp = e;
		z.remainder = r;
	}

	return z;
}

// routine 4
floating accumulation(vector<floating> vect, bool round)
{
	floating s = vect[0];
	int k = 1;
	while (k < vect.size())
	{
		s = addition(s, vect[k], round);
		k = k + 1;
	}
	return s;
}


int main()
{

	default_random_engine generator;
	uniform_real_distribution<double> distribution(-1, 1);


	// input precision
	precision = 5;

	ofstream file;

	// input size of vector as i
	i = 100;

	// correctness 1: validate the error of translate routine, rounding strategy is round towards 0
	vector<double> vect = generate(-99999, 99999);

	file.open("correctness1.txt", ios::out);

	for (int j = 0; j < i; j++)
	{
		// check for the condition of the problem
		if (vect[j] > pow(10, precision - 9 - 1))
		{
			floating z = translate(vect[j], false);
			double error = abs((z.mantissa * pow(10, z.exp) - vect[j]) / vect[j]);

			file << error << "\n";
		}
	}

	// export the bound
	file << "bound is " << pow(10, 1 - precision);
	file.close();

	i = 100;
	// correctness 2: validate the bound of addition model, rounding strategy is round towards 0
	vector<double> vect1 = generate(-99999, 99999);
	vector<double> vect2 = generate(-99999, 99999);

	file.open("correctness2.txt", ios::out);

	for (int j = 0; j < i; j++)
	{
		double x = vect1[j] + vect2[j];
		floating y = addition(translate(vect1[j], false), translate(vect2[j], false), false);
		double error = abs((x - y.mantissa * pow(10, y.exp)) / x);

		// "error" is the relative error
		file << error << "\n";
	}

	// export the bound
	file << "bound is " << pow(10, 1 - precision);
	file.close();

	i = 100;
	// correctness 3: Investigate the error distribution of two rounding strategies

	// regenerate two vectors from uniform distribution [1.0, 2.0]
	vect1 = generate(1, 2);
	vect2 = generate(1, 2);

	// Round to nearest, the boolean variable input for translate and addition routines is true
	file.open("correctness3_1.txt", ios::out);

	for (int j = 0; j < i; j++)
	{
		double x = vect1[j] + vect2[j];
		floating y = addition(translate(vect1[j], true), translate(vect2[j], true), true);
		double error = (x - y.mantissa * pow(10, y.exp)) / x;

		// "error" is the relative error
		file << error << "\n";
	}

	// export the bound
	file << "bound is " << pow(10, 1 - precision);
	file.close();

	// Round towards 0, the boolean variable input for translate and addition routines is false

	file.open("correctness3_2.txt", ios::out);

	for (int j = 0; j < i; j++)
	{
		double x = vect1[j] + vect2[j];
		floating y = addition(translate(vect1[j], false), translate(vect2[j], false), false);
		double error = (x - y.mantissa * pow(10, y.exp)) / x;

		// "error" is the relative error
		file << error << "\n";
	}

	// export the bound
	file << "bound is " << pow(10, 1 - precision);
	file.close();


	i = 100;
	// Accumulation 1

	file.open("accumulation1.txt", ios::out);

	// run the program countmax times, you could increase countmax
	int count = 0;
	int countmax = 10;

	while (count < countmax)
	{
		vect = generate(-99999, 99999);
		double sum = 0;

		vector<floating> vect3;

		// make a vector "vect3" of floating point representations of numbers in vector "vect"
		for (int j = 0; j < i; j++)
		{
			sum = sum + vect[j];
			vect3.push_back(translate(vect[j], false));
		}

		floating floating_sum = accumulation(vect3, false);

		double error = abs(floating_sum.mantissa * pow(10, floating_sum.exp) - sum);

		// calculate the bound
		double bound = 0;

		// first calculate the norm of vector "vect"
		for (int j = 0; j < i; j++)
		{
			bound = bound + abs(vect[j]);
		}

		// obtain the bound by multiplying the norm by (n-1) and unit roundoff
		bound = bound * (i - 1) * pow(10, 1 - precision);

		// output the ratio error/bound to make comparison
		file << error / bound << "\n";
		count = count + 1;
	}
	
	file.close();

	// Accumulation 2

	file.open("accumulation2.txt", ios::out);
	
	// run the program countmax times, you could increase countmax
	count = 0;
	countmax = 10;
	i = 100;

	// increase precision to match IEEE format
	precision = 10;

	// these variables are used to calculate the relative error of the sequence of 100 numbers
	double exact_sum = 0;
	floating float_sum;
	float_sum.mantissa = 0;
	float_sum.exp = 0;
	float_sum.remainder = 0;

	
	double k_rel = 0;
	vect = generate(-99999, 99999);

	// save the sequence, calculate the exact sum and the sum of the floating point representations
	for (int j = 0; j < i; j++)
	{
		file << vect[j] << "\n";
		exact_sum = exact_sum + vect[j];
		float_sum = addition(float_sum, translate(vect[j], false), false);
	}

	while (count < countmax)
	{
		// these variables are used to calculate c_rel
		double normx = 0;
		double normp = 0;
		double sumx = 0;
		double sump = 0;

		// generate a vector with similar size as vector "vect" for perturbations
		vector<double> perturb = generate(-1, 1);

		for (int j = 0; j < i; j++)
		{
			perturb[j] = perturb[j] * pow(10, 1 - precision) * abs(vect[j]);
			normx = normx + abs(vect[j]);
			normp = normp + abs(perturb[j]);
			sumx = sumx + vect[j];
			sump = sump + perturb[j];
		}

		// calculating relative conditioning number by taking the maximum over all perturbations
		double c_rel = abs(sump) * normx / (abs(sumx) * normp);
		k_rel = max(k_rel, c_rel);
		count = count + 1;
	}
	
	// export k_rel, calculate the relative error in the exact sum of the sequence of 100 numbers
	file << "k_rel is " << k_rel << "\n";
	file << "error is " << abs((exact_sum - float_sum.mantissa * pow(10, float_sum.exp)) / exact_sum);
	file.close();

	// Accumulation 3

	precision = 5;

	ofstream file1;

	file.open("accumulation3_1.txt", ios::out);
	file1.open("accumulation3_2.txt", ios::out);
	
	vector<floating> vect3;
	
	count = 0;
	countmax = 10;
	i = 8;

	// This vector is taken from the problem
	vect = { 5, 5, 5, 5, -5, -5, -5, -5 };

	// run the program countmax times, you could increase countmax
	while (count < countmax)
	{
		count = count + 1;
		double sum = 0;
		vector<double> perturb = generate(-1, 1);	

		for (int j = 0; j < i; j++)
		{
			perturb[j] = perturb[j] * 50 * pow(10, 1 - precision) + vect[j];
			sum = sum + perturb[j];

			// translate new perturbation vector to floating point representations in vect3
			// notice that the rounding strategy here does not matter because of randomness
			vect3.push_back(translate(perturb[j], true));
		}

		// round to nearest absolute error
		floating sum1 = accumulation(vect3, true);
		double error1 = abs(sum1.mantissa * pow(10, sum1.exp) - sum);

		// round to 0 absolute error
		floating sum2 = accumulation(vect3, false);
		double error2 = abs(sum2.mantissa * pow(10, sum2.exp) - sum);

		// "error1" is the error for round to nearest, "error2" is the error for round towards 0
		file << error1 << "\n";
		file1 << error2 << "\n";
	
	}
	
	file.close();
	file1.close();

	// Accumulation 4

	file.open("accumulation4_1.txt", ios::out);
	file1.open("accumulation4_2.txt", ios::out);
	
	double sum1, sum2, sumpositive, sumnegative;
	floating intermediate_sum1, intermediate_sum2;
	count = 0;
	countmax = 5;

	// i = 8 since the vectors have 8 elements
	i = 8;

	// these vectors are taken from the problem
	vect1 = { 5, 5, 5, 5, -5, -5, -5, -5 };
	vect2 = { 5, -5, 5, -5, 5, -5, 5, -5 };

	// run the program countmax times for vect1, you could increase countmax
	while (count < countmax)
	{

		// generate new perturbations and reset all variables at the beginning of the run
		vector<double> perturb = generate(-1, 1);

		intermediate_sum1.mantissa = 0;
		intermediate_sum1.exp = 0;
		intermediate_sum1.remainder = 0;
		intermediate_sum2.mantissa = 0;
		intermediate_sum2.exp = 0;
		intermediate_sum2.remainder = 0;

		sum1 = 0;
		sum2 = 0;
		sumpositive = 0;
		sumnegative = 0;

		for (int j = 0; j < i; j++)
		{
			// multiply the perturbations generated with the bound to achieve suitable perturbations
			perturb[j] = perturb[j] * 50 * pow(10, 1 - precision) + vect1[j];
			sum1 = sum1 + perturb[j];
			floating p = translate(perturb[j], true);
	
			// calculating round to nearest intermediate sum, sum2 is used to calculate round to nearest bound
			intermediate_sum1 = addition(intermediate_sum1, p, true);
			sum2 = sum2 + abs(intermediate_sum1.mantissa * pow(10, intermediate_sum1.exp));

			// calculating round to 0 intermediate sum, sumpositive and sumnegative are used to calculate round to 0 bound
			intermediate_sum2 = addition(intermediate_sum2, p, false);
			if (intermediate_sum2.mantissa >= 0)
			{
				sumpositive = sumpositive + intermediate_sum2.mantissa * pow(10, intermediate_sum2.exp);
			}
			else
			{
				sumnegative = sumnegative - intermediate_sum2.mantissa * pow(10, intermediate_sum2.exp);
			}
		}

		// calculating round to nearest bound as bound1, round to 0 bound as bound2;
		double bound1 = sum2 * pow(10, 1 - precision) / 2;
		double bound2 = pow(10, 1 - precision);
		if (sumpositive > sumnegative)
		{
			bound2 = bound2 * sumpositive;
		}
		else
		{
			bound2 = bound2 * sumnegative;
		}
		
		// calculating running error for round to nearest as error1 and round to 0 as error2
		double error1 = abs(sum1 - intermediate_sum1.mantissa * pow(10, intermediate_sum1.exp));
		double error2 = abs(sum1 - intermediate_sum2.mantissa * pow(10, intermediate_sum2.exp));

		// write running error ratio of round to nearest to "accumulation4_1.txt" for vect1
		file << error1 / bound1 << "\n";

		// write running error ratio of round to 0 to "accumulation4_2.txt" for vect1
		file1 << error2 / bound2 << "\n";
		
		count = count + 1;
	}
	
	file.close();
	file1.close();

	file.open("accumulation4_3.txt", ios::out);
	file1.open("accumulation4_4.txt", ios::out);
	count = 0;
	countmax = 5;

	// run the program countmax times for vect2, you could increase countmax
	// the operations are exactly the same as for vect1
	while (count < countmax)
	{
		// create new perturbations and reset all variables
		vector<double> perturb = generate(-1, 1);

		intermediate_sum1.mantissa = 0;
		intermediate_sum1.exp = 0;
		intermediate_sum1.remainder = 0;
		intermediate_sum2.mantissa = 0;
		intermediate_sum2.exp = 0;
		intermediate_sum2.remainder = 0;

		sum1 = 0;
		sum2 = 0;
		sumpositive = 0;
		sumnegative = 0;

		for (int j = 0; j < i; j++)
		{
			perturb[j] = perturb[j] * 50 * pow(10, 1 - precision) + vect2[j];
			sum1 = sum1 + perturb[j];
			floating p = translate(perturb[j], true);

			// calculating round to nearest intermediate sum, sum2 is used to calculate round to nearest bound
			intermediate_sum1 = addition(intermediate_sum1, p, true);
			sum2 = sum2 + abs(intermediate_sum1.mantissa * pow(10, intermediate_sum1.exp));

			// calculating round to 0 intermediate sum, sumpositive and sumnegative are used to calculate round to 0 bound
			intermediate_sum2 = addition(intermediate_sum2, p, false);
			if (intermediate_sum2.mantissa >= 0)
			{
				sumpositive = sumpositive + intermediate_sum2.mantissa * pow(10, intermediate_sum2.exp);
			}
			else
			{
				sumnegative = sumnegative - intermediate_sum2.mantissa * pow(10, intermediate_sum2.exp);
			}
		}

		// calculating round to nearest bound as bound1, round to 0 bound as bound2;
		double bound1 = sum2 * pow(10, 1 - precision) / 2;
		double bound2 = pow(10, 1 - precision);
		if (sumpositive > sumnegative)
		{
			bound2 = bound2 * sumpositive;
		}
		else
		{
			bound2 = bound2 * sumnegative;
		}

		// calculating running error for round to nearest as error1 and round to 0 as error2
		double error1 = abs(sum1 - intermediate_sum1.mantissa * pow(10, intermediate_sum1.exp));
		double error2 = abs(sum1 - intermediate_sum2.mantissa * pow(10, intermediate_sum2.exp));

		// write running error ratio of round to nearest to "accumulation4_3.txt" for vect2
		file << error1 / bound1 << "\n";

		// write running error ratio of round to 0 to "accumulation4_4.txt" for vect2
		file1 << error2 / bound2 << "\n";

		count = count + 1;
	}

	file.close();
	file1.close();

	return 0;
}
