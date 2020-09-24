// This program is used to perform LU factorization of an input matrix
// algorithms for no pivoting, partial pivoting and complete pivoting are included

#include <iostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>
#include <chrono>
#include <time.h>
#include <sstream>

using namespace std;

int n;
double threshold = 0.00001;
// the following vector is used to save the permutations in complete pivoting for vector x
vector<int> position = {};

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

// generate matrix A, get_A(true) if randomly generated, get_A(false) if input from user's data
vector<vector<double>> get_A(bool m)
{
	vector<vector<double>> A = {};
	if (m == true)
	{
		ofstream file;
		file.open("A.csv", ios::out);

		// randomly generate A
		for (int i = 0; i < n; i++)
		{
			A.push_back(create(-10, 10, n));
		}

		// run this code if you want A to be diagonal-dominant
		// make the diagonal element in each row equal to 100 * the element with the max magnitude of that row
		/*for (int i = 0; i < n; i++)
		{
			double max = abs(A[i][0]);
			int track = 0;
			for (int j = 1; j < n; j++)
			{
				if (abs(A[i][j]) > max)
				{
					max = abs(A[i][j]);
					track = j;
				}
			}
			A[i][i] = 100 * A[i][track];
		}

		// run the following code to permutate a diagonally-dominant matrix
		// we simply interchange row 0 and row n-1 of the matrix
		vector<double> temp = {};
		temp = A[0];
		A[0] = A[n - 1];
		A[n - 1] = temp;*/

		// output A to a file
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				file << A[i][j] << ", ";
			}
			file << "\n";
		}
		file.close();
	}
	else
	{
		ifstream file1("A.txt");
		// read A from the file
		for (int i = 0; i < n; i++)
		{
			vector<double> vec = {};
			for (int j = 0; j < n; j++)
			{
				string temp1;
				getline(file1, temp1, ',');				
				vec.push_back(stod(temp1));
			}
			A.push_back(vec);
		}
		file1.close();
	}
	return A;
}

// generate vector b, get_b(true) if randomly generated, get_b(false) if input from user's data
vector<double> get_b(bool m)
{
	vector<double> b = {};
	if (m == true)
	{
		ofstream file;
		file.open("b.csv", ios::out);
		b = create(-10, 10, n);

		for (int i = 0; i < n; i++)
		{
			file << b[i] << "\n";
		}
		file.close();
	}
	else
	{
		ifstream file("b.csv");

		for (int i = 0; i < n; i++)
		{
			double temp = 0;
			file >> temp;
			b.push_back(temp);
		}
		file.close();
	}
	return b;
}

// LU factorization without pivoting
int no_pivoting(vector<vector<double>> A, vector<double> b)
{
	ofstream file1;
	ofstream file2;
	ofstream file3;
	file1.open("U.csv", ios::out);
	file2.open("L.csv", ios::out);
	file3.open("factorized_b.csv", ios::out);

	// variable i runs through the columns to calculate l 
	// variable j runs through the rows to adjust the new values
	// variable k runs through the elements of row j to adjust the new values
	for (int i = 0; i < n - 1; i++)
	{
		if (abs(A[i][i]) < threshold)
		{
			return 1;
		}

		for (int j = i + 1; j < n; j++)
		{
			// calculate lambda and save it in the respective column i
			A[j][i] = A[j][i] / A[i][i];
			// make changes in A and b
			for (int k = i + 1; k < n; k++)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
			b[j] = b[j] - A[j][i] * b[i];
		}
	}

	// output the result of U in file1, L in file2 and b in file3
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j)
			{
				file1 << 0 << ", ";
				file2 << A[i][j] << ", ";
			}
			else if (i < j)
			{
				file1 << A[i][j] << ", ";
				file2 << 0 << ", ";
			}
			else
			{
				file1 << A[i][j] << ", ";
				file2 << 1 << ", ";
			}
		}
		file1 << "\n";
		file2 << "\n";
		file3 << b[i] << "\n";
	}
	file1.close();
	file2.close();
	file3.close();
	return 0;
}

// LU factorization with partial pivoting
int partial_pivoting(vector<vector<double>> A, vector<double> b)
{
	ofstream file1;
	ofstream file2;
	ofstream file3;
	ofstream file4;
	file1.open("U.csv", ios::out);
	file2.open("L.csv", ios::out);
	file3.open("permutated_A.csv", ios::out);
	file4.open("factorized_b.csv", ios::out);
	
	// B will be the value of permutated A along with the pivoting process, B = PA
	vector<vector<double>> B = A;

	// variable i runs through the columns to calculate l 
	// variable j runs through the rows to adjust the new values
	// variable k runs through the elements of row j to adjust the new values
	for (int i = 0; i < n - 1; i++)
	{
		// perform partial pivoting
		double max = abs(A[i][i]);
		int track = i;
		// check the i column for the element with maximum magnitude
		for (int j = i + 1; j < n; j++)
		{
			if (abs(A[j][i]) > max)
			{
				max = abs(A[j][i]);
				track = j;
			}
		}
		
		// interchange rows if needed, interchanging rows also causes interchange of rows in b
		if (track != i)
		{
			vector<double> temp = A[i];
			A[i] = A[track];
			A[track] = temp;
			
			temp = B[i];
			B[i] = B[track];
			B[track] = temp;

			double temp1 = b[i];
			b[i] = b[track];
			b[track] = temp1;
		}

		if (abs(A[i][i]) < threshold)
		{
			return 1;
		}

		for (int j = i + 1; j < n; j++)
		{
			// calculate lambda and save it in the respective column i
			A[j][i] = A[j][i] / A[i][i];

			// make change in A and b
			for (int k = i + 1; k < n; k++)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
			b[j] = b[j] - A[j][i] * b[i];
		}
	}

	// output the result of U in file1, L in file2, permutated A (as B) in file3 and b in file4
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j)
			{
				file1 << 0 << ", ";
				file2 << A[i][j] << ", ";
			}
			else if (i < j)
			{
				file1 << A[i][j] << ", ";
				file2 << 0 << ", ";
			}
			else
			{
				file1 << A[i][j] << ", ";
				file2 << 1 << ", ";
			}
			file3 << B[i][j] << ", ";
		}
		file1 << "\n";
		file2 << "\n";
		file3 << "\n";
		file4 << b[i] << "\n";
	}
	file1.close();
	file2.close();
	file3.close();
	file4.close();
	return 0;
}

// LU factorization with complete pivoting
int complete_pivoting(vector<vector<double>> A, vector<double> b)
{
	ofstream file1;
	ofstream file2;
	ofstream file3;
	ofstream file4;
	file1.open("U.csv", ios::out);
	file2.open("L.csv", ios::out);
	file3.open("permutated_A.csv", ios::out);
	file4.open("factorized_b.csv", ios::out);

	// B will be the value of permutated A along with the pivoting process, B = PA
	vector<vector<double>> B = A;
	
	// variable i runs through the columns to calculate l 
	// variable j runs through the rows to adjust the new values
	// variable k runs through the elements of row j to adjust the new values
	for (int i = 0; i < n - 1; i++)
	{
		// perform complete pivoting
		double max = abs(A[i][i]);
		int row = i;
		int column = i;
		// check the active matrix for the element with maximum magnitude
		for (int j = i; j < n; j++)
		{
			for (int k = i; k < n; k++)
			{
				if (abs(A[j][k]) > max)
				{
					max = abs(A[j][k]);
					row = j;
					column = k;
				}
			}
		}

		// interchange rows and columns if needed
		// interchanging rows causes interchange of rows in b, interchanging column causes interchange of rows in x
		if (row != i)
		{
			vector<double> temp = A[i];
			A[i] = A[row];
			A[row] = temp;

			temp = B[i];
			B[i] = B[row];
			B[row] = temp;

			double temp1 = b[i];
			b[i] = b[row];
			b[row] = temp1;
		}

		if (column != i)
		{

			for (int j = 0; j < n; j++)
			{
				double temp = A[j][i];
				A[j][i] = A[j][column];
				A[j][column] = temp;

				temp = B[j][i];
				B[j][i] = B[j][column];
				B[j][column] = temp;
			}
			// this vector reflects the position of x
			double temp2 = 0;
			temp2 = position[i];
			position[i] = position[column];
			position[column] = temp2;
		}

		if (abs(A[i][i]) < threshold)
		{
			return 1;
		}

		for (int j = i + 1; j < n; j++)
		{
			// calculate lambda and save it in the respective column i
			A[j][i] = A[j][i] / A[i][i];
			// make change in A and b
			for (int k = i + 1; k < n; k++)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
			b[j] = b[j] - A[j][i] * b[i];
		}
	}

	// 	// output the result of U in file1, L in file2, permutated A (as B) in file3 and b in file4
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j)
			{
				file1 << 0 << ", ";
				file2 << A[i][j] << ", ";
			}
			else if (i < j)
			{
				file1 << A[i][j] << ", ";
				file2 << 0 << ", ";
			}
			else
			{
				file1 << A[i][j] << ", ";
				file2 << 1 << ", ";
			}
			file3 << B[i][j] << ", ";
		}
		file1 << "\n";
		file2 << "\n";
		file3 << "\n";
		file4 << b[i] << "\n";
	}
	file1.close();
	file2.close();
	file3.close();
	file4.close();
	return 0;
}

// this function solves for x in Ax = b
void solve(int m)
{
	// read the lower triangular L from LU factorization to solve for x
	ifstream file1("U.csv");
	ifstream file2("factorized_b.csv");
	ofstream file3;
	file3.open("computed_x.csv", ios::out);
	vector<vector<double>> U = {};
	vector<double> b = {};
	for (int i = 0; i < n; i++)
	{
		vector<double> vec = {};
		for (int j = 0; j < n; j++)
		{
			string temp1;
			getline(file1, temp1, ',');
			vec.push_back(stod(temp1));
		}
		U.push_back(vec);
		double temp2 = 0;
		file2 >> temp2;
		b.push_back(temp2);
	}
	file1.close();
	file2.close();

	// perform solution for x in the same vector b
	b[n - 1] = b[n - 1] / U[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		for (int j = i + 1; j < n; j++)
		{
			b[i] = b[i] - U[i][j] * b[j];
		}
		b[i] = b[i] / U[i][i];
	}

	// for complete pivoting, since the position of x might have been changed
	// we need to permutate them back to original position for comparison

	vector<double> temp = b;
	for (int i = 0; i < n; i++)
	{
		b[position[i]] = temp[i];
	}

	// output the computed x
	for (int i = 0; i < n; i++)
	{
		file3 << b[i] << "\n";
	}
	file3.close();
}

int main()
{
	// change n as desire
	n = 20;

	// generate matrix A, get_A(true) if randomly generated, get_A(false) if input from user's data
	// generate vector b, get_b(true) if randomly generated, get_b(false) if input from user's data
	vector<vector<double>> A = get_A(true);
	vector<double> b = get_b(true);

	// the following "for" loop is only useful for complete pivoting
	for (int i = 0; i < n; i++)
	{
		position.push_back(i);
	}

	// perform LU factorization
	int result = complete_pivoting(A, b);

	if (result == 1)
	{
		cout << "Cannot perform factorization due to singularity or values being too small";
	}

	// solve for x
	solve(0);
}
