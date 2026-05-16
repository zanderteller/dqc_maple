#include "matrix.h"
#include "array3d.h"
#include <complex.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cctype>
#include "distributed.h"
#include "net.h"

char *after(char *haystack, char const *needle)
{
	char*ret = strstr(haystack, needle);
	if (ret) {
		ret += strlen(needle);
	}
	return ret;
}

Matrix<double> *read_real_matrix(char *cur, char **out)
{
	int rows = strtol(cur, &cur, 10);
	cur++;
	int cols = strtol(cur, &cur, 10);
	auto ret = new Matrix<double>(cols, rows);
	for (int y = 0; y < rows; y++)
	{
		for (int x = 0; x < cols; x++)
		{
			while (*cur != '-' && *cur != '.' && !isdigit(*cur))
			{
				cur++;
			}
			ret->set(x, y, strtod(cur, &cur));
		}
	}
	if (out)
	{
		*out = cur;
	}
	return ret;
}

int main(int argc, char **argv)
{
	int num_frames = 2;
	if(argc < 2) {
		return 1;
	}
	FILE *f = fopen(argv[1], "rb");
	
	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	fseek(f, 0, SEEK_SET);
	
	char *text = new char[fsize+1];
	if (fsize != fread(text, 1, fsize, f)) {
		return 1;
	}
	fclose(f);
	
	bool split(false);
	for (int i = 2; i < argc; i++)
	{
		if (!strcmp("--split", argv[i])) {
			split = true;
		} else {
			AddSecondary(argv[i], 1.0f);
		}
	}
	
	text[fsize] = 0;
	char *expHText = after(text, "expH := Matrix(");
	
	char *cur;
	int rows = strtol(expHText, &cur, 10);
	cur++;
	int cols = strtol(cur, &cur, 10);
	Matrix<complex<double>> expH(cols, rows);
	for (int y = 0; y < rows; y++)
	{
		for (int x = 0; x < cols; x++)
		{
			cur = after(cur, "complex<double>(");
			double r,i;
			r = strtod(cur, &cur);
			while (*cur != '-' && *cur != '.' && !isdigit(*cur))
			{
				cur++;
			}
			i = strtod(cur, &cur);
			expH.set(x, y, complex<double>(r, i));
		}
	}
	char *xopsText = after(text, "xops := Array(");
	int first = strtol(xopsText, &cur, 10);
	cur += 2;
	int size = 1 + strtol(cur, &cur, 10) - first;
	vector<Matrix<double> > xops;
	for (int i = 0; i < size; i++)
	{
		
		cur = after(cur, "Matrix(");
		auto el = read_real_matrix(cur, &cur);
		xops.push_back(*el);
		delete el;
	}
	
	char *SimText = after(text, "Sim := Matrix(");
	auto Sim = read_real_matrix(SimText, 0);
	
	char *SimdagText = after(text, "Simdag := Matrix(");
	auto Simdag = read_real_matrix(SimdagText, 0);
	
	char *Traj1Text = after(text, "Traj1 := Array(");
	first = strtol(Traj1Text, &cur, 10);
	cur += 2;
	rows = 1 + strtol(cur, &cur, 10) - first;
	cur++;
	first = strtol(cur, &cur, 10);
	cur += 2;
	cols = 1 + strtol(cur, &cur, 10) - first;
	cur++;
	first = strtol(cur, &cur, 10);
	cur += 2;
	int zsize = 1 + strtol(cur, &cur, 10) - first;
	Array3D Traj1(cols, rows, zsize + num_frames);
	for (int y = 0; y < rows; y++)
	{
		for (int x = 0; x < cols; x++)
		{
			for (int z = 0; z < zsize; z++)
			{
				while (*cur != '-' && *cur != '.' && !isdigit(*cur))
				{
					cur++;
				}
				Traj1.set(x, y, z, strtod(cur, &cur));
			}
		}
	}
	if (split) {
		DistributedAddTrajectories(0, num_frames/2, Traj1, expH, xops, *Simdag, expH.rows(), xops.size(), 0.09, 0.6, 0.2);
		DistributedAddTrajectories(num_frames/2, num_frames, Traj1, expH, xops, *Simdag, expH.rows(), xops.size(), 0.09, 0.6, 0.2);
	} else {
		DistributedAddTrajectories(0, num_frames, Traj1, expH, xops, *Simdag, expH.rows(), xops.size(), 0.09, 0.6, 0.2);
	}
	//Array3D output(RunIter(0, Traj1, expH, xops, *Simdag, expH.rows(), xops.size(), 0.09, 0.6, 0.2));
	for (int z = 0; z < Traj1.depth(); z++)
	{
		for (int y = 0; y < Traj1.rows(); y++)
		{
			for (int x = 0; x < Traj1.cols(); x++)
			{
				printf("%f\n", Traj1.get(x,y,z));
			}
		}
	}
	
	delete[] text;
	return 0;
}