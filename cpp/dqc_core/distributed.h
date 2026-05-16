#ifndef DISTRIBUTED_H_
#define DISTRIBUTED_H_

#include <complex>
#include <tuple>
#include "matrix.h"
#include "array3d.h"
#include "operators.h"

using std::complex;
using std::tuple;

//Commands that can be sent to a secondary
enum DistCommand {
	CMD_MATRIXMULT,
	CMD_ADDTRAJECTORIES,
	CMD_BUILDFRAMESAUTO,
	CMD_AGGREGATEPOTENTIALCONTRIBUTIONS,
	CMD_PING
};

enum ValType {
	TYPE_DOUBLE,
	TYPE_COMPLEX
};

//writes a 32-bit integer to a buffer in little-endian format
inline char *write32(char *buf, int num)
{
	*(buf++) = num;
	*(buf++) = num >> 8;
	*(buf++) = num >> 16;
	*(buf++) = num >> 24;
	return buf;
}

//writes a double to a buffer in the host endian format
//NOTE: depending on the host's endianness is technically bad, but big-endian systems are very uncommon now
inline char *writeDouble(char *buf, double val)
{
	memcpy(buf, &val, sizeof(val));
	return buf + sizeof(val);
}

// writes a boolean to a buffer
inline char* writeBoolean(char* buf, bool val)
{
	*(buf++) = val;
	return buf;
}

// identical in function to the normal AddTrajectories, but the results are distributed across all the current secondaries
// note that all nodes need to process nHambasis rows, so you'll only get a good speedup from using secondaries when the
// hambasis is much smaller than the whole result set
void DistributedAddTrajectories(int startFrame, int endFrame, Array3D& tFinal, Matrix<complex<double>>& expH, vector<Matrix<double> >& xops, Matrix<double>& simDag, int nHambasis, int nPotbasis, double sigma, double step, double mass);

// identical in function to BuildFramesAuto, but the calculations are distributed across all current secondary/remote
// machines (and the local machine)
// note that this setup does *not* suffer from the problem described above for DistributedAddTrajectories with nHambasis
void DistributedBuildFramesAuto(Array3D& all_frames, Matrix<double>& basis_rows, Matrix<complex<double>>& exph, vector<Matrix<double> >& xops,
								Matrix<double>& simdag, double sigma, double stopping_threshold, bool rescale);

// identical in function to AggregatePotentialContributions, but the calculations are distributed across all current secondary/remote
// machines (and the local machine)
tuple<Matrix<double>, Matrix<double>> DistributedAggregatePotentialContributions(Matrix<double>& mat, Matrix<double>& mat_ham_sub, double sigma);

#endif //DISTRIBUTED_H_
