#ifndef TRAJECTORIES_H_
#define TRAJECTORIES_H_
#include <complex>

using std::vector;
using std::complex;

void AddTrajectories(int startFrame, int endFrame, Array3D &tFinal, Matrix<complex<double>> &expH, vector<Matrix<double> > &xops, Matrix<double> &simDag,
					int nHambasis, int nPotbasis, double sigma, double step, double mass);

void BuildFramesAuto(Array3D &all_frames, Matrix<double> &basis_rows, Matrix<complex<double>> &exph, vector<Matrix<double> > &xops, Matrix<double> &simdag,
					double sigma, double stopping_threshold, bool rescale);

#endif //TRAJECTORIES_H_
