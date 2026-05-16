#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <cstdio>
#include <ctime>
#include <thread>
#include <functional>
#include <tuple>
#include "net.h"
#include "matrix.h"
#include <complex>
#include "array3d.h"
#include "operators.h"
#include "trajectories.h"
#include "distributed.h"

using std::tuple;
using std::make_tuple;
using std::get;

void receiver_thread(FILE *log, CommScope &comms, vector<char *> &results, vector<int> &work_units, size_t result_size, bool size_is_per_unit)
{
	time_t start = time(NULL);
	if (log) { fprintf(log, "Receiver thread: started at %s", asctime(localtime(&start))); }
	comms.ReceiveResults(results, work_units, result_size, size_is_per_unit);
	time_t end = time(NULL);
	if (log) { fprintf(log, "Receiver thread: all remote results received after %d seconds -- done at %s", int(end - start), asctime(localtime(&end))); }
}

// call AddTrajectories on multiple machines (via the 'secondary' executable) and aggregate the results
void DistributedAddTrajectories(int startFrame, int endFrame, Array3D &tFinal, Matrix<complex<double>> &expH, vector<Matrix<double> > &xops, Matrix<double> &simDag, int nHambasis, int nPotbasis, double sigma, double step, double mass)
{
	if (tFinal.rows() == nHambasis)
	{
		//work can't be distributed, just run the operation locally
		AddTrajectories(startFrame, endFrame, tFinal, expH, xops, simDag, nHambasis, nPotbasis, sigma, step, mass);
		return;
	}
	vector<int> work_units;
	DistributeUnits(tFinal.rows() - nHambasis, work_units);
	char base_buffer[1 + 11 * 4 + 3 * sizeof(double)];
	base_buffer[0] = CMD_ADDTRAJECTORIES;
	char * cur = base_buffer + 1;
	cur = write32(cur, startFrame ? 1 : 0);
	cur = write32(cur, endFrame - startFrame + (startFrame ? 1 : 0));
	cur = write32(cur, nHambasis);
	cur = write32(cur, nPotbasis);
	cur = write32(cur, expH.cols());
	cur = write32(cur, expH.rows());
	cur = write32(cur, simDag.cols());
	cur = write32(cur, simDag.rows());
	cur = write32(cur, xops.size());
	cur = write32(cur, tFinal.cols());
	cur = writeDouble(cur, sigma);
	cur = writeDouble(cur, step);
	char *end_base = writeDouble(cur, mass);
	int cur_row = nHambasis + work_units[0];
	vector <char *> results;
	FILE *f = NULL; // fopen("dist_addtrajectories_log.txt", "w");
	time_t start = time(NULL);
	if (f) { fprintf(f, "Started at %s", asctime(localtime(&start))); }
	
	{
		//Handles acquiring the comms mutex and setting up keepalive thread
		CommScope comms;
		for (int i = 1; i < work_units.size(); i++)
		{
			cur = write32(end_base, work_units[i]);
			comms.SendData(i-1, base_buffer, cur - base_buffer);
			if (startFrame) {
				Array3D hambasisFrame0(tFinal.view(0, tFinal.cols(), 0, nHambasis, 0, 1));
				comms.SendData(i-1,  (char *)hambasisFrame0.raw(), hambasisFrame0.cols() * hambasisFrame0.rows() * sizeof(double));
				Array3D workFrame0(tFinal.view(0, tFinal.cols(), cur_row, work_units[i], 0, 1));
				comms.SendData(i-1, (char *)workFrame0.raw(), workFrame0.cols() * workFrame0.rows() * sizeof(double));
			}
			
			Array3D hambasis(tFinal.view(0, tFinal.cols(), 0, nHambasis, startFrame, 1));
			comms.SendData(i-1, (char *)hambasis.raw(), hambasis.cols() * hambasis.rows() * sizeof(double));
			Array3D work(tFinal.view(0, tFinal.cols(), cur_row, work_units[i], startFrame, 1));
			comms.SendData(i-1, (char *)work.raw(), work.cols() * work.rows() * sizeof(double));
			comms.SendData(i-1, (char *)expH.raw(), expH.cols() * expH.rows() * sizeof(complex<double>));
			comms.SendData(i-1, (char *)simDag.raw(), simDag.cols() * simDag.rows() * sizeof(double));
			for (int j = 0; j < xops.size(); j++)
			{
				char buffer[8];
				cur = write32(buffer, xops[j].cols());
				cur = write32(cur, xops[j].rows());
				comms.SendData(i-1, buffer, sizeof(buffer));
				comms.SendData(i-1, (char *)xops[j].raw(), xops[j].rows() * xops[j].cols() * sizeof(double));
			}
			results.push_back(new char [tFinal.cols() * work_units[i] * (endFrame - startFrame) * sizeof(double)]);
			cur_row += work_units[i];
		}
		time_t sent_time = time(NULL);
		if (f) { fprintf(f, "Finished sending data to remotes after %d seconds -- done at %s", int(sent_time - start), asctime(localtime(&sent_time))); }
		Array3D local(tFinal.view(0, tFinal.cols(), 0, nHambasis + work_units[0], 0, tFinal.depth()));
		bool size_is_per_unit = true; // total result size is proportional to number of work units
		std::thread receiver(receiver_thread, f, std::ref(comms), std::ref(results), std::ref(work_units), tFinal.cols() * (endFrame - startFrame) * sizeof(double), size_is_per_unit);
		if (f) { fprintf(f, "Calling AddTrajectories to calculate %d frames with data size %d cols, %d rows\n", endFrame - startFrame, local.cols(), local.rows()); }
		AddTrajectories(startFrame, endFrame, local, expH, xops, simDag, nHambasis, nPotbasis, sigma, step, mass);
		time_t local_time = time(NULL);
		if (f) { fprintf(f, "Finished local part of calculation in %d seconds -- done at %s", int(local_time - start), asctime(localtime(&local_time))); }
		receiver.join();
	}

	cur_row = nHambasis + work_units[0];
	for (int i = 0; i < results.size(); i++)
	{
		Array3D result(tFinal.cols(), work_units[i+1], endFrame - startFrame, (double *)results[i]);
		for (int z = 0; z < result.depth(); z++)
		{
			for (int y = 0; y < result.rows(); y++)
			{
				for (int x = 0; x < result.cols(); x++)
				{
					tFinal.set(x, y+cur_row, startFrame+z+1, result.get(x, y, z));
				}
			}
		}
		cur_row += result.rows();
		delete[] results[i];
	}
	time_t end = time(NULL);
	if (f)
	{
		fprintf(f, "Finished combining results after %d seconds -- done at %s", int(end - start), asctime(localtime(&end)));
		fclose(f);
	}
} // end function DistributedAddTrajectories

// call BuildFramesAuto on multiple machines (via the 'secondary' executable) and aggregate the results
// note that this signature is identical to the signature of BuildFramesAuto
void DistributedBuildFramesAuto(Array3D& all_frames, Matrix<double>& basis_rows, Matrix<complex<double>>& exph, vector<Matrix<double> >& xops,
								Matrix<double>& simdag, double sigma, double stopping_threshold, bool rescale)
{
	FILE* f = NULL; // fopen("dist_buildframesauto_log.txt", "w");
	time_t start_time = time(NULL);
	if (f) { fprintf(f, "Started at %s", asctime(localtime(&start_time))); }

	int num_cols = all_frames.cols();

	// note: this includes the initial frame (passed in) and all subsequent frames (which we build here)
	int num_frames = all_frames.depth();

	int num_basis_rows = basis_rows.rows();

	assert(num_cols == basis_rows.cols());
	assert(num_cols == xops.size());

	// divide up the rows between the machines
	vector<int> work_units;
	DistributeUnits(all_frames.rows(), work_units);

	if (work_units.size() == 1)
	{
		// no remote machines (or very small job) -- just do the work locally

		if (f) { fprintf(f, "No remote machines or very small job -- doing all work locally\n"); };

		BuildFramesAuto(all_frames, basis_rows, exph, xops, simdag, sigma, stopping_threshold, rescale);

		if (f)
		{
			time_t end_time = time(NULL);
			fprintf(f, "Calculated local results in %d seconds -- done at %s", int(end_time - start_time), asctime(localtime(&end_time)));
			fclose(f);
		}
		return;
	} // end if only running locally

	// send scalar data -- need room for the command, 2 doubles, 8 ints, and 1 boolean
	// note: the last parameter (work units -- i.e., rows) is added separately for each remote in the for loop below
	char base_buffer[1 + 2 * sizeof(double) + 8 * 4 + 1];
	base_buffer[0] = CMD_BUILDFRAMESAUTO; // only need 1 byte for this
	char* cur = base_buffer + 1;
	// add double values
	cur = writeDouble(cur, sigma);
	cur = writeDouble(cur, stopping_threshold);
	// add int values
	cur = write32(cur, num_frames);
	cur = write32(cur, num_cols);
	cur = write32(cur, num_basis_rows);
	cur = write32(cur, exph.rows());
	cur = write32(cur, exph.cols());
	cur = write32(cur, simdag.rows());
	cur = write32(cur, simdag.cols());
	char* end_base = writeBoolean(cur, rescale);

	int cur_row = work_units[0]; // work_units[0] will be done by this machine, so start after that
	vector <char*> results;

	// 2FIX: is this code block here purely so that the CommScope destructor gets called at the end of it?
	// (seems like there should be a cleaner way to have that organized -- maybe just move the code out of the destructor)
	{
		CommScope comms; // handles acquiring the comms mutex and setting up keepalive thread

		// index 0 in work_units is the local machine -- start at 1 for the remotes
		for (int i = 1; i < work_units.size(); i++)
		{
			// for each remote machine...
			// note: i is index of remote machine in 'work_units', and i - 1 is index of remote machine in 'remotes'

			// add number or work units (number of rows) and send all scalar parameters
			cur = write32(end_base, work_units[i]);
			comms.SendData(i - 1, base_buffer, cur - base_buffer);

			// send basis rows
			comms.SendData(i - 1, (char*)basis_rows.raw(), basis_rows.rows() * basis_rows.cols() * sizeof(double));

			// send rows to evolve
			Matrix<double> work(all_frames.view2D(0, num_cols, cur_row, work_units[i], 0));
			comms.SendData(i - 1, (char*)work.raw(), work.rows() * work.cols() * sizeof(double));
			cur_row += work_units[i]; // increment for the next iteration of the loop

			// send operators
			comms.SendData(i - 1, (char*)exph.raw(), exph.rows() * exph.cols() * sizeof(complex<double>));
			comms.SendData(i - 1, (char*)simdag.raw(), simdag.rows() * simdag.cols() * sizeof(double));

			for (int j = 0; j < xops.size(); j++)
			{
				char buffer[8]; // room for 2 ints
				cur = write32(buffer, xops[j].rows());
				cur = write32(cur, xops[j].cols());
				comms.SendData(i - 1, buffer, sizeof(buffer));
				comms.SendData(i - 1, (char*)xops[j].raw(), xops[j].rows() * xops[j].cols() * sizeof(double));
			}

			// allocate memory for what we're expecting to receive back from the remote
			// note: remotes do *not* send back the initial frame we sent to them
			results.push_back(new char[work_units[i] * num_cols * (num_frames - 1) * sizeof(double)]);
		} // end for each remote

		time_t sent_time = time(NULL);
		if (f) { fprintf(f, "Finished sending data to remotes after %d seconds -- done at %s", int(sent_time - start_time), asctime(localtime(&sent_time))); }

		// start receiver thread
		bool size_is_per_unit = true; // total result size is proportional to number of work units
		std::thread receiver(receiver_thread, f, std::ref(comms), std::ref(results), std::ref(work_units), num_cols * (num_frames - 1) * sizeof(double), size_is_per_unit);

		// set up and run the work assigned to this machine
		// note: num_frames currently includes the initial frame
		Array3D local(all_frames.view(0, num_cols, 0, work_units[0], 0, num_frames));
		// 2FIX: ALSO FIX HERE TO MATCH NUM FRAMES (MEANS NUM NEW FRAMES BUILT OR INCLUDES INITIAL FRAME?)
		if (f) { fprintf(f, "Calling BuildFramesAuto to calculate %d new frames with data size %d rows, %d cols\n",
						num_frames - 1, local.rows(), local.cols()); }
		BuildFramesAuto(local, basis_rows, exph, xops, simdag, sigma, stopping_threshold, rescale);
		time_t local_time = time(NULL);
		if (f) { fprintf(f, "Finished local part of calculations in %d seconds -- done at %s", int(local_time - start_time), asctime(localtime(&local_time))); }

		// wait for remotes to finish their work
		receiver.join();

	} // end code block indented for mysterious reasons (for scoping purposes?  and/or having to do with threading?)

	// put results together into final array
	cur_row = work_units[0]; // local work was done directly in a view of all_frames, so local results are already in place
	for (int i = 0; i < results.size(); i++)
	{
		// note: initial frame is not returned by the remote
		Array3D result(num_cols, work_units[i + 1], num_frames - 1, (double*)results[i]);
		for (int z = 0; z < result.depth(); z++)
		{
			for (int y = 0; y < result.rows(); y++)
			{
				for (int x = 0; x < result.cols(); x++)
				{
					all_frames.set(x, cur_row + y, z + 1, result.get(x, y, z));
				}
			}
		}
		cur_row += result.rows();
		delete[] results[i];
	}

	time_t end_time = time(NULL);
	if (f)
	{
		fprintf(f, "Finished combining results after %d seconds -- done at %s", int(end_time - start_time), asctime(localtime(&end_time)));
		fprintf(f, "Done\n");
		fclose(f);
	}
} // end function DistributedBuildFramesAuto

// call AggregatePotentialConbtributions on multiple machines (via the 'secondary' executable) and aggregate the results
// note that this signature is identical to the signature of AggregatePotentialContributions
tuple<Matrix<double>, Matrix<double>> DistributedAggregatePotentialContributions(Matrix<double>& mat, Matrix<double>& mat_ham_sub, double sigma)
{
	FILE* f = NULL; // fopen("dist_aggregatepotentialcontributions_log.txt", "w");
	time_t start_time = time(NULL);
	if (f) { fprintf(f, "Started at %s", asctime(localtime(&start_time))); }

	int num_basis_rows = mat_ham_sub.rows();
	int num_cols = mat.cols();
	assert(num_cols == mat_ham_sub.cols()); // mat and mat_ham_sub must have the same number of columns

	// divide up the rows between the machines
	vector<int> work_units;
	DistributeUnits(mat.rows(), work_units);

	if (work_units.size() == 1)
	{
		// no remote machines (or very small job) -- just do the work locally

		if (f) { fprintf(f, "No remote machines or very small job -- doing all work locally\n"); };

		tuple<Matrix<double>, Matrix<double>> results = AggregatePotentialContributions(mat, mat_ham_sub, sigma);

		if (f)
		{
			time_t end_time = time(NULL);
			fprintf(f, "Calculated local results in %d seconds -- done at %s", int(end_time - start_time), asctime(localtime(&end_time)));
			fclose(f);
		}

		return results;
	} // end if only running locally

	// send scalar data -- need room for the command, 1 double (sigma), and 2 ints (num_basis_rows and num_cols)
	// note: the last parameter (work units -- i.e., rows) is added separately for each remote in the for loop below
	char base_buffer[1 + sizeof(double) + 2 * 4];
	base_buffer[0] = CMD_AGGREGATEPOTENTIALCONTRIBUTIONS; // only need 1 byte for this
	char* cur = base_buffer + 1;
	// add double values (just sigma)
	cur = writeDouble(cur, sigma);
	// add int values (num_basis_rows and num_cols)
	cur = write32(cur, num_basis_rows);
	char* end_base = write32(cur, num_cols);

	int cur_row = work_units[0]; // work_units[0] will be done by this machine, so start after that
	vector <char*> results;

	// returned results from each remote will be 2 double matrices, each of size <num_basis_rows x num_basis_rows>
	size_t result_size = 2 * num_basis_rows * num_basis_rows * sizeof(double);

	Matrix<double> pot_mat(num_basis_rows, num_basis_rows);
	Matrix<double> psi_mat(num_basis_rows, num_basis_rows);

	// 2FIX: is this code block here purely so that the CommScope destructor gets called at the end of it?
	// (seems like there should be a cleaner way to have that organized -- maybe just move the code out of the destructor)
	{
		CommScope comms; // handles acquiring the comms mutex and setting up keepalive thread

		// index 0 in work_units is the local machine -- start at 1 for the remotes
		for (int i = 1; i < work_units.size(); i++)
		{
			// for each remote machine...
			// note: i is index of remote machine in 'work_units', and i - 1 is index of remote machine in 'remotes'

			// add work units for this remote and send scalar parameters
			cur = write32(end_base, work_units[i]);
			comms.SendData(i - 1, base_buffer, cur - base_buffer);

			// send basis rows
			comms.SendData(i - 1, (char*)mat_ham_sub.raw(), mat_ham_sub.rows() * mat_ham_sub.cols() * sizeof(double));

			// send rows to process
			Matrix<double> work(mat.view(0, num_cols, cur_row, work_units[i]));
			comms.SendData(i - 1, (char*)work.raw(), work.rows() * work.cols() * sizeof(double));
			cur_row += work_units[i]; // increment for the next iteration of the loop

			// allocate memory for what we're expecting to receive back from the remote
			results.push_back(new char[result_size]);
		} // end for each remote

		time_t sent_time = time(NULL);
		if (f) { fprintf(f, "Finished sending data to remotes after %d seconds -- done at %s", int(sent_time - start_time), asctime(localtime(&sent_time))); }

		// start receiver thread
		bool size_is_per_unit = false; // total result size is not proportional to number of work units
		std::thread receiver(receiver_thread, f, std::ref(comms), std::ref(results), std::ref(work_units), result_size, size_is_per_unit);

		// set up and run the work assigned to this machine
		Matrix<double> mat_local(mat.view(0, num_cols, 0, work_units[0]));
		if (f) { fprintf(f, "Calling AggregatePotentialContributions to calculate contributions to potential for %d rows (with %d cols)\n",
						mat_local.rows(), mat_local.cols()); }
		tuple<Matrix<double>, Matrix<double>> local_results = AggregatePotentialContributions(mat_local, mat_ham_sub, sigma);
		Matrix<double> pot_mat_local(get<0>(local_results));
		Matrix<double> psi_mat_local(get<1>(local_results));
		memcpy(pot_mat.raw(), pot_mat_local.raw(), num_basis_rows * num_basis_rows * sizeof(double));
		memcpy(psi_mat.raw(), psi_mat_local.raw(), num_basis_rows * num_basis_rows * sizeof(double));

		time_t local_time = time(NULL);
		if (f) { fprintf(f, "Finished local part of calculations in %d seconds -- done at %s", int(local_time - start_time), asctime(localtime(&local_time))); }

		// wait for remotes to finish their work
		receiver.join();
	} // end code block indented for mysterious reasons (for scoping purposes?  and/or having to do with threading?)

	// aggregate all results into pot_mat and psi_mat
	for (int i = 0; i < results.size(); i++)
	{
		Matrix<double> pot_mat_temp(num_basis_rows, num_basis_rows, (double*)results[i]);
		Matrix<double> psi_mat_temp(num_basis_rows, num_basis_rows, ((double*)results[i]) + num_basis_rows * num_basis_rows);

		#pragma omp parallel for shared(num_basis_rows, pot_mat, psi_mat, pot_mat_temp, psi_mat_temp)
		for (int idx = 0; idx < num_basis_rows * num_basis_rows; idx++)
		{
			int row_idx = idx / num_basis_rows;
			int col_idx = idx % num_basis_rows;
			pot_mat.set(col_idx, row_idx, pot_mat.get(col_idx, row_idx) + pot_mat_temp.get(col_idx, row_idx));
			psi_mat.set(col_idx, row_idx, psi_mat.get(col_idx, row_idx) + psi_mat_temp.get(col_idx, row_idx));			
		} // end for each entry int matrices

		delete[] results[i];
	} // end for each remote result

	time_t end_time = time(NULL);
	if (f)
	{
		fprintf(f, "Finished combining results after %d seconds -- done at %s", int(end_time - start_time), asctime(localtime(&end_time)));
		fprintf(f, "Done\n");
		fclose(f);
	}

	return make_tuple(pot_mat, psi_mat);
} // end function DistributedAggregatePotentialContributions
