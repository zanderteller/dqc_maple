#include <vector>
#include <ctime>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <tuple>
#include "net.h"
#include <complex>
#include "matrix.h"
#include "array3d.h"
#include "trajectories.h"
#include "operators.h"
#include "distributed.h"

#include <cstdio>
#include <cstdint>
#include <stdexcept>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <Windef.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#define close closesocket
#else
extern "C" {
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
}
#endif

using std::complex;
using std::runtime_error;
using std::chrono::system_clock;
using std::chrono::minutes;
using std::tuple;
using std::get;

// this is a shortcut: enum represents each value with an integer, starting at zero -- we make use of this
// fact below
// NOTE: the order matters here -- it must be the same as the order in which the parameters are assembled
// in the Distributed<Foo> function that actually sends the data
enum ADDTRAJ_PARAMS {
	ATP_START_FRAME,
	ATP_FRAMES,
	ATP_NHAMBASIS,
	ATP_NPOTBASIS,
	ATP_EXPH_COLS,
	ATP_EXPH_ROWS,
	ATP_SIMDAG_COLS,
	ATP_SIMDAG_ROWS,
	ATP_XOPS_SIZE,
	ATP_COLS,
	ATP_NUM_INT_PARAMS
};

// this is a shortcut: enum represents each value with an integer, starting at zero -- we make use of this
// fact below
// NOTE: the order matters here -- it must be the same as the order in which the parameters are assembled
// in the Distributed<Foo> function that actually sends the data
enum BUILDAUTO_PARAMS {
	BAP_NUM_FRAMES,
	BAP_NUM_COLS,
	BAP_NUM_BASIS_ROWS,
	BAP_EXPH_ROWS,
	BAP_EXPH_COLS,
	BAP_SIMDAG_ROWS,
	BAP_SIMDAG_COLS,
	BAP_NUM_INT_PARAMS
};

// this is a shortcut: enum represents each value with an integer, starting at zero -- we make use of this
// fact below
// NOTE: the order matters here -- it must be the same as the order in which the parameters are assembled
// in the Distributed<Foo> function that actually sends the data
enum AGGPOT_PARAMS {
	APP_NUM_BASIS_ROWS,
	APP_NUM_COLS,
	APP_NUM_INT_PARAMS
};

template <class T, class U>
void secondaryMatrixMult(int sockfd, Matrix<T> &a, Matrix<U> &b)
{
	RecvAll(sockfd, (char *)a.raw(), sizeof(T) * a.cols() * a.rows());
	RecvAll(sockfd, (char *)b.raw(), sizeof(U) * b.cols() * b.rows());
	auto result(a * b);
	SendAll(sockfd, (char *)result.raw(), sizeof(T) * result.rows() * result.cols());
}

bool done;
std::mutex m;
std::condition_variable done_notify;

void keepalive(int sockfd)
{
	std::unique_lock<std::mutex> lock(m);
	while (!done_notify.wait_until(lock, system_clock::now() + minutes(1), [](){ return done; }))
	{
		char cmd = 0;
		send(sockfd, &cmd, 1, 0);
		size_t bytes = recv(sockfd, &cmd, 1, 0);
		if (bytes != 1 || cmd != 0)
		{
			printf("Error receiving keepalive ping -  bytes received %d, command: %d\n", (int)bytes, cmd);
			return;
		}
		time_t t = time(NULL);
		printf("Received ping at %s", asctime(localtime(&t)));
	}
}


// Designed to run as a standalone executable on a seconary machine.
// The program will run forever, listening for a job from the primary server.
int main(int argc, char ** argv)
{
	for (;;) // waiting for connection
	{
		printf("Secondary listening on port %s\n", SOCKET_PORT);
		int sockfd = ListenOnce(SOCKET_PORT);
		puts("Received connection");
		try
		{
			for(;;) // waiting for command
			{
				char cmd;
				RecvAll(sockfd, &cmd, sizeof(cmd));
				printf("Received command: %d\n", cmd);
				switch (cmd)
				{
				case CMD_PING: {
					send(sockfd, &cmd, 1, 0);
					printf("PING\n\n");
					break;
				}
				case CMD_MATRIXMULT: {
					char eltypes[2];
					RecvAll(sockfd, eltypes, sizeof(eltypes));
					int32_t dims[4];
					RecvAll(sockfd, (char *)&dims, sizeof(dims));
					if (eltypes[0] == TYPE_DOUBLE) {
						//TODO: assert eltypes[1] also TYPE_DOUBLE
						Matrix<double> a(dims[0], dims[1]), b(dims[2], dims[3]);
						secondaryMatrixMult(sockfd, a, b);
					} else if (eltypes[1] == TYPE_DOUBLE) {
						Matrix<complex<double>> a(dims[0], dims[1]);
						Matrix<double> b(dims[2], dims[3]);
						secondaryMatrixMult(sockfd, a, b);
					} else {
						//TODO: assert we actually got TYPE_COMPLEX
						Matrix<complex<double>> a(dims[0], dims[1]), b(dims[2], dims[3]);
						secondaryMatrixMult(sockfd, a, b);
					}
					break;
				} // end case CMD_MATRIXMULT
				case CMD_ADDTRAJECTORIES: {
					int32_t int_params[ATP_NUM_INT_PARAMS];
					double double_params[3];
					RecvAll(sockfd, (char *)int_params, sizeof(int_params));
					RecvAll(sockfd, (char *)double_params, sizeof(double_params));
					int32_t rows;
					RecvAll(sockfd, (char *)&rows, sizeof(rows));
					
					Array3D tFinal(int_params[ATP_COLS], int_params[ATP_NHAMBASIS] + rows, int_params[ATP_FRAMES] + 1);
					RecvAll(sockfd, (char *)tFinal.raw(), tFinal.rows() * tFinal.cols() * (int_params[ATP_START_FRAME] + 1) * sizeof(double));
					
					Matrix<complex<double>> expH(int_params[ATP_EXPH_COLS], int_params[ATP_EXPH_ROWS]);
					RecvAll(sockfd, (char *)expH.raw(), expH.cols() * expH.rows() * sizeof(complex<double>));
					
					Matrix<double> simDag(int_params[ATP_SIMDAG_COLS], int_params[ATP_SIMDAG_ROWS]);
					RecvAll(sockfd, (char *)simDag.raw(), simDag.cols() * simDag.rows() * sizeof(double));
					time_t t = time(NULL);
					printf("CMD_ADDTRAJECTORIES: Received parameters at %s", asctime(localtime(&t)));
					
					vector<Matrix<double> > xops;
					xops.reserve(int_params[ATP_XOPS_SIZE]);
					for (int i = 0; i < int_params[ATP_XOPS_SIZE]; i++)
					{
						int32_t dims[2];
						RecvAll(sockfd, (char *)dims, sizeof(dims));
						Matrix<double> xop(dims[0], dims[1]);
						RecvAll(sockfd, (char *)xop.raw(), xop.cols() * xop.rows() * sizeof(double));
						xops.push_back(xop);
					}
					time_t start = time(NULL);
					done = false;
					std::thread pinger(keepalive, sockfd);
					printf("CMD_ADDTRAJECTORIES: Calling AddTrajectories for %d frames on data of size %d cols, %d rows at %d\n", int_params[ATP_FRAMES], tFinal.cols(), tFinal.rows(), (int)start);
					AddTrajectories(
						int_params[ATP_START_FRAME], int_params[ATP_FRAMES], tFinal, expH, xops, simDag, int_params[ATP_NHAMBASIS],
						int_params[ATP_NPOTBASIS], double_params[0], double_params[1], double_params[2]
					);
					time_t calc_end = time(NULL);
					printf("CMD_ADDTRAJECTORIES: Calculated results took %d seconds -- done at %s", (int)(calc_end - start), asctime(localtime(&calc_end)));
					{
						std::lock_guard<std::mutex> lock(m);
						done = true;
					}
					done_notify.notify_one();
					pinger.join();
					cmd = 1;
					send(sockfd, &cmd, 1, 0);
					for (int i = 1 + int_params[ATP_START_FRAME]; i < tFinal.depth(); i++)
					{
						Matrix<double> frame(tFinal.view2D(0, tFinal.cols(), int_params[ATP_NHAMBASIS], rows, i));
						SendAll(sockfd, (char *)frame.raw(), frame.cols() * frame.rows() * sizeof(double));
					}
					time_t sent = time(NULL);
					printf("CMD_ADDTRAJECTORIES: Results sent after %d seconds -- done at %s", (int)(sent-start), asctime(localtime(&sent)));
					break;
				} // end case CMD_ADDTRAJECTORIES
				case CMD_BUILDFRAMESAUTO: {
					// see DistributedBuildFramesAuto for the order in which parameter values are added to the data stream

					// receive double parameters
					double double_params[2];
					RecvAll(sockfd, (char*)double_params, sizeof(double_params));

					// receive int parameters
					int32_t int_params[BAP_NUM_INT_PARAMS];
					RecvAll(sockfd, (char*)int_params, sizeof(int_params));

					// receive boolean parameter
					bool rescale = false;
					RecvAll(sockfd, (char*)rescale, sizeof(rescale));

					// receive num_rows parameter
					int32_t num_rows;
					RecvAll(sockfd, (char*)&num_rows, sizeof(num_rows));

					// receive basis rows
					Matrix<double> basis_rows(int_params[BAP_NUM_COLS], int_params[BAP_NUM_BASIS_ROWS]);
					RecvAll(sockfd, (char*)basis_rows.raw(), basis_rows.rows() * basis_rows.cols() * sizeof(double));

					// receive rows to evolve and put them into an Array3D with memory allocated for all frames to be built
					// note that we're not filling all of the memory we allocate here -- only the first slice in dimension 3
					Array3D all_frames(int_params[BAP_NUM_COLS], num_rows, int_params[BAP_NUM_FRAMES]);
					RecvAll(sockfd, (char*)all_frames.raw(), all_frames.rows() * all_frames.cols() * sizeof(double));

					// receive operators: exph, then simdag, then xops
					Matrix<complex<double>> exph(int_params[BAP_EXPH_COLS], int_params[BAP_EXPH_ROWS]);
					RecvAll(sockfd, (char*)exph.raw(), exph.rows() * exph.cols() * sizeof(complex<double>));
					Matrix<double> simdag(int_params[BAP_SIMDAG_COLS], int_params[BAP_SIMDAG_ROWS]);
					RecvAll(sockfd, (char*)simdag.raw(), simdag.rows() * simdag.cols() * sizeof(double));
					vector<Matrix<double>> xops;
					xops.reserve(int_params[BAP_NUM_COLS]); // xops has one operator matrix for each column/dimension
					for (int i = 0; i < int_params[BAP_NUM_COLS]; i++)
					{
						int32_t dims[2];
						RecvAll(sockfd, (char*)dims, sizeof(dims));
						Matrix<double> xop(dims[1], dims[0]);
						RecvAll(sockfd, (char*)xop.raw(), xop.rows() * xop.cols() * sizeof(double));
						xops.push_back(xop);
					}
					time_t t = time(NULL);
					printf("CMD_BUILDFRAMESAUTO: Received parameters at %s", asctime(localtime(&t)));

					time_t start = time(NULL);
					done = false;
					std::thread pinger(keepalive, sockfd);
					printf("CMD_BUILDFRAMESAUTO: Calling BuildFramesAuto to build %d frames on data of size %d rows, %d cols at %d\n",
							int_params[BAP_NUM_FRAMES] - 1, all_frames.rows(), all_frames.cols(), (int)start);
					// double_params has [sigma, stopping_threshold] in it
					BuildFramesAuto(all_frames, basis_rows, exph, xops, simdag, double_params[0], double_params[1], rescale);
					time_t calc_end = time(NULL);
					printf("CMD_BUILDFRAMESAUTO: Calculating results took %d seconds -- done at %s", (int)(calc_end - start), asctime(localtime(&calc_end)));
					{
						std::lock_guard<std::mutex> lock(m);
						done = true;
					}
					done_notify.notify_one();
					pinger.join();

					// ping the primary machine
					cmd = 1;
					send(sockfd, &cmd, 1, 0);

					// send results
					// note: we send back only the new frames, *not* the initial frame we received
					for (int frame_idx = 1; frame_idx < all_frames.depth(); frame_idx++)
					{
						Matrix<double> frame(all_frames.view2D(0, all_frames.cols(), 0, all_frames.rows(), frame_idx));
						SendAll(sockfd, (char*)frame.raw(), frame.rows() * frame.cols() * sizeof(double));
					}
					time_t sent = time(NULL);
					printf("CMD_BUILDFRAMESAUTO: Results sent after %d seconds -- done at %s", (int)(sent - start), asctime(localtime(&sent)));
					printf("CMD_BUILDFRAMESAUTO: Done\n\n");
					break;
				} // end case CMD_BUILDFRAMESAUTO
				case CMD_AGGREGATEPOTENTIALCONTRIBUTIONS: {
					// see DistributedAggregatePotentialContributions for the order in which parameter values are added to the data stream

					// receive the only double parameter, sigma
					double sigma;
					RecvAll(sockfd, (char*)&sigma, sizeof(sigma));

					// receive int parameters
					int32_t int_params[APP_NUM_INT_PARAMS];
					RecvAll(sockfd, (char*)int_params, sizeof(int_params));

					// receive num_rows parameter
					int32_t num_rows;
					RecvAll(sockfd, (char*)&num_rows, sizeof(num_rows));

					// receive basis rows
					Matrix<double> basis_rows(int_params[APP_NUM_COLS], int_params[APP_NUM_BASIS_ROWS]);
					RecvAll(sockfd, (char*)basis_rows.raw(), basis_rows.rows() * basis_rows.cols() * sizeof(double));

					// receive data rows to process
					Matrix<double> data_rows(int_params[APP_NUM_COLS], num_rows);
					RecvAll(sockfd, (char*)data_rows.raw(), data_rows.rows() * data_rows.cols() * sizeof(double));

					time_t t = time(NULL);
					printf("CMD_AGGREGATEPOTENTIALCONTRIBUTIONS: Received parameters at %s", asctime(localtime(&t)));

					time_t start = time(NULL);
					done = false;
					std::thread pinger(keepalive, sockfd);
					printf("CMD_AGGREGATEPOTENTIALCONTRIBUTIONS: Calling AggregatePotentialContributions for data of size %d rows, %d cols at %d\n",
							data_rows.rows(), data_rows.cols(), (int)start);

					tuple<Matrix<double>, Matrix<double>> results = AggregatePotentialContributions(data_rows, basis_rows, sigma);
					Matrix<double> pot_mat(get<0>(results));
					Matrix<double> psi_mat(get<1>(results));

					time_t calc_end = time(NULL);
					printf("CMD_AGGREGATEPOTENTIALCONTRIBUTIONS: Calculating results took %d seconds -- done at %s", (int)(calc_end - start), asctime(localtime(&calc_end)));
					{
						std::lock_guard<std::mutex> lock(m);
						done = true;
					}
					done_notify.notify_one();
					pinger.join();

					// ping the primary machine
					cmd = 1;
					send(sockfd, &cmd, 1, 0);

					// send results (note: the order matters -- send pot_mat back first)
					SendAll(sockfd, (char*)pot_mat.raw(), pot_mat.rows() * pot_mat.cols() * sizeof(double));
					SendAll(sockfd, (char*)psi_mat.raw(), psi_mat.rows() * psi_mat.cols() * sizeof(double));

					time_t sent = time(NULL);
					printf("CMD_AGGREGATEPOTENTIALCONTRIBUTIONS: Results sent after %d seconds -- done at %s", (int)(sent - start), asctime(localtime(&sent)));
					printf("CMD_AGGREGATEPOTENTIALCONTRIBUTIONS: Done\n\n");
					break;
				} // end case CMD_AGGREGATEPOTENTIALCONTRIBUTIONS
				default:
					throw runtime_error("Received unrecognized command");
				}
			} // end for (ever -- waiting for command)
		}
		catch(runtime_error e)
		{
			close(sockfd);
			fprintf(stderr, "Exception: %s\n\n", e.what());
		}
	} // end for (ever -- trying to start socket listener and wait for command)
	return 0;
} // end function main
