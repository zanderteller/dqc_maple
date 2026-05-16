#define  _CRT_SECURE_NO_WARNINGS

#include "maplec.h"
#include <vector>
#include <cstdio>
#include <stdexcept>
#include "matrix.h"
#include "array3d.h"
#include <complex>  // we were using our homegrown "complex.h"
#include "trajectories.h"
#include "distributed.h"
#include "net.h"
#include "operators.h"

#ifdef _WIN32
#include <Windows.h>
BOOL APIENTRY DllMain(HMODULE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved)
{
	return TRUE;
}
#define DLLEXPORT __declspec(dllexport)
#else
#define __declspec(A)
#define DLLEXPORT
#endif

using std::complex;
using std::tuple;
using std::get;

extern "C"
{
	/*
	print list of current secondary machines to Maple standard output
	
	note: this is a "Maple-aware" function

	Maple inputs
	* (none)
	
	Maple outputs
	* (none)
	*/
	DLLEXPORT ALGEB ListSecondariesC(MKernelVector kv, ALGEB* args)
	{
		vector<Secondary> rem = GetRemotes();
		if (rem.size() == 0)
		{
			MaplePrintf(kv, "no secondary machines\n");
		}
		else
		{
			MaplePrintf(kv, "%d secondary machines:\n", rem.size());
			for (int i = 0; i < rem.size(); i++)
				MaplePrintf(kv, "\t%s (%0.2f)\n", rem[i].address.c_str(), rem[i].weight);
		}

		return args[0];
	} // end function ListSecondariesC

	// return a count of currently available secondary machines
	DLLEXPORT int SecondaryCountC()
	{
		vector<Secondary> rem = GetRemotes();
		return rem.size();
	} // end function SecondaryCountC

	/*
	create and return DQC operators

	note: this is a "Maple-aware" function that deals with parameter checking and conversion of values
	to and from Maple's internal format.

	Maple inputs
	* mat: data matrix of doubles
	* nhambasis: use this many rows (starting from first row) as the basis
	* npotbasis: use this many rows (starting from first row) to build the potential
	* sigma: value of sigma
	* step: value of time step
	* mass: value of mass

	Maple outputs
	* eth: complex matrix -- the Hamiltonian time-evolution exponentiation operator
	* xops: vector of double matrices -- the position-expectation operators (one for each dimension)
	* sim_t: transpose of 'similarity' matrix -- used to to convert from representation in the 'raw'
	    basis of basis rows to representation in the orthonormal basis of eigenstates
	*/
	DLLEXPORT ALGEB MakeOperatorsC(MKernelVector kv, ALGEB* args)
	{
		try
		{
			if (7 != MapleNumArgs(kv, (ALGEB)args))
			{
				MapleRaiseError(kv, (char*)"invalid argument count");
			}

			ALGEB mat = args[1], nhambasis = args[2], npotbasis = args[3], sigma = args[4], step = args[5], mass = args[6], newham = args[7];
			if (!IsMapleRTable(kv, mat))
			{
				MapleRaiseError(kv, (char*)"mat must be an rtable");
			}
			if (!IsMapleInteger(kv, nhambasis))
			{
				MapleRaiseError(kv, (char*)"nhambasis must be an integer");
			}
			if (!IsMapleInteger(kv, npotbasis))
			{
				MapleRaiseError(kv, (char*)"npotbasis must be an integer");
			}
			if (!IsMapleNumeric(kv, sigma))
			{
				MapleRaiseError(kv, (char*)"sigma must be numeric");
			}
			if (!IsMapleNumeric(kv, step))
			{
				MapleRaiseError(kv, (char*)"step must be numeric");
			}
			if (!IsMapleNumeric(kv, mass))
			{
				MapleRaiseError(kv, (char*)"mass must be numeric");
			}
			if (!IsMapleInteger(kv, newham))
			{
				MapleRaiseError(kv, (char*)"newham must be an integer (to be converted to boolean)");
			}

			// convert integer/boolean arguments
			int nhambasisI, npotbasisI;
			nhambasisI = MapleToInteger32(kv, nhambasis);
			npotbasisI = MapleToInteger32(kv, npotbasis);
			bool newhamB = (bool)MapleToInteger32(kv, newham);

			// convert floating-point arguments
			double sigmaD, stepD, massD;
			sigmaD = MapleToFloat64(kv, sigma);
			stepD = MapleToFloat64(kv, step);
			massD = MapleToFloat64(kv, mass);

			// convert matrix argument
			RTableSettings rtsettings;
			RTableGetSettings(kv, &rtsettings, mat);
			if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
				rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
			{
				MapleRaiseError(kv, (char*)"mat must be a rectangular float64 matrix");
			}
			if (rtsettings.order != RTABLE_C)
			{
				rtsettings.order = RTABLE_C;
				rtsettings.foreign = FALSE;
				mat = RTableCopy(kv, &rtsettings, mat);
			}
			M_INT num_rows = RTableUpperBound(kv, mat, 1);
			M_INT num_cols = RTableUpperBound(kv, mat, 2);
			Matrix<double> matM(num_cols, num_rows, (double *)RTableDataBlock(kv, mat));

			// call MakeOperators
			tuple<Matrix<complex<double>>, vector<Matrix<double>>, Matrix<double>> results = MakeOperators(matM, nhambasisI, npotbasisI, sigmaD, stepD, massD, newhamB);
			Matrix<complex<double>> eth(get<0>(results));
			vector<Matrix<double>> xopvec(get<1>(results));
			Matrix<double> sim_t = (get<2>(results));

			// return eth as Maple output
			RTableGetDefaults(kv, &rtsettings);
			rtsettings.subtype = RTABLE_MATRIX;
			rtsettings.data_type = RTABLE_COMPLEX;
			rtsettings.storage = RTABLE_RECT;
			rtsettings.index_functions = ToMapleNULL(kv);
			rtsettings.order = RTABLE_FORTRAN;
			rtsettings.foreign = FALSE;
			rtsettings.num_dimensions = 2;
			M_INT eth_bounds[4] = { 1, eth.rows(), 1, eth.cols()};
			ALGEB eth_table = RTableCreate(kv, &rtsettings, NULL, eth_bounds);
			complex<double>* eth_data = (complex<double> *)RTableDataBlock(kv, eth_table);
			for (int col_idx = 0; col_idx < eth.cols(); col_idx++)
			{
				for (int row_idx = 0; row_idx < eth.rows(); row_idx++)
				{
					*(eth_data++) = eth.get(col_idx, row_idx);
				}
			}

			// return xopvec as Maple output
			RTableGetDefaults(kv, &rtsettings);
			rtsettings.subtype = RTABLE_ARRAY;
			rtsettings.data_type = RTABLE_DAG;
			rtsettings.storage = RTABLE_RECT;
			rtsettings.index_functions = ToMapleNULL(kv);
			rtsettings.order = RTABLE_FORTRAN;
			rtsettings.foreign = FALSE;
			rtsettings.num_dimensions = 1;
			M_INT bounds[4] = { 1, (M_INT)xopvec.size()};
			ALGEB xop_table = RTableCreate(kv, &rtsettings, NULL, bounds);

			// build table for each operator matrix and add it to the output array
			for (int idx = 0; idx < xopvec.size(); idx++)
			{
				Matrix<double> op_mat(xopvec[idx]);
				RTableGetDefaults(kv, &rtsettings);
				rtsettings.subtype = RTABLE_MATRIX;
				rtsettings.data_type = RTABLE_FLOAT64;
				rtsettings.storage = RTABLE_RECT;
				rtsettings.index_functions = ToMapleNULL(kv);
				rtsettings.order = RTABLE_FORTRAN;
				rtsettings.foreign = FALSE;
				rtsettings.num_dimensions = 2;
				M_INT bounds[4] = { 1, op_mat.rows(), 1, op_mat.cols() };
				ALGEB op_mat_table = RTableCreate(kv, &rtsettings, NULL, bounds);
				double* op_mat_data = (double*)RTableDataBlock(kv, op_mat_table);
				for (int col_idx = 0; col_idx < op_mat.cols(); col_idx++)
				{
					for (int row_idx = 0; row_idx < op_mat.rows(); row_idx++)
					{
						*(op_mat_data++) = op_mat.get(col_idx, row_idx);
					}
				}
                M_INT table_idx = (M_INT)idx + 1;
				RTableData val;
				val.dag = op_mat_table;
				RTableAssign(kv, xop_table, &table_idx, val);
			}

			// return sim_t as Maple output
			RTableGetDefaults(kv, &rtsettings);
			rtsettings.subtype = RTABLE_MATRIX;
			rtsettings.data_type = RTABLE_FLOAT64;
			rtsettings.storage = RTABLE_RECT;
			rtsettings.index_functions = ToMapleNULL(kv);
			rtsettings.order = RTABLE_FORTRAN;
			rtsettings.foreign = FALSE;
			rtsettings.num_dimensions = 2;
			M_INT sim_t_bounds[4] = { 1, sim_t.rows(), 1, sim_t.cols() };
			ALGEB sim_t_table = RTableCreate(kv, &rtsettings, NULL, sim_t_bounds);
			double* sim_t_data = (double*)RTableDataBlock(kv, sim_t_table);
			for (int col_idx = 0; col_idx < sim_t.cols(); col_idx++)
			{
				for (int row_idx = 0; row_idx < sim_t.rows(); row_idx++)
				{
					*(sim_t_data++) = sim_t.get(col_idx, row_idx);
				}
			}

			return ToMapleExpressionSequence(kv, 3, eth_table, xop_table, sim_t_table);

		}
		catch (std::runtime_error e)
		{
			MapleRaiseError(kv, (char*)e.what());
			return NewMapleExpressionSequence(kv, 3);
		}
	} // end function MakeOperatorsC

	/*
	build and return new frames added on to the frames passed in.

	note: this is a "Maple-aware" function that deals with parameter checking and conversion of values
	to and from Maple's internal format.

	Maple inputs
	* nframes: number of frames to add (note: we always add 1 fewer frames than requested)
	* traj1: 3-D array of existing frames
	* exph: complex matrix -- the Hamiltonian time-evolution exponentiation operator
	* xops: vector of double matrices -- the position-expectation operators (one for each dimension)
	* simdag: transpose ( or 'dag', for 'dagger') of 'similarity' matrix -- used to to convert from representation
		in the 'raw' basis of basis rows to representation in the orthonormal basis of eigenstates
	* nhambasis: use this many rows (starting from first row) as the basis
	* npotbasis: number of rows to use to build the potential [NEVER USED HERE -- ALREADY BAKED INTO THE OPERATORS]
	* sigma: value of sigma
	* step: value of time step [NEVER USED HERE -- ALREADY BAKED INTO THE OPERATORS]
	* mass: value of mass [NEVER USED HERE -- ALREADY BAKED INTO THE OPERATORS]

	Maple outputs
	* out: 3-D array with new frames added on to existing frames passed in
	*/
	DLLEXPORT ALGEB MakeTrajectoriesP(MKernelVector kv, ALGEB *args)
	{
		if (10 != MapleNumArgs(kv, (ALGEB)args))
		{
			MapleRaiseError(kv, (char *)"invalid argument count");
		}
		ALGEB nframes = args[1], traj1 = args[2], exph = args[3], xops = args[4],
			simdag = args[5], nhambasis = args[6], npotbasis = args[7], sigma = args[8],
			step = args[9], mass = args[10];
		if (!IsMapleInteger(kv, nframes))
		{
			MapleRaiseError(kv, (char *)"nframes must be an integer");
		}
		if (!IsMapleRTable(kv, traj1))
		{
			MapleRaiseError(kv, (char *)"traj1 must be an rtable");
		}
		if (!IsMapleRTable(kv, exph))
		{
			MapleRaiseError(kv, (char *)"exph must be an rtable");
		}
		if (!IsMapleRTable(kv, xops))
		{
			MapleRaiseError(kv, (char *)"xops must be an rtable");
		}
		if (!IsMapleRTable(kv, simdag))
		{
			MapleRaiseError(kv, (char *)"simdag must be an rtable");
		}
		if (!IsMapleInteger(kv, nhambasis))
		{
			MapleRaiseError(kv, (char *)"nhambasis must be an integer");
		}
		if (!IsMapleInteger(kv, npotbasis))
		{
			MapleRaiseError(kv, (char *)"npotbasis must be an integer");
		}
		if (!IsMapleNumeric(kv, sigma))
		{
			MapleRaiseError(kv, (char *)"sigma must be numeric");
		}
		if (!IsMapleNumeric(kv, step))
		{
			MapleRaiseError(kv, (char *)"step must be numeric");
		}
		if (!IsMapleNumeric(kv, mass))
		{
			MapleRaiseError(kv, (char *)"mass must be numeric");
		}
		int nframesI, nhambasisI, npotbasisI;
		nframesI = MapleToInteger32(kv, nframes);
		nhambasisI = MapleToInteger32(kv, nhambasis);
		npotbasisI = MapleToInteger32(kv, npotbasis);
		double sigmaD, stepD, massD;
		sigmaD = MapleToFloat64(kv, sigma);
		stepD = MapleToFloat64(kv, step);
		massD = MapleToFloat64(kv, mass);
		RTableSettings rtsettings;
		RTableGetSettings(kv, &rtsettings, traj1);
		if (rtsettings.subtype != RTABLE_ARRAY || rtsettings.data_type != RTABLE_FLOAT64 ||
			rtsettings.storage != RTABLE_RECT || rtsettings.num_dimensions != 3 ||
			!IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char *)"traj1 must be a rectangular 3-dimensional float64 array");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			traj1 = RTableCopy(kv, &rtsettings, traj1);
		}
		M_INT rows = RTableUpperBound(kv, traj1, 1) - RTableLowerBound(kv, traj1, 1) + 1;
		M_INT cols = RTableUpperBound(kv, traj1, 2) - RTableLowerBound(kv, traj1, 2) + 1;
		M_INT frames = RTableUpperBound(kv, traj1, 3) - RTableLowerBound(kv, traj1, 3) + 1;
		Array3D tFinal(cols, rows, frames + nframesI - 1);
		double *data = (double *)RTableDataBlock(kv, traj1);
		for (int z = 0; z < frames; z++)
		{
			for (int y = 0; y < rows; y++)
			{
				for (int x = 0; x < cols; x++)
				{
					tFinal.set(x, y, z, *(data++));
				}
			}
		}

		M_INT startFrame = frames - 1;
		RTableGetSettings(kv, &rtsettings, exph);
		if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_COMPLEX ||
			rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char *)"exph must be a rectangular complex matrix");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			exph = RTableCopy(kv, &rtsettings, exph);
		}
		rows = RTableUpperBound(kv, exph, 1);
		cols = RTableUpperBound(kv, exph, 2);
		Matrix<complex<double>> exphMat(cols, rows, (complex<double> *)RTableDataBlock(kv, exph));

		RTableGetSettings(kv, &rtsettings, xops);
		if (rtsettings.subtype != RTABLE_ARRAY || rtsettings.data_type != RTABLE_DAG ||
			rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char *)"xops must be a rectangular array of maple objects");
		}
		frames = RTableUpperBound(kv, xops, 1) - RTableLowerBound(kv, xops, 1) + 1;
		vector<Matrix<double> > xopsVec;
		xopsVec.reserve(frames);
		ALGEB *entries = (ALGEB *)RTableDataBlock(kv, xops);
		for (int i = 0; i < frames; i++)
		{
			ALGEB entry = entries[i];
			if (!IsMapleRTable(kv, entry))
			{
				MapleRaiseError(kv, (char *)"xops must only contain rtables");
			}
			RTableGetSettings(kv, &rtsettings, entry);
			
			if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
				rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions) ||
				rtsettings.order != RTABLE_C)
			{
				rtsettings.subtype = RTABLE_MATRIX;
				rtsettings.data_type = RTABLE_FLOAT64;
				rtsettings.storage = RTABLE_RECT;
				rtsettings.index_functions = ToMapleNULL(kv);
				rtsettings.order = RTABLE_C;
				rtsettings.foreign = FALSE;
				entry = RTableCopy(kv, &rtsettings, entry);
			}
			rows = RTableUpperBound(kv, entry, 1);
			cols = RTableUpperBound(kv, entry, 2);
			xopsVec.push_back(Matrix<double>(cols, rows, (double *)RTableDataBlock(kv, entry)));
		}

		RTableGetSettings(kv, &rtsettings, simdag);
		if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
			rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char *)"simdag must be a rectangular float64 matrix");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			simdag = RTableCopy(kv, &rtsettings, simdag);
		}
		rows = RTableUpperBound(kv, simdag, 1);
		cols = RTableUpperBound(kv, simdag, 2);
		Matrix<double> simdagMat(cols, rows, (double *)RTableDataBlock(kv, simdag));
		try
		{
			// note that DistributedAddTrajectories defaults to using just the primary machine when no secondaries are active,
			// so we can just call it here in either case (secondaries or no secondaries)
			DistributedAddTrajectories(startFrame, tFinal.depth() - 1, tFinal, exphMat, xopsVec, simdagMat, nhambasisI, npotbasisI, sigmaD, stepD, massD);
		}
		catch (std::runtime_error e)
		{
			MapleRaiseError(kv, (char *)e.what());
		}
		
		RTableGetDefaults(kv, &rtsettings);
		rtsettings.subtype = RTABLE_ARRAY;
		rtsettings.data_type = RTABLE_FLOAT64;
		rtsettings.storage = RTABLE_RECT;
		rtsettings.index_functions = ToMapleNULL(kv);
		rtsettings.order = RTABLE_FORTRAN;
		rtsettings.foreign = FALSE;
		rtsettings.num_dimensions = 3;
		M_INT bounds[6] = { 1, tFinal.rows(), 1, tFinal.cols(), 1, tFinal.depth() };
		ALGEB out = RTableCreate(kv, &rtsettings, NULL, bounds);
		data = (double *)RTableDataBlock(kv, out);
		for (int z = 0; z < tFinal.depth(); z++)
		{
			for (int x = 0; x < tFinal.cols(); x++)
			{
				for (int y = 0; y < tFinal.rows(); y++)
				{
					*(data++) = tFinal.get(x, y, z);
				}
			}
		}
		return out;
	} // end function MakeTrajectoriesP

	/*
	build and return new frames added on to the frames passed in.

	'auto' here means the code stops evolving individual rows as they stop moving.  (Stopping is defined by
	stopping_threshold -- see BuildFramesAuto for details.)

	note: this is a "Maple-aware" function that deals with parameter checking and conversion of values
	to and from Maple's internal format.

	Maple inputs
	* numframestobuild: number of new frames to build
	* currentframe: 3-D array of existing frames.  (note: does *not* need to include the actual first frame,
	    before any DQC evolution has taken place)
	* basisrows: matrix of basis rows from first frame (before any DQC evolution has taken place)
	* exph: complex matrix -- the Hamiltonian time-evolution exponentiation operator
	* xops: vector of double matrices -- the position-expectation operators (one for each dimension)
	* simdag: transpose ( or 'dag', for 'dagger') of 'similarity' matrix -- used to to convert from representation
		in the 'raw' basis of basis rows to representation in the orthonormal basis of eigenstates
	* sigma: value of sigma (for coherence, must be the same value that was used to build the operators)
	* stoppingthreshold: threshold for deciding when a row has stopped moving (see BuildFramesAuto for details)

	Maple outputs
	* frames_out: 3-D array of frames, with new frames added on to the existing frames passed in
	*/
	DLLEXPORT ALGEB BuildFramesAutoC(MKernelVector kv, ALGEB *args)
	{
		if (9 != MapleNumArgs(kv, (ALGEB)args))
		{
			MapleRaiseError(kv, (char *)"invalid argument count");
		}
		ALGEB numframestobuild = args[1], currentframe = args[2], basisrows = args[3], exph = args[4],
			xops = args[5], simdag = args[6], sigma = args[7], stoppingthreshold = args[8], rescale = args[9];

		if (!IsMapleInteger(kv, numframestobuild))
		{
			MapleRaiseError(kv, (char*)"numframestobuild must be an integer");
		}
		if (!IsMapleRTable(kv, currentframe))
		{
			MapleRaiseError(kv, (char *)"currentframe must be an rtable");
		}
		if (!IsMapleRTable(kv, basisrows))
		{
			MapleRaiseError(kv, (char*)"basisrows must be an rtable");
		}
		if (!IsMapleRTable(kv, exph))
		{
			MapleRaiseError(kv, (char *)"exph must be an rtable");
		}
		if (!IsMapleRTable(kv, xops))
		{
			MapleRaiseError(kv, (char *)"xops must be an rtable");
		}
		if (!IsMapleRTable(kv, simdag))
		{
			MapleRaiseError(kv, (char *)"simdag must be an rtable");
		}
		if (!IsMapleNumeric(kv, sigma))
		{
			MapleRaiseError(kv, (char *)"sigma must be numeric");
		}
		if (!IsMapleNumeric(kv, stoppingthreshold))
		{
			MapleRaiseError(kv, (char *)"stoppingthreshold must be numeric");
		}
		if (!IsMapleInteger(kv, rescale))
		{
			MapleRaiseError(kv, (char *)"rescale must be an integer (to be converted to boolean)");
		}

		// set up scalar arguments
		double sigmaD = MapleToFloat64(kv, sigma);
		double stoppingthresholdD = MapleToFloat64(kv, stoppingthreshold);
		int numframestobuildI = MapleToInteger32(kv, numframestobuild);
		bool rescaleB = (bool)MapleToInteger32(kv, rescale);

		// set up allframesArr: allocate memory for all frames and put input ('current') frame in the first frame
		RTableSettings rtsettings;
		RTableGetSettings(kv, &rtsettings, currentframe);
		if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
			rtsettings.storage != RTABLE_RECT || rtsettings.num_dimensions != 2 ||
			!IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char*)"currentframe must be a rectangular float64 matrix");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			currentframe = RTableCopy(kv, &rtsettings, currentframe);
		}
		M_INT num_rows = RTableUpperBound(kv, currentframe, 1);
		M_INT num_cols = RTableUpperBound(kv, currentframe, 2);
		if (num_cols == 1)
		{
			// corner case: a single row may be reinterpreted along the way as a single column
			// ZT 2022-03-09: TURNING THIS OFF TO BE ABLE TO RUN WITH A SINGLE RAW-DATA COLUMN
			// MapleRaiseError(kv, (char*)"currentframe must have multiple columns");
		}
		Array3D allframesArr(num_cols, num_rows, numframestobuildI + 1);
		// put current frame into first dim-3 slice of allframesArr
		double* data = (double *)RTableDataBlock(kv, currentframe);
		for (int y = 0; y < num_rows; y++)
		{
			for (int x = 0; x < num_cols; x++)
			{
				allframesArr.set(x, y, 0, *(data++));
			}
		}

		// set up basisrowsMat
		RTableGetSettings(kv, &rtsettings, basisrows);
		if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
			rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char*)"basisrows must be a rectangular float64 matrix");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			basisrows = RTableCopy(kv, &rtsettings, basisrows);
		}
		num_rows = RTableUpperBound(kv, basisrows, 1);
		num_cols = RTableUpperBound(kv, basisrows, 2);
		Matrix<double> basisrowsMat(num_cols, num_rows, (double *)RTableDataBlock(kv, basisrows));

		// set up exphMat
		RTableGetSettings(kv, &rtsettings, exph);
		if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_COMPLEX ||
			rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char*)"exph must be a rectangular complex matrix");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			exph = RTableCopy(kv, &rtsettings, exph);
		}
		num_rows = RTableUpperBound(kv, exph, 1);
		num_cols = RTableUpperBound(kv, exph, 2);
		Matrix<complex<double>> exphMat(num_cols, num_rows, (complex<double> *)RTableDataBlock(kv, exph));

		// set up xops vector of matrices
		RTableGetSettings(kv, &rtsettings, xops);
		if (rtsettings.subtype != RTABLE_ARRAY || rtsettings.data_type != RTABLE_DAG ||
			rtsettings.storage != RTABLE_RECT || rtsettings.num_dimensions != 1 ||
			!IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char*)"xops must be a rectangular 1-dimensional array of maple objects");
		}
		M_INT num_xops = RTableUpperBound(kv, xops, 1) - RTableLowerBound(kv, xops, 1) + 1;
		vector<Matrix<double>> xopsVec;
		xopsVec.reserve(num_xops); // note: we can't do this in the vector constructor because Matrix doesn't have a default constructor
		ALGEB* entries = (ALGEB*)RTableDataBlock(kv, xops);
		for (int idx = 0; idx < num_xops; idx++)
		{
			ALGEB entry = entries[idx];
			if (!IsMapleRTable(kv, entry))
			{
				MapleRaiseError(kv, (char*)"xops must only contain rtables");
			}
			RTableGetSettings(kv, &rtsettings, entry);
			if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
				rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
			{
				MapleRaiseError(kv, (char*)"xops must only contain rectangular float64 matrices");
			}
			if (rtsettings.order != RTABLE_C)
			{
				rtsettings.order = RTABLE_C;
				rtsettings.foreign = FALSE;
				entry = RTableCopy(kv, &rtsettings, entry);
			}
			num_rows = RTableUpperBound(kv, entry, 1);
			num_cols = RTableUpperBound(kv, entry, 2);
			xopsVec.push_back(Matrix<double>(num_cols, num_rows, (double *)RTableDataBlock(kv, entry)));
		} // end for (setting up each matrix entry in xops)

		// set up simdagMat
		RTableGetSettings(kv, &rtsettings, simdag);
		if (rtsettings.subtype != RTABLE_MATRIX || rtsettings.data_type != RTABLE_FLOAT64 ||
			rtsettings.storage != RTABLE_RECT || !IsMapleNULL(kv, rtsettings.index_functions))
		{
			MapleRaiseError(kv, (char *)"simdag must be a rectangular float64 matrix");
		}
		if (rtsettings.order != RTABLE_C)
		{
			rtsettings.order = RTABLE_C;
			rtsettings.foreign = FALSE;
			simdag = RTableCopy(kv, &rtsettings, simdag);
		}
		num_rows = RTableUpperBound(kv, simdag, 1);
		num_cols = RTableUpperBound(kv, simdag, 2);
		Matrix<double> simdagMat(num_cols, num_rows, (double *)RTableDataBlock(kv, simdag));

		try
		{
			// this function will distribute work to remote machines if any are available (via AddSecondary) -- otherwise,
			// the work is done locally.
			DistributedBuildFramesAuto(allframesArr, basisrowsMat, exphMat, xopsVec, simdagMat, sigmaD, stoppingthresholdD, rescaleB);
		}
		catch (std::runtime_error e)
		{
			// MapleRaiseError(kv, (char *)e.what());
			MaplePrintf(kv, "BuildFramesAutoC Failed: %s\n", e.what());
			return currentframe;
		}

		// set up frames_out
		RTableGetDefaults(kv, &rtsettings);
		rtsettings.subtype = RTABLE_ARRAY;
		rtsettings.data_type = RTABLE_FLOAT64;
		rtsettings.storage = RTABLE_RECT;
		rtsettings.index_functions = ToMapleNULL(kv);
		rtsettings.order = RTABLE_FORTRAN;
		rtsettings.foreign = FALSE;
		rtsettings.num_dimensions = 3;
		M_INT bounds[6] = { 1, allframesArr.rows(), 1, allframesArr.cols(), 1, allframesArr.depth() };
		ALGEB frames_out = RTableCreate(kv, &rtsettings, NULL, bounds);
		data = (double *)RTableDataBlock(kv, frames_out);
		// NOTE: we've created frames_out in Fortran order, so loop over columns before rows
		for (int z = 0; z < allframesArr.depth(); z++)
		{
			for (int x = 0; x < allframesArr.cols(); x++)
			{
				for (int y = 0; y < allframesArr.rows(); y++)
				{
					*(data++) = allframesArr.get(x, y, z);
				}
			}
		}

		return frames_out;
	} // end function BuildFramesAutoC

} // end extern "C"
