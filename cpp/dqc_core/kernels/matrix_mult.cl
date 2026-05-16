//OpenCL doesn't support single precision by default, so we need to declare our
//use of the extension with a pragma in all kernels that need it
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
//the __kernel prefix indicates this function is a kernel (i.e. a unit of computation
//that will be executed a bunch of times with different parameters)
//the __global prefix on some of the parameters indicate that they are pointers to objects
//in the "global" address space. The "global" space is shared across all instances of a kernel
__kernel void matrix_mult(__global double *left, __global double *right, __global double *out, int size, int left_stride, int right_stride, int out_stride)
{
	//When executing a kernel, you tell OpenCL the number of dimmensions (from 1-3) and the size of those dimmensions
	//OpenCL will execute the kernel once for each unique combination of values within the ranges specified
	//get_global_id gets the value for a given dimmension (first is 0, second is 1 and third is 2)
	//For matrix multiplication, 2 dimmensions makes sense with the two values representing a coordinate
	//in the output matrix
	const int dstx = get_global_id(0);
	const int dsty = get_global_id(1);
	//Here we calculate the multiplication result for a single cell of the matrix
	double accum = 0.0;
	for (int src = 0; src < size; src++)
	{
		//the performance of this could be improved by using vectorized types
		//unfortunately the optimal vector size is GPU dependent
		//A library like CLBlast can generate optimal OpenCL kernels
		//based on GPU parameters at runtime and is probably the sensible tool for
		//doing basic BLAS operations in OpenCL
		accum += left[src + dsty * left_stride] * right[dstx + src * right_stride];
	}
	out[dstx + dsty * out_stride] = accum;
}
