#include "matrix.h"
#include "complex.h"
#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include <stdexcept>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
extern "C" {
#include <CL/cl.h>
#include <time.h>
}

#define MAX_PLATFORMS 16
#define MAX_DEVICES 16

using std::runtime_error;

void checkResult(const char *name, cl_int res)
{
	if (res != CL_SUCCESS) {
		fprintf(stderr, "OpenCL call %s returned %d\n", name, res);
		throw runtime_error("OpenCL call failed");
	}
}

#define checkResultCall(func, ...) checkResult(#func, func(__VA_ARGS__))

void context_error_callback(const char *errinfo, const void *private_info, size_t private_size, void *user_data)
{
	fprintf(stderr, "OpenCL Context Error: %s\n", errinfo);
	fflush(stderr);
}

cl_program load_cl_program(cl_context context, const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "Failed to open kernel program file: %s\n", path);
        throw runtime_error("Failed to open kernel program");
    }
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    char * buffer = new char[fsize+1];
    if (fread(buffer, 1, fsize, f) != fsize) {
        fprintf(stderr, "Failed to read kernel program file: %s\n", path);
        delete[] buffer;
        throw runtime_error("Failed to read kernel program");
    }
    buffer[fsize] = 0;
    cl_int result;
    cl_program prog = clCreateProgramWithSource(context, 1, (const char **)&buffer, NULL, &result);
    delete[] buffer;
    checkResult("clCreateProgramWithSource", result);
    return prog;
}

int main(int argc, char **argv)
{
	cl_platform_id platform_ids[MAX_PLATFORMS];
	cl_uint num_platforms;
    //The first step in setting up an OpenCL context is selecting a platform (essentially an OpenCL implementation)
    //We get a list of aavailable platform as here and then print them
	checkResultCall(clGetPlatformIDs, MAX_PLATFORMS, platform_ids, &num_platforms);
	
	for (cl_uint i = 0; i < num_platforms; i++)
	{
		char name_buf[128];
		char vendor_buf[128];
		char extensions_buf[512];
		checkResultCall(clGetPlatformInfo, platform_ids[i], CL_PLATFORM_NAME, sizeof(name_buf), (void *)name_buf, NULL);
		checkResultCall(clGetPlatformInfo, platform_ids[i], CL_PLATFORM_VENDOR, sizeof(vendor_buf), (void *)vendor_buf, NULL);
		checkResultCall(clGetPlatformInfo, platform_ids[i], CL_PLATFORM_EXTENSIONS, sizeof(extensions_buf), (void *)extensions_buf, NULL);
		printf("%d: %s - %s - %s\n", i, name_buf, vendor_buf, extensions_buf);
	}
	
	//Currently we just select the first platform available
    //It would be better for this to be selectable with some kind of config file, environment variable or command line parameter
    //but most of the time there is only one platform available
	cl_platform_id selected_platform = platform_ids[0];
	
	cl_device_id device_ids[MAX_DEVICES];
	cl_uint num_devices;
    //Each plaform can have multiple devices (i.e. a system with 4 AMD GPUs would probably have one OpenCL platform, with 4 devices)
    //We get a list of available devices for our selected platform and print some info about them here.
	checkResultCall(clGetDeviceIDs, selected_platform, CL_DEVICE_TYPE_GPU|CL_DEVICE_TYPE_ACCELERATOR, MAX_DEVICES, device_ids, &num_devices);
	for (cl_uint i = 0; i < num_devices; i++)
	{
		char name_buf[128];
		checkResultCall(clGetDeviceInfo, device_ids[i], CL_DEVICE_NAME, sizeof(name_buf), (void *)name_buf, NULL);
		printf("%d: %s\n", i, name_buf);
	}
	
	const cl_context_properties context_props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform_ids[0], 0};
	cl_int result;
    //Now that we've selected a platform and have a list of devices we want to use (we are using all of them)
    //we can create our context
	cl_context context = clCreateContext(context_props, num_devices, device_ids, context_error_callback, NULL, &result);
	checkResult("clCreateContext", result);
	
    //Now we need to load and build our OpenCL code, we only have a single program in this little sample
    //but in a large project we might have multiple ones
	cl_program prog = load_cl_program(context, "kernels/matrix_mult.cl");
	
	result = clBuildProgram(prog, num_devices, device_ids, NULL, NULL, NULL);
	if (result == CL_BUILD_PROGRAM_FAILURE) {
		size_t log_size;
		clGetProgramBuildInfo(prog, device_ids[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		char *log = new char[log_size];
		clGetProgramBuildInfo(prog, device_ids[0], CL_PROGRAM_BUILD_LOG, log_size, log, &log_size);
		fprintf(stderr, "OpenCL build error: %s\n", log);
		delete[] log;
		return 1;
	} else {
		checkResult("clBuildProgram", result);
	}
	
    //now that our code is compiled, we can create a kernel from one of the kernel
    //functions defined in the code we loaded compiled
    //most realistic projects will have multiple kernels, but this demo only has one
	cl_kernel matrix_mult = clCreateKernel(prog, "matrix_mult", &result);
	checkResult("clCreateKernel", result);
	
	Matrix<double> a(1024,1024), b(1024,1024);
	Matrix<double> out(b.cols(), a.rows()), out_cl(b.cols(), a.rows());
	
    //Here we initialize some matrices with random data within a reasonable range
	for (int y = 0; y < a.rows(); y++)
	{
		for (int x = 0; x < a.cols(); x++)
		{
			a.set(x,y, rand() * 5.0 / RAND_MAX);
			b.set(x,y, rand() * 5.0 / RAND_MAX);
		}
	}
	struct timespec start, end;
	clock_gettime(CLOCK_REALTIME, &start);
    //Here we execute the multiply on the CPU so we have a baseline result to compare to
	a.mult(b, out);
	clock_gettime(CLOCK_REALTIME, &end);
	double elapsed = end.tv_sec - start.tv_sec;
	elapsed += ((double)end.tv_nsec)/ 1000000000.0 - ((double)start.tv_nsec)/ 1000000000.0;

	clock_gettime(CLOCK_REALTIME, &start);
    //Here we create some buffers for transferring data to and from the GPU
	cl_mem left = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(double) * a.cols() * a.rows(), (void *)a.raw(), &result);
	checkResult("clCreateBuffer", result);
	
	cl_mem right = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(double) * b.cols() * b.rows(), (void *)b.raw(), &result);
	checkResult("clCreateBuffer", result);
	
	cl_mem out_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(double) * out_cl.cols() * out_cl.rows(), NULL, &result);
	checkResult("clCreateBuffer", result);
	
    //Work is submitted to the GPU via a command queue, we create that here
	cl_command_queue command_queue = clCreateCommandQueue(context, device_ids[0], 0, &result);
	checkResult("clCreateBuffer", result);
	
    //Here we set the parameters that will be passed to our kernel function
    //these values will be used for each invocation of the kernel
	clSetKernelArg(matrix_mult, 0, sizeof(cl_mem), &left);
	clSetKernelArg(matrix_mult, 1, sizeof(cl_mem), &right);
	clSetKernelArg(matrix_mult, 2, sizeof(cl_mem), &out_buf);
	int mm_size = a.cols();
	int mm_left_stride = a.cols();
	int mm_right_stride = b.cols();
	int mm_out_stride = out_cl.cols();
	clSetKernelArg(matrix_mult, 3, sizeof(int), &mm_size);
	clSetKernelArg(matrix_mult, 4, sizeof(int), &mm_left_stride);
	clSetKernelArg(matrix_mult, 5, sizeof(int), &mm_right_stride);
	clSetKernelArg(matrix_mult, 6, sizeof(int), &mm_out_stride);
	
    //Here we define the number of dimmensions, the range of those dimensions and enqueue the execution of our kernel
	const size_t work_size[] = {(size_t)out_cl.cols(), (size_t)out_cl.rows()};
	cl_event mult_event;
	checkResultCall(clEnqueueNDRangeKernel, command_queue, matrix_mult, 2, NULL, work_size, NULL, 0, NULL, &mult_event);
    //Work in OpenCL normally executes asynchronously. mult_event now has a handle
    //that can be used to wait on the completion of the execution of our kernel on this data
	
    //While our input buffers were implicitly copied to the GPU on creation due to the flags used, the same is not true of our result buffer
    //we need to manually read that back
	checkResultCall(clEnqueueReadBuffer, command_queue, out_buf, CL_TRUE, 0, sizeof(double) * out_cl.cols() * out_cl.rows(), (void *)out_cl.raw(), 1, &mult_event, NULL);
    //Because we passed CL_TRUE for the blocking_read paramater above, this call will not return until the buffer read has been completeted
    //because we passed mult_event in as the event_wait_list parameter, the buffer copy will not start until the kernel is finished executing
	clock_gettime(CLOCK_REALTIME, &end);
	double elapsed2 = end.tv_sec - start.tv_sec;
	elapsed2 += ((double)end.tv_nsec)/ 1000000000.0 - ((double)start.tv_nsec)/ 1000000000.0;
	printf("CPU time: %f, GPU time: %f\n", elapsed, elapsed2);
    //Here we compare the results from the CPU and GPU
	for (int y = 0; y < a.rows(); y++)	
	{
		for (int x = 0; x < a.cols(); x++)
		{
			if (fabs(out.get(x, y) - out_cl.get(x, y)) > out.get(x, y) * 4.0 * DBL_EPSILON) {
				printf("%d, %d - difference %f\n", x, y, out.get(x, y) - out_cl.get(x, y));
			}
		}
	}
	return 0;
}