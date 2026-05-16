#define _CRT_SECURE_NO_WARNINGS
#include <stdexcept>
#include <algorithm>
#include <thread>
#include <chrono>
#include <condition_variable>
#include "matrix.h"
#include "complex.h"
#include "array3d.h"
#include "distributed.h"
#include "net.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <Windef.h>
#include <winsock2.h>
#include <ws2tcpip.h>
typedef long long ssize_t;
#define close closesocket

#include <cstdlib>

//Windows requires explicit initialization and cleanup of its socket library
void SocketCleanup()
{
	WSACleanup();
}

void SocketInit()
{
	static bool did_init;
	static WSAData wsa_junk;
	if (!did_init)
	{
		did_init = true;
		WSAStartup(MAKEWORD(2,2), &wsa_junk);
		atexit(SocketCleanup);
	}
}

#else
extern "C" {
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/select.h>
#include <netdb.h>
#include <unistd.h>
#include <netinet/tcp.h>
#define SocketInit()
}
#endif

using std::vector;
using std::runtime_error;
using std::chrono::system_clock;
using std::chrono::minutes;

static vector<Secondary> remotes;
static float total_weight = 1.0;
static std::thread *pinger;
//this mutex is used to keep the keepalive thread from sending pings during a pending remote command
static std::mutex comm_mutex;
static bool keepalive_exit;
static std::condition_variable exit_notify;

// returns a copy of the remotes vector
vector<Secondary> GetRemotes()
{
	vector<Secondary> r = remotes; // make a copy
	return r;
} // end function GetRemotes

//this function is run in a separate thread and sends periodic keepalive messages to all connected secondaries
//ping failures will cause the secondary to be marked as disconnected. An attempt to reconnect will be made on
//the next remote request
//NOTE: It does not handle keepalive during long running remote commands which is less than ideal. Fixing that
//requires more sophisticated framing/command handling. So instead a separate keepalive runs on the secondary
//during a distributed call
void keepalive()
{
	std::unique_lock<std::mutex> lock(comm_mutex);
	while (!exit_notify.wait_until(lock, system_clock::now() + minutes(1), [](){ return keepalive_exit; }))
	{
		for (int i = 0; i < remotes.size(); i++)
		{
			if (remotes[i].sockfd != -1)
			{
				char cmd = CMD_PING;
				send(remotes[i].sockfd, &cmd, 1, 0);
				size_t bytes = recv(remotes[i].sockfd, &cmd, 1, 0);
				if (bytes != 1 || cmd != CMD_PING)
				{
					//ping failed, close connection
					close(remotes[i].sockfd);
					remotes[i].sockfd = -1;
				}
			}
		}
	}
}

void AddSecondary(char const *address, float weight)
{
	remotes.push_back(Secondary(address, weight));
	total_weight += weight;
}

extern "C" {
	DLLEXPORT void AddSecondaryC(char const *address, float weight)
	{
		AddSecondary(address, weight);
	}

	DLLEXPORT void ClearSecondariesC()
	{
		if (pinger)
		{
			{
				std::unique_lock<std::mutex> lock(comm_mutex);
				keepalive_exit = true;
				exit_notify.notify_one();
			}
			pinger->join();
			delete pinger;
			pinger = nullptr;
		}
		for (int i = 0; i < remotes.size(); i++)
		{
			if (remotes[i].sockfd != -1)
			{
				close(remotes[i].sockfd);
			}
		}
		remotes.clear();
		total_weight = 1.0;
	}
}

void DistributeUnits(int units, vector<int> &out)
{
	// corner case: if number of work units is small, just do all the work locally
	if (units < 10 * total_weight)
	{
		out.push_back(units);
		return;
	}

	int cur_processor = int(units / total_weight);
	out.push_back(cur_processor);
	int remaining = units - cur_processor;
	float max_weight = 1.0f;
	auto best(out.begin());
	for (auto cur(remotes.begin()), end(remotes.end()); cur != end; ++cur)
	{
		cur_processor = int(cur->weight * units / total_weight);
		out.push_back(cur_processor);
		remaining -= cur_processor;
		if (cur->weight >= max_weight)
		{
			max_weight = cur->weight;
			best = out.end();
			--best;
		}
	}
	if (remaining > 0)
	{
		(*best) += remaining;
	}
}

void SendAll(int sockfd, char const *data, size_t length)
{
	while (length > 0)
	{
		ssize_t sent = send(sockfd, data, length, 0);
		if (sent <= 0) {
			fprintf(stderr, "send failed: %s (%d)\n", strerror(errno), errno);
			throw runtime_error("Failed to send data to remote");
		}
		length -= sent;
		data += sent;
	}
}

void RecvAll(int sockfd, char *data, size_t length)
{
	while (length > 0)
	{
		ssize_t recvd = recv(sockfd, data, length, 0);
		if (recvd <= 0) {
			
			throw runtime_error("Failed to receive data from remote");
		}
		length -= recvd;
		data += recvd;
	}
}

CommScope::CommScope() :
lock(comm_mutex)
{
}

CommScope::~CommScope()
{
	if (!pinger)
	{
		pinger = new std::thread(keepalive);
	}
}

void CommScope::SendData(int remote, char const *data, size_t length)
{
	if (remotes[remote].sockfd == -1) {
		SocketInit();
		struct addrinfo hints = {};
		struct addrinfo *result;
		hints.ai_family = AF_UNSPEC;
		hints.ai_socktype = SOCK_STREAM;
		hints.ai_flags = 0;
		hints.ai_protocol = 0;
		int status = getaddrinfo(remotes[remote].address.c_str(), SOCKET_PORT, &hints, &result);
		if (status != 0) {
			throw runtime_error("Failed to get address of remote");
		}
		for (struct addrinfo *cur = result; cur; cur = cur->ai_next)
		{
			remotes[remote].sockfd = socket(cur->ai_family, cur->ai_socktype, cur->ai_protocol);
			if (remotes[remote].sockfd == -1) {
				continue;
			}
			
			if (!connect(remotes[remote].sockfd, cur->ai_addr, cur->ai_addrlen)) {
				break;
			}
			close(remotes[remote].sockfd);
			remotes[remote].sockfd = -1;
		}
		if (remotes[remote].sockfd == -1) {
			throw runtime_error("Failed to connect to remote");
		}
	}
	SendAll(remotes[remote].sockfd, data, length);
}

/*
receive results from the remote machines as they come in

inputs:
* buffers: vector of pre-allocated char arrays to store received results
* work_units: vector of number of work units (i.e., data rows) for each remote machine
* size: size of each result
* size_is_per_unit: boolean -- if true, total size of result from each remote machine is size * work_units[i].
    if false, total size of result from each remote machine is just size.

output:
* none
*/
void CommScope::ReceiveResults(vector<char *> buffers, vector<int> work_units, size_t size, bool size_is_per_unit)
{
	vector<size_t> received(remotes.size(), 0);
	fd_set pending;
	FD_ZERO(&pending);
	int num_pending =  0, highest = -1;
	for (int i = 0; i < remotes.size(); i++)
	{
		FD_SET(remotes[i].sockfd, &pending);
		int val = 1;
		setsockopt(remotes[i].sockfd, IPPROTO_TCP, TCP_NODELAY, (const char *)&val, sizeof(val));
		highest = std::max(remotes[i].sockfd, highest);
		num_pending++;
	}
	while  (num_pending)
	{
		int status = select(highest + 1, &pending, 0, 0, 0);
		if (status <= 0)
		{
			throw runtime_error("Error while waiting on data from remotes");
		}
		highest = -1;
		num_pending = 0;
		
		// if size_is_per_unit is true, full_size will be overwritten in each iteration of the loop
		size_t full_size = size;
		
		for (int i = 0; i < remotes.size(); i++)
		{
			if (size_is_per_unit)
				full_size = work_units[i+1] * size;
			if (FD_ISSET(remotes[i].sockfd, &pending)) {
				char cmd;
				ssize_t read;
				if (received[i]) {
					cmd = 1;
				} else {
					read = recv(remotes[i].sockfd, &cmd, 1, 0);
					if (read != 1) {
						throw runtime_error("Error while reading from remote socket");
					}
				}
				if (cmd) {
					read = recv(remotes[i].sockfd, buffers[i] + received[i], full_size - received[i], 0);
					if (read <= 0) {
						throw runtime_error("Error while reading from remote socket");
					}
					received[i] += read;
				} else {
					//keep alive message, just echo it back
					send(remotes[i].sockfd, &cmd, 1, 0);
				}
				FD_CLR(remotes[i].sockfd, &pending);
			}
			if (received[i] < full_size) {
				FD_SET(remotes[i].sockfd, &pending);
				highest = std::max(remotes[i].sockfd, highest);
				num_pending++;
			}
		}
	}
} // end method ReceiveResults

int ListenOnce(char const *port)
{
	SocketInit();
	struct addrinfo hints = {}, *result;
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;
	
	getaddrinfo(NULL, port, &hints, &result);
	struct addrinfo *pref = NULL;
	for (struct addrinfo *cur=result; cur; cur=cur->ai_next)
	{
		if (cur->ai_family == AF_INET) {
			pref = cur;
			break;
		} else if (!pref) {
			pref = cur;
		}
	}
	//TODO: error checking
	int listen_sock = socket(pref->ai_family, pref->ai_socktype, pref->ai_protocol);
	bind(listen_sock, pref->ai_addr, pref->ai_addrlen);
	listen(listen_sock, 1);
	int client_sock = accept(listen_sock, 0, 0);
	close(listen_sock);
	return client_sock;
}
