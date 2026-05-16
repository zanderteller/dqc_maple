#ifndef NET_H_
#define NET_H_

#ifdef _WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#include <string>
#include <cstddef>
#include <mutex>

using std::vector;
using std::string;

struct Secondary
{
	string const address;
	float  weight;
	int    sockfd;

	Secondary(char const* address, float weight) : address(address), weight(weight), sockfd(-1) {};
};

static const char* SOCKET_PORT = "6881"; // both primary and secondary machines communicate on this port

//Adds a secondary server to the secondary list with the given weight
//A weight of 1.0 means the secondary will receive the same number of work units as the primary node
//higher weights will receive more work units and lower weights will receive fewer units
void AddSecondary(char const *address, float weight);

// calls AddSecondary (for use by external code)
extern "C" { DLLEXPORT void AddSecondaryC(char const *address, float weight); }

// Returns a copy of the remotes vector
vector<Secondary> GetRemotes();

//given a total number of work units, produces a set of per-node unit counts in out
//note that the local node is always returned in index 0, the first remote node starts at 1
void DistributeUnits(int units, vector<int> &out);

//repeatedly calls send on a socket until length bytes are sent or the connection drops
void SendAll(int sockfd, char const *data, size_t length);

//repeatedly calls recv on a socket until length bytes are received or the connection drops
void RecvAll(int sockfd, char *data, size_t length);

//listens on the given port and then stops listening once a connection received
//returns the newly connected socket file descriptor
int ListenOnce(char const *port);

/**
 * This class primarily exists to automatically acquire the mutex that keeps the keepalive
 * thread from trying to send a ping in the middle of a remote operation
 * A CommScope object should be created right before sending a command and should be kept
 * in scope until results have been received
 */
class CommScope
{
	std::unique_lock<std::mutex> lock;
	
public:
	CommScope();
	~CommScope();
	//Sends data to the secondary indicated by the remote parameter
	void SendData(int remote, char const *data, size_t length);
	//Receives results from all secondaries. The same work_units vector that was used
	//to distribute work to the individual units should be passed to this method.
	//if size_is_per_unit is true, size specifies the size in bytes of a single work unit.
	//is size_is_per_unit is false, size is the total size of the result from each secondary.
	void ReceiveResults(vector<char *> buffers, vector<int> work_units, size_t size, bool size_is_per_unit);
};

#endif //NET_H_
