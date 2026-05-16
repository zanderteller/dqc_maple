#include <windows.h>
#include <vector>
#include <atomic>
#include <curl.h>

static atomic_bool logger_running;
static CRITICAL_SECTION log_lock;
static CONDITION_VARIABLE log_cond;
static vector<string *> entries;
static size_type read_pos;


static DWORD WINAPI log_thread(LPVOID param)
{
	CURL *curl = curl_easy_init();
	curl_easy_setopt(curl, CURLOPT_URL, LOG_URL);
	struct curl_slist *headers = curl_slist_append(NULL, "Content-Type: application/json");
	curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
	string payload;
	for (;;)
	{
		EnterCriticalSection(&lock_lock);
			while (read_pos == entries.size())
			{
				SleepConditionVariableCS(&log_cond, &lock_lock, INVINITE);
			}
			payload = "[";
			for (size_t n = entries.size(); read_pos < n; read_pos++)
			{
				if (n) {
					payload.push_back(',');
				}
				payload.append(entries[read_pos]);
				delete entries[read_pos];
			}
			payload.push_back('}');
		LeaveCriticalSection(&lock_lock);
		curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, payload.size());
		curl_easy_setopt(curl, CURLOPT_POSTFIELDS, payload.c_str());
		curl_easy_perform(curl);
	}
	return 0;
}

void log(vector<std::pair<string, string>> const &properties)
{
	//TODO compare_and_exchange
	bool expected = false;
	if (logger_running.compare_exchange_strong(expected, true)) {
		InitializeCriticalSection(&log_lock);
		InitializeConditionVariable(&log_cond);
		CreateThread(NULL, 0, log_thread, NULL, 0, NULL);
	}
	string *entry = new string();
	entry->push_back('{');
	for (size_type i = 0, n = properties.size(); i < n; i++)
	{
		auto &kv(properties[i]);
		if (i) {
			entry->push_back(',');
		}
		entry->push_back('"');
		entry->append(kv.first);
		entry->append("\":\"");
		entry->append(kv.second);
		entry->push_back('"');
	}
	entry->push_back('}');
	EnterCriticalSection(&log_lock);
		entries.push_back(entry);
		if (start_index > entries.size() / 2) {
			size_type i, n;
			for (i = 0, n = entries.size(); read_pos < n; i++, read_pos++)
			{
				entries[i] = entries[read_pos];
			}
			entries.resize(i);
			read_pos = 0;
		}
	LeaveCriticalSection(&log_lock);
	WakeConditionVariable(&log_cond);
}
