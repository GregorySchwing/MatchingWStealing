
#ifndef THREAD_FACTORY_H
#define THREAD_FACTORY_H
#include <vector>
#include <thread>
#include <mutex>
#include <iostream>
#include <cassert>
#include <iostream>
#include <cstdint>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "Graph.h"
#include "Frontier.h"
#include "Matcher.h"

struct BenchResult {
  size_t num_readers;
  size_t written_messages;
  std::vector<size_t> read_messages;
  std::chrono::milliseconds duration;
};

void print_results(const BenchResult &results) {
  size_t tot_read_messages = 0;
  int ni = 0;
  for (size_t n : results.read_messages){
    printf("reader: \t\t %d read_messages %zu \n", ni++,n);
    tot_read_messages += n;
  }

  printf("duration: \t\t %zu millseconds\n", results.duration.count());
  printf("num_readers: \t\t %zu reader\n", results.num_readers);
  printf("written_msgs: \t\t %f message/millseconds\n",
         (float)results.written_messages / (results.duration.count()));
  printf("avg_read_msgs: \t\t %f message/millseconds\n",
         (float)(tot_read_messages / results.num_readers) /
             results.duration.count());
  printf("\n");
}


#include <iostream>
#include <vector>
#include <thread>

class ThreadFactory {
public:
    template <typename IT, typename VT>
    static void create_threads_concurrentqueue_baseline(std::vector<std::thread> &threads,
                                                        unsigned num_threads,
                                                        std::vector<size_t> & read_messages,
                                                        Graph<IT, VT>& graph,
                                                        std::vector<Frontier<IT>*> & frontiers,
                                                        volatile bool &foundPath,
                                                        volatile bool &finished);
};

template <typename IT, typename VT>
void ThreadFactory::create_threads_concurrentqueue_baseline(std::vector<std::thread> &threads,
                                                                   unsigned num_threads,
                                                                   std::vector<size_t> & read_messages,
                                                                   Graph<IT, VT>& graph,
                                                                   std::vector<Frontier<IT>*> & frontiers,
                                                                   volatile bool &foundPath,
                                                                   volatile bool &finished) {

    for (unsigned i = 1; i < num_threads+1; ++i) {
        threads[i-1] = std::thread(&Matcher::hello_world, i);
        /*
        threads[i-1] = std::thread(&Matcher::search_persistent_baseline<IT,VT>,
                                    std::ref(graph),
                                    std::ref(frontiers),
                                    std::ref(read_messages),
                                    std::ref(foundPath),
                                    std::ref(finished),
                                    i);*/

        // Create a cpu_set_t object representing a set of CPUs. Clear it and mark
        // only CPU i as set.
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i, &cpuset);
        int rc = pthread_setaffinity_np(threads[i-1].native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
            std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
        }
    }
}


#endif // THREAD_FACTORY_H