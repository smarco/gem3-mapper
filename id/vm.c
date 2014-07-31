#include <sys/resource.h>
#include <sys/time.h>
#include <stdint.h>

#include <stdio.h>
#include <unistd.h>

void get_total_free_mem()
{
    uint64_t vm_size = 0;
    FILE *statm = fopen("/proc/self/statm", "r");
    fscanf(statm,"%ld", &vm_size);
    vm_size = (vm_size + 1) * 1024;

    printf("VirtualMemSize = %lu GB (used by self => /proc/self/statm)\n", vm_size/1024/1024);

    struct rlimit lim;
    getrlimit(RLIMIT_AS, &lim);
    printf("VirtualMemLimit = %lu GB\n", lim.rlim_cur/1024/1024/1024);

    const uint64_t total_virtual_memory = lim.rlim_cur - vm_size; // Bytes
    printf("VirtualMemAvail = %lu GB\n", total_virtual_memory/1024/1024/1024);

    // Physical
    const uint64_t page_size = sysconf(_SC_PAGESIZE)/1024; // KB
    const uint64_t num_pages = sysconf(_SC_AVPHYS_PAGES);
    const double available_memory = (double)(page_size*num_pages)/1024.0/1024.0; // GB
    printf("PageSize = %lu KB '_SC_PAGESIZE'\n", page_size);
    printf("NumPages = %lu    '_SC_AVPHYS_PAGES'\n", num_pages);
    printf("AvailableMemory = %f GB \n", available_memory);
}

int main() {
//  get_total_free_mem();

  // Open /proc/meminfo
  FILE* const meminfo = fopen("/proc/meminfo", "r");
  // Parse /proc/meminfo
  char *line = NULL;
  uint64_t size=0,line_length=0,chars_read=0;
  while ((chars_read=getline(&line,&line_length,meminfo))!=-1) {
    if (strncmp(line,"Cached:",7)==0) {
      sscanf(line,"%*s %lu",&size);
      printf("Cached: ---> %lu\n",size);
    }
  }

  return 0;
}


//#include <windows.h>
//
//size_t getTotalSystemMemory()
//{
//    MEMORYSTATUSEX status;
//    status.dwLength = sizeof(status);
//    GlobalMemoryStatusEx(&status);
//    return status.ullTotalPhys;
//}




//#include <ios>
//#include <iostream>
//#include <fstream>
//#include <string>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

//void process_mem_usage(double& vm_usage, double& resident_set)
//{
//   using std::ios_base;
//   using std::ifstream;
//   using std::string;
//
//   vm_usage     = 0.0;
//   resident_set = 0.0;
//
//   // 'file' stat seems to give the most reliable results
//   //
//   ifstream stat_stream("/proc/self/stat",ios_base::in);
//
//   // dummy vars for leading entries in stat that we don't care about
//   //
//   string pid, comm, state, ppid, pgrp, session, tty_nr;
//   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
//   string utime, stime, cutime, cstime, priority, nice;
//   string O, itrealvalue, starttime;
//
//   // the two fields we want
//   //
//   unsigned long vsize;
//   long rss;
//
//   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
//               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
//               >> utime >> stime >> cutime >> cstime >> priority >> nice
//               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
//
//   stat_stream.close();
//
//   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
//   vm_usage     = vsize / 1024.0;
//   resident_set = rss * page_size_kb;
//}
//
//int main()
//{
//   using std::cout;
//   using std::endl;
//
//   double vm, rss;
//   process_mem_usage(vm, rss);
//   cout << "VM: " << vm << "; RSS: " << rss << endl;
//}

