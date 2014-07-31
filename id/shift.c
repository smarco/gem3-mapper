#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int main(int argc,char** argv) {
  uint64_t i = 1ull;
  uint64_t pi = 1ull<<63;
  uint64_t mod = atoll(argv[1]);
  printf("mod=%lu left %lx right %lx\n",mod,(uint64_t)i<<mod,(uint64_t)pi>>mod);
}
