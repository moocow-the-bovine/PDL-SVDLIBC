#include <stdio.h>

int main(void) {
  printf("sizeof(short) = %d\n", sizeof(short));
  printf("sizeof(int) = %d\n", sizeof(int));
  printf("sizeof(long) = %d\n", sizeof(long));
  printf("sizeof(long long) = %d\n", sizeof(long long));
  printf("sizeof(void *) = %d\n", sizeof(void*));
  printf("sizeof(float) = %d\n", sizeof(float));
  printf("sizeof(double) = %d\n", sizeof(double));
  /* not valid c: bummer
#if sizeof(long)==sizeof(long long)
  printf(": long==long long\n");
#else
  printf(": long!=long long\n");
#endif
  */
  return 0;
}
