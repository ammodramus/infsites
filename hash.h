#ifndef HASH_HEADER
#define HASH_HEADER

#include <stdlib.h>
#include <stdint.h>

#define DEFAULT_HASH_TABLE_SIZE 10007
// (used to be 10000, but figured a prime number might reduce collisons)
// (no performance difference at all in a test run)

uint32_t SuperFastHash (const char * data, int len);
uint32_t jenkins_one_at_a_time_hash(char * key, size_t len);

#endif
