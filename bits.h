#ifndef BITS_UTILITY_H
#define BITS_UTILITY_H

#include <cstdint>

#define R2(n)    n,     n + 2*64,     n + 1*64,     n + 3*64
#define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )

const uint32_t lookuptable[256] = { R6(0), R6(2), R6(1), R6(3) };

const uint32_t MASKS[21] = {1,3,7,15,31,63,127,255,511,1023,2047,4095,8191,16383,32767,65535,131071,262143,524287,1048575,2097151};

inline uint32_t bits_reverse(int len, uint32_t val)
{
  return (lookuptable[(val      ) & 0xff] << 24 |
          lookuptable[(val >>  8) & 0xff] << 16 |
          lookuptable[(val >> 16) & 0xff] <<  8 |
          lookuptable[(val >> 24) & 0xff]
         ) >> (sizeof(uint32_t)*8 - len);
}

inline int hamming_weight(unsigned int n) {
    return __builtin_popcount(n);
}

inline int hamming_weight(unsigned long long n) {
    return __builtin_popcountll(n);
}

inline int hamming_distance(unsigned int x, unsigned int y)
{
    return __builtin_popcount(x ^ y);
}

inline int hamming_distance(unsigned long long x, unsigned long long y)
{
    return __builtin_popcountll(x ^ y);
}

#endif //BITS_UTILITY_H
