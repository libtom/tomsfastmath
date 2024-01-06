#include <tfm_private.h>

#ifdef TFM_NO_STDLIB
/* rely on compiler to provide optimized version */
void fp_memset(void *restrict dst, unsigned char c, size_t n) {
  uint8_t *d = dst;
  for (size_t i = 0; i < n; ++i) d[i] = c;
}
#endif
