#include <tfm_private.h>

#ifdef TFM_NO_STDLIB
/* rely on compiler to provide optimized version */
void fp_memcpy(void *restrict dst, const void *restrict src, size_t n) {
  uint8_t *d = dst;
  const uint8_t *s = src;
  for (size_t i = 0; i < n; ++i) d[i] = s[i];
}
#endif
