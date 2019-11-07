# TomsFastMath, a fast ISO C bignum library. -- Tom St Denis
# SPDX-License-Identifier: Unlicense

CFLAGS += -Wall -W -O3 -Wno-unused-result

default: mtest

.PHONY: mtest
mtest: mtest.o
	$(CC) $(CFLAGS) mtest.o -ltommath -o mtest

clean:
	rm -f *.o mtest *~ mtest.exe
