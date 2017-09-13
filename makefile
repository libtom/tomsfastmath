#makefile for TomsFastMath
#
#
ifndef LIBNAME
	LIBNAME=libtfm.a
endif

INSTALL_CMD = install
UNINSTALL_CMD = rm


include makefile_include.mk

ifeq ($V,1)
silent=
else
silent=@
endif

ifeq ($(COVERAGE),1)
CFLAGS += -fprofile-arcs -ftest-coverage
LDFLAGS += -lgcov
LIB_PRE = -Wl,--whole-archive
LIB_POST = -Wl,--no-whole-archive
endif

%.o: %.c
ifneq ($V,1)
	@echo "   * ${CC} $@"
endif
	${silent} ${CC} ${CFLAGS} -c $< -o $@

#START_INS
OBJECTS=src/addsub/fp_add.o src/addsub/fp_add_d.o src/addsub/fp_addmod.o src/addsub/fp_cmp.o \
src/addsub/fp_cmp_d.o src/addsub/fp_cmp_mag.o src/addsub/fp_sub.o src/addsub/fp_sub_d.o \
src/addsub/fp_submod.o src/addsub/s_fp_add.o src/addsub/s_fp_sub.o src/bin/fp_radix_size.o \
src/bin/fp_read_radix.o src/bin/fp_read_signed_bin.o src/bin/fp_read_unsigned_bin.o \
src/bin/fp_reverse.o src/bin/fp_signed_bin_size.o src/bin/fp_s_rmap.o src/bin/fp_toradix.o \
src/bin/fp_toradix_n.o src/bin/fp_to_signed_bin.o src/bin/fp_to_unsigned_bin.o \
src/bin/fp_unsigned_bin_size.o src/bit/fp_cnt_lsb.o src/bit/fp_count_bits.o src/bit/fp_div_2.o \
src/bit/fp_div_2d.o src/bit/fp_lshd.o src/bit/fp_mod_2d.o src/bit/fp_rshd.o src/divide/fp_div.o \
src/divide/fp_div_d.o src/divide/fp_mod.o src/divide/fp_mod_d.o src/exptmod/fp_2expt.o \
src/exptmod/fp_exptmod.o src/misc/fp_ident.o src/misc/fp_rand.o src/misc/fp_set.o \
src/mont/fp_montgomery_calc_normalization.o src/mont/fp_montgomery_reduce.o \
src/mont/fp_montgomery_setup.o src/mul/fp_mul_2.o src/mul/fp_mul_2d.o src/mul/fp_mul.o \
src/mul/fp_mul_comba_12.o src/mul/fp_mul_comba_17.o src/mul/fp_mul_comba_20.o src/mul/fp_mul_comba_24.o \
src/mul/fp_mul_comba_28.o src/mul/fp_mul_comba_32.o src/mul/fp_mul_comba_3.o src/mul/fp_mul_comba_48.o \
src/mul/fp_mul_comba_4.o src/mul/fp_mul_comba_64.o src/mul/fp_mul_comba_6.o src/mul/fp_mul_comba_7.o \
src/mul/fp_mul_comba_8.o src/mul/fp_mul_comba_9.o src/mul/fp_mul_comba.o \
src/mul/fp_mul_comba_small_set.o src/mul/fp_mul_d.o src/mul/fp_mulmod.o src/numtheory/fp_gcd.o \
src/numtheory/fp_invmod.o src/numtheory/fp_isprime.o src/numtheory/fp_isprime_ex.o \
src/numtheory/fp_lcm.o src/numtheory/fp_prime_miller_rabin.o src/numtheory/fp_prime_random_ex.o \
src/sqr/fp_sqr.o src/sqr/fp_sqr_comba_12.o src/sqr/fp_sqr_comba_17.o src/sqr/fp_sqr_comba_20.o \
src/sqr/fp_sqr_comba_24.o src/sqr/fp_sqr_comba_28.o src/sqr/fp_sqr_comba_32.o src/sqr/fp_sqr_comba_3.o \
src/sqr/fp_sqr_comba_48.o src/sqr/fp_sqr_comba_4.o src/sqr/fp_sqr_comba_64.o src/sqr/fp_sqr_comba_6.o \
src/sqr/fp_sqr_comba_7.o src/sqr/fp_sqr_comba_8.o src/sqr/fp_sqr_comba_9.o src/sqr/fp_sqr_comba.o \
src/sqr/fp_sqr_comba_generic.o src/sqr/fp_sqr_comba_small_set.o src/sqr/fp_sqrmod.o

HEADERS_PUB:=src/headers/tfm.h
HEADERS=src/headers/tfm_private.h $(HEADERS_PUB)

#END_INS

$(OBJECTS): $(HEADERS)

$(LIBNAME): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@

install: .common_install

uninstall: .common_uninstall

.PHONY: test
test: $(LIBNAME) demo/test.o
	$(CC) $(CFLAGS) demo/test.o $(LIB_PRE) $(LIBNAME) $(LIB_POST) $(PROF) -o test

test_standalone: CFLAGS+=-DTFM_DEMO_TEST_VS_MTEST=0

.PHONY: test_standalone
test_standalone: $(LIBNAME) demo/test.o
	$(CC) $(CFLAGS) demo/test.o $(LIB_PRE) $(LIBNAME) $(LIB_POST) $(PROF) -o test

testme: test mtest
	./mtest/mtest -15 | ./test

timing: $(LIBNAME) demo/timing.o
	$(CC) $(CFLAGS) demo/timing.o $(LIBNAME) $(PROF) -o timing

profiled:
	CC="$(CC)" CROSS_COMPILE="${CROSS_COMPILE} CFLAGS="${CFLAGS} -fprofile-generate" MAKE=${MAKE} ${MAKE} timing
	./test
	rm -f `find . -type f -name "*.o" | xargs`
	rm -f `find . -type f -name "*.a" | xargs`
	rm -f test
	CC=$(CC) CROSS_COMPILE="${CROSS_COMPILE} CFLAGS="${CFLAGS} -fprofile-use" MAKE=${MAKE} ${MAKE} timing

# target that pre-processes all coverage data
lcov-single-create:
	lcov --capture --no-external --directory src -q --output-file coverage_std.info

# target that removes all coverage output
cleancov-clean:
	rm -f `find . -type f -name "*.info" | xargs`
	rm -rf coverage/

# generates html output from all coverage_*.info files
lcov:
	lcov `find -name 'coverage_*.info' -exec echo -n " -a {}" \;` -o coverage.info -q 2>/dev/null
	genhtml coverage.info --output-directory coverage -q

# combines all necessary steps to create the coverage from a single testrun with e.g.
lcov-single:
	$(MAKE) cleancov-clean
	$(MAKE) lcov-single-create
	$(MAKE) lcov


#make the code coverage of the library
coverage: CFLAGS += -fprofile-arcs -ftest-coverage
coverage: LDFLAGS += -lgcov
coverage: LIB_PRE = -Wl,--whole-archive
coverage: LIB_POST = -Wl,--no-whole-archive

coverage: testme
	$(MAKE) lcov-single

stest: $(LIBNAME) demo/stest.o
	$(CC) $(CFLAGS) demo/stest.o $(LIBNAME) -o stest

rsatest: $(LIBNAME) demo/rsa.o
	$(CC) $(CFLAGS) demo/rsa.o $(LIBNAME) -o rsatest

# ref:         $Format:%D$
# git commit:  $Format:%H$
# commit time: $Format:%ai$
