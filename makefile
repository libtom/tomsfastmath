#makefile for TomsFastMath
#
#
VERSION=0.14-develop

CFLAGS += -Wall -W -Wshadow -Isrc/headers

# Compiler and Linker Names
ifndef PREFIX
  PREFIX=
endif

ifeq ($(CC),cc)
  CC = $(PREFIX)gcc
endif
LD=$(PREFIX)ld
AR=$(PREFIX)ar
RANLIB=$(PREFIX)ranlib

ifndef MAKE
   MAKE=make
endif

ifeq ($V,1)
silent=
else
silent=@
endif

%.o: %.c
ifneq ($V,1)
	@echo "   * ${CC} $@"
endif
	${silent} ${CC} ${CFLAGS} -c $< -o $@

ifdef COMPILE_DEBUG
#debug
CFLAGS += -g3
else
ifndef IGNORE_SPEED

CFLAGS += -O3 -funroll-loops

#profiling
#PROF=-pg -g
#CFLAGS += $(PROF)

#speed
CFLAGS += -fomit-frame-pointer

endif
endif

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

ifndef LIBPATH
   LIBPATH=/usr/lib
endif

ifndef INCPATH
   INCPATH=/usr/include
endif

ifndef INSTALL_GROUP
   GROUP=wheel
else
   GROUP=$(INSTALL_GROUP)
endif

ifndef INSTALL_USER
   USER=root
else
   USER=$(INSTALL_USER)
endif

ifndef LIBNAME
	LIBNAME=libtfm.a
endif

default: $(LIBNAME)

$(OBJECTS): $(HEADERS)

$(LIBNAME): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@

install: $(LIBNAME)
	install -d -g $(GROUP) -o $(USER) $(DESTDIR)$(LIBPATH)
	install -d -g $(GROUP) -o $(USER) $(DESTDIR)$(INCPATH)
	install -g $(GROUP) -o $(USER) $(LIBNAME) $(DESTDIR)$(LIBPATH)
	install -g $(GROUP) -o $(USER) $(HEADERS_PUB) $(DESTDIR)$(INCPATH)

HEADER_FILES=$(notdir $(HEADERS_PUB))
uninstall:
	rm $(DESTDIR)$(LIBPATH)/$(LIBNAME)
	rm $(HEADER_FILES:%=$(DESTDIR)$(INCPATH)/%)

.PHONY: mtest
mtest: $(LIBNAME)
	cd mtest; CC="$(CC)" CFLAGS="$(CFLAGS) -I../" MAKE=${MAKE} ${MAKE} mtest

demo/test.o: CFLAGS+=-Wno-unused-result

.PHONY: test
test: $(LIBNAME) demo/test.o
	$(CC) $(CFLAGS) demo/test.o $(LIBNAME) $(PROF) -o test

test_standalone: CFLAGS+=-DTFM_DEMO_TEST_VS_MTEST=0

.PHONY: test_standalone
test_standalone: $(LIBNAME) demo/test.o
	$(CC) $(CFLAGS) demo/test.o $(LIBNAME) $(PROF) -o test

timing: $(LIBNAME) demo/timing.o
	$(CC) $(CFLAGS) demo/timing.o $(LIBNAME) $(PROF) -o timing

profiled:
	CC="$(CC)" PREFIX="${PREFIX} CFLAGS="${CFLAGS} -fprofile-generate" MAKE=${MAKE} ${MAKE} timing
	./test
	rm -f `find . -type f -name "*.o" | xargs`
	rm -f `find . -type f -name "*.a" | xargs`
	rm -f test
	CC=$(CC) PREFIX="${PREFIX} CFLAGS="${CFLAGS} -fprofile-use" MAKE=${MAKE} ${MAKE} timing
	
stest: $(LIBNAME) demo/stest.o
	$(CC) $(CFLAGS) demo/stest.o $(LIBNAME) -o stest

rsatest: $(LIBNAME) demo/rsa.o
	$(CC) $(CFLAGS) demo/rsa.o $(LIBNAME) -o rsatest

docdvi: tfm.tex
	cp tfm.tex tfm.bak
	touch --reference=tfm.tex tfm.bak
	(printf "%s" "\def\fixedpdfdate{"; date +'D:%Y%m%d%H%M%S%:z' -d @$$(stat --format=%Y tfm.tex) | sed "s/:\([0-9][0-9]\)$$/'\1'}/g") > tfm-deterministic.tex
	printf "%s\n" "\pdfinfo{" >> tfm-deterministic.tex
	printf "%s\n" "  /CreationDate (\fixedpdfdate)" >> tfm-deterministic.tex
	printf "%s\n}\n" "  /ModDate (\fixedpdfdate)" >> tfm-deterministic.tex
	cat tfm.tex >> tfm-deterministic.tex
	mv tfm-deterministic.tex tfm.tex
	touch --reference=tfm.bak tfm.tex
	touch tfm.ind
	latex tfm >/dev/null
	latex tfm >/dev/null
	makeindex tfm
	latex tfm >/dev/null

docs: docdvi
	latex tfm >/dev/null
	pdflatex tfm >/dev/null
	sed -b -i 's,^/ID \[.*\]$$,/ID [<0> <0>],g' tfm.pdf
	mv tfm.bak tfm.tex
	mv -f tfm.pdf doc

#This rule cleans the source tree of all compiled code, not including the pdf
#documentation.
clean:
	rm -f `find . -type f -name "*.o" | xargs`
	rm -f `find . -type f -name "*.lo"  | xargs`
	rm -f `find . -type f -name "*.a" | xargs`
	rm -f `find . -type f -name "*.la"  | xargs`
	rm -f `find . -type f -name "*.obj" | xargs`
	rm -f `find . -type f -name "*.lib" | xargs`
	rm -f `find . -type f -name "*.exe" | xargs`
	rm -f `find . -type f -name "*.gcov" | xargs`
	rm -f `find . -type f -name "*.gcda" | xargs`
	rm -f `find . -type f -name "*.gcno" | xargs`
	rm -f `find . -type f -name "*.il" | xargs`
	rm -f `find . -type f -name "*.dyn" | xargs`
	rm -f `find . -type f -name "*.dpi" | xargs`
	rm -rf `find . -type d -name "*.libs" | xargs`
	rm -f tfm.aux  tfm.dvi  tfm.idx  tfm.ilg  tfm.ind  tfm.lof  tfm.log  tfm.out  tfm.toc  test  test.exe
	cd mtest; MAKE=${MAKE} ${MAKE} clean

.PHONY: pre_gen
pre_gen:
	perl gen.pl
	sed -e 's/[[:blank:]]*$$//' mpi.c > pre_gen/mpi.c
	rm mpi.c

zipup:
	rm -rf ../tomsfastmath-$(VERSION) && rm -f ../tfm-$(VERSION).zip ../tfm-$(VERSION).tar.bz2 && \
	expsrc.sh -i . -o ../tomsfastmath-$(VERSION) --svntags --no-fetch -p '*.c' -p '*.h' && \
	MAKE=${MAKE} ${MAKE} -C ../tomsfastmath-$(VERSION) docs && \
	tar -c ../tomsfastmath-$(VERSION)/* | xz -cz > ../tfm-$(VERSION).tar.xz && \
	find ../tomsfastmath-$(VERSION)/ -type f -exec unix2dos -q {} \; && \
	zip -9 -r ../tfm-$(VERSION).zip ../tomsfastmath-$(VERSION)/* && \
	gpg -b -a ../tfm-$(VERSION).tar.xz && gpg -b -a ../tfm-$(VERSION).zip

new_file:
	bash updatemakes.sh

# $Source$
# $Revision$
# $Date$
