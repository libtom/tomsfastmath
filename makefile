#makefile for TomsFastMath
#
#
CFLAGS += -Wall -W -Wshadow -I./ -O3 -funroll-all-loops

#profiling
#PROF=-pg -g
#CFLAGS += $(PROF)

#speed
CFLAGS += -fomit-frame-pointer

VERSION=0.03

default: libtfm.a

OBJECTS = \
fp_set.o \
\
fp_rshd.o fp_lshd.o fp_div_2d.o fp_mod_2d.o fp_mul_2d.o fp_2expt.o \
fp_mul_2.o fp_div_2.o  \
\
fp_cnt_lsb.o \
\
fp_add.o fp_sub.o fp_mul.o fp_sqr.o fp_div.o fp_mod.o \
s_fp_add.o s_fp_sub.o \
\
fp_cmp_d.o fp_add_d.o fp_sub_d.o fp_mul_d.o fp_div_d.o fp_mod_d.o \
fp_addmod.o fp_submod.o fp_mulmod.o fp_sqrmod.o fp_invmod.o \
fp_gcd.o fp_lcm.o fp_prime_miller_rabin.o fp_isprime.o \
fp_prime_random_ex.o fp_mul_comba.o fp_sqr_comba.o \
\
fp_montgomery_setup.o fp_montgomery_calc_normalization.o fp_montgomery_reduce.o \
\
fp_exptmod.o \
\
fp_cmp.o fp_cmp_mag.o \
\
fp_unsigned_bin_size.o fp_read_unsigned_bin.o fp_to_unsigned_bin.o \
fp_signed_bin_size.o fp_read_signed_bin.o fp_to_signed_bin.o \
fp_read_radix.o fp_toradix.o fp_radix_size.o fp_count_bits.o fp_reverse.o fp_s_rmap.o \
\
fp_ident.o 

libtfm.a: $(OBJECTS)
	$(AR) $(ARFLAGS) libtfm.a $(OBJECTS)
	ranlib libtfm.a

mtest/mtest: mtest/mtest.c
	cd mtest ; make mtest

test: libtfm.a demo/test.o mtest/mtest
	$(CC) $(CFLAGS) demo/test.o libtfm.a $(PROF) -o test

stest: libtfm.a demo/stest.o 
	$(CC) demo/stest.o libtfm.a -o stest

docdvi: tfm.tex
	touch tfm.ind
	latex tfm >/dev/null
	latex tfm >/dev/null
	makeindex tfm
	latex tfm >/dev/null

docs: docdvi
	latex tfm >/dev/null
	dvipdf tfm
	mv -f tfm.pdf doc

clean:
	rm -f $(OBJECTS) *.a demo/*.o test tfm.aux  tfm.dvi  tfm.idx  tfm.ilg  tfm.ind  tfm.lof  tfm.log  tfm.toc stest *~
	cd mtest ; make clean

zipup: docs clean
	perl gen.pl ; mv mpi.c pre_gen/ ; \
	cd .. ; rm -rf tfm* tomsfastmath-$(VERSION) ; mkdir tomsfastmath-$(VERSION) ; \
	cp -R ./tomsfastmath/* ./tomsfastmath-$(VERSION)/ ; \
	tar -c tomsfastmath-$(VERSION)/* | bzip2 -9vvc > tfm-$(VERSION).tar.bz2 ; \
	zip -9r tfm-$(VERSION).zip tomsfastmath-$(VERSION)/*
