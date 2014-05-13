LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE := tfm

LOCAL_SRC_FILES := \
	src/addsub/fp_add.c \
	src/addsub/fp_add_d.c \
	src/addsub/fp_addmod.c \
	src/addsub/fp_cmp.c \
	src/addsub/fp_cmp_d.c \
	src/addsub/fp_cmp_mag.c \
	src/addsub/fp_sub.c \
	src/addsub/fp_sub_d.c \
	src/addsub/fp_submod.c \
	src/addsub/s_fp_add.c \
	src/addsub/s_fp_sub.c \
	src/bin/fp_radix_size.c \
	src/bin/fp_read_radix.c \
	src/bin/fp_read_signed_bin.c \
	src/bin/fp_read_unsigned_bin.c \
	src/bin/fp_reverse.c \
	src/bin/fp_signed_bin_size.c \
	src/bin/fp_s_rmap.c \
	src/bin/fp_toradix.c \
	src/bin/fp_to_signed_bin.c \
	src/bin/fp_to_unsigned_bin.c \
	src/bin/fp_unsigned_bin_size.c \
	src/bit/fp_cnt_lsb.c \
	src/bit/fp_count_bits.c \
	src/bit/fp_div_2.c \
	src/bit/fp_div_2d.c \
	src/bit/fp_lshd.c \
	src/bit/fp_mod_2d.c \
	src/bit/fp_rshd.c \
	src/divide/fp_div.c \
	src/divide/fp_div_d.c \
	src/divide/fp_mod.c \
	src/divide/fp_mod_d.c \
	src/exptmod/fp_2expt.c \
	src/exptmod/fp_exptmod.c \
	src/misc/fp_ident.c \
	src/misc/fp_set.c \
	src/mont/fp_montgomery_calc_normalization.c \
	src/mont/fp_montgomery_reduce.c \
	src/mont/fp_montgomery_setup.c \
	src/mul/fp_mul_2.c \
	src/mul/fp_mul_2d.c \
	src/mul/fp_mul.c \
	src/mul/fp_mul_comba_12.c \
	src/mul/fp_mul_comba_17.c \
	src/mul/fp_mul_comba_20.c \
	src/mul/fp_mul_comba_24.c \
	src/mul/fp_mul_comba_28.c \
	src/mul/fp_mul_comba_32.c \
	src/mul/fp_mul_comba_3.c \
	src/mul/fp_mul_comba_48.c \
	src/mul/fp_mul_comba_4.c \
	src/mul/fp_mul_comba_64.c \
	src/mul/fp_mul_comba_6.c \
	src/mul/fp_mul_comba_7.c \
	src/mul/fp_mul_comba_8.c \
	src/mul/fp_mul_comba_9.c \
	src/mul/fp_mul_comba.c \
	src/mul/fp_mul_comba_small_set.c \
	src/mul/fp_mul_d.c \
	src/mul/fp_mulmod.c \
	src/numtheory/fp_gcd.c \
	src/numtheory/fp_invmod.c \
	src/numtheory/fp_isprime.c \
	src/numtheory/fp_lcm.c \
	src/numtheory/fp_prime_miller_rabin.c \
	src/numtheory/fp_prime_random_ex.c \
	src/sqr/fp_sqr.c \
	src/sqr/fp_sqr_comba_12.c \
	src/sqr/fp_sqr_comba_17.c \
	src/sqr/fp_sqr_comba_20.c \
	src/sqr/fp_sqr_comba_24.c \
	src/sqr/fp_sqr_comba_28.c \
	src/sqr/fp_sqr_comba_32.c \
	src/sqr/fp_sqr_comba_3.c \
	src/sqr/fp_sqr_comba_48.c \
	src/sqr/fp_sqr_comba_4.c \
	src/sqr/fp_sqr_comba_64.c \
	src/sqr/fp_sqr_comba_6.c \
	src/sqr/fp_sqr_comba_7.c \
	src/sqr/fp_sqr_comba_8.c \
	src/sqr/fp_sqr_comba_9.c \
	src/sqr/fp_sqr_comba.c \
	src/sqr/fp_sqr_comba_generic.c \
	src/sqr/fp_sqr_comba_small_set.c \
	src/sqr/fp_sqrmod.c

LOCAL_C_INCLUDES := $(LOCAL_PATH)/src/headers

LOCAL_CFLAGS += -DTFM_ARM

ifeq ($(TARGET_ARCH_ABI),armeabi-v7a)
# Possible optimizations:
#  -ftree-vectorize: have GCC attempt to automatically vectorize loops
#  -ftree-vectorizer-verbose=2: verbose output during compile
# Note: not all V7-a targets support NEON!
LOCAL_ARM_NEON := true
LOCAL_CFLAGS += -DTFM_ARM_V7A -ftree-vectorize
else ifeq ($(TARGET_ARCH_ABI),armeabi)
LOCAL_CFLAGS += -DTFM_ARM_V5TE
else
LOCAL_CFLAGS += -DTFM_ARM_V4M
endif

include $(BUILD_STATIC_LIBRARY)
