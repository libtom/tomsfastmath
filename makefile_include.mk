#
# Include makefile used by makefile + makefile.shared
#  (GNU make only)

ifndef INSTALL_CMD
$(error your makefile must define INSTALL_CMD)
endif
ifndef UNINSTALL_CMD
$(error your makefile must define UNINSTALL_CMD)
endif

#
# The Version of the library
#
VERSION=0.13.1-next
VERSION_LIB=1:0:0
VERSION_PC=0.13.1

# Compiler and Linker Names
ifndef CROSS_COMPILE
  CROSS_COMPILE:=
endif

ifeq ($(CC),cc)
  CC := $(CROSS_COMPILE)gcc
endif
LD:=$(CROSS_COMPILE)ld
AR:=$(CROSS_COMPILE)ar
RANLIB=$(CROSS_COMPILE)ranlib

ifndef MAKE
   MAKE=make
endif

#
# Compilation flags
#
# Note that we're extending the environments' CFLAGS.
# If you think that our CFLAGS are not nice you can easily override them
# by giving them as a parameter to make:
#  make CFLAGS="-I./src/headers/ -DLTC_SOURCE ..." ...
#

CFLAGS += -Wall -W -Wshadow -Isrc/headers

ifdef COMPILE_DEBUG
#debug
CFLAGS += -g3
else
ifndef IGNORE_SPEED

CFLAGS += -O3

PLATFORM := $(shell uname | sed -e 's/_.*//')
ifneq ($(PLATFORM), Darwin)
CFLAGS += -funroll-loops
endif

#profiling
#PROF=-pg -g
#CFLAGS += $(PROF)

#speed
CFLAGS += -fomit-frame-pointer

endif
endif

#
# (Un)Install related variables
#
DESTDIR  ?=
PREFIX   ?= /usr/local
LIBPATH  ?= $(PREFIX)/lib
INCPATH  ?= $(PREFIX)/include


#
# Build targets
#

default: $(LIBNAME)


demo/test.o: CFLAGS+=-Wno-unused-result

.PHONY: mtest
mtest: $(LIBNAME)
	CC="$(CC)" CFLAGS="$(CFLAGS) -I../" MAKE=${MAKE} ${MAKE} -C mtest/ mtest

.common_install: $(LIBNAME)
	install -d $(DESTDIR)$(LIBPATH)
	$(INSTALL_CMD) $(LIBNAME) $(DESTDIR)$(LIBPATH)/$(LIBNAME)
	install -d $(DESTDIR)$(INCPATH)
	install $(HEADERS_PUB) $(DESTDIR)$(INCPATH)


HEADER_FILES=$(notdir $(HEADERS_PUB))
.common_uninstall:
	$(UNINSTALL_CMD) $(DESTDIR)$(LIBPATH)/$(LIBNAME)
	rm $(HEADER_FILES:%=$(DESTDIR)$(INCPATH)/%)


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

docs:
	$(MAKE) -C doc/ $@ V=$(V)

doc/tfm.pdf:
	$(MAKE) -C doc/ tfm.pdf V=$(V)

.PHONY: pre_gen
pre_gen:
	perl gen.pl
	sed -e 's/[[:blank:]]*$$//' mpi.c > pre_gen/mpi.c
	rm mpi.c

zipup: pre_gen doc/tfm.pdf
	@# Update the index, so diff-index won't fail in case the pdf has been created.
	@#   As the pdf creation modifies tfm.tex, git sometimes detects the
	@#   modified file, but misses that it's put back to its original version.
	@git update-index --refresh
	@git diff-index --quiet HEAD -- || ( echo "FAILURE: uncommited changes or not a git" && exit 1 )
	rm -rf tomsfastmath-$(VERSION) tfm-$(VERSION).*
	@# files/dirs excluded from "git archive" are defined in .gitattributes
	git archive --format=tar --prefix=tomsfastmath-$(VERSION)/ HEAD | tar x
	mkdir -p tomsfastmath-$(VERSION)/doc
	cp doc/tfm.pdf tomsfastmath-$(VERSION)/doc/tfm.pdf
	tar -c tomsfastmath-$(VERSION)/ | xz -6e -c - > tfm-$(VERSION).tar.xz
	zip -9rq tfm-$(VERSION).zip tomsfastmath-$(VERSION)
	rm -rf tomsfastmath-$(VERSION)
	gpg -b -a tfm-$(VERSION).tar.xz
	gpg -b -a tfm-$(VERSION).zip

new_file:
	bash updatemakes.sh
