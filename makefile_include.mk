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

CFLAGS += -O3 -funroll-loops

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
DESTDIR  ?= /usr/local
LIBPATH  ?= $(DESTDIR)/lib
INCPATH  ?= $(DESTDIR)/include


#
# Build targets
#

default: $(LIBNAME)


.common_install: $(LIBNAME)
	install -d $(LIBPATH)
	$(INSTALL_CMD) $(LIBNAME) $(LIBPATH)/$(LIBNAME)
	install -d $(INCPATH)
	install $(HEADERS_PUB) $(INCPATH)


HEADER_FILES=$(notdir $(HEADERS_PUB))
.common_uninstall:
	$(UNINSTALL_CMD) $(LIBPATH)/$(LIBNAME)
	rm $(HEADER_FILES:%=$(INCPATH)/%)


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
