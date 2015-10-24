#!/bin/bash -e
#
# Can be run with e.g. ./testme.sh "gcc-4.8 gcc-4.9", defaults to ./testme.sh "gcc"

_runtest()
{
  echo -n "Run test $1 $2"
  trap 'echo " - build not successful, errors are:" && cat test_gcc_errors.txt' INT TERM EXIT
  make clean > /dev/null
  CC="${1}" make test_standalone -j9 $2 > /dev/null 2>test_gcc_errors.txt
  trap - INT TERM EXIT
  local outfile="test_$(echo ${1}${2} | tr -d '\"' | tr ' ' '_').txt"
  trap 'echo " - tests not successful, failed at:" && tail ${outfile}' INT TERM EXIT
  ./test > ${outfile}
  echo " successful"
  trap - INT TERM EXIT
}

gccopt="-m32 -m64 -mx32"
if [ $# -ge 1 ]
then
  gccver=$1
else
  gccver="gcc"
fi

for gopt in ${gccopt};
do
	for gccv in ${gccver};
	do
		_runtest "${gccv} ${gopt}" "-f makefile.shared"
		_runtest "${gccv} ${gopt}" ""
	done
done

