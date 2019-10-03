#!/bin/bash
# SPDX-License-Identifier: Unlicense
export a=`find . -type f | sort | grep "[.]/src" | grep "[.]c" | grep -v generators | sed -e 'sE\./EE' | sed -e 's/\.c/\.o/' | xargs`
perl ./parsenames.pl OBJECTS "$a"
export a=`find . -type f | grep [.]/src | grep [.]h | sed -e 'se\./ee' | xargs`
perl ./parsenames.pl HEADERS "$a"
