#!/bin/bash

bash genlist.sh > tmplist

perl filter.pl makefile tmplist
sed -e 's/ *$//' < tmp.delme > makefile
rm -f tmp.delme

perl filter.pl makefile.shared tmplist
sed -e 's/ *$//' < tmp.delme > makefile.shared
rm -f tmp.delme

rm -f tmplist
rm -f tmp.delme
