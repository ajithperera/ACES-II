#!/bin/ksh

test -f make.error || exit

egrep -v \
   -e "ANSI extension: IMPLICIT NONE statement" \
   -e "ANSI extension: nonStandard form of DO statement" \
   -e "ANSI extension: input contains lower case letters" \
   -e "ANSI extension: dollar sign or underscore in identifier" \
   -e "ANSI extension: symbolic name.* longer than 6 characters" \
   -e "ANSI extension: nonStandard intrinsic function \"rshift\"" \
   -e "Warning: -xarch=native has been explicitly specified" \
   make.error

#(implicit none|form of do|arch=native|lower case|underscore|longer than) make.error
