
# This script should behave just like $PAGER.

# This script runs $PAGER (or 'more' if unset) on a file, directory, or stdin.
# If the file does not exist, it searches for the highest source authority.
# For example, if a directory contains the files: main.F, main.f, and main.o,
# then `m.sh main.` will run $PAGER on main.F while ignoring main.f and main.o.
# - Anthony Yau

tmp_pager=${PAGER:-more}
tmp_ls="/bin/ls -F"

if [ "x$1" = "x" ]
then $tmp_pager
else for arg in "$@"
     do
        if [ -d $arg ]
        then printf "\nListing contents of $arg\n\n"; $tmp_ls $arg
        else if   [ -f ${arg}    ]; then $tmp_pager ${arg}
             elif [ -f ${arg}cpp ]; then $tmp_pager ${arg}cpp
             elif [ -f ${arg}c   ]; then $tmp_pager ${arg}c
             elif [ -f ${arg}F95 ]; then $tmp_pager ${arg}F95
             elif [ -f ${arg}F90 ]; then $tmp_pager ${arg}F90
             elif [ -f ${arg}F   ]; then $tmp_pager ${arg}F
             elif [ -f ${arg}fpp ]; then $tmp_pager ${arg}fpp
             elif [ -f ${arg}f95 ]; then $tmp_pager ${arg}f95
             elif [ -f ${arg}f90 ]; then $tmp_pager ${arg}f90
             elif [ -f ${arg}f   ]; then $tmp_pager ${arg}f
             else echo "No such primary source file exists: $arg"
             fi # [ -f ${arg} ]
        fi # [ -d $arg ]
     done
fi # [ ! $1 ]

