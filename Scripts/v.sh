
# This script runs $EDITOR (or 'vi' if unset) on a file, directory, or stdin.
# If the file does not exist, it searches for the highest source authority.
# For example, if a directory contains the files: main.F, main.f, and main.o,
# then `v.sh main.` will run $EDITOR on main.F while ignoring main.f and main.o.
# - Anthony Yau

tmp_editor=${EDITOR:-vi}

if [ ! $1 ]
then $tmp_editor
else for arg in "$@"
     do
        if [ -d $arg ]
        then printf "\nv: cannot edit \"$arg\" since the directory exists\n\n"
        else if   [ -f ${arg}    ]; then $tmp_editor ${arg}
             elif [ -f ${arg}cpp ]; then $tmp_editor ${arg}cpp
             elif [ -f ${arg}c   ]; then $tmp_editor ${arg}c
             elif [ -f ${arg}F95 ]; then $tmp_editor ${arg}F95
             elif [ -f ${arg}F90 ]; then $tmp_editor ${arg}F90
             elif [ -f ${arg}F   ]; then $tmp_editor ${arg}F
             elif [ -f ${arg}fpp ]; then $tmp_editor ${arg}fpp
             elif [ -f ${arg}f95 ]; then $tmp_editor ${arg}f95
             elif [ -f ${arg}f90 ]; then $tmp_editor ${arg}f90
             elif [ -f ${arg}f   ]; then $tmp_editor ${arg}f
             else $tmp_editor $arg
             fi # [ -f ${arg} ]
        fi # [ -d $arg ]
     done
fi # [ ! $1 ]

