
# Usage: mvfF.sh [-r] [old extension] [new extension]

# This script scans files of the type *.${old extension} for strings that
# are generally used by the C pre-processor. If any files are found, they
# are renamed %.${new extension}. If the flag "-r" is passed in as argument 1,
# then mvfF.sh will recurse the entire directory's subdirectories.
# - Anthony Yau

old="f"
new="F"

if [ "x$1" = "x-r" ]
then
   test "x$2" = "x" || old="$2"
   test "x$3" = "x" || new="$3"
else
   test "x$1" = "x" || old="$1"
   test "x$2" = "x" || new="$2"
fi

# rename current files
files=$(egrep -l -e "(^\#|__FILE__|__LINE__)" *.${old} 2> /dev/null)
if [ "x$files" != "x" ]
then
   for i in $files
   do
      j=$(echo $i | sed -e "s/${old}$/${new}/")
      if [ -f $j ]
      then
         printf "ERROR: Cannot rename $PWD/$i since $PWD/$j already exists\n"
         exit 1
      fi
      echo Moving $i to $j
      mv $i $j
   done
fi

# enter all subdirectories
if [ "x$1" = "x-r" ]
then
   files=$(ls)
   if [ "x$files" != "x" ]
   then
      for i in *
      do
         if [ -d $i ]
         then
            cd $i
            printf "$(pwd):\n"
            $0 -r $old $new || exit 1
            cd ..
         fi
      done
   fi
fi

