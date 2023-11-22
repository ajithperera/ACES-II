#!/bin/ksh

# force the run
test "x$1" != "x-f" && printf "
Usage: $0 -f

" && exit 1

# limit new records to crisp
if test $(uname -n) != "crisp.qtp.ufl.edu"
then printf "\nNew test data should be appended on crisp.qtp.ufl.edu only!\n\n"
     exit 1
fi

alias g="egrep -ai"
MAXERR="1.E-9"

for recf in rec.NULL rec.*
do file=${recf#rec.}

 # cycle if there are no records to process
   test -s $recf || continue

 # attempt to remove the COORD record from tests that do not optimize geometries
   outf=out.$file
   if test -f $outf
   then i=$(g -hc 'JODA beginning optimization cycle' ${outf})
        if test $i -eq 0
        then echo "removing COORD record from $recf"
             sed '/d COORD/,/^d/{/d COORD/d;/^ *[0-9-]/d;}' $recf > tmp \
              && mv tmp $recf
        fi
   fi

 # replace any existing test data from the test file and remove the record file
   test -f $file || continue
   echo "adding $recf to $file"
   sed "
# remove trailing blanks from the ZMAT file
s/ *$//

# replace existing test data in the ZMAT file with the record file
/TEST\.DAT/,\${
               \${
                  s/.*/TEST\.DAT/
                  r $recf
                  b
                 }
               d
               b
              }

# add new test data
\${
   a\\
TEST.DAT
   r $recf
  }
" $file > tmp \
 && fail=0 || fail=1
if test $fail -eq 0
then if true
     then # record files already have error bars
          mv tmp $file
     else # record files need default error bars
          sed "/TEST\.DAT/,\${/^ *[0-9-]/ s/\$/ $MAXERR/;}" tmp > $file
     fi
     rm -f $recf tmp
fi

done

