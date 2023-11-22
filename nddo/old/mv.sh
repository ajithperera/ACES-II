test $# -eq 0 && exit
for sub in $*
do files=$(egrep -l -i "^ *subroutine *$sub *($|\()" unused/*.f*)
   c=0; for f in junk $files; do let c+=1; done
   test $c -eq 1 && echo $sub NOT FOUND && continue
   test $c -eq 2 && mv -i $files . && continue
   c=0; file=""
   for f in $files
   do sed -n "/^ *interface/,/^ *end *interface/d;
              /^ *subroutine *$sub *$/p;
              /^ *subroutine *$sub *[^0-9A-Za-z]/p;" $f | wc -l | read i
      test $i -eq 1 && file=$f
      let c+=$i
   done
   test $c -eq 1 && mv -i $file . || echo $sub MULTIPLY DEFINED
done
