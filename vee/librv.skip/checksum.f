      
       subroutine checksum(name,a,n)
C
        implicit double precision(a-h,o-z)
        character*(*) name
        dimension a(n)
C
        return
        sum=0.0d0
        sum1=0.0d0
        sum2=0.0d0
        do 1 i=1,n
        sum=sum+a(i)*a(i)
        sum1=sum1+abs(a(i))
        sum2=sum2+a(i)
1       continue
        write(6,300)name,sum,sum1,sum2
300     format(t3,a,1x,3(d20.10,1x))
        return
        end
