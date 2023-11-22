      
       subroutine checksm2(name,a,n1,n2,sum,sum1,sum2)
C
        implicit double precision(a-h,o-z)
        character*8 name
        dimension a(n1*n2)
c        return
C   
        ind=0
        do 2 ij=1,n1
        sum=0.0
        sum1=0.0
        sum2=0.0
        do 1 i=1,n2
        ind=ind+1
        sum=sum+a(ind)*a(ind)
        sum1=sum1+abs(a(ind))
        sum2=sum2+a(ind)
1       continue
        write(*,*) 'checksum',name,'ij ',sum,sum1,sum2
2       continue
        return
        end
