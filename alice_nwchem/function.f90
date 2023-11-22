module mod_fun
contains

subroutine make_angindex(ind2,ind3,ind4,str2,bastype)
   implicit none
   integer i,ind2,ind3,ind4,bastype
   character(4) str2
   character(1) orbp(3)
   character(2) orbd(5)
   data (orbp(i),i=1,3) /'x','y','z'/
   data (orbd(i),i=1,5) /' 0','-2',' 1',' 2','-1'/

   if (str2(1:1)=='s') then
      ind3=1 
      ind4=1
   elseif (str2(1:1)=='p') then
      ind3=2 
      do i = 1, 3
         if (str2(2:2)==orbp(i)) ind4=i
      enddo
   elseif (str2(1:1)=='d') then
      ind3=3 
      if (bastype==1) then
!        cartesian
      elseif (bastype==2) then
!        spherical
         do i = 1, 5
            if (str2(3:4)==orbd(i)) ind4=i
         enddo
      endif
   elseif (str2(1:1)=='f') then
      ind3=4 
   endif

end subroutine make_angindex


subroutine rot_aoind(ntot,aoind,ind1,ind2)
   implicit none
   integer ntot,aoind(ntot,5),ind1,ind2
   integer k,aotmp(5)
   do k = 1, 5 
      aotmp(k)      = aoind(ind1,k) 
      aoind(ind1,k) = aoind(ind2,k)
      aoind(ind2,k) = aotmp(k)
   enddo
end subroutine rot_aoind


subroutine cart2sph(mat,irow,icol,str)
   implicit none
   character(len=*) str
   integer irow,icol,i,j,ind
   double precision mat(irow,icol) 
   double precision mat_d(36) 
!  sph(d0,d-2,d1,d2,d-1,xx) = mat_d * cart(xx,xy,xz,yy,yz,zz)
   data (mat_d(i),i = 1,36) &
         !  xx    xy   xz   yy   yz   zz
         / '-1', '0', '0','-1', '0', '2', &  !d0
           ' 0', '1', '0',' 0', '0', '0', &  !d-2
           ' 0', '0', '1',' 0', '0', '0', &  !d+1
           ' 1', '0', '0','-1', '0', '0', &  !d+2
           ' 0', '0', '0',' 0', '1', '0', &  !d-1
           ' 1', '0', '0',' 1', '0', '1' /   !xx

   ind = 0
   do i = 1, irow
      do j = 1, icol
         ind = ind+1
         if (str=="d") mat(i,j)=mat_d(ind) 
      enddo
   enddo
end subroutine cart2sph


subroutine sph2cart(mat,irow,icol,str)
   implicit none
   character(len=*) str
   integer irow,icol,i,j,ind
   double precision mat(irow,icol) 
   double precision mat_d(36) 
   data (mat_d(i),i = 1,36) &
         !  2      5    1    3      4    6
         !  d0     d-2  d+1  d+2    d-1  xx
         / '-1/6', '0', '0',' 1/2', '0', '0', &  !1
           '   0', '1', '0','   0', '0', '0', &  !2
           '   0', '0', '1','   0', '0', '0', &  !3
           '-1/6', '0', '0','-1/2', '0', '0', &  !4
           '   0', '0', '0','   0', '1', '0', &  !5
           ' 1/3', '0', '0','   0', '0', '0' /   !6

   ind = 0
   do i = 1, irow
      do j = 1, icol
         ind = ind+1
         if (str=="d") mat(i,j)=mat_d(ind)
      enddo
   enddo
end subroutine sph2cart


subroutine fun_ang_index1(val,strlen,str)
   implicit none
   integer val,strlen
   character(len=*) str
   integer i,j
   val = 0
   do i = 1, strlen
      if (str(i:i)=="S" .or. str(i:i)=="s") then
         val = 0
      elseif (str(i:i)=="X" .or. str(i:i)=="x") then
         val = val + 100 
      elseif (str(i:i)=="Y" .or. str(i:i)=="y") then
         val = val + 10 
      elseif (str(i:i)=="Z" .or. str(i:i)=="z") then
         val = val + 1 
      else
         read(str(2:4),*) val
         return
      endif
   enddo
end subroutine fun_ang_index1


subroutine fun_ang_index2(ntot,i,ao_ind,ao_ang)
   implicit none
   integer ntot,i
   integer ao_ind(ntot,3),ao_ang(ntot,3)
   ao_ang(i,3) = 1
   if (i==1) return
   if (ao_ind(i,1)/=ao_ind(i-1,1)) return
   if (ao_ang(i,1)/=ao_ang(i-1,1)) return
   if (ao_ang(i,2)/=ao_ang(i-1,2)) return
   ao_ang(i,3) = ao_ang(i-1,3)+1 
   return
end subroutine fun_ang_index2


subroutine fun_exchint(num1,num2)
   implicit none
   integer num1,num2,itemp
   itemp = num1
   num1 = num2
   num2 = itemp
end subroutine fun_exchint

subroutine fun_exchstr(str1,str2)
   implicit none
   character(3) str1,str2,stemp
   stemp = str1
   str1 = str2
   str2 = stemp
end subroutine fun_exchstr


FUNCTION Upper(s1)  RESULT (s2)
   CHARACTER(*)       :: s1
   CHARACTER(LEN(s1)) :: s2
   CHARACTER          :: ch
   INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
   INTEGER            :: i
  
   DO i = 1,LEN(s1)
      ch = s1(i:i)
      IF (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
      s2(i:i) = ch
   END DO
END FUNCTION Upper


FUNCTION Lower(s1)  RESULT (s2)
   CHARACTER(*)       :: s1
   CHARACTER(LEN(s1)) :: s2
   CHARACTER          :: ch
   INTEGER,PARAMETER  :: DUC = ICHAR('A') - ICHAR('a')
   INTEGER            :: i
  
   DO i = 1,LEN(s1)
      ch = s1(i:i)
      IF (ch >= 'A'.AND.ch <= 'Z') ch = CHAR(ICHAR(ch)-DUC)
      s2(i:i) = ch
   END DO
END FUNCTION Lower

end module mod_fun
