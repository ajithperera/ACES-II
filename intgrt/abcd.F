         subroutine abcd(itotwt,iradpt,iangpt,
     &   n,m,totwt,natoms,icntr)

       implicit none

      double precision itotwt,totwt
      integer  n,m,iradpt,iangpt,natoms,icntr

      dimension itotwt(natoms,n,m)

      itotwt(icntr,iradpt,iangpt)=totwt

      end

