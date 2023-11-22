      subroutine ehttran(c2z,z2c,iz,nbasxeht,natoms)
      implicit none
      double precision c2z,z2c
      integer iz,nbasxeht,natoms
c
      integer iatom,len,ioff,i
c
      dimension c2z(nbasxeht,nbasxeht),z2c(nbasxeht,nbasxeht)
      dimension iz(natoms)
c
      dimension len(50),ioff(50)
c
c-----------------------------------------------------------------------
c     Temporary routine for C1 minimal basis set transformation matrices
c-----------------------------------------------------------------------
c
      call zero(c2z,nbasxeht*nbasxeht)
      call zero(z2c,nbasxeht*nbasxeht)
c
      do 10 iatom=1,natoms
      if(iz(iatom).le. 2)                       len(iatom) = 1
      if(iz(iatom).ge. 3 .and. iz(iatom).le.10) len(iatom) = 5
      if(iz(iatom).ge.11 .and. iz(iatom).le.18) len(iatom) = 9
   10 continue
c
      do 20 iatom=1,natoms
      if(iatom.eq.1)then
       ioff(iatom) = 0
      else
       ioff(iatom) = ioff(iatom-1) + len(iatom-1)
      endif
   20 continue
c
      do  40 iatom=1,natoms
c
      if(len(iatom).eq.1 .or. len(iatom).eq. 5)then
c
       do  30 i=1,len(iatom)
c
       c2z(ioff(iatom)+i,ioff(iatom)+i) = 1.0D+00
       z2c(ioff(iatom)+i,ioff(iatom)+i) = 1.0D+00
   30  continue
      else
c
       c2z(ioff(iatom)+1,ioff(iatom)+1) = 1.0D+00
       c2z(ioff(iatom)+2,ioff(iatom)+2) = 1.0D+00
       c2z(ioff(iatom)+3,ioff(iatom)+3) = 1.0D+00
       c2z(ioff(iatom)+4,ioff(iatom)+4) = 1.0D+00
       c2z(ioff(iatom)+7,ioff(iatom)+5) = 1.0D+00
       c2z(ioff(iatom)+5,ioff(iatom)+6) = 1.0D+00
       c2z(ioff(iatom)+8,ioff(iatom)+7) = 1.0D+00
       c2z(ioff(iatom)+6,ioff(iatom)+8) = 1.0D+00
       c2z(ioff(iatom)+9,ioff(iatom)+9) = 1.0D+00
c
       z2c(ioff(iatom)+1,ioff(iatom)+1) = 1.0D+00
       z2c(ioff(iatom)+2,ioff(iatom)+2) = 1.0D+00
       z2c(ioff(iatom)+3,ioff(iatom)+3) = 1.0D+00
       z2c(ioff(iatom)+4,ioff(iatom)+4) = 1.0D+00
       z2c(ioff(iatom)+5,ioff(iatom)+7) = 1.0D+00
       z2c(ioff(iatom)+6,ioff(iatom)+5) = 1.0D+00
       z2c(ioff(iatom)+7,ioff(iatom)+8) = 1.0D+00
       z2c(ioff(iatom)+8,ioff(iatom)+6) = 1.0D+00
       z2c(ioff(iatom)+9,ioff(iatom)+9) = 1.0D+00
c
      endif
c
   40 continue
      return
      end
