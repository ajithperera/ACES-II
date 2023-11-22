      SUBROUTINE ATOMIC_GUESS(DENS,LDIM1,IUHF)
      IMPLICIT NONE
      INTEGER LDIM1,IUHF,II,JJ,KK,starts,ldimtmp,nbastmp,
     &  lili,lilj,bigi,bigj,indx,maxlen
      DOUBLE PRECISION DENS((IUHF+1)*LDIM1)
      CHARACTER BAS*30,BASTMP*30
#include "machsp.com"

C this routine only works if symmetry is off, so need to check this
     
      
C read the ZMAT.BAS file to get each atom and basis function
C then look for them in ATOMDENS and put them in the density matrix

      maxlen=1000000
      bigi=0
      bigj=0
      open(66,file='ZMAT.BAS')
      DO JJ=1,maxlen
        read(66,*,err=100)BAS
        starts=INDEX(BAS,":")
        If ((starts.ne.0)) Then
          open(67,file='ATOMDENS')
          rewind 67
          Do II=1,maxlen
            If (II.eq.maxlen) GOTO 666
            read(67,*,err=101)BASTMP
            If (BAS(1:5).eq.BASTMP(1:5)) Then
              read(67,*)nbastmp
              read(67,*)ldimtmp
              Do KK=1,ldimtmp
         lili=1+(-1+aint(dsqrt(8.D0*KK+0.999D0)))/2
         lilj=KK-(lili*(lili-1))/2
         indx=bigj+lilj+((bigi+lili)*((bigi+lili)-1))/2
         read(67,111)DENS(indx)
              ENDDO
              bigi=bigi+nbastmp
              bigj=bigj+nbastmp
              GOTO 100
            EndIf
 101   Continue
          EndDo
          close(67)
        EndIf
 100  Continue
      ENDDO
      close(66)

      If (IUHF.eq.1) Then
        Do II=1,LDIM1
          DENS(LDIM1+II)=DENS(II)
        EndDo
      Else
        Do II=1,LDIM1
          DENS(II)=DENS(II)*2.D0
        EndDo
      EndIf

      RETURN

C report an error if the basis set to found in ATOMDENS
 666  Continue
         WRITE(6,'(a,/)')
     &            "@ATOMIC_DENSITY: "//BAS(1:10)//"is not in ATOMDENS."
         CALL ERREX
         STOP
 111     FORMAT(50(F20.15))
      END
