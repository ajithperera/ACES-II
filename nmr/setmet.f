  
      
      SUBROUTINE SETMET(IUHF)
      IMPLICIT INTEGER(A-Z)
      CHARACTER*7 NAME
      CHARACTER*11 NAME1
      CHARACTER*4 NAME2
      CHARACTER*5 NAME3
      CHARACTER*6 NAME4
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL SEMI,ROHF,QRHF
      LOGICAL NOABCD,YESNO,GIAO,JSO
C
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/REF/ROHF,QRHF,SEMI
      COMMON/ABCD/NOABCD
      COMMON/GAUGE/GIAO
      COMMON/JNMR/JSO
C
      EQUIVALENCE(METHOD,IFLAGS(2))
C
3000  FORMAT(' Chemical shifts are calculated at the ',A7, 
     &       ' level.')       
3001  FORMAT(' The paramagnetic SO contribution to J is',
     &       ' calculated at the ',A7,' level.')       
c     &       ' for ',A7,' second derivative calculation.')
c3000  FORMAT(' Dipole derivatives are calculated at the ',A7, 
c     &       ' level.')       
c     &       ' for ',A7,' second derivative calculation.')
3005  FORMAT(' The reference state is a QRHF wave functions.')
3006  FORMAT(' The reference state is a non HF wave functions.')
3007  FORMAT(' The reference state is a ROHF wave function.')
3008  FORMAT(' Semi-canonical orbitals are used.')
C
      NOABCD=.TRUE.
      MBPT2=.FALSE.
      MBPT3=.FALSE.
      M4DQ=.FALSE.
      M4SDQ=.FALSE.
      M4SDTQ=.FALSE.
      CCD=.FALSE.
      QCISD=.FALSE.
      CCSD=.FALSE.
      ROHF=.FALSE.
      SEMI=.FALSE.
      GIAO=.TRUE. 
      JSO=.FALSE.
C
      IF(IFLAGS(18).EQ.5) GIAO=.FALSE.
      IF(IFLAGS(18).EQ.8) JSO=.TRUE.
C
      IF(METHOD.EQ.1) THEN
       MBPT2=.TRUE.
       NAME='MBPT(2)'
       IF(JSO) THEN
        write(6,3001)NAME
       ELSE
        write(6,3000)NAME
       ENDIF
      ELSE IF(METHOD.EQ.2) THEN
       MBPT3=.TRUE.
       NAME='MBPT(3)'
       IF(JSO) THEN
        write(6,3001)NAME
       ELSE
        write(6,3000)NAME
       ENDIF
       call errex
      ELSE IF(METHOD.EQ.3) THEN
       M4SDQ=.TRUE.
c       NAME1='SDQ-MBPT(4)'
c       write(6,3001) NAME1
       call errex   
      ELSE IF(METHOD.EQ.4) THEN
       M4SDTQ=.TRUE.
c       NAME='MBPT(4)'
c       write(6,3000)NAME
       call errex
      ELSE IF(METHOD.EQ.8) THEN
       CCD=.TRUE.
c       NAME2='CCD'
c       write(6,3002)NAME2
       call errex
      ELSE IF(METHOD.EQ.10) THEN
       CCSD=.TRUE.
c       NAME3='CCSD'
c       write(6,3003)NAME3
       call errex
      ELSE IF(METHOD.EQ.23) THEN
       QCISD=.TRUE.
c       NAME4='QCISD'
c       write(6,3004)NAME4
       call errex
      ENDIF
C
C CHECK IF WE ARE DEALING WITH A NON-HF CASE 
C
      IF(IFLAGS(32).NE.0.OR.IFLAGS(33).NE.0) THEN
C  
       QRHF=.TRUE.
C
       WRITE(6,3005)
C
      ENDIF
      IF(IFLAGS(11).EQ.2) THEN
C
       ROHF=.TRUE.
       WRITE(6,3007)
       IUHF=1
       IF(IFLAGS(39).EQ.1) THEN
        SEMI=.TRUE.
        WRITE(6,3008)
       ENDIF
      ENDIF
C
C CHECK IF THE MOABCD FILE EXISTS, IF THEN SET NOABCD TO .FALSE.
C
c YAU : old
c      INQUIRE(FILE='MOABCD',EXIST=YESNO)
c      IF(YESNO) NOABCD=.FALSE.
c YAU : new
c This approach will have to be redesigned for CC nmr calculations.
      call aces_io_remove(52,'MOABCD')
      NOABCD=.true.
c YAU : end
C
      RETURN
      END
