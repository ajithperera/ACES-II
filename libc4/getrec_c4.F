CLIC
C
C...   Copyright (c) 1992 - 2005 by the authors of Aces II (see below).
C...   All Rights Reserved.
C...
C...         The source code in this file is part of 
C...
C...    Aces II Mainz-Austin-Budapest-Version 2005, 
C...    a quantum chemical program package written by 
C...    J.F. Stanton, J. Gauss, J .D. Watts, P.G. Szalay, R.J. Bartlett
C...    with contribution from 
C...    A.A. Auer, D.B. Bernholdt, O. Christiansen, M.E. Harding, M. Heckert,
C...    O. Heun, C. Huber, D. Jonsson, J. Juselius, W.J. Lauderdale, 
C...    T. Metzroth, K. Ruud and the integral packages
C...    MOLECULE (J. Almlöf and P.R. Taylor),Props (P.R. Taylor), and 
C...    ABACUS (T. Helgaker, H.AA. Jensen, P.Jørgensen, and J. Olsen).
C...    See also J.F. Stanton, J. Gauss, J. D. Watts, W.J. Lauderdale, 
C...    R.J. Bartlett, Int. J. Quantum Chem. Symp. 26, 879 (1992)
C...
C...    as well as : http://www.aces2.de for the current version.
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of AcesII,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...                      license@aces2.de
C...
C...   For information on how to get a licence see:
C...                       http://www.aces2.de 
C
CLIC

      SUBROUTINE GETREC_C4(LUJNK,FILNAM,LABEL,LENGTH,ISTUFF)
      IMPLICIT INTEGER(A-Z)
      PARAMETER (LUFIL = 75)
      PARAMETER (IBUFLN = 128)
      CHARACTER*(*) LABEL
      CHARACTER*8 MARKER
      CHARACTER*(*) FILNAM
      DIMENSION ISTUFF(LENGTH),IBUF(IBUFLN)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /JOBARC/ MARKER(1000),LOC(1000),SIZE(1000),NRECS,
     &                IRECWD,IRECLN
      ILOC=LOCCHR_C4(1000,MARKER,1,LABEL)
      IF (ILOC.EQ.1001) THEN
         IF (LUJNK.GT.0) THEN
            WRITE(LUOUT,50)LUFIL,LABEL
 50         FORMAT(T3,'@GETREC-F, File attached to unit ',I3,
     $           ' does not contain record ',A8,'.')
            CALL ERREX
         else if ( lujnk .eq. 0 ) then
            length = -1
            return
         ELSE
            CALL IZERO(ISTUFF,MAX(1,LENGTH))
            RETURN
         ENDIF 
      ENDIF
      if ( lujnk .eq. 0 ) then
         length = size(iloc)
         return
      endif
      if ( length .gt. size(iloc) ) then
         write (6,*) '@GETREC: requested length of vector ', label,
     $        ' is longer than actual'
         write (6,*) 'Actual -- ', size(iloc), ' Requested -- ', length
      endif         
      IWRDST=LOC(ILOC)
      IRECST=1+IWRDST/IRECWD
      IOFFST=MOD(IWRDST,IRECWD)
      IF(IOFFST.EQ.0)THEN
       IRECST=IRECST-1
       IOFFST=IRECWD
      ENDIF
      NWRDST=IRECWD-IOFFST+1
      IWRDEN=LOC(ILOC)+LENGTH-1
      IRECEN=1+IWRDEN/IRECWD
      IOFFEN=MOD(IWRDEN,IRECWD)
      IF(IOFFEN.EQ.0)THEN
       IRECEN=IRECEN-1
       IOFFEN=IRECWD
      ENDIF
      NWRDEN=IOFFEN
      IF(IRECEN.EQ.IRECST)THEN
       READ(LUFIL,REC=IRECST,ERR=999,IOSTAT=IOS)IBUF
CSSS       CALL MYICOPY(IBUF(IOFFST),ISTUFF,LENGTH)
       CALL ICOPY(LENGTH,IBUF(IOFFST),1,ISTUFF,1)
       RETURN
      ELSE
       IOFF=1
       READ(LUFIL,REC=IRECST,ERR=999,IOSTAT=IOS)IBUF
CSSS       CALL MYICOPY(IBUF(IOFFST),ISTUFF(IOFF),NWRDST)
       CALL ICOPY(NWRDST,IBUF(IOFFST),1,ISTUFF(IOFF),1)
       IOFF=IOFF+NWRDST
       IF(IRECEN.EQ.IRECST)RETURN
       DO 10 I=IRECST+1,IRECEN-1
        READ(LUFIL,REC=I,ERR=999,IOSTAT=IOS)IBUF
CSSS        CALL MYICOPY(IBUF,ISTUFF(IOFF),IRECWD)
        CALL ICOPY(IRECWD,IBUF,1,ISTUFF(IOFF),1)
        IOFF=IOFF+IRECWD
10     CONTINUE
       READ(LUFIL,REC=IRECEN,ERR=999,IOSTAT=IOS)IBUF
CSSS       CALL ICOPY(IBUF,ISTUFF(IOFF),NWRDEN)
       CALL ICOPY(NWRDEN,IBUF,1,ISTUFF(IOFF),1)
       RETURN
      ENDIF
999   CONTINUE
      WRITE(LUOUT,51)LUFIL,LABEL,IOS
51    FORMAT(T3,'@GETREC_C4-F, I/O error on unit ',I3,' searching for ',
     &A10,'.',/,T3,' System IOSTAT error number is ',I3,'.')
      CALL ERREX
      END
