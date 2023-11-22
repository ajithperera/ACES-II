CIC
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

      SUBROUTINE PUTREC_C4(LUJNK,FILNAM,LABEL,LENGTH,ISTUFF)
      IMPLICIT INTEGER(A-Z)
      PARAMETER (IBUFLN = 128)
      PARAMETER (LUFIL = 75)
      CHARACTER*(*) LABEL
      CHARACTER*8 MARKER
      CHARACTER*(*) FILNAM
      DIMENSION ISTUFF(LENGTH),IBUF(IBUFLN)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /JOBARC/ MARKER(1000),LOC(1000),SIZE(1000),NRECS,
     &                IRECWD,IRECLN
      IF(LENGTH.EQ.0)RETURN
      ILOC=LOCCHR_C4(1000,MARKER,1,LABEL)
      IF(ILOC.EQ.1001)THEN
       ILOC=LOCCHR_C4(1000,MARKER,1,'OPENSLOT')
       IF(ILOC.GE.1001)THEN
        WRITE(LUOUT,1000)
1000    FORMAT(T3,'@PUTREC_C4-F, Maximum allowable number of records ',
     &            'exceeded in JOBARC.',/,T3,'@PUTREC-F, Redimension ',
     &            'arrays in common block JOBARC.')
       ENDIF
       MARKER(ILOC)=LABEL
       SIZE(ILOC)=LENGTH
       IF(ILOC.NE.1)THEN
        LOC(ILOC)=LOC(ILOC-1)+SIZE(ILOC-1)
       ELSE
        LOC(ILOC)=1
       ENDIF
      ELSE
       IF(LENGTH.GT.SIZE(ILOC))THEN
        WRITE(LUOUT,1001) MARKER(iloc), LENGTH,SIZE(ILOC)
1001    FORMAT(T3,'@PUTREC-F, Size of record ',A8,'[',I8,'words ] ',
     &           'exceeds original length of ',I8,'.')
        CALL ERREX
       ENDIF
      ENDIF
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
       IF(IRECST.LE.NRECS)READ(LUFIL,REC=IRECST,ERR=999,
     &                         IOSTAT=IOS)IBUF
CSSS       CALL MYICOPY(ISTUFF,IBUF(IOFFST),LENGTH)
       CALL ICOPY(LENGTH,ISTUFF,1,IBUF(IOFFST),1)
       WRITE(LUFIL,REC=IRECST,ERR=999,IOSTAT=IOS)IBUF
       NRECS=MAX(NRECS,IRECST)
       RETURN
      ELSE
       IOFF=1
       IF(IRECST.LE.NRECS)READ(LUFIL,REC=IRECST,ERR=999,
     &                         IOSTAT=IOS)IBUF
CSSS       CALL MYICOPY(ISTUFF(IOFF),IBUF(IOFFST),NWRDST)
       CALL ICOPY(NWRDST,ISTUFF(IOFF),1,IBUF(IOFFST),1)
       WRITE(LUFIL,REC=IRECST,ERR=999,IOSTAT=IOS)IBUF
       NRECS=MAX(NRECS,IRECST)
       IOFF=IOFF+NWRDST
       DO 10 I=IRECST+1,IRECEN-1
        IF(I.LE.NRECS)READ(LUFIL,REC=I,ERR=999,
     &                     IOSTAT=IOS)IBUF
CSSS        CALL MYICOPY(ISTUFF(IOFF),IBUF,IRECWD)
        CALL ICOPY(IRECWD,ISTUFF(IOFF),1,IBUF,1)
        WRITE(LUFIL,REC=I,ERR=999,IOSTAT=IOS)IBUF
        NRECS=MAX(NRECS,I)
        IOFF=IOFF+IRECWD
10     CONTINUE
       IF(IRECEN.LE.NRECS)READ(LUFIL,REC=IRECEN,ERR=999,
     &                         IOSTAT=IOS)IBUF
CSSS       CALL MYICOPY(ISTUFF(IOFF),IBUF,NWRDEN)
       CALL ICOPY(NWRDEN,ISTUFF(IOFF),1,IBUF,1)
       WRITE(LUFIL,REC=IRECEN,ERR=999,IOSTAT=IOS)IBUF
       NRECS=MAX(NRECS,IRECEN)
       RETURN
      ENDIF
999   CONTINUE
      WRITE(LUOUT,51)LUFIL,LABEL,IOS
      CALL ERREX
51    FORMAT(T3,'@PUTREC_C4-F, I/O error on unit ',I3,' while trying ',
     &'to write ',A10,'.',/,T3,' System IOSTAT error number is ',
     &I3,'.')
      END
