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

      SUBROUTINE GETLST_C4(ZLIST,NFIRST,NBLOCK,USECXN,ISPIN,ITYPE)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ZLIST(1)
      COMMON / / ICORE(1)
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &               MOIODS(10,500),MOIOFL(10,500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FILSPC/ ILNBUF,IPRCLN,IPRCWD
      COMMON /IOPOS/ ICRSIZ,ICHCSZ,IOFFX(2),LENREC
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /WHERE/ QUIKGET(10,500)
      IF(QUIKGET(ISPIN,ITYPE).NE.0)THEN
       CALL GETFAST_C4(ZLIST,NFIRST,NBLOCK,ISPIN,ITYPE)
       RETURN
      ENDIF
      IF(NFIRST+NBLOCK-1.GT.MOIODS(ISPIN,ITYPE))THEN
       WRITE(LUOUT,2000)ISPIN,ITYPE,NFIRST+NBLOCK-1,MOIODS(ISPIN,ITYPE)
2000   FORMAT(T3,'@GETLST-F, Error reading list [',I2,',',I3,']',/,
     &        T3,'           Distribution ',I9,' requested, but ',
     &           'only ',I9,' are available.')
       CALL ERREX
      ENDIF
      FILE=MOIOFL(ISPIN,ITYPE)
      NWORDS=NBLOCK*MOIOSZ(ISPIN,ITYPE)*IINTFP
      NRECST=MOIO(ISPIN,ITYPE)
      NWRDST=MOIOWD(ISPIN,ITYPE)
      NRECF=NRECST
      NOFFF=(NFIRST-1)*MOIOSZ(ISPIN,ITYPE)*IINTFP
      NWRDF=NWRDST+NOFFF
      IF(NWRDF.GT.IPRCWD)THEN
       INCREM=(NWRDF-1)/IPRCWD
       NWRDF=NWRDF-IPRCWD*INCREM
       NRECF=NRECST+INCREM
      ENDIF
      IREC=NRECF
      IOFF=1
1     NGET=MIN(NWORDS,IPRCWD-NWRDF+1)
      CALL FETCH_C4(ZLIST(IOFF),FILE,IREC,NWRDF,NGET)
      IREC=IREC+1
      NWORDS=NWORDS-NGET
      NWRDF=1
      IOFF=IOFF+NGET
      IF(NWORDS.NE.0)GOTO 1
      RETURN
      END
