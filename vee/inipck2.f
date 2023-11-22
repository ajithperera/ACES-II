      SUBROUTINE INIPCK2(IRREPX,SYTYPL,SYTYPR,LIST,IARG1,IARGX2,ISET)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: inipck2.f,v 1.1.1.1 2003/04/02 19:22:05 aces Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     inipck2 -- Create a complete MOIO list
C
C SYNOPSIS
      Integer IrrepX, SyTypL, SyTypR, List, IArg1, IArgX2, ISet
C
C ARGUMENTS
C     IrrepX  Overall symmetry irrep of list (IN)
C     SyTypL  Symmetry type of LHS indices (IN)
C     SyTypR  Symmetry type of RHS indices (IN)
C     List    MOIO list number (IN)
C     IArg1   Control for subsidiary UPDMOI call. (IN)
C             = 0  Normal
C             = 1  If MOIO file containing List does not yet exist
C             = 2  Flush TOTRECMO to JOBARC
C             = 3  Load TOTRECMO from JOBARC
C     IArgX2  Control for subsidiary UPDMOI call. (IN)
C             = -1       Begin list at physical record boundary
C             otherwise  Use next available space.
C     ISet    Set ISYTYP entry for ths list (IN)
C             = 0        Don't set ISYTYP (not recommended)
C             otherwise  Set ISYTYP.
C
C DESCRIPTION
C     Initialize all irreps of LIST using a series of calls to UPDMOI.
C     Differs from INIPCK in that this routine guarantees the list
C     will be large enough to handle any possible IRREPX (via later
C     NEWTYP2 calls).
C
C ROUTINES REQUIRED
      External IrMxSymSz, INIPCK, NewTyp2
      Integer IrMxSymSz
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     IrrepY  Overall irrep of largest size
C
      Integer IrrepY
C
      IF(SYTYPL.EQ.0.OR.SYTYPR.EQ.0)RETURN
C     Irrep for which list will be largest
C
      IrrepY = IrMxSymSz(SyTypL, SyTypR)
C
C     Create list for largest irrep
C
      Call INIPCK(IrrepY, SyTypL, SyTypR, List, IArg1, IArgX2, ISet)
C
C     Now switch the list to the one they really wanted
C
      Call NewTyp2( IrrepX, List, SyTypL, SyTypR, .TRUE.)
C
c      Write (6, 9000) List, IrrepY, IrrepX
c 9000 Format(1X,'@INIPCK2-I Creating list ', I3, ' sized for irrep ',
c     $   I1, ' initialized for irrep ', I1, '.')
C
      Return
      End
