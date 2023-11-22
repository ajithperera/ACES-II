      Integer Function IrMxSymSz(ITypL, ITypR)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C $Id: irmxsymsz.f,v 1.1 2008/03/02 22:01:52 taube Exp $
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     IrMxSymSz -- Return the irrep of the largest ISymSz of a list
C
C SYNOPSIS
      Integer ITypL, ITypR
C
C ARGUMENTS
C     ITypL   Symmetry type of left hand indices (IN)
C     ITypR   Symmetry type of right hand indices (IN)
C
C DESCRIPTION
C     Returns the overall irrep for which a list composed of ITypL and
C     ITypR data would be largest.  Useful for determining the irrep
C     a list should be created with to insure that there is adequate
C     space to handle all possible irreps.
C
C ROUTINES REQUIRED
      External ISymSz2
      Integer ISymSz2
C
C COMMON BLOCKS
      Integer NStart, NIrrep, IrrepA, IrrepB, DIrPrd
      Common /SymInf/ NStart, NIrrep, IrrepA(255), IrrepB(255),
     $   DirPrd(8, 8)
      Integer IrpDPD, ISyTyp, NTot
      Common /SymPop/ IrpDPD(8, 22), ISyTyp(2, 500), NTot(18)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
C     IrrepT  Overall symmetry of list
C     Size    Size of current IrrepT
C     SizeMx  Current largest size
C
      Integer IrrepT, Size, SizeMx
C
      SizeMx = 0
      IrMxSymSz = 0
C
      Do 1000 IrrepT = 1, NIrrep
C
C        Determine the size for this overall symmetry
C
         Size = ISymSz2(IrrepT, ITypL, ITypR)
C
c         Write (6, 9000) IrrepT, Size
c 9000    Format(1X, '@IrMxSymSz-I, Irrep ', I1, ' has size ', I10, '.')
C
C        See how this size compares to the other symmetries
C
         If (Size .gt. SizeMx) then
            IrMxSymSz = IrrepT
            SizeMx = Size
         EndIf
 1000 Continue
C
      Return
      End
