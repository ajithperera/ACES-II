
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
      SUBROUTINE DAXPBY (N,DZ,DA,DX,DB,DY)
C
C     Purpose:
C     This subroutine computes DZ = DA * DX + DB * DY.  Several special
C     cases are handled separately:
C        DA =  0.0, DB =  0.0 => DZ = 0.0
C        DA =  0.0, DB =  1.0 => DZ = DY  (this is COPY)
C        DA =  0.0, DB = -1.0 => DZ = -DY
C        DA =  0.0, DB =   DB => DZ = DB * DY  (this is SCAL)
C        DA =  1.0, DB =  0.0 => DZ = DX  (this is COPY)
C        DA =  1.0, DB =  1.0 => DZ = DX + DY
C        DA =  1.0, DB = -1.0 => DZ = DX - DY
C        DA =  1.0, DB =   DB => DZ = DX + DB * DY (this is AXPY)
C        DA = -1.0, DB =  0.0 => DZ = -DX
C        DA = -1.0, DB =  1.0 => DZ = -DX + DY
C        DA = -1.0, DB = -1.0 => DZ = -DX - DY
C        DA = -1.0, DB =   DB => DZ = -DX + DB * DY
C        DA =   DA, DB =  0.0 => DZ = DA * DX  (this is SCAL)
C        DA =   DA, DB =  1.0 => DZ = DA * DX + DY  (this is AXPY)
C        DA =   DA, DB = -1.0 => DZ = DA * DX - DY
C        DA =   DA, DB =   DB => DZ = DA * DX + DB * DY
C     DZ may be the same as DX or DY.
C
C     Parameters:
C     N  = the dimension of the vectors (input).
C     DZ = the vector result (output).
C     DA = scalar multiplier for DX (input).
C     DX = one of the vectors (input).
C     DB = scalar multiplier for DY (input).
C     DY = the other vector (input).
C
C     Noel M. Nachtigal
C     March 23, 1993
C
C**********************************************************************
C
      INTEGER N
      DOUBLE PRECISION DA, DB, DX(N), DY(N), DZ(N)
C
C     Local variables.
C
      INTEGER I
C
      IF (N.LE.0) RETURN
C
      IF (DA.EQ.0.0D0) THEN
         IF (DB.EQ.0.0D0) THEN
C           DA = 0.0, DB = 0.0 => DZ = 0.0.
            DO 10 I = 1, N
               DZ(I) = 0.0D0
 10         CONTINUE
         ELSE IF (DB.EQ.1.0D0) THEN
C           DA = 0.0, DB = 1.0 => DZ = DY (this is COPY).
            DO 20 I = 1, N
               DZ(I) = DY(I)
 20         CONTINUE
         ELSE IF (DB.EQ.-1.0D0) THEN
C           DA = 0.0, DB = -1.0 => DZ = -DY.
            DO 30 I = 1, N
               DZ(I) = -DY(I)
 30         CONTINUE
         ELSE
C           DA = 0.0, DB = DB => DZ = DB * DY (this is SCAL).
            DO 40 I = 1, N
               DZ(I) = DB * DY(I)
 40         CONTINUE
         END IF
      ELSE IF (DA.EQ.1.0D0) THEN
         IF (DB.EQ.0.0D0) THEN
C           DA = 1.0, DB = 0.0 => DZ = DX (this is COPY).
            DO 50 I = 1, N
               DZ(I) = DX(I)
 50         CONTINUE
         ELSE IF (DB.EQ.1.0D0) THEN
C           DA = 1.0, DB = 1.0 => DZ = DX + DY.
            DO 60 I = 1, N
               DZ(I) = DX(I) + DY(I)
 60         CONTINUE
         ELSE IF (DB.EQ.-1.0D0) THEN
C           DA = 1.0, DB = -1.0 => DZ = DX - DY.
            DO 70 I = 1, N
               DZ(I) = DX(I) - DY(I)
 70         CONTINUE
         ELSE
C           DA = 1.0, DB = DB => DZ = DX + DB * DY (this is AXPY).
            DO 80 I = 1, N
               DZ(I) = DX(I) + DB * DY(I)
 80         CONTINUE
         END IF
      ELSE IF (DA.EQ.-1.0D0) THEN
         IF (DB.EQ.0.0D0) THEN
C           DA = -1.0, DB = 0.0 => DZ = -DX
            DO 90 I = 1, N
               DZ(I) = -DX(I)
 90         CONTINUE
         ELSE IF (DB.EQ.1.0D0) THEN
C           DA = -1.0, DB = 1.0 => DZ = -DX + DY
            DO 100 I = 1, N
               DZ(I) = -DX(I) + DY(I)
 100        CONTINUE
         ELSE IF (DB.EQ.-1.0D0) THEN
C           DA = -1.0, DB = -1.0 => DZ = -DX - DY.
            DO 110 I = 1, N
               DZ(I) = -DX(I) - DY(I)
 110        CONTINUE
         ELSE
C           DA = -1.0, DB = DB => DZ = -DX + DB * DY
            DO 120 I = 1, N
               DZ(I) = -DX(I) + DB * DY(I)
 120        CONTINUE
         END IF
      ELSE
         IF (DB.EQ.0.0D0) THEN
C           DA = DA, DB = 0.0 => DZ = DA * DX (this is SCAL).
            DO 130 I = 1, N
               DZ(I) = DA * DX(I)
 130        CONTINUE
         ELSE IF (DB.EQ.1.0D0) THEN
C           DA = DA, DB = 1.0 => DZ = DA * DX + DY (this is AXPY)
            DO 140 I = 1, N
               DZ(I) = DA * DX(I) + DY(I)
 140        CONTINUE
         ELSE IF (DB.EQ.-1.0D0) THEN
C           DA = DA, DB = -1.0 => DZ = DA * DX - DY.
            DO 150 I = 1, N
               DZ(I) = DA * DX(I) - DY(I)
 150        CONTINUE
         ELSE
C           DA = DA, DB = DB => DZ = DA * DX + DB * DY.
            DO 160 I = 1, N
               DZ(I) = DA * DX(I) + DB * DY(I)
 160        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C
C**********************************************************************
C**********************************************************************
C
      SUBROUTINE DHPSRT (N,X,Y,Z)
C
C     Purpose:
C     This subroutine sorts the  array X using HeapSort.  It rearranges
C     the elements of Y and Z at the same time.
C
C     Parameters:
C     N = the length of the arrays (input).
C     X = the primary array to be used in sorting (input/output).
C     Y = another array, sorted in the same order as X (input/output).
C     Z = another array, sorted in the same order as X (input/output).
C
C     Noel M. Nachtigal
C     October 4, 1990
C
C**********************************************************************
C
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), Z(N)
C
C     Local variables.
C
      INTEGER I, J, K, L
      DOUBLE PRECISION TMPX, TMPY, TMPZ
C
      IF (N.LE.1) RETURN
C
      L = N / 2 + 1
      K = N
 10   IF (L.GT.1) THEN
         L = L - 1
         TMPX = X(L)
         TMPY = Y(L)
         TMPZ = Z(L)
      ELSE
         TMPX = X(K)
         TMPY = Y(K)
         TMPZ = Z(K)
         X(K) = X(1)
         Y(K) = Y(1)
         Z(K) = Z(1)
         K = K - 1
         IF (K.LE.1) THEN
            X(1) = TMPX
            Y(1) = TMPY
            Z(1) = TMPZ
            RETURN
         END IF
      END IF
      I = L
      J = L + L
 20   IF (J.LE.K) THEN
         IF (J.LT.K) THEN
            IF (X(J).GT.X(J+1)) J = J + 1
         END IF
         IF (TMPX.GT.X(J)) THEN
            X(I) = X(J)
            Y(I) = Y(J)
            Z(I) = Z(J)
            I    = J
            J    = J + J
         ELSE
            J = K + 1
         END IF
         GO TO 20
      END IF
      X(I) = TMPX
      Y(I) = TMPY
      Z(I) = TMPZ
      GO TO 10
      END
C
C**********************************************************************
C**********************************************************************
C
      SUBROUTINE DRANDN (N,DX,SEED)
C
C     Purpose:
C     Fills the vector DX with random numbers  between 0 and 1.  If the
C     SEED is given, it should be odd and positive.  The generator is a
C     fairly unsophisticated one, from Pearson's  "Numerical methods in
C     engineering and science" book.
C
C     Parameters:
C     N    = the dimension of the vector (input).
C     DX   = the vector to fill with random numbers (output).
C     SEED = the seed for the generator (input).
C
C     Noel M. Nachtigal
C     April 23, 1993
C
C**********************************************************************
C
      INTRINSIC DBLE, IABS, MOD
C
      INTEGER N, SEED
      DOUBLE PRECISION DX(N)
C
C     Local variables.
C
      INTEGER I, J
C
C     Local variables that are saved from one call to the next.
C
      DOUBLE PRECISION DMAX
      INTEGER IM, IMAX, IS
      SAVE DMAX, IM, IMAX, IS
      DATA IM/0/
C
C     Initialize the generator data.
C
      IF (IM.EQ.0) THEN
         J  = 0
         IM = 1
         DO 10 I = 1, 31
            J = J + 1
            IF (IM*2.LE.IM) GO TO 20
            IM = IM * 2
 10      CONTINUE
 20      IMAX = (IM-1) * 2 + 1
         DMAX = DBLE(IMAX)
         DO 30 I = 1, MOD(J,3)
            J = J - 1
            IM = IM / 2
 30      CONTINUE
         IM = IM + 5
         IS = IABS(MOD(IM*30107,IMAX))
      END IF
C
C     Check whether we have a new seed.
C
      IF (SEED.GT.0) IS = (SEED / 2) * 2 + 1
C
C     Here goes the rest.
C
      DO 40 I = 1, N
         DX(I) = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
 40   CONTINUE
C
      RETURN
      END
C
C**********************************************************************
C**********************************************************************
C
C     Copyright (C) 1991-1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file "cpyrit.doc" in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
C     This file contains the routines for  the eigenvalue solver, using
C     the  three-term  recurrence  variant  of  the look-ahead  Lanczos
C     algorithm.
C
C**********************************************************************
C
      SUBROUTINE DULAL (NDIM,NLEN,NLIM,MAXN,MAXVW,M,NORM,DWK,IDX,IWK,
     $                  HWK,VECS,INFO)
C
C     Purpose:
C     This  subroutine uses  the Lanczos  algorithm  with look-ahead to
C     compute eigenvalue  estimates.  The algorithm was first described
C     in  the  RIACS Technical  Report 90.45, `An Implementation of the
C     Look-Ahead Lanczos Algorithm for Non-Hermitian Matrices, Part I`,
C     by R.W. Freund, M.H. Gutknecht and N.M. Nachtigal, November 1990.
C     The routine will first run the look-ahead Lanczos algorithm for a
C     number of steps.  It will then prompt the user  for the number of
C     eigenvalues to compute.  In addition, the routine will  also find
C     the common eigenvalues in two successive runs.
C
C     Parameters:
C     For a description of  the parameters, see the file "dulal.doc" in
C     the current directory.
C
C     External routines used:
C     subroutine daxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     double precision ddot(n,dx,incx,dy,incy)
C        BLAS-1 routine, computes y^H * x.
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C     double precision dnrm2(n,dx,incx)
C        BLAS-1 routine, computes the 2-norm of x.
C     subroutine drandn(n,x,seed)
C        Library routine, fills x with random numbers.
C     subroutine dulal1(ndim,nlen,m,maxvw,numchk,n,l,lstar,nl,nlstar,vf,
C        ierr,adjust,norm1,norma,tmax,tnrm,dwk,iwk,idx,vecs)
C        Low-level routine, builds the vectors v_{n+1} and w_{n+1}.
C     subroutine dulalc(maxn,neig,nold,hwk,tf,vf)
C        Library routine, computes common eigenvalues.
C     subroutine dulale(maxn,neig,hwk,tf,vf)
C        Library routine, computes Lanczos eigenvalue approximations.
C
C     Noel M. Nachtigal
C     October 24, 1990
C
C**********************************************************************
C
      INTRINSIC DABS, DSQRT, MAX0, MOD
      EXTERNAL DAXPBY, DDOT, DLAMCH, DNRM2, DULAL1, DULALC, DULALE
      DOUBLE PRECISION DDOT, DLAMCH, DNRM2
C
      INTEGER INFO(4), M, MAXN, MAXVW, NDIM, NLEN, NLIM
      INTEGER IDX(3,NLIM+2), IWK(M,4)
      DOUBLE PRECISION DWK(M,5*M+7), HWK(MAXN,2*MAXN+8), NORM
      DOUBLE PRECISION VECS(NDIM,2*MAXVW+4+2)
C
C     Common block variables.
C
C
C     Common block DULALX.
C
      DOUBLE PRECISION NORMA
      COMMON /DULALX/NORMA
C
C     Miscellaneous parameters.
C
      DOUBLE PRECISION DONE, DTEN, DZERO
      PARAMETER (DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, L, LSTAR, MVWBLT, N, NL, NLSTAR, NMAX, NUMCHK
      SAVE    IERR, L, LSTAR, MVWBLT, N, NL, NLSTAR, NMAX, NUMCHK
      INTEGER RETLBL, TF, VF
      SAVE    RETLBL, TF, VF
      DOUBLE PRECISION ADJUST, NORM1, TMAX, TMIN, TNRM
      SAVE             ADJUST, NORM1, TMAX, TMIN, TNRM
C
C     Local variables, transient.
C
      INTEGER I, J, K, NEIG, REVCOM
      DOUBLE PRECISION DTMP, DTMP1, DTMP2
C
C     Local variables.
C
      INTEGER NLAL, NOLD
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C        REVCOM   RETLBL      Comment
C           0        0    first call, go to label 10
C           1       60    returning from AXB, go to label 60
C           2       70    returning from ATXB, go to label 70
C

#ifdef _DEBUG_LVL0
        Write(6,*) 
        Write(6,"(a)") "---Entered dulal ---"
        Write(6,*)
#endif

      Write(6,*) 
      REVCOM  = INFO(2)
      INFO(2) = 0

#ifdef _DEBUG_LVL0
      Write(6,"(a,2i4)") " @-dulal revcom", revcom,RETLBL
      Write(6,*) 
#endif 

      IF (REVCOM.EQ.0) THEN
         N      = 0
         NLAL   = 0
         MVWBLT = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.60) THEN
            GO TO 60
         END IF
      ELSE IF (REVCOM.EQ.2) THEN
         IF (RETLBL.EQ.70) GO TO 70
      END IF
      IERR = 1
      GO TO 110
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)      IERR = 2
      IF (NLEN.LT.1)      IERR = 2
      IF (NLIM.LT.1)      IERR = 2
      IF (MAXVW.LT.1)     IERR = 2
      IF (NLEN.GT.NDIM)   IERR = 2
      IF (NLIM.GT.MAXN-2) IERR = 2
      IF (M.LT.2*MAXVW+2) IERR = 2
      Write(6,*) "@-Dulal Ierr", Ierr
      IF (IERR.NE.0) GO TO 110

C
C     Extract from INFO the output units TF and VF,  and the starting
C     vector flags J and K.
C
      VF    = MAX0(INFO(1),0)
      J     = VF / 100000
      VF    = VF - J * 100000
      K     = VF / 10000
      VF    = VF - K * 10000
      TF    = VF / 100
      VF    = VF - TF * 100

      Write(6,*) J,K
C
C     Extract the norms.
C
      NORMA = DZERO
      NORM1 = DABS(NORM)
C
C     Set the adjustment parameters.
C
      NUMCHK = 25
      ADJUST = DTEN
C
C     Extract and check the various tolerances and norm estimates.
C
      TNRM = DLAMCH('E') * DTEN
      TMIN = DSQRT(DSQRT(DLAMCH('S')))
      TMAX = DONE / TMIN
C
C     Set up wrapped indices.  The following indices are used:
C        IDX(1,I) = indices used for IWK and DWK, row-dimensioned M;
C        IDX(2,I) = indices used for v_i;
C        IDX(3,I) = indices used for w_i.
C
      DO 20 I = 1, NLIM+2
         IDX(1,I) = MOD(I-1,M) + 1
         IDX(2,I) = 2 + 1 + 0*(MAXVW+2) + MOD(I-1,MAXVW+2)
         IDX(3,I) = 2 + 1 + 1*(MAXVW+2) + MOD(I-1,MAXVW+2)
 20   CONTINUE
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal;Nlim+2", Nlim+2
      Write(6,"(12(1x,i3))") ((idx(i,j),j=1,nlim+2),i=1,3)
      Write(6,*)
#endif 
C
C     Check whether any starting vectors must be supplied.
C
      IF (K.EQ.0) CALL DRANDN (NLEN,VECS(1,1),1)
      CALL DAXPBY (NLEN,VECS(1,IDX(2,1)),DONE,VECS(1,1),DZERO,VECS(1,
     $IDX(2,1)))
      IF (J.EQ.0) CALL DRANDN (NLEN,VECS(1,2),1)
      CALL DAXPBY (NLEN,VECS(1,IDX(3,1)),DONE,VECS(1,2),DZERO,VECS(1,
     $IDX(3,1)))
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal;the starting R and L vecs"
      call checksum("@-Dulal-recieve", VECS(1,IDX(2,1)),NLEN,s)
      call checksum("@-Dulal-recieve", VECS(1,IDX(3,1)),NLEN,s)
#endif 
C
C     Zero out the Lanczos matrix.
C
      DO 40 I = 1, MAXN
         DO 30 J = 1, MAXN
            HWK(I,J) = DZERO
 30      CONTINUE
 40   CONTINUE
C
C     Scale the first pair of Lanczos vectors and check for invariant
C     subspaces.
C
      DTMP1 = DNRM2(NLEN,VECS(1,IDX(2,1)),1)
      DTMP2 = DNRM2(NLEN,VECS(1,IDX(3,1)),1)
#ifdef _DEBUG_LVL0
      Write(6,*) "@-dulal; the norms of input", DTMP1,DTMP2
#endif 
      IF (DTMP1.LT.TNRM) IERR = IERR + 16
      IF (DTMP2.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 110
      DWK(IDX(1,1),5*M+4)   = DONE
      DWK(IDX(1,1),IDX(1,1)) = DDOT(NLEN,VECS(1,IDX(2,1)),1,VECS(1,IDX(3
     $,1)),1) / ( DTMP1 * DTMP2 )
      IF ((DTMP1.GE.TMAX).OR.(DTMP1.LE.TMIN)) THEN
         DTMP = DONE / DTMP1
         CALL DAXPBY (NLEN,VECS(1,IDX(2,1)),DTMP,VECS(1,IDX(2,1)),DZERO,
     $VECS(1,IDX(2,1)))
         DTMP1 = DONE
      END IF
      IF ((DTMP2.GE.TMAX).OR.(DTMP2.LE.TMIN)) THEN
         DTMP = DONE / DTMP2
         CALL DAXPBY (NLEN,VECS(1,IDX(3,1)),DTMP,VECS(1,IDX(3,1)),DZERO,
     $VECS(1,IDX(3,1)))
         DTMP2 = DONE
      END IF
      DWK(IDX(1,1),5*M+5) = DONE / DTMP1
      DWK(IDX(1,1),5*M+6) = DONE / DTMP2
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal;the starting R and L vecs after processing"
      call checksum("@-Dulal-recieve_1", VECS(1,IDX(2,1)),NLEN,s)
      call checksum("@-Dulal-recieve_1", VECS(1,IDX(3,1)),NLEN,s)
      Write(6,*) 
#endif 
C
C     Initialize the counters.
C
      L      = 1
      N      = 1
      NLAL   = 0
      NMAX   = 0
      IWK(IDX(1,0+1),1) = 1
      IWK(IDX(1,1+1),1) = 1
      MVWBLT = 1
      NL     = IWK(IDX(1,L+1),1)
C
C     This is one step of the coupled two-term Lanczos algorithm.
C     V_n, W_n, D_{n-1}, F_{n-1}, H_{n-1}, and D_{nn} are given.
C
C     Note that the previous step was  not necessarily step N-1,  as it
C     could have been a restart.
C     Except at the first step, the following hold:
C       NL.LE.N, N-NL.LE.MAXVW-1.
C
      IF (VF.NE.0) WRITE (VF,'(A29)') 'Running look-ahead Lanczos...'
 50   IWK(IDX(1,N),4) = L
      IWK(IDX(1,N),3) = N

C     Build v_{n+1} and w_{n+1}.
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,IDX(2,N)),VECS(1,1))
C
      INFO(2) = 1
      INFO(3) = IDX(2,N)
      INFO(4) = 1
      RETLBL  = 60
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal Sending R,IDX(2,N)", IDX(2,N)
      call checksum("@-Dulal-send", VECS(1,IDX(2,N)),NLEN,s)
#endif 
      RETURN
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal Recieving R"
      call checksum("@-Dulal-recieve", VECS(1,1),NLEN,s)
      call checksum("@-Dulal-recieve", VECS(1,IDX(2,N)),NLEN,s)
#endif 
C
C     Have the caller carry out ATXB, then return here.
C        CALL ATXB (VECS(1,IDX(3,N)),VECS(1,2))
C
 60   INFO(2) = 2
      INFO(3) = IDX(3,N)
      INFO(4) = 2
      RETLBL  = 70
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal Sending L,IDX(2,N)", IDX(2,N)
      call checksum("@-Dulal-send", VECS(1,IDX(3,N)),NLEN,s)
#endif 
      RETURN
#ifdef _DEBUG_LVL0
      Write(6,*) "@-Dulal Recieving L"
      call checksum("@-Dulal-recieve", VECS(1,2),NLEN,s)
      call checksum("@-Dulal-recieve", VECS(1,IDX(2,N)),NLEN,s)
#endif 
#ifdef _DEBUG_LVL0
      Write(6,"(1x,i8,1x,i8)") N, Nmax
#endif 

 70   IF ((VF.NE.0).AND.(N.LE.NMAX)) WRITE (VF,'(A11,I8)') 'Rebuilding:'
     $, N+1

      CALL DULAL1 (NDIM,NLEN,M,MAXVW,NUMCHK,N,L,LSTAR,NL,NLSTAR,VF,IERR,
     $            ADJUST,NORM1,NORMA,TMAX,TMIN,TNRM,DWK,IWK,IDX,VECS)
      IF ((IERR.NE.0).AND.(IERR.NE.16).AND.(IERR.NE.32)
     $    .AND.(IERR.NE.48)) GO TO 110
      MVWBLT = MAX0(MVWBLT,NL-IWK(IDX(1,L-1+1),1))
      IWK(IDX(1,N),2) = NLSTAR
C
C     Update the counter for steps taken.
C
      NLAL = N
      NMAX = MAX0(NMAX,NLAL)
C
C     Extract the column of H.
C
      CALL DAXPBY (N+1,HWK(1,N),DZERO,HWK(1,N),DZERO,HWK(1,N))
      DO 80 I = NLSTAR, N+1
         HWK(I,N) = DWK(IDX(1,I),2*M+IDX(1,N))
 80   CONTINUE
C
C     Do the next step.
C
      IF ((IERR.NE.0).OR.(N.GE.NLIM)) GO TO 90
      N = N + 1
      GO TO 50
C
C     The eigenvalue code starts here.
C
 90   NOLD = 0
      RETURN 
 100  NEIG = 0
      WRITE (6,'()')
      WRITE (6,'(A43,I4)') 'Number of Lanczos steps completed        : '
     $, NLAL
      WRITE (6,'(A43,$)') 'Number of eigenvalues to compute (0=done): '
      RETURN
      READ (5,*) NEIG
      IF ((NEIG.GT.0).AND.(NEIG.LE.NLAL)) THEN
         CALL DULALE (MAXN,NEIG,HWK,TF,VF)
         IF (NEIG.EQ.0) GO TO 100
         IF (NOLD.NE.0) CALL DULALC (MAXN,NEIG,NOLD,HWK,TF,VF)
         NOLD = NEIG
         CALL DAXPBY (NOLD,HWK(1,2*MAXN+7),DONE,HWK(1,2*MAXN+4),DZERO,
     $HWK(1,2*MAXN+7))
         CALL DAXPBY (NOLD,HWK(1,2*MAXN+8),DONE,HWK(1,2*MAXN+5),DZERO,
     $HWK(1,2*MAXN+8))
      ELSE IF (NEIG.EQ.0) THEN
         GO TO 110
      END IF
      GO TO 100
C
C     Done.
C
 110  RETLBL  = 0
      NLIM    = NLAL
      INFO(1) = IERR
      NORM    = NORM1
      MAXVW   = MVWBLT

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Exit dulal ---"
        Write(6,*)
#endif
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE DULALC (MAXN,NEIG,NOLD,HWK,TF,VF)
C
C     Purpose:
C     This subroutine finds eigenvalues common  to the last two runs of
C     the eigenvalue routine.
C
C     Parameters:
C     MAXN  = the dimensioned size of the array HWK (input).
C     NEIG  = the number of eigs in the current run (input).
C     NOLD  = on entry, the number of eigs in the previous run; on exit
C             the number of common eigs (input/output).
C     HWK   = double precision work array (input/output).
C     TF    = if not 0, then the routine will output the common eigs to
C             this unit (input).
C     VF    = if not 0, then the routine will produce  verbose messages
C             on this unit (input).
C
C     External routines used:
C     subroutine dhpsrt(n,wr,wi,wm)
C        Library routine, sorts wr and wi and wm.
C
C     Noel M. Nachtigal
C     August 6, 1993
C
C**********************************************************************
C
      INTRINSIC DABS
      EXTERNAL DHPSRT
C
      INTEGER MAXN, NEIG, NOLD, TF, VF
      DOUBLE PRECISION HWK(MAXN,2*MAXN+8)
C
C     Local variables.
C
      INTEGER I, J, K
      DOUBLE PRECISION DISTM, EIMAG, EMAGN, EREAL, OIMAG, OREAL, OMAGN
      DOUBLE PRECISION TOL, TOLM
C
C     Get the number of eigenvalues to find.
C
 10   WRITE (6,'()')
      WRITE (6,'(A43,$)') 'Find the common eigenvalues (1=Yes/0=No) ? '
      READ (5,*) I
      IF (I.NE.1) RETURN
      WRITE (6,'(A43,$)') 'Enter separation tolerance               : '
      READ (5,*) TOL
C
C     Check for common eigenvalues.
C
      K = 0
      DO 30 I = 1, NEIG
         EREAL = HWK(I,2*MAXN+4)
         EIMAG = HWK(I,2*MAXN+5)
         EMAGN = HWK(I,2*MAXN+6)
         DO 20 J = 1, NOLD
            OREAL = HWK(J,2*MAXN+7)
            OIMAG = HWK(J,2*MAXN+8)
            OMAGN = DSQRT(OREAL**2 + OIMAG**2)
            TOLM  = 0.5D0 * TOL * (EMAGN + OMAGN)
            IF (DABS(EMAGN-OMAGN).GE.TOLM) GO TO 20
            DISTM = DSQRT((EREAL-OREAL)**2 + (EIMAG-OIMAG)**2)
            IF (DISTM.GE.TOLM) GO TO 20
            K        = K + 1
            HWK(K,2*MAXN+1) = EREAL
            HWK(K,2*MAXN+2) = EIMAG
            HWK(K,2*MAXN+3) = EMAGN
            GO TO 30
 20      CONTINUE
 30   CONTINUE
C
C     Sort the final eigenvalues by decreasing real part.
C
      CALL DHPSRT (K,HWK(1,2*MAXN+1),HWK(1,2*MAXN+2),HWK(1,2*MAXN+3))
C
C     Output the results.
C
      IF (VF.NE.0) WRITE (VF,'(A43,I4)')
     $             'Number of common eigenvalues found       : ', K
      IF (TF.NE.0) THEN
         WRITE (TF,'(A30,E12.7)') 'Common eigenvalues, tolerance: ', TOL
         WRITE (TF,'(2E25.17)') (HWK(I,2*MAXN+1),HWK(I,2*MAXN+2),I=1,K)
         WRITE (TF,'()')
      END IF
      GO TO 10
C
      END
C
C**********************************************************************
C
      SUBROUTINE DULALE (MAXN,NEIG,HWK,TF,VF)
C
C     Purpose:
C     This subroutine computes the  eigenvalues estimates obtained from
C     the NEIG x NEIG upper left corner of the Lanczos matrix.  It uses
C     the heuristic proposed  by Cullum and Willoughby  to separate the
C     genuine eigenvalues from the spurious ones.
C
C     Parameters:
C     MAXN  = the dimensioned size of the array HWK (input).
C     NEIG  = the size of  the submatrix  whose eigs  are computed.  On
C             exit, it is the number of eigs found (input/output).
C     HWK   = double precision work array (input/output).
C     TF    = if not 0, then the routine  will output  the Lanczos eigs
C             to this unit (input).
C     VF    = if not 0, then the routine will produce  verbose messages
C             on this unit (input).
C
C     External routines used:
C     subroutine balanc(nm,n,a,low,igh,scale)
C        EISPACK routine, balances a matrix.
C     subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
C        EISPACK routine, computes eigenvalues of an upper Hessenberg
C        matrix.
C     subroutine daxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     subroutine dhpsrt(n,wr,wi,wm)
C        Library routine, sorts wr and wi and wm.
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C
C     Noel M. Nachtigal
C     August 6, 1993
C
C**********************************************************************
C
      INTRINSIC DABS, DSQRT
      EXTERNAL BALANC, DAXPBY, DHPSRT, DLAMCH, HQR
      DOUBLE PRECISION DLAMCH
C
      INTEGER MAXN, NEIG, TF, VF
      DOUBLE PRECISION HWK(MAXN,2*MAXN+8)
C
C     Miscellaneous parameters.
C
      DOUBLE PRECISION DONE, DZERO
      PARAMETER (DONE = 1.0D0,DZERO = 0.0D0)
C
C     Local variables.
C
      INTEGER I, IERR, J, K, NBEG, NCHK, NEND, NLAL
      DOUBLE PRECISION DISTM, EIMAG, EMAGN, EREAL, TOL, TOLM
      LOGICAL MULTIP
C
C     Compute the separation tolerance.
C
      TOL = DSQRT(DLAMCH('E'))
C
C     Extract the NEIG x NEIG matrix.
C
      DO 10 I = 1, NEIG
         CALL DAXPBY (I+1,HWK(1,MAXN+I),DONE,HWK(1,I),DZERO,HWK(1,MAXN+I
     $))
 10   CONTINUE
C
C     Balance it and compute its eigenvalues.
C
      IF (VF.NE.0) WRITE (VF,'(A24)') 'Computing eigenvalues...'
      CALL BALANC (MAXN,NEIG,HWK(1,MAXN+1),I,J,HWK(1,2*MAXN+1))
      CALL HQR (MAXN,NEIG,I,J,HWK(1,MAXN+1),HWK(1,2*MAXN+1),HWK(1,2*MAXN
     $+2),IERR)
      IF (IERR.NE.0) THEN
         WRITE (6,'(A23,I2)') 'Error return from HQR: ', IERR
         RETURN
      END IF
C
C     Compute the magnitude of the eigenvalues.
C
      DO 20 I = 1, NEIG
         HWK(I,2*MAXN+3) = DSQRT(HWK(I,2*MAXN+1)**2 + HWK(I,2*MAXN+2)**2
     $)
 20   CONTINUE
C
C     Sort the eigenvalues in decreasing order of magnitude.
C
      CALL DHPSRT (NEIG,HWK(1,2*MAXN+3),HWK(1,2*MAXN+1),HWK(1,2*MAXN+2))
      IF (TF.NE.0) THEN
         WRITE (TF,'(A33,I4)') 'All Lanczos eigenvalues at step: ', NEIG
     $ 
         WRITE (TF,'(2E25.17)') (HWK(I,2*MAXN+1),HWK(I,2*MAXN+2),I=1,
     $NEIG)
         WRITE (TF,'()')
      END IF
C
C     Separate the repeated eigenvalues from the simple ones.
C
      I    = 1
      NEND = 0
      NLAL = NEIG
      NBEG = NEIG + 1
 30   IF (I.GT.NLAL) GO TO 70
C
C     Outer loop, over all eigenvalues.
C
      J      = I + 1
      EREAL  = HWK(I,2*MAXN+1)
      EIMAG  = HWK(I,2*MAXN+2)
      EMAGN  = HWK(I,2*MAXN+3)
      MULTIP = .FALSE.
 40   IF (J.GT.NLAL) GO TO 60
C
C     Inner loop, find eigenvalues J close to eigenvalue I.  "Close" is
C     defined relative to the average of the two points,  so as to be a
C     symmetric relationship.
C
      TOLM = 0.5D0 * TOL * (EMAGN + HWK(J,2*MAXN+3))
C
C     If the difference  in magnitude  between I and J  is greater than
C     the tolerance, skip the remaining eigenvalues, they will  be even
C     further away (the eigenvalues are assumed sorted in some order of
C     magnitude).
C
      IF (EMAGN-HWK(J,2*MAXN+3).GE.TOLM) GO TO 60
C
C     Compute the distance between the two eigenvalues.
C
      DISTM = DSQRT((EREAL-HWK(J,2*MAXN+1))**2 + (EIMAG-HWK(J,2*MAXN+2))
     $**2)
C
C     Check the eigenvalue J.
C
      IF (DISTM.GE.TOLM) THEN
         J = J + 1
      ELSE
         MULTIP = .TRUE.
         DO 50 K = J, NLAL-1
            HWK(K,2*MAXN+1) = HWK(K+1,2*MAXN+1)
            HWK(K,2*MAXN+2) = HWK(K+1,2*MAXN+2)
            HWK(K,2*MAXN+3) = HWK(K+1,2*MAXN+3)
 50      CONTINUE
         NLAL = NLAL - 1
      END IF
      GO TO 40
C
C     Collect the multiple eigenvalues in ELAL(1:NEND),  and the simple
C     ones in ELAL(NBEG:NEIG).
C
 60   IF (MULTIP) THEN
         NEND        = NEND + 1
         HWK(NEND,2*MAXN+4) = EREAL
         HWK(NEND,2*MAXN+5) = EIMAG
         HWK(NEND,2*MAXN+6) = EMAGN
      ELSE
         NBEG = NBEG - 1
         HWK(NBEG,2*MAXN+4) = EREAL
         HWK(NBEG,2*MAXN+5) = EIMAG
         HWK(NBEG,2*MAXN+6) = EMAGN
      END IF
      I = I + 1
      GO TO 30
 70   CONTINUE
C
C     Check whether they are all considered genuine.
C
      IF (NEND.EQ.NLAL) GO TO 120
C
C     Extract the (NEIG-1) x (NEIG-1) matrix.
C
      NCHK = NEIG - 1
      DO 80 I = 1, NCHK
         CALL DAXPBY (I+1,HWK(1,MAXN+I),DONE,HWK(2,I+1),DZERO,HWK(1,MAXN
     $+I))
 80   CONTINUE
C
C     Balance it and compute its eigenvalues.
C
      IF (VF.NE.0) WRITE (VF,'(A30)') 'Computing check eigenvalues...'
      CALL BALANC (MAXN,NCHK,HWK(1,MAXN+1),I,J,HWK(1,2*MAXN+1))
      CALL HQR (MAXN,NCHK,I,J,HWK(1,MAXN+1),HWK(1,2*MAXN+1),HWK(1,2*MAXN
     $+2),IERR)
      IF (IERR.NE.0) THEN
         WRITE (6,'(A23,I2)') 'Error return from HQR: ', IERR
         RETURN
      END IF
C
C     Compute the magnitude of the check eigenvalues.
C
      DO 90 I = 1, NCHK
         HWK(I,2*MAXN+3) = DSQRT(HWK(I,2*MAXN+1)**2 + HWK(I,2*MAXN+2)**2
     $)
 90   CONTINUE
C
C     Sort the eigenvalues in decreasing order of magnitude.
C
      CALL DHPSRT (NCHK,HWK(1,2*MAXN+3),HWK(1,2*MAXN+1),HWK(1,2*MAXN+2))
C
C     Check the  Lanczos eigenvalues that were simple;  keep those that
C     appear only in the Lanczos matrix.
C     The approach used is brute force; it could be more sophisticated.
C
      DO 110 I = NBEG, NEIG
         EREAL = HWK(I,2*MAXN+4)
         EIMAG = HWK(I,2*MAXN+5)
         EMAGN = HWK(I,2*MAXN+6)
         DO 100 J = 1, NCHK
            TOLM = 0.5D0 * TOL * (EMAGN + HWK(J,2*MAXN+3))
            IF (DABS(EMAGN-HWK(J,2*MAXN+3)).GE.TOLM) GO TO 100
            DISTM = DSQRT((EREAL-HWK(J,2*MAXN+1))**2 + (EIMAG-HWK(J,
     $2*MAXN+2))**2)
            IF (DISTM.LT.TOLM) GO TO 110
 100     CONTINUE
C
C     It appears only in the Lanczos matrix, so it is good.
C
         NEND        = NEND + 1
         HWK(NEND,2*MAXN+4) = EREAL
         HWK(NEND,2*MAXN+5) = EIMAG
         HWK(NEND,2*MAXN+6) = EMAGN
 110  CONTINUE
C
C     Sort the final eigenvalues by decreasing real part.
C
 120  CALL DHPSRT (NEND,HWK(1,2*MAXN+4),HWK(1,2*MAXN+5),HWK(1,2*MAXN+6))
C
C     Output the results.
C
      IF (VF.NE.0) WRITE (VF,'(A43,I4)')
     $             'Number of Lanczos eigenvalues found      : ', NEND
      IF (TF.NE.0) THEN
         WRITE (TF,'(A30,I4)') 'Filtered eigenvalues at step: ', NEIG
         WRITE (TF,'(2E25.17)') (HWK(I,2*MAXN+4),HWK(I,2*MAXN+5),I=1,
     $NEND)
         WRITE (TF,'()')
      END IF
      NEIG = NEND
C
C     Done.
C
      RETURN
      END
C
C**********************************************************************
C
      DOUBLE PRECISION FUNCTION DULALH(I,J)
C
C     Purpose:
C     Returns the recurrence coefficients for inner vectors.
C
C     Parameters:
C     I = row index of the coefficient (input).
C     J = column index of the coefficient (input).
C
C     Noel M. Nachtigal
C     July 9, 1993
C
C**********************************************************************
C
C
C     Common block DULALX.
C
      DOUBLE PRECISION NORMA
      COMMON /DULALX/NORMA
C
C
      INTEGER I, J
C
      IF ((I.LT.1).OR.(J.LT.1)) THEN
         DULALH = 0.0D0
      ELSE IF (I.EQ.J) THEN
         DULALH = NORMA
      ELSE
         DULALH = 0.0D0
      END IF
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE DULAL1 (NDIM,NLEN,M,MAXVW,NUMCHK,N,L,LSTAR,NL,NLSTAR,
     $                  VF,IERR,ADJUST,NORM1,NORMA,TMAX,TMIN,TNRM,DWK,
     $                  IWK,IDX,VECS)
C
C     Purpose:
C     This subroutine builds a new pair of vectors v_{n+1} and w_{n+1}.
C     It is called internally by  the Lanczos code, and is not meant to
C     to be called directly by the user code.
C
C     Parameters:
C        See descriptions in "dulal.doc" or "duqmr.doc".
C
C     External routines used:
C     subroutine daxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     double precision ddot(n,x,incx,y,incy)
C        BLAS-1 routine, computes y^H * x.
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C     double precision dnrm2(n,x,incx)
C        BLAS-1 routine, computes the 2-norm of x.
C     subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
C        LINPACK routine, computes the QR factorization of x.
C     subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
C        LINPACK routine, applies the QR factorization of x.
C     subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
C        LINPACK routine, computes the SVD of x.
C     subroutine dulal2 (ndim,nlen,m,dwk,iwk,idx,vecs)
C        Forces closure of an inner block.
C     double precision dulalh(i,n)
C        User-supplied routine, computes inner recurrence coefficients.
C
C     Noel M. Nachtigal
C     July 9, 1993
C
C**********************************************************************
C
      INTRINSIC DABS, DMAX1, MAX0
      EXTERNAL DAXPBY, DDOT, DLAMCH, DNRM2, DQRDC, DQRSL, DSVDC
      EXTERNAL DULAL2, DULALH
      DOUBLE PRECISION DDOT, DLAMCH, DNRM2, DULALH
C
      INTEGER IERR, L, LSTAR, M, MAXVW, N, NDIM, NL, NLEN, NLSTAR
      INTEGER IDX(3,3), IWK(M,4), NUMCHK, VF
      DOUBLE PRECISION ADJUST, NORM1, NORMA, TMAX, TMIN, TNRM
      DOUBLE PRECISION DWK(M,5*M+7), VECS(NDIM,2*MAXVW+4+2)
C
C     Miscellaneous parameters.
C
      DOUBLE PRECISION DONE, DZERO
      PARAMETER (DONE = 1.0D0,DZERO = 0.0D0)
C
C     Local variables.
C
      INTEGER I, IJ, J, LBLKSZ, NP1
      DOUBLE PRECISION DTMP, DTMP1, DTMP2, DTMP3, DTMP4, DTMP5
      DOUBLE PRECISION IDTMP1, IDTMP2, IDTMP3, IDTMP4
      LOGICAL IBUILT, INNER, RERUN
C
C     Initialize local variables.
C
      NP1        = N + 1
      DWK(IDX(1,N),5*M+7) = DZERO
      RERUN      = .FALSE.
      IBUILT     = .FALSE.
      IERR       = 0
      LBLKSZ     = N - NL + 1
C
C     Update the norm estimates.
C
      IF (N.LE.NUMCHK) THEN
         DTMP1 = DNRM2(NLEN,VECS(1,1),1) * DWK(IDX(1,N),5*M+5)
         DTMP2 = DNRM2(NLEN,VECS(1,2),1) * DWK(IDX(1,N),5*M+6)
         NORMA = DMAX1(DTMP1,DTMP2,NORMA)
      END IF
      NORM1 = DMAX1(ADJUST*NORMA,NORM1)
C
C     Clear current column of H.
C
#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,i4,1x,i4)")" @-dulal1-The values of M and N",
     +                            M,N
#endif 

      DO 100 I = 1, M
         DWK(IDX(1,I),2*M+IDX(1,N)) = DZERO
 100  CONTINUE
C
C     Update D^{(n-1)} to D^{(n)}.
C
      DO 120 I = NL, N-1
         DTMP = DZERO
         DO 110 J = NL, N-1
            DTMP = DTMP + DWK(IDX(1,I),IDX(1,J)) * DWK(IDX(1,J),2*M+
     $IDX(1,N-1))
 110     CONTINUE
         DTMP     = ( DWK(IDX(1,I),M+IDX(1,N-1)) - DTMP ) / DWK(IDX(1,N)
     $,2*M+IDX(1,N-1))
         DWK(IDX(1,I),IDX(1,N)) = DTMP
         DWK(IDX(1,N),IDX(1,I)) = DTMP * DWK(IDX(1,N),5*M+4) / DWK(IDX(1
     $,I),5*M+4)
 120  CONTINUE
C
C     Compute l^\star.
C
      LSTAR  = MAX0(1,L-1)
      NLSTAR = IWK(IDX(1,LSTAR+1),1)
C
C     Compute F_{n,1:n-1} and F_{1:n-1,n}.  Delay computing F_{nn}.
C
      DO 140 I = MAX0(1,NL-1), N-1
         DTMP = DZERO
         DO 130 J = NL, I+1
            DTMP = DTMP + DWK(IDX(1,N),IDX(1,J)) * DWK(IDX(1,J),2*M+
     $IDX(1,I))
 130     CONTINUE
         DWK(IDX(1,N),M+IDX(1,I)) = DTMP
         DWK(IDX(1,I),M+IDX(1,N)) = DTMP * DWK(IDX(1,I),5*M+4) / 
     $DWK(IDX(1,N),5*M+4)
         NORMA    = DMAX1(NORMA,DABS(DWK(IDX(1,N),M+IDX(1,I))),
     $DABS(DWK(IDX(1,I),M+IDX(1,N))))
         NORM1    = DMAX1(ADJUST*NORMA,NORM1)
 140  CONTINUE
C
C     Check whether D_l is nonsingular.
C
      IF (LBLKSZ.EQ.1) THEN
         INNER = DABS(DWK(IDX(1,NL),IDX(1,NL))).EQ.DZERO
      ELSE
         DO 160 I = NL, N
            DO 150 J = NL, N
               DWK(I-NL+1,4*M+J-NL+1) = DWK(IDX(1,I),IDX(1,J))
 150        CONTINUE
 160     CONTINUE
         CALL DSVDC (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,5*M+1),DWK(1,5*M
     $+3),DZERO,0,DZERO,0,DWK(1,5*M+2),0,I)
         Write(6,*) "What is the value of I after dsvdc", i
         IF (I.NE.0) THEN
            IERR = -I
            GOTO 370
         END IF
         DTMP = DZERO
         IF (DWK(1,5*M+1).NE.DZERO) DTMP = DWK(LBLKSZ,5*M+1) / DWK(1,5*M
     $+1)
         INNER = DTMP.LT.(DBLE(NLEN) * DLAMCH('E'))
      END IF
      IF ((VF.NE.0).AND.INNER)
     $    WRITE (VF,'(A31)') '... moment matrix D is singular'
C
C     Compute H_{n_{l-1}:n_l-1,n}, for n_l > 1.
C     Also, build the part common to both inner and regular vectors.
C     Assume that VECS(1,1) = A v_n and VECS(1,2) = A^T w_n.
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The starting vecs at start dulal1 processing"
      call checksum("@-Dulal1",VECS(1,1),Nlen,s)
      call checksum("@-Dulal1",VECS(1,2),Nlen,s)
      write(6,*)
#endif 

      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DZERO
      IF (NL.EQ.1) GO TO 180
      IF (NL.EQ.NLSTAR+1) THEN
         DWK(IDX(1,NLSTAR),2*M+IDX(1,N)) = DWK(IDX(1,NLSTAR),M+IDX(1,N))
     $ / DWK(IDX(1,NLSTAR),IDX(1,NLSTAR))
         DTMP = DWK(IDX(1,NLSTAR),2*M+IDX(1,N)) * DWK(IDX(1,NLSTAR),5*M+
     $5) / DWK(IDX(1,N),5*M+5)
         CALL DAXPBY (NLEN,VECS(1,1),DONE,VECS(1,1),-DTMP,VECS(1,IDX(2,
     $NLSTAR)))
         DTMP = DWK(IDX(1,NLSTAR),2*M+IDX(1,N)) * DWK(IDX(1,NLSTAR),5*M+
     $6) / DWK(IDX(1,N),5*M+6) * DWK(IDX(1,N),5*M+4) / DWK(IDX(1,NLSTAR)
     $,5*M+4)
         CALL DAXPBY (NLEN,VECS(1,2),DONE,VECS(1,2),-DTMP,VECS(1,IDX(3,
     $NLSTAR)))
      ELSE IF (NL.GT.NLSTAR+1) THEN
C
C     If a block of size larger than 1 was just closed, combine all its
C     vectors into one.
C
         DTMP2 = DZERO
         DO 170 I = NL-1, NLSTAR, -1
            DWK(IDX(1,I),2*M+IDX(1,N)) = DWK(IDX(1,NL-1),M+IDX(1,N)) * 
     $DWK(IDX(1,I),3*M+IDX(1,L-1))
            IF (N.NE.NL) GO TO 170
            DTMP1 = DWK(IDX(1,I),3*M+IDX(1,L-1)) * DWK(IDX(1,I),5*M+5)
            CALL DAXPBY (NLEN,VECS(1,IDX(2,NL-1)),DTMP2,VECS(1,IDX(2,NL-
     $1)),DTMP1,VECS(1,IDX(2,I)))
            DTMP1 = DWK(IDX(1,I),3*M+IDX(1,L-1)) * DWK(IDX(1,I),5*M+6) /
     $ DWK(IDX(1,I),5*M+4)
            CALL DAXPBY (NLEN,VECS(1,IDX(3,NL-1)),DTMP2,VECS(1,IDX(3,NL-
     $1)),DTMP1,VECS(1,IDX(3,I)))
            DTMP2 = DONE
 170     CONTINUE
         DTMP1 = DWK(IDX(1,NL-1),M+IDX(1,N)) / DWK(IDX(1,N),5*M+5)
         CALL DAXPBY (NLEN,VECS(1,1),DONE,VECS(1,1),-DTMP1,VECS(1,IDX(2,
     $NL-1)))
         DTMP1 = DWK(IDX(1,NL-1),M+IDX(1,N)) * DWK(IDX(1,N),5*M+4) / 
     $DWK(IDX(1,N),5*M+6)
         CALL DAXPBY (NLEN,VECS(1,2),DONE,VECS(1,2),-DTMP1,VECS(1,IDX(3,
     $NL-1)))
      END IF
C
C     Compute F_{nn}.
C
 180  DWK(IDX(1,N),M+IDX(1,N)) = DDOT(NLEN,VECS(1,IDX(3,N)),1,VECS(1,1),
     $1) * DWK(IDX(1,N),5*M+5) * DWK(IDX(1,N),5*M+6)
      NORMA    = DMAX1(NORMA,DABS(DWK(IDX(1,N),M+IDX(1,N))))
      NORM1    = DMAX1(ADJUST*NORMA,NORM1)
      IF (INNER) GO TO 320
C
C     Compute H_{n_l:n,n}.
C
      IWK(IDX(1,L+1+1),1) = NP1
      IF (LBLKSZ.EQ.1) THEN
         DWK(IDX(1,NL),2*M+IDX(1,N)) = DWK(IDX(1,NL),M+IDX(1,N)) / 
     $DWK(IDX(1,NL),IDX(1,NL))
      ELSE
         DO 200 J = NL, N
            DWK(J-NL+1,5*M+1) = DWK(IDX(1,J),M+IDX(1,N))
            DO 190 IJ = NL, N
               DWK(J-NL+1,4*M+IJ-NL+1) = DWK(IDX(1,J),IDX(1,IJ))
 190        CONTINUE
 200     CONTINUE
         CALL DQRDC (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,5*M+3),0,DZERO,0
     $)
         CALL DQRSL (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,5*M+3),DWK(1,5*M
     $+1),
     $               DZERO,DWK(1,5*M+1),DWK(1,5*M+2),DWK(1,5*M+1),DZERO,
     $100,J)
         DO 210 J = NL, N
            DWK(IDX(1,J),2*M+IDX(1,N)) = DWK(J-NL+1,5*M+2)
 210     CONTINUE
      END IF
C
C     Either D_l is nonsingular, or VW is being rerun.  For the latter,
C     just finish building the vectors.  For the former,  the check for
C     H_{n_{l-1}:n,n} could be done.  However, the look-ahead  strategy
C     requires  the smaller of  the norm checks for H_{n_{l-1}:n,n} and
C     H_{n_l:n,n+1}, so both norms are computed first and then checked.
C
C     Build regular vectors.
C
      DTMP = DWK(IDX(1,N),2*M+IDX(1,N))
      CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,1),-DTMP,VECS(1,
     $IDX(2,N)))
      CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,2),-DTMP,VECS(1,
     $IDX(3,N)))
      DO 220 I = NL, N-1
         DTMP = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),5*M+5) / 
     $DWK(IDX(1,N),5*M+5)
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1)),-
     $DTMP,VECS(1,IDX(2,I)))
         DTMP = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),5*M+6) / 
     $DWK(IDX(1,N),5*M+6) * DWK(IDX(1,N),5*M+4) / DWK(IDX(1,I),5*M+4)
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1)),-
     $DTMP,VECS(1,IDX(3,I)))
 220  CONTINUE
C
C     Compute scale factors for the new vectors.
C
 230  DTMP3      = DNRM2(NLEN,VECS(1,IDX(2,NP1)),1)
      DTMP4      = DNRM2(NLEN,VECS(1,IDX(3,NP1)),1)
      DTMP1      = DWK(IDX(1,N),5*M+5) * DTMP3
      DTMP2      = DWK(IDX(1,N),5*M+6) * DTMP4
      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DTMP1
      IF (DTMP1.LT.TNRM) IERR = IERR + 16
      IF (DTMP2.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GOTO 360
      DWK(IDX(1,NP1),5*M+4)     = DWK(IDX(1,N),5*M+4) * DTMP1 / DTMP2
      DWK(IDX(1,NP1),IDX(1,NP1)) = DDOT(NLEN,VECS(1,IDX(3,NP1)),1,VECS(1
     $,IDX(2,NP1)),1) / ( DTMP3 * DTMP4 )
      IF ((DTMP3.GE.TMAX).OR.(DTMP3.LE.TMIN)) THEN
         DTMP = DONE / DTMP3
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DTMP,VECS(1,IDX(2,NP1)),
     $DZERO,VECS(1,IDX(2,NP1)))
         DTMP3 = DONE
      END IF
      IF ((DTMP4.GE.TMAX).OR.(DTMP4.LE.TMIN)) THEN
         DTMP = DONE / DTMP4
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DTMP,VECS(1,IDX(3,NP1)),
     $DZERO,VECS(1,IDX(3,NP1)))
         DTMP4 = DONE
      END IF
      DWK(IDX(1,NP1),5*M+5) = DONE / DTMP3
      DWK(IDX(1,NP1),5*M+6) = DONE / DTMP4
C
C     Compute the last column of D_l^{-1}.
C
      IF (LBLKSZ.EQ.1) THEN
         DWK(IDX(1,NL),3*M+IDX(1,L)) = DONE / DWK(IDX(1,NL),IDX(1,NL))
      ELSE
         DO 240 I = 1, LBLKSZ-1
            DWK(I,5*M+1) = DZERO
 240     CONTINUE
         DWK(LBLKSZ,5*M+1) = DONE
         CALL DQRSL (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,5*M+3),DWK(1,5*M
     $+1),
     $               DZERO,DWK(1,5*M+1),DWK(1,5*M+2),DWK(1,5*M+1),DZERO,
     $100,J)
         DO 250 I = NL, N
            DWK(IDX(1,I),3*M+IDX(1,L)) = DWK(I-NL+1,5*M+2)
 250     CONTINUE
      END IF
C
C     If VW is being rerun, then skip to the end.
C
      IF (RERUN) GO TO 360
C
C     Compute the norm of H_{n_{l-1}:n,n}.
C
      DTMP1 = DZERO
      DTMP2 = DZERO
      DO 260 I = IWK(IDX(1,L-1+1),1), N
         DTMP  = DABS(DWK(IDX(1,I),2*M+IDX(1,N)))
         DTMP1 = DTMP1 + DTMP
         DTMP2 = DTMP2 + DTMP / DWK(IDX(1,I),5*M+4)
 260  CONTINUE
      DTMP2 = DTMP2 * DWK(IDX(1,N),5*M+4)
C
C     Build the 2nd term for the next step, H_{n_l:n,n+1}, regular.
C     Compute the norm of H_{n_l:n,n+1}.
C
      DTMP3 = DZERO
      DTMP4 = DZERO
      DTMP5 = DWK(IDX(1,NP1),IDX(1,NP1)) * DWK(IDX(1,NP1),2*M+IDX(1,N)) 
     $* DWK(IDX(1,N),5*M+4) / DWK(IDX(1,NP1),5*M+4)
      DO 270 I = NL, N
         DTMP  = DABS(DTMP5 * DWK(IDX(1,I),3*M+IDX(1,L)))
         DTMP3 = DTMP3 + DTMP
         DTMP4 = DTMP4 + DTMP / DWK(IDX(1,I),5*M+4)
 270  CONTINUE
      DTMP4 = DTMP4 * DWK(IDX(1,NP1),5*M+4)
C
C     Check H_{n_{l-1}:n,n} and H_{n_l:n-1,n+1}.
C
      DTMP  = DMAX1(DTMP1,DTMP2,DTMP3,DTMP4)
      INNER = DTMP.GT.NORM1
      IF (.NOT.INNER) GO TO 360
      DWK(IDX(1,N),5*M+7) = DTMP
C
C     If H_{n_{l-1}:n,n} is bad, build inner vectors.
C
      IF (DMAX1(DTMP1,DTMP2).GT.NORM1) GO TO 320
C
C     If H_{n_l:n-1,n+1} is bad, check the inner vectors.
C     This only applies if n_l > 1.
C
      IF (NL.LE.1) GO TO 320
C
C     Build the inner vectors to compute the 2nd term at the next step.
C     Get the coefficients H_{n_l:n+1,n} in a temporary location.
C     Build inner vectors in VECS(1,1) and VECS(1,2).
C
      IBUILT = .TRUE.
      DO 280 I = NL, N
         DWK(IDX(1,I),2*M+IDX(1,NP1)) = DULALH(I,N)
         DTMP       = DWK(IDX(1,I),2*M+IDX(1,NP1)) * DWK(IDX(1,I),5*M+5)
     $ / DWK(IDX(1,N),5*M+5)
         CALL DAXPBY (NLEN,VECS(1,1),DONE,VECS(1,1),-DTMP,VECS(1,IDX(2,I
     $)))
         DTMP       = DWK(IDX(1,I),2*M+IDX(1,NP1)) * DWK(IDX(1,I),5*M+6)
     $ / DWK(IDX(1,N),5*M+6) * DWK(IDX(1,N),5*M+4) / DWK(IDX(1,I),5*M+4)
         CALL DAXPBY (NLEN,VECS(1,2),DONE,VECS(1,2),-DTMP,VECS(1,IDX(3,I
     $)))
 280  CONTINUE
      IDTMP3       = DNRM2(NLEN,VECS(1,1),1)
      IDTMP4       = DNRM2(NLEN,VECS(1,2),1)
      IDTMP1       = DWK(IDX(1,N),5*M+5) * IDTMP3
      IDTMP2       = DWK(IDX(1,N),5*M+6) * IDTMP4
      DWK(IDX(1,NP1),2*M+IDX(1,NP1)) = IDTMP1
      IF (IDTMP1.LT.TNRM) IERR = IERR + 16
      IF (IDTMP2.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GOTO 330
C
C     Build the 2nd term for the next step, H_{n_{l-1}:n_l-1,n+1},
C     inner.  Compute the norm of H_{n_{l-1}:n_l-1,n+1}.
C
      DTMP = DZERO
      DO 290 J = NLSTAR, NL-1
         DTMP = DTMP + DWK(IDX(1,NL),IDX(1,J)) * DWK(IDX(1,I),2*M+IDX(1,
     $N))
 290  CONTINUE
      DO 300 J = NL, N
         DTMP = DTMP + DWK(IDX(1,NL),IDX(1,J)) * DWK(IDX(1,I),2*M+IDX(1,
     $NP1))
 300  CONTINUE
      DTMP1 = DZERO
      DTMP2 = DZERO
      DTMP3 = DMAX1(DTMP3,DTMP4)
      DTMP  = ( DWK(IDX(1,NL),M+IDX(1,N)) - DTMP ) / DWK(IDX(1,NP1),2*M+
     $IDX(1,NP1))
      DTMP5 = DTMP * DWK(IDX(1,NL),2*M+IDX(1,NL-1)) * DWK(IDX(1,NL-1),
     $5*M+4) / DWK(IDX(1,NL),5*M+4)
      DO 310 I = IWK(IDX(1,L-1+1),1), NL-1
         DTMP  = DABS(DWK(IDX(1,I),3*M+IDX(1,L-1)) * DTMP5)
         DTMP1 = DTMP1 + DTMP
         DTMP2 = DTMP2 + DTMP / DWK(IDX(1,I),5*M+4)
 310  CONTINUE
      DTMP2 = DTMP2 * DWK(IDX(1,N),5*M+4)
C
C     Compare the inner and regular versions of the 2nd term at the next
C     step.  Build the vector corresponding to the smaller term.
C
      INNER = DTMP3.GT.DMAX1(DTMP1,DTMP2)
      IF (.NOT.INNER) GO TO 360
C
C     Build inner vectors.
C     Check whether the block has to be forced to close.
C
 320  IF (VF.NE.0) WRITE (VF,'(A7,I5,A9)') 'Vector ',NP1,' is inner'
      IF (NP1-NL.EQ.MAXVW) THEN
         CALL DULAL2 (NDIM,NLEN,M,N,L,LSTAR,NL,NLSTAR,VF,IERR,ADJUST,
     $               NORM1,DWK,IWK,IDX,VECS)
         IF (IERR.NE.0) GO TO 370
         LBLKSZ = N - NL + 1
         INNER  = .FALSE.
         RERUN  = .TRUE.
         NP1    = N + 1
         GO TO 230
      END IF
C
C     The temporary vectors  contain either just partial inner vectors,
C     or the  completed ones,  depending on whether  IBUILT is  TRUE or
C     FALSE.  In either case, replace the regular vectors.
C
 330  CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,1),DZERO,VECS(1,
     $IDX(2,NP1)))
      CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,2),DZERO,VECS(1,
     $IDX(3,NP1)))
C
C     Get the coefficients H_{n_l:n+1,n} and build inner vectors.
C
      IF (IBUILT) THEN
         DO 340 I = NL, NP1
            DWK(IDX(1,I),2*M+IDX(1,N)) = DWK(IDX(1,I),2*M+IDX(1,NP1))
 340     CONTINUE
         DTMP1 = IDTMP1
         DTMP2 = IDTMP2
         DTMP3 = IDTMP3
         DTMP4 = IDTMP4
      ELSE
         DO 350 I = NL, N
            DWK(IDX(1,I),2*M+IDX(1,N)) = DULALH(I,N)
            DTMP     = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),5*M+5) 
     $/ DWK(IDX(1,N),5*M+5)
            CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1))
     $,-DTMP,VECS(1,IDX(2,I)))
            DTMP     = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),5*M+6) 
     $/ DWK(IDX(1,N),5*M+6) * DWK(IDX(1,N),5*M+4) / DWK(IDX(1,I),5*M+4)
            CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1))
     $,-DTMP,VECS(1,IDX(3,I)))
 350     CONTINUE
         DTMP3 = DNRM2(NLEN,VECS(1,IDX(2,NP1)),1)
         DTMP4 = DNRM2(NLEN,VECS(1,IDX(3,NP1)),1)
         DTMP1 = DWK(IDX(1,N),5*M+5) * DTMP3
         DTMP2 = DWK(IDX(1,N),5*M+6) * DTMP4
      END IF
      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DTMP1
      IF (DTMP1.LT.TNRM) IERR = IERR + 16
      IF (DTMP2.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GOTO 360
      DWK(IDX(1,NP1),5*M+4)     = DWK(IDX(1,N),5*M+4) * DTMP1 / DTMP2
      DWK(IDX(1,NP1),IDX(1,NP1)) = DDOT(NLEN,VECS(1,IDX(3,NP1)),1,VECS(1
     $,IDX(2,NP1)),1) / ( DTMP3 * DTMP4 )
      IF ((DTMP3.GE.TMAX).OR.(DTMP3.LE.TMIN)) THEN
         DTMP = DONE / DTMP3
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DTMP,VECS(1,IDX(2,NP1)),
     $DZERO,VECS(1,IDX(2,NP1)))
         DTMP3 = DONE
      END IF
      IF ((DTMP4.GE.TMAX).OR.(DTMP4.LE.TMIN)) THEN
         DTMP = DONE / DTMP4
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DTMP,VECS(1,IDX(3,NP1)),
     $DZERO,VECS(1,IDX(3,NP1)))
         DTMP4 = DONE
      END IF
      DWK(IDX(1,NP1),5*M+5) = DONE / DTMP3
      DWK(IDX(1,NP1),5*M+6) = DONE / DTMP4
C
C     If regular vectors were built, update the counters.
C
 360  IF (.NOT.INNER) THEN
         L  = L + 1
         NL = IWK(IDX(1,L+1),1)
      END IF
C
 370  RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE DULAL2 (NDIM,NLEN,M,N,L,LSTAR,NL,NLSTAR,VF,IERR,
     $                  ADJUST,NORM1,DWK,IWK,IDX,VECS)
C
C     Purpose:
C     This subroutine rebuilds the data for vectors v_{n+1} and w_{n+1}
C     at the point where a block  was forced  to close.  The routine is
C     called internally by the look-ahead Lanczos code and is not meant
C     to be called directly by the user code.
C
C     Parameters:
C        See descriptions in "dulal.doc" or "duqmr.doc".
C
C     External routines used:
C     subroutine daxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
C        LINPACK routine, computes the QR factorization of x.
C     subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
C        LINPACK routine, applies the QR factorization of x.
C
C     Noel M. Nachtigal
C     May 31, 1993
C
C**********************************************************************
C
      EXTERNAL DAXPBY, DQRDC, DQRSL
C
      INTEGER IERR, L, LSTAR, M, N, NDIM, NL, NLEN, NLSTAR
      INTEGER IDX(3,3), IWK(M,4), VF
      DOUBLE PRECISION ADJUST, DWK(M,5*M+7), NORM1, VECS(NDIM,*)
C
C     Miscellaneous parameters.
C
      DOUBLE PRECISION DONE, DZERO
      PARAMETER (DONE = 1.0D0,DZERO = 0.0D0)
C
C     Local variables.
C
      INTEGER I, IJ, J, LBLKSZ, NP1
      DOUBLE PRECISION DTMP1, DTMP2, SCALV, SCALW
C
C     Find the index of the vector pair with the smallest pass value.
C
      IERR = 0
      IF (VF.NE.0) WRITE (VF,'(A20)') 'block did not close:'
      J     = NL
      DTMP1 = DWK(IDX(1,J),5*M+7)
      DO 100 I = NL+1, N
         DTMP2 = DWK(IDX(1,I),5*M+7)
         IF (DTMP2.GT.DZERO) THEN
            IF ((DTMP1.EQ.DZERO).OR.(DTMP2.LT.DTMP1)) THEN
               J     = I
               DTMP1 = DTMP2
            END IF
         END IF
 100  CONTINUE
      IF (DTMP1.EQ.DZERO) THEN
         IF (VF.NE.0) WRITE (VF,'(A47)') 
     $'... no new norm estimates available (aborting).'
         IERR = 8
         RETURN
      END IF
      NORM1 = ADJUST * DTMP1
      IF (VF.NE.0) WRITE (VF,'(A40,I5,E11.4)')
     $      '... updated norms, restarting from step:', IWK(IDX(1,J),3),
     $ NORM1
      IF (IWK(IDX(1,J),3).EQ.N) RETURN
      N  = IWK(IDX(1,J),3)
      L  = IWK(IDX(1,N),4)
      NL = IWK(IDX(1,L+1),1)
C
C     Initialize local variables.
C
      NP1        = N + 1
      DWK(IDX(1,N),5*M+7) = DZERO
      LBLKSZ     = N - NL + 1
C
C     Compute l^\star.
C
      LSTAR  = L - 1
      NLSTAR = IWK(IDX(1,LSTAR+1),1)
C
C     Compute H_{n_l:n,n}.  Save the old coefficients.
C
      IWK(IDX(1,L+1+1),1) = NP1
      IF (LBLKSZ.EQ.1) THEN
         DWK(IDX(1,NL),2*M+IDX(1,NP1)) = DWK(IDX(1,NL),2*M+IDX(1,N))
         DWK(IDX(1,NL),2*M+IDX(1,N))   = DWK(IDX(1,NL),M+IDX(1,N)) / 
     $DWK(IDX(1,NL),IDX(1,NL))
      ELSE
         DO 120 J = NL, N
            DWK(J-NL+1,5*M+1) = DWK(IDX(1,J),M+IDX(1,N))
            DO 110 IJ = NL, N
               DWK(J-NL+1,4*M+IJ-NL+1) = DWK(IDX(1,J),IDX(1,IJ))
 110        CONTINUE
 120     CONTINUE
         CALL DQRDC (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,5*M+3),0,DZERO,0
     $)
         CALL DQRSL (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,5*M+3),DWK(1,5*M
     $+1),
     $               DZERO,DWK(1,5*M+1),DWK(1,5*M+2),DWK(1,5*M+1),DZERO,
     $100,J)
         DO 130 J = NL, N
            DWK(IDX(1,J),2*M+IDX(1,NP1)) = DWK(IDX(1,J),2*M+IDX(1,N))
            DWK(IDX(1,J),2*M+IDX(1,N)) = DWK(J-NL+1,5*M+2)
 130     CONTINUE
      END IF
C
C     Convert inner vectors to regular vectors.
C
      SCALV = DWK(IDX(1,NP1),2*M+IDX(1,N)) * DWK(IDX(1,NP1),5*M+5)
      SCALW = DWK(IDX(1,NP1),2*M+IDX(1,N)) * DWK(IDX(1,NP1),5*M+6) * 
     $DWK(IDX(1,N),5*M+4) / DWK(IDX(1,NP1),5*M+4)
      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DZERO
      DO 140 I = NL, N
         DTMP1 = DWK(IDX(1,I),2*M+IDX(1,N)) - DWK(IDX(1,I),2*M+IDX(1,NP1
     $))
         DTMP2 = DTMP1 * DWK(IDX(1,I),5*M+5) / SCALV
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1)),-
     $DTMP2,VECS(1,IDX(2,I)))
         DTMP2 = DTMP1 * DWK(IDX(1,I),5*M+6) / SCALW * DWK(IDX(1,N),5*M+
     $4) / DWK(IDX(1,I),5*M+4)
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1)),-
     $DTMP2,VECS(1,IDX(3,I)))
 140  CONTINUE
C
      RETURN
      END
C
C**************** 
