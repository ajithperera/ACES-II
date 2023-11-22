      subroutine fixfbar(icore, maxcor,iuhf)
C
C   The Lambda intermediates F(ea) and F(mi) are missing the terms
C
C - 1/2 SUM M F(E,M) T(A,M)            
C   1/2 SUM E F(E,M) T(E,I)
C
C Respectively.  This subroutine is here to fix that.
C   
C Finally, the lambda intermediates have the diagonal fock elements subtracted
C  out and stored in the JOBARC file.  They are added back in here.
c 
C  Constructed mostly out of pieces of the FMICONT and FEACONT routines by
C   Renee P Mattie, Feb. 1991
C  Don Comeau for remembering that the diagonal fock elements needed
C     to be added here.  
CEND
C imp.inc should contain the appropriate implicit "undefined", "none" or
c  "logical" statement for this machine
C      include 'imp.inc'
      integer maxcor, iuhf, icore(maxcor)
C
      integer iintln, ifltln, iintfp, ialone, ibitwd
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD  
      integer NT, NF1, NF2, pop, vrt
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),
     $     NF1(2),NF2(2)
      integer nstart, nirrep, irreps, dirprd
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2), DIRPRD(8,8)
      integer irpdpd, isytyp, ntot
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      integer iflags
      COMMON /FLAGS/IFLAGS(100)
C
      integer newcor, iofft1(8,2),
     $     i0za, i0aa, itop, ioffea, ioffem, irrep,nvrt, nocc,
     $     i0zb, i0bb, i0e,  ioffmi,    leneig, 
     $     noe, nve, luja
      double precision one, half, halfm
      parameter ( one = 1.0d0,  half = .5d+00, halfm = -half)
C
      luja = -1
C
C  Load t1's into top of core.  Newcor is the 'new' amount of core remaining 
c   below the t1's
C
      call gett1(icore,maxcor,newcor,iuhf,iofft1)
      if (newcor .lt. 0) then
         write (*,*) 'Error in fixfbar: need at least ',-newcor,
     $        'more icore to call gett1'
         call errex
      endif
C
C   Calculate sizes, starting locations for F(AI), F(MI), and F(EA)
C     F(AI) starts at I0ZA and has the same length as T1AA.
C     F(MI) or F(EA) will start at I0AA and have lengths 
C      NF1(1)   and NF2(1)    respectively
C
      noe = 0
      nve = 0
      do 10, irrep = 1, nirrep
         noe = noe + pop(irrep,1)
         nve = nve + vrt(irrep,1)
 10   continue
      noe    = noe * iintfp
      leneig = noe + nve* iintfp
C
      I0ZA = 1
      I0AA = I0ZA + NT(1)*IINTFP
      I0E  = I0AA + max(NF1(1),NF2(1))*IINTFP
      itop = I0E + leneig
      if (itop .gt. newcor ) then
         write (*,*) 'fixfbar needs at least ', itop - newcor,
     $        ' more icore to construct F(AI)'
         call errex
      endif
C
      CALL GETLST(ICORE(I0ZA),1,1,1,1,93)
      call getrec(luja, 'JOBARC', 'SCFEVALA', leneig, ICORE(I0E) )
C                                   _
C - 1/2 SUM M F(E,M) T(A,M)   => <a|f|e>
C
      CALL GETLST(ICORE(I0AA),1,1,1,1,92)
      IOFFEA=I0AA
      IOFFEM=I0ZA
      DO 100 IRREP=1,NIRREP
         NOCC=POP(IRREP,1)
         NVRT=VRT(IRREP,1)
         if (min(NOCC,NVRT) .ne. 0)
     $        CALL XGEMM('N','T',NVRT,NVRT,NOCC,HALFM,
     $        ICORE(IOFFEM),NVRT,
     $        ICORE(iofft1(irrep,1)),NVRT,ONE,ICORE(IOFFEA),NVRT)
         IOFFEA=IOFFEA+NVRT*NVRT*IINTFP 
         IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
 100  CONTINUE
c      call addfdiag(icore(i0e+noe), icore(i0aa),nirrep,vrt(1,1) )
      CALL PUTLST(ICORE(I0AA),1,1,1,1,92)
C                                 _
C  1/2 SUM E  F(E,M) T(E,J) => <J|f|M>
C
      CALL GETLST(ICORE(I0AA),1,1,1,1,91)
      IOFFMI=I0AA
      IOFFEM=I0ZA
      DO 200 IRREP=1,NIRREP
         NOCC=POP(IRREP,1)
         NVRT=VRT(IRREP,1)
         IF(MIN(NVRT,NOCC).NE.0)THEN
            CALL XGEMM('T','N',NOCC,NOCC,NVRT,HALF,ICORE(IOFFEM),NVRT,
     $           ICORE(iofft1(irrep,1)),NVRT,ONE,ICORE(IOFFMI),NOCC)
         ENDIF
         IOFFMI=IOFFMI+NOCC*NOCC*IINTFP
         IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
 200  CONTINUE
c      call addfdiag(icore(i0e), icore(i0aa),nirrep, pop(1,1) )
      CALL PUTLST(ICORE(I0AA),1,1,1,1,91)
C     
C
C   NOW the beta part
C   See comments above on storage locations and lengths 
C     Note that the same space used for the alpha part is
C     Being reused for the Beta part.

        IF(IUHF.EQ.1) THEN
           noe = 0
           nve = 0
           do 210, irrep = 1, nirrep
              noe = noe + pop(irrep,2)
              nve = nve + vrt(irrep,2)
 210       continue
           noe    = noe * iintfp
           leneig = noe + nve* iintfp
C
           I0ZB = 1
           I0BB = I0ZB + NT(2)*IINTFP
           i0e  = I0BB + max(NF1(2),NF2(2))*IINTFP
           itop = i0e + leneig
           if (itop .gt. newcor ) then
              write (*,*) 'fixfbar needs at least ', itop - newcor,
     $             ' more icore to construct F(ai)'
              call errex
           endif
C
           CALL GETLST(ICORE(I0ZB),1,1,1,2,93)
           call getrec(luja, 'JOBARC', 'SCFEVALA', leneig, ICORE(I0E) )
C                                    _
C - 1/2 SUM M F(e,m) T(a,m)   =>  <a|f|e>
C
           CALL GETLST(ICORE(I0BB),1,1,1,2,92)
           IOFFEA=I0BB
           IOFFEM=I0ZB
           DO 300 IRREP=1,NIRREP
              NOCC=POP(IRREP,2)
              NVRT=VRT(IRREP,2)
              if (min(NOCC,NVRT).NE. 0)
     $             CALL XGEMM('N','T',NVRT,NVRT,NOCC,HALFM,
     $             ICORE(IOFFEM), NVRT,
     $             ICORE(iofft1(irrep,2)),NVRT,ONE,ICORE(IOFFEA),NVRT)
              IOFFEA=IOFFEA+NVRT*NVRT*IINTFP 
              IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
 300       CONTINUE
c           call addfdiag(icore(i0e+noe), icore(i0bb),nirrep,vrt(1,2) )
           CALL PUTLST(ICORE(I0BB),1,1,1,2,92)
C
C                                   _
C  1/2 SUM E  F(e,m) T(e,j)   =. <j|f|m>
C
           CALL GETLST(ICORE(I0BB),1,1,1,2,91)
           IOFFMI=I0BB
           IOFFEM=I0ZB
           DO 400 IRREP=1,NIRREP
              NOCC=POP(IRREP,2)
              NVRT=VRT(IRREP,2)
              IF (MIN(NOCC,NVRT) .NE. 0)
     $             CALL XGEMM('T','N',NOCC,NOCC,NVRT,HALF,
     $             ICORE(IOFFEM),NVRT,
     $             ICORE(iofft1(irrep,2)),NVRT,ONE,ICORE(IOFFMI),NOCC)
              IOFFMI=IOFFMI+NOCC*NOCC*IINTFP
              IOFFEM=IOFFEM+NOCC*NVRT*IINTFP
 400       CONTINUE
c           call addfdiag(icore(i0e), icore(i0bb),nirrep,pop(1,2) )
           CALL PUTLST(ICORE(I0BB),1,1,1,2,91)
           
        ENDIF
C
      return
      end
