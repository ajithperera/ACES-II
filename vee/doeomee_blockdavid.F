      SUBROUTINE DOEOMEE_BLOCKDAVID(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE DRIVES THE CALCULATION OF EXCITATION ENERGIES
C BY THE CC EQUATION OF MOTION METHOD.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION EIGVAL,EIGVAL_T,OSCSTR,R
      CHARACTER*5 HAND(2)
      CHARACTER*1 NATURE(100,8)
      LOGICAL CONVRG
      LOGICAL ESPROP
      LOGICAL NONSTD
      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      LOGICAL EMINFOL,EVECFOL
      INTEGER END,BGN,BGN_IRP,END_IRP
      DIMENSION ICORE(MAXCOR)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/EXTINF/NDIMR,IOLDEST
      COMMON/EXTINF3/IROOT,LOCROOT,ITROOT
      COMMON/RMAT/ R(10000)
      COMMON /CNVRGE/ EMINFOL,EVECFOL
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/STATSYM/IRREPX
      COMMON/GUESS/DOUBLE,NONSTD
      COMMON/CALCINFO/NROOT(8)
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC 
      COMMON/PROPGRAD/ESPROP,IOPTROOT,IOPTSYM
      COMMON/EIGPROB/ISIDE
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
C
      DATA LISTT2 /444/
      DATA HAND /'right',' left'/

#ifdef _DEBUG_VEE
      write(*,*) '@VEE: Entered ', __FILE__
#endif

      IF(CC .AND. (.NOT. CCD) )WRITE(6,6000)
      IF(MBPT2)WRITE(6,6001)
      IF(CCD)WRITE(6,6002)
      IF(RCCD)WRITE(6,6003)
      IF(DRCCD)WRITE(6,6004)
      IF(LCCD)WRITE(6,6005)
      IF(LCCSD)WRITE(6,6006)
      IF(CC2)WRITE(6,6007)

      IF (CCD .OR. RCCD .OR. DRCCD) THEN
         DO ISPIN=1,IUHF+1 
            CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
         END DO
C
         IF (IUHF .EQ. 0) THEN
             CALL ZERO(ICORE,NT(1))
             CALL PUTLST(ICORE,1,1,1,1,190)
         ELSE 
             DO ISPIN = 1, 1+IUHF
                CALL ZERO(ICORE,NT(ISPIN))
                CALL PUTLST(ICORE,1,1,1,ISPIN,190)
             ENDDO
         ENDIF
      ENDIF 
      CORE_WINDOW = 0
      DO IRREP = 1, NIRREP
         CORE_WINDOW=CORE_WINDOW + IWINDOW(IRREP)
      ENDDO
C     
C  put eigval and oscstr to zero
C
      CALL ZERO(EIGVAL,800)
      CALL ZERO(OSCSTR,800)
C
C
C READ IN THE NUMBER OF ROOTS REQUESTED IN EACH SYMMETRY BLOCK
C
      ITROOT=0
      DO 10 IRREPX=1,NIRREP
       TOTITER = 0
       ISIDE=1
       IF(NROOT(IRREPX).EQ.0)GOTO 10
       WRITE(6,7000)HAND(ISIDE)
       CALL ZERO(R,10000)
       WRITE(6,2000)IRREPX,NROOT(IRREPX)
       CALL NEWLST(IRREPX,ICORE,MAXCOR,IUHF)
       IF(NSIZEC.EQ.0)GOTO 10
C
C PUT DIAGONAL PART OF HBAR ON LIST 472 (SINGLES ONLY RIGHT NOW)
C
       CALL HBARDIAG(IRREPX,ICORE,MAXCOR/IINTFP,IUHF)
C
C  DETERMINE EXCITATION PATTERN 
C
       DO IOPT = 1,2
          CALL CALCEXCP(IUHF, ICORE, MAXCOR/IINTFP, IRREPX,
     &                 .TRUE., NONSTD, NSIZEC,IOPT)
          CALL PUTEXCP(IUHF, ICORE, MAXCOR/IINTFP, IRREPX, 
     &                 IPROJECT, LISTH0, ICOLPR1, NSIZEC,IOPT)
       ENDDO

C IFLAGS2(119) = 1 turns on the core EE search. With this option and %EXCITE
C core EEs are obtained in state specic manner. The requested core EE must be
C specified using %EXCITE* namelist in ZMAT.  When %EXCITE is set NONSTD 
C is true. The %EXCITE format is as follows (see EOM_EXCITE file from a 
C small EE calculation. 
C %EXCITE*
C 1  No. of root.
C 1  No. of guess vectors
C 1 1 0 10 0 1.0 (spin-type, I, A, J, B, C, where spin type AA=1,BB=2,
C AAAA=1,BBBB=2,ABAB=3, I=origin,A=destination,j=origin,B=destination.
C The IFLAGS2(120) is EOM_PRJCT and option 4 indicates core.
C
      
      IF (IFLAGS2(119) .EQ. 1) THEN

         IF (IFLAGS2(120) .EQ. 4) THEN

         IF (NONSTD) THEN
            Write(6,*) 
            If (CORE_WINDOW .GT.0) Write(6,"(a,a)") " %EXCITE",
     &         " specification supersedes core-window specification."
             Write(6,*)
             Write(6,"(a,a)") "  Core excitations via core-projection",
     &                        " using %EXCITE*."
             Write(6,*)
                           
             CALL DRIVE_CORE_EE_STATE_SP(ICORE,MAXCOR/IINTFP,IRREPX,
     &                          IUHF,IROOT+1) 
         ELSE 
C 
C The Iflags2(172) is core-window. It can be set to the number of
C core orbitals that need to be included. 
C
             IF (CORE_WINDOW .GT. 0) THEN
             Write(6,*)
             Write(6,"(a,a)") "  Core excitations via",
     &                        " core-projection using a core-window."
             Write(6,*)
                 CALL DRIVE_CORE_EE_CORE_WINDOW(ICORE,MAXCOR/IINTFP,
     &                                          IRREPX,IUHF) 
             ENDIF
         ENDIF 
         ENDIF 

         Write(6,"(a)") " ------Summary of projection schemes------"
         Write(6,*)

         IF (NONSTD) Then

             Write(6,"(a,a)") " The excitation mask is built",
     &                        " for a specific core using %EXCITE*"
             Write(6,"(a,a)") " If core-window is also used its", 
     &                        " specifications are ignored." 
             Write(6,*)
             Write(6,"(a,a)") " Warning! If limited to a one core", 
     &                        " orbital this can lead to wrong"
             Write(6,"(a)")   " core EEs."

         Elseif (Iflags2(172) .GT. 0) Then

             Write(6,"(a,a)") " The excitation mask is built",
     &                        " for a core-window." 
             Write(6,"(a,a)") " If %EXCITE* is also used window will", 
     &                        " narrows to the specified orbital." 
         Else 
             Write(6,"(a,a)") " Marcel Nooijen's original projection",
     &                        " scheme is used. Conceptually this is "
             Write(6,"(a,a)") " identical to core-valence separation",
     &                        "(unpublished work from ~1996)."
             Write(6,*) 
             If (Evecfol) THEN
                 Write(6,"(a,a)") " Root following procedure: overlap",
     &                            " guess vector."
             Else
                 Write(6,"(a,a)") " Root following procedure: Minimum",
     &                            " energy."
             Endif 

         Endif 

         Write(6,*)
         Write(6,"(a)") "------End Summary of projection schemes------"
         Write(6,*)

         DO ITOP = 1,2 
            CALL PUTEXCP(IUHF, ICORE, MAXCOR/IINTFP, IRREPX, 
     &                IPROJECT, LISTH0, ICOLPR1, NSIZEC,ITOP)
         ENDDO 
      ENDIF 
C
C Handle the EOM_PRJCT=NTO; Y. Park and A. Perera, 02/2017,
C
      IF (IFLAGS2(120) .EQ. 5) THEN
         Write(6,"(a)") "Excitations energies via NTO-projection."

          CALL DRIVE_NTO_EE_STATE_PRJCT(ICORE,MAXCOR/IINTFP,IRREPX,
     &                                  IUHF,IROOT+1)
          DO ITOP = 1,2
             CALL PUTEXCP(IUHF, ICORE, MAXCOR/IINTFP, IRREPX, 
     &                   IPROJECT, LISTH0, ICOLPR1, NSIZEC,ITOP)   
          ENDDO 

      ENDIF 

      MAXITER =  IFLAGS2(105)
      TOL=10.0D0**(-IFLAGS(98))

       CALL BLOCK_DAVID_DRIVER(IRREPX,ICORE,MAXCOR,NSIZEC,IUHF,ISIDE,
     &                         NTOL,MAXITER)
C
       WRITE(6,*)
       WRITE(6,"(2a)") " Comments: The message printed above only",
     &                 " state the number of roots requested"
       WRITE(6,"(4a)") " per irrep. and the total number",
     &                 " of iterations took to converge or for",
     &                 " conver."
       WRITE(6,"(2a)") " criteria to fall below 0.01 for those roots",
     &                 " that do not strictly converge."

       WRITE(6,"(3a)") " As a result If there are roots that are not",
     &                 " converged then this count does not"
       WRITE(6,"(2a)") " reflect the actual total number of",
     &                 " iterations!!"
       WRITE(6,*)
10    CONTINUE
C
#ifdef _DEBUG_LVLM
       call modf_amps(Icore, maxcor, iuhf, 400, .true., "R")
#endif
      TOTROOT = 0
      DO I = 1, NIRREP
        TOTROOT = TOTROOT + NROOT(I)
      ENDDO
      IF (TOTROOT .GT. 1) CALL TABLEE(EIGVAL,EIGVAL_T,OSCSTR,NATURE,
     &                                BGN,BGN_IRP,END,END_IRP,
     &                                IFLAGS(87))
C
      RETURN
1007  FORMAT(T2,'Solution of ',I3,' roots required ',I4,
     &       ' iterations.')
1008  FORMAT(T3,'Too many iterations.  Giving up on this symmetry.')
2000  FORMAT(/,T3,'Beginning symmetry block ',I3,'.',I4,
     &       ' roots requested.')
6000  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'EOM-CCSD method.')
6001  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'EOM-MBPT(2) method.')
6002  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'CCD method.')
6003  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'rCCD method.')
6004  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'drCCD method.')
6005  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'LCCD method.')
6006  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'LCCSD method.')
6007  FORMAT(T3,'@DOEOMEE-I, Excitation energies computed by the ',
     &          'CC2 method.')
7000  FORMAT(T3,'@DOEOMEE-I, ',A,'-hand eigenvectors will be computed.')
      END
