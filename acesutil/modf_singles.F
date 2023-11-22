      SUBROUTINE MODF_SINGLES(ISPIN,ILIST,NT1SIZ,NOCC,NVRT,NSMSZ,
     &                        T,SYVEC,I,A,TYPE,ACT,IUNIT,MODE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER A(NT1SIZ), SYVEC(NSMSZ), POP, VRT, A_IND
      INTEGER NSTART, NIRREP, IRREPS, DIRPRD, NDRPOP(8), NDRVRT(8),
     &        IPOP_OFF(8), IVRT_OFF(8)
      CHARACTER*1 TYPE
      CHARACTER*2 SPCASE(3)
      CHARACTER*4 ACT, MODE
C
      DIMENSION T(NT1SIZ),I(NT1SIZ)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), ID(18)
      COMMON /SYMINF/ NSTART, NIRREP, IRREPS(255,2), DIRPRD(8,8)
C
      DATA SPCASE /'AA','BB','AB'/
C
      IGETI(IX,NUMA)=1+(IX-1)/NUMA
      IGETA(IX,NUMA)=IX-(IGETI(IX,NUMA)-1)*NUMA
      NNM1O2(IX)=(IX*(IX-1))/2
C
      IF (ACT .EQ. "SAVE") Write(IUNIT, "(I3)") NT1SIZ
      IF(NT1SIZ.EQ.0)RETURN
C
      IF (MODE .EQ. "AMPS") THEN
         CALL GETLST(T,1,1,1,ISPIN,ILIST)
      ELSE IF (MODE .EQ. "RESD") THEN
         CALL GETLST(T,1,1,1,ISPIN+2,ILIST)
      ENDIF 
  
      IF (ACT .EQ. "ZERO") THEN
         IF(MODE .EQ. "AMPS") CALL PUTLST(T,1,1,1,ISPIN,ILIST)
         IF(MODE .EQ. "RESD") CALL PUTLST(T,1,1,1,ISPIN+2,ILIST)
         RETURN
      ENDIF
C
#ifdef _DEBUG_LVLM1
      Write(6,*) 
      Write(6, "(a)") "T1 amplitudes"
      Write(6, "(6F13.8)") (T(it1), it1=1, NT1SIZ)
      Write(6,*)
#endif

      DO IRANK=1,NT1SIZ
         IUNPK=SYVEC(IRANK)
         A(IRANK)=IGETA(IUNPK,NVRT)+NOCCO(ISPIN)
         I(IRANK)=IGETI(IUNPK,NVRT)
      ENDDO
C
#ifdef _DEBUG_LVLM1
      Write(6,"(a,a,a)") type,"1 amplitudes for spin case ",
     &                   spcase(ispin)
      write(6,*)
      do irank = 1, nt1siz
         Write(6,"(i3,1x,i3,1x,F13.8)") I(Irank), A(irank), T(irank)
      enddo
      write(6,*)
#endif
C
      CALL GETREC(-1,'JOBARC','NDROPPOP',NIRREP,NDRPOP)
      CALL GETREC(-1,'JOBARC','NDROPVRT',NIRREP,NDRVRT)
      CALL GETREC(-1,'JOBARC','NDROTPOP',1,NDRPOPT)
      CALL GETREC(-1,'JOBARC','NDROTVRT',1,NDRVRTT)
C
#ifdef _DEBUG_LVLM1 
       Write(6,"(a,8I4)") "Occ-Drop", (NDRPOP(k), k=1, Nirrep)
       Write(6,"(a,8I4)") "Vrt-Drop", (NDRVRT(k), k=1, Nirrep)
#endif
      IPOP_OFF(1) = 0
      IVRT_OFF(1) = NDRPOPT
      DO IRREP = 2, NIRREP
         IPOP_OFF(IRREP) = IPOP_OFF(IRREP-1) + NDRPOP(IRREP-1)
         IVRT_OFF(IRREP) = IVRT_OFF(IRREP-1) + NDRVRT(IRREP-1)
      ENDDO

      IRREPX = 1
      Ioff   = 0
C
#ifdef _DEBUG_LVLM1
      Write(6,*) 
      Write(6,"(6x,a,3x,a)") 'i  [i_SYM]  ', 'a  [a_SYM]  '
      Write(6,*)
#endif
c
      DO IRREPR= 1, NIRREP
C
         IRREPL=DIRPRD(IRREPR, IRREPX)
C
         DO NOCC = 1, POP(IRREPR, ISPIN)
C
            DO NVRT = 1, VRT(IRREPL, ISPIN)

               Ioff  = Ioff + 1
               I_IND = I(Ioff) + NDRPOP(IRREPR) + IPOP_OFF(IRREPR)
               A_IND = A(IOff) + IVRT_OFF(IRREPL)
C
               IF (IOFF .LE. NT1SIZ) THEN
C
#ifdef _DEBUG_LVLM1
                   WRITE(6,"(3X, I4, 3X, ' [',I1,'],  ',1X, I4,3X, 
     &                   ' [',I1,']   ;',10X, F13.9)")
     &                      I_IND, IRREPL, A_IND, IRREPR, T(Ioff)
#endif
C
                  IF (ACT .EQ. "SAVE") WRITE(IUNIT, "(3X,I4,3X,I1,
     &                                       3X,I4,3X, I1, 3X, 
     &                                       F15.12)") I_IND, IRREPL,
     &                                       A_IND, IRREPR, T(Ioff)
C
                  IF (ACT .EQ. "READ") CALL T1MATCH(I_IND, IRREPL, 
     &                                              A_IND, IRREPR,
     &                                              T(Ioff))
C
               ENDIF
C
            ENDDO
               
          ENDDO 
      ENDDO
#ifdef _DEBUG_LVLM1
      Write(6,*)
      Write(6, "(a)") "Read T1 amplitudes"
      Write(6, "(6F13.8)") (T(it1), it1=1, NT1SIZ)
      Write(6,*)
#endif
    
      IF (ACT .EQ. "READ") THEN
         IF(MODE .EQ. "AMPS") CALL PUTLST(T,1,1,1,ISPIN,ILIST)
         IF(MODE .EQ. "RESD") CALL PUTLST(T,1,1,1,ISPIN+2,ILIST)
      ENDIF
C
      RETURN
      END
