      SUBROUTINE CHECK_AB2AA(S2AB,N,WORK,MAXCOR,ISPIN)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION WORK(MAXCOR),S2AB(N)

      INTEGER DISSYS,XIRREP,DIRPRD,VRT,SIRREP,POP

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SINFO/ NS(8), SIRREP
      COMMON/SLISTS/LS1IN, LS1OUT, LS2IN(2,2), LS2OUT(2,2)

      LISTS2EX = LS2OUT(ISPIN,ISPIN+1)

      IOFF = 1 + POP(SIRREP,ISPIN)

      DO XIRREP = 1, NIRREP

         IRREPA = DIRPRD(XIRREP,SIRREP)

         NROW=IRPDPD(XIRREP,20+ISPIN)
         NCOL=VRT(IRREPA,ISPIN) * NS(SIRREP)
         NDIM= IRPDPD(XIRREP,2+ISPIN)

         I000 = 1
         I010 = I000 + NROW*NCOL
         IEND = I010 + NDIM*NCOL
     
          CALL DCOPY(N,S2AB(IOFF),1,WORK(I000),1)
          CALL ASSYM2A(XIRREP,POP(1,1),NROW,NCOL,WORK(I000),
     +                WORK(I010),WORK(I010))
         
          CALL SQSYM(XIRREP,POP(1,1),NDIM,NROW,NCOL,WORK(I010),
     +               WORK(I000))
#ifdef _DEBUG_LVL0
          write(6,*) 
          write(6,"(a)") "S2(ABAB) --> S2(AAAA)"
          Call output(Work(I010),1,Ndim,1,Ncol,Ndim,Ncol,1)
#endif 
      ENDDO 

      RETURN
      END 
