
#include "disttype.h"

      SUBROUTINE NEWLST(IRREPX,ICORE,MAXCOR,IUHF)
C
C RESETS LISTS FOR A PARTICULAR IRREP
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL ESPROP
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /PROPGRAD/ ESPROP,IOPTROOT,IOPTSYM
      COMMON /EXTRAP/ MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON /FLAGS/ IFLAGS(100)
C
C SINGLES LISTS
C
      DO 5 ISPIN=1,1+IUHF
       T1SIZE=IRPDPD(IRREPX,8+ISPIN)
       ABSIZE=IRPDPD(IRREPX,18+ISPIN)
       IJSIZE=IRPDPD(IRREPX,20+ISPIN)
       call aces_list_resize(ISPIN,  490,      T1SIZE)
       call aces_list_resize(ISPIN+2,490,      T1SIZE)
       call aces_list_resize(9,      447+ISPIN,T1SIZE)
       call aces_list_resize(ISPIN,  491,      IJSIZE)
       call aces_list_resize(ISPIN+2,491,      IJSIZE)
       call aces_list_resize(ISPIN,  492,      ABSIZE)
       call aces_list_resize(ISPIN+2,492,      ABSIZE)
       call aces_list_resize(ISPIN,  493,      T1SIZE)
       call aces_list_resize(ISPIN+2,493,      T1SIZE)
5     CONTINUE
C
C LISTS USED BY DAVIDSON EXTRAPOLATOR
C
      LEN=0
      DO 10 ISPIN=1,1+IUHF
       LEN=LEN+IRPDPD(IRREPX,8+ISPIN)
10    CONTINUE
      LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,446),ISYTYP(2,446))
      IF(IUHF.NE.0)THEN
       LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,444),ISYTYP(2,444))
       LEN=LEN+IDSYMSZ(IRREPX,ISYTYP(1,445),ISYTYP(2,445))
      ENDIF
      call aces_list_resize(1,470,LEN)
      call aces_list_resize(1,471,LEN)
      call aces_list_resize(1,472,LEN)
c     call aces_list_resize(1,473,LEN)
      call aces_list_resize(2,470,LEN)
      call aces_list_resize(2,471,LEN)
      call aces_list_resize(2,472,LEN)
c     call aces_list_resize(2,473,LEN)
      NSIZEC=LEN
C
      DO 11 ISPIN=3,3-2*IUHF,-1
C
C T2-SIZE LISTS
C
       TTYPEL=ISYTYP(1,43+ISPIN)
       TTYPER=ISYTYP(2,43+ISPIN)
       CALL NEWTYP2(IRREPX,443+ISPIN,TTYPEL,TTYPER,.TRUE.)
       CALL NEWTYP2(IRREPX,447+ISPIN,TTYPEL,TTYPER,.TRUE.)
       CALL NEWTYP2(IRREPX,460+ISPIN,TTYPEL,TTYPER,.TRUE.)
       CALL NEWTYP2(IRREPX,453+ISPIN,TTYPEL,TTYPER,.TRUE.)
       IF(IFLAGS(93).EQ.2)THEN
        TTYPEL=15
        CALL NEWTYP2(IRREPX,413+ISPIN,TTYPEL,TTYPER,.TRUE.)
        CALL NEWTYP2(IRREPX,213+ISPIN,TTYPEL,TTYPER,.TRUE.)
        CALL NEWTYP2(IRREPX,463+ISPIN,TTYPEL,TTYPER,.TRUE.)
       ENDIF
C
11    CONTINUE
C
C RESORTED T2 LISTS - reduce to those needed for rhf later
C
c YAU : old
c      DO 20 ILIST=434,443
c       TTYPEL=ISYTYP(1,ILIST)
c       TTYPER=ISYTYP(2,ILIST)
c       CALL NEWTYP2(IRREPX,ILIST,TTYPEL,TTYPER,.TRUE.)
c20    CONTINUE
c YAU : new
      CALL NEWTYP2(IRREPX,434,h_DistType_AI,h_DistType_AI,.TRUE.)
      CALL NEWTYP2(IRREPX,435,h_DistType_ai,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,436,h_DistType_ai,h_DistType_AI,.TRUE.)
      CALL NEWTYP2(IRREPX,437,h_DistType_AI,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,438,h_DistType_aI,h_DistType_Ai,.TRUE.)
      CALL NEWTYP2(IRREPX,439,h_DistType_Ai,h_DistType_aI,.TRUE.)
      CALL NEWTYP2(IRREPX,440,h_DistType_AI,h_DistType_AI,.TRUE.)
      CALL NEWTYP2(IRREPX,441,h_DistType_ai,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,442,h_DistType_AI,h_DistType_ai,.TRUE.)
      CALL NEWTYP2(IRREPX,443,h_DistType_Ai,h_DistType_aI,.TRUE.)
c YAU : end
C
C TWO-PARTICLE DIPOLE MOMENT LISTS
C
C _
C M(Ab,Ij)
C
      DO 30 ISPIN=3,3-2*IUHF,-1
       ISIZE=IDSYMSZ(IRREPX,ISYTYP(1,13+ISPIN),ISYTYP(2,13+ISPIN))
       call aces_list_resize(ISPIN,482,ISIZE)
30    CONTINUE
      RETURN
      END
