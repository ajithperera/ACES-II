         
      SUBROUTINE DFAI(IRREPX
     &                ,FAIA,dfaia)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      CHARACTER*1 SPCASE(2)
      LOGICAL GEOM,FIELD,THIRD,MAGN,STERM,SPIN

      DIMENSION DFAIA(100),UAIA(100),SIJA(100),SABA(100),
     &          EVAL(100),FAIA(100),FIJA(100),FABA(100),SCR(100)

      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/PERT/NTPERT,NPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/PERT2/FIELD,GEOM,THIRD,MAGN,SPIN
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)

      DATA SPCASE /'A','B'/
      DATA ONE,ONEM,HALF,HALFM,TWOM /1.D0,-1.D0,0.5D0,-0.5D0,-2.D0/
      ispin=1
                         
c      DO 10 ISPIN=1,2
         IOFFD=1
         IOFFS=1
         DO 30 IRREPR=1,NIRREP
          IRREPL=DIRPRD(IRREPR,IRREPX)
          NOCCR=POP(IRREPR,ISPIN)
          NVRTL=VRT(IRREPL,ISPIN)
          NOCCL=POP(IRREPL,ISPIN)
          IOFFF=1
          write(*,*)'vrt=',irrepr,vrt(irrepr,1)
         write(*,*) noccr,nvrtl,irrepl,irrepR,irrepx
           do 89 irrep=1,nirrep
           if(irrep .eq. irrepl) write(*,*) irrep,irrepl
c          DO 89 IRREP=1,IRREPL-1
           IRREP1=DIRPRD(IRREP,IRREPX)
           IOFFF=IOFFF+VRT(IRREP,ISPIN)*POP(IRREP,ISPIN)
89        CONTINUE
           write(*,*)'cp' 
          call dcopy(noccr*nvrtl,FAIA(iofff),1,DFAIA(ioffd),1)
          IOFFD=IOFFD+NVRTL*NOCCR
30       CONTINUE
        write(*,*) 'faia'
        call kkk(nt(1),faia)
        call kkk(irpdpd(irrepx,9),dfaia)
       call updmoi(nirrep,irpdpd(irrepx,9),irrepx,183,0,0)
       CALL PUTLST(DFAIA,1,1,1,IRREPX,183)
c10    CONTINUE
C
      RETURN
      END 
