      

      SUBROUTINE FAOPOP(POP,VRT,AOPOP,NIRREP)
      IMPLICIT INTEGER (A-Z)
      DIMENSION POP(8,2),VRT(8,2),AOPOP(8)
      DO 20 IRREP=1,NIRREP
       AOPOP(IRREP)=POP(IRREP,1)+VRT(IRREP,1)
20    CONTINUE
      RETURN
      END
