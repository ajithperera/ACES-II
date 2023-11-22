


      SUBROUTINE PUTBLK(TARGET,APACK,TYPE,SCR,ISYM,NBAS,PCKLEN,
     &                  FULLEN,ISPIN)
C
C THIS ROUTINE ACCEPTS A PARTICULAR BLOCK (VIRTUAL-VIRTUAL,
C  OCCUPIED-OCCUPIED, VIRTUAL-OCCUPIED OR OCCUPIED-VIRTUAL)
C  FOR A PARTICULAR SPIN CASE (ISPIN=1 OR ISPIN=2) OF A TWO
C  INDEX QUANTITY WHICH IS SYMMETRY PACKED A,B I,J A,I OR I,A
C  AND STUFFS IT INTO THE FULL UNPACKED MATRIX.
C
C INPUT:
C       APACK - THE SYMMETRY PACKED BLOCK TO BE PLACED IN THE
C                TARGET MATRIX. (LENGTH: NBAS*NBAS)
C       TYPE  - 'VV' FOR VIRT-VIRT, 'OO' FOR OCC-OCC, 'VO' FOR
C                VIRT-OCC, 'OV' FOR OCC-VIRT.
C       NBAS  - THE NUMBER OF BASIS FUNCTIONS (DIMENSION OF TARGET).
C       PCKLEN- THE LENGTH OF APACK.
C       FULLEN- LENGTH OF BLOCK IN UNPACKED FORM
C                (NVRT*NVRT FOR 'VV', NOCC*NOCC FOR 'OO', ETC.)
C
C SCRATCH:
C       SCR   - USED TO HOLD THE UNPACKED BLOCK. (LENGTH: FULLEN)
C       ISYM  - USED TO HOLD THE SYMMETRY VECTORS (LENGTH: FULLEN)
C
C OUTPUT:
C       TARGET- THE FULL UNPACKED TWO-INDEX MATRIX.
C
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*2 TYPE
      CHARACTER*8 LABEL
      DOUBLE PRECISION TARGET(NBAS*NBAS),APACK(PCKLEN)
      DOUBLE PRECISION SCR(1)
      DIMENSION ISYM(FULLEN)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /dropgeo/ ndrgeo
      CALL ZERO(SCR,FULLEN)
      IF(TYPE.EQ.'OO')THEN
       IF(ISPIN.EQ.1)LABEL='SOAOA2  '
       IF(ISPIN.EQ.2)LABEL='SOBOB2  '
       IROWOFF=1
       ICOLOFF=1
       NROW=NOCCO(ISPIN)
       NCOL=NOCCO(ISPIN)
      ELSEIF(TYPE.EQ.'VV')THEN
       IF(ISPIN.EQ.1)LABEL='SVAVA2  '
       IF(ISPIN.EQ.2)LABEL='SVBVB2  '
       IROWOFF=NOCCO(ISPIN)+1
       ICOLOFF=NOCCO(ISPIN)+1
       NROW=NVRTO(ISPIN)
       NCOL=NVRTO(ISPIN)
      ELSEIF(TYPE.EQ.'OV')THEN
       IF(ISPIN.EQ.1)LABEL='SOAVA2  '
       IF(ISPIN.EQ.2)LABEL='SOBVB2  '
       IROWOFF=1
       ICOLOFF=NOCCO(ISPIN)+1
       NROW=NOCCO(ISPIN)
       NCOL=NVRTO(ISPIN)
      ELSEIF(TYPE.EQ.'VO')THEN
       IF(ISPIN.EQ.1)LABEL='SVAOA2  '
       IF(ISPIN.EQ.2)LABEL='SVBOB2  '
       IROWOFF=NOCCO(ISPIN)+1
       ICOLOFF=1
       NROW=NVRTO(ISPIN)
       NCOL=NOCCO(ISPIN)
      ENDIF
c------------------------------------------------------
      if (ndrgeo.eq.0) then
        CALL GETREC(20,'JOBARC',LABEL,FULLEN,ISYM)
      else
       call aces_ja_fin
       istate = ishell('mv JOBARC JOBARC_DM')
       istate = ishell('mv JAINDX JAINDX_DM')
       istate = ishell('mv JOBARC_AM JOBARC')
       istate = ishell('mv JAINDX_AM JAINDX')
       call aces_ja_init
        CALL GETREC(20,'JOBARC',LABEL,FULLEN,ISYM)
       call aces_ja_fin
       istate = ishell('mv JOBARC JOBARC_AM')
       istate = ishell('mv JAINDX JAINDX_AM')
       istate = ishell('mv JOBARC_DM JOBARC')
       istate = ishell('mv JAINDX_DM JAINDX')
       call aces_ja_init
      endif
      CALL SCATTER(PCKLEN,SCR,ISYM,APACK)
      CALL BLKCPY(SCR,NROW,NCOL,TARGET,NBAS,NBAS,IROWOFF,
     &            ICOLOFF)
      RETURN
      END
