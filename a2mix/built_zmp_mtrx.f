      SUBROUTINE BUILT_ZMP_MATRX(ZMPMAT, WORK, ILEFT, NORBS)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION ZMPMAT(NORBS, NORBS), WORK(ILEFT), IBUF(600)
C
      MATDIM = NORBS*NORBS
C
      I000 = 1
      I010 = I000 + matdim
      I020 = I010 + matdim
      I030 = I020 + matdim*matdim
      I040 = I030 + matdim
      I050 = I040 + matdim
      I060 = I050 + matdim
      I070 = I060 + matdim
      I080 = I070 + 600

      IF (I080 .GE. ILEFT) CALL INSMEM("@built_zmp_mtrx", I080, ILEFT)
C
C get overlap and onehao unpacked matrices and 2e integrals
C stored in list format
C ******* Must run with symmetry=off***********
C onehao: I000->I010, ovrlp: I010->I020, 2e ints: I020-I030
C trnfor: I030->I040,tmpmat:I040->I050,pmat:I050->I060
C unitar: I060->I070, buf: I070->I080,ibuf: I080->I090

      Call getint(WORK(I000),WORK(I010),WORK(I020),WORK(I030),
     &            WORK(I040),WORK(I050),WORK(I060),WORK(I070),
     &            ibuf, norbs)
C
C perform ZMP procedure
C onehao: I000->I010, twoint: I020-I030, trnfor: I030->I040,
C norbs,fock: I050->I060,
C dens: I060->I070,tmp_1: I110->I120,tmp_2: I120->I130,
C tmp_3: I010->I020, tmp_4: I040->I050
C
      Call getzmp(WORK(I000),WORK(I020),WORK(I030),
     &            norbs,WORK(I050),WORK(I060),
     &            WORK(I110),WORK(I120),WORK(I010),WORK(I040))

      RETURN
      END


