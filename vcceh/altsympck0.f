C
      SUBROUTINE ALTSYMPCK0(WIN, WOUT, NSIZIN, NSIZOT, ISCR, 
     &                      IRREPX, SPTYPE)
C
C This routine accepts a symmetry packed four-index list and returns
C the same list but with an alternative scheme for symmetry packing.
C
C The List (AI,BJ) is presumed to be packed (Ai,Bj). This routine return
C the list packed (Aj,Bi). Useful for rings and also for the T1 contribution
C to the ring intermediate.  
C
C INPUT: 
C
C WIN    - The symmetry packed Ai-Bj list.
C NSIZIN - The total size of the sym. packed input vector.
C NSIZOT - The total size of the sym. packed output vector.
C SPTYPE - The spin type for the input list
C          'Aaaa' for (Ai,Bj)  (Aj,Bi returned)
C          'Bbbb' for (Ai,Bj)  (Aj,Bi returned)
C          'Abab' for (Ai,Bj)  (Aj,Bi returned)
C          'Baba' for (Ai,Bj)  (Aj,Bi returned)
C          'Abba' for (Ai,Bj)  (Aj,Bi returned (Type Aabb)
C          'Baab' for (Ai,Bj)  (Aj,Bi returned (Type Bbaa)
C          'Aabb' for (Ai,Bj)  (Aj,Bi returned (Type Abba)
C          'Bbaa' for (Ai,Bj)  (Aj,Bi returned (Type Baab)
C OUTPUT: 
C
C  WOUT  - The Symmetry Packed Aj-Bi List.
C       
C SCRATCH:
C
C  ISCR  - Scratch area to hold the symmetry vectors and inverse
C          symmetry vectors which are needed.
C
C          (Size: 2*Nocci*Nvrta If spcase is Aaaa, Bbbb, Abab or
C          Baba;  Nvrta*(Nocci+Noccj)+Nvrtb(Nocci+Noccj) otherwise.
C         
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION WIN(NSIZIN),WOUT(NSIZOT)
      CHARACTER*4 SPTYPE
      DIMENSION ISCR(1),ISPIN(4)
C
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      DO 10 I = 1, 4
         IF(SPTYPE(I:I).EQ.'A') ISPIN(I) = 1
         IF(SPTYPE(I:I).EQ.'B') ISPIN(I) = 2
 10   CONTINUE
C     
      CALL SSTGEN(WIN, WOUT, NSIZIN, VRT(1,ISPIN(1)), POP(1,ISPIN(2)),
     &            VRT(1,ISPIN(3)), POP(1,ISPIN(4)), ISCR, IRREPX, 
     &            '1432')
C
      RETURN
      END
