C
      SUBROUTINE OPRBAR1R(MIJ,MAB,MAI,MIJBAR,MABBAR,MAIBAR,
     &                    MIABAR,SCR,MAXCOR,IUHF,IRREPX,ISIDE,
     &                    MBAR00,IEOMPROP)
C
C THIS ROUTINE CALCULATES ALL ONE-PARTICLE PARTS OF THE TRANSFORMED
C ONE-PARTICLE OPERATOR:
C
C                _
C                M = exp(-T) M exp(T)
C
C TO WIT:
C    _ 
C    M (ij) = M(ij) + T(je) M(ei) 
C    _
C    M (ab) = M(ab) - T(ma) M(mb)
C    _                                       _             _
C    M (ai) = M(ai) + T(im,ae) M(me) + T(ie) M(ae) - T(ma) M(mi)
C    _                                        
C    M (ia) = M(ia) (If antisymmetric M(ia) = -M(ai))
C
C NOTE THAT THE M(ae) ELEMENTS USED IN THE THIRD EQUATION CORRESPOND
C TO THE *TRANSFORMED* DIPOLE OPERATOR, WHICH IS HOW THE T1*T1 TERMS
C ARE INCLUDED. 
C
      IMPLICIT INTEGER (A-Z)
      LOGICAL PRINT
      DOUBLE PRECISION ONE,ONEM,ZILCH,SNRM2,TWO,SCR,SDOT,
     &                 SIGN1, SIGN2
      DOUBLE PRECISION MIJ(*),MAB(*),MAI(*),MEXPECT,MBAR00
      DOUBLE PRECISION MIJBAR(*),MABBAR(*),MAIBAR(*),MIABAR(*)
      DIMENSION SCR(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
C
      DATA ONE,ONEM,ZILCH,TWO /1.0D0,-1.0D0,0.0D0,2.0D0/
C
      PRINT=IFLAGS(1).GE.10 
C
C 08/94 Ajith. See the comments below.
C Extended to include antisymmetric perturbations and documented
C the changes made to take care the antisymmetry. When is not
C specificaly mentioned contractions are identical for both symmetric
C and anti-symmetric perturbations.
C
      IF (IFLAGS(18) .EQ. 8) THEN
         IANTI = 1
      ELSE
         IANTI = 0
      ENDIF
C
      LISTT1=90+100*(ISIDE-1)
C
      SIGN1 = ONE
      SIGN2 = ONE
      IF (ISIDE .EQ. 1 .AND. IANTI .NE. 0) SIGN1 = ONEM
      IF (ISIDE .EQ. 2 .AND. IANTI .NE. 0) SIGN2 = ONEM
C
C DO OCCUPIED-OCCUPIED PART
C      
C Iside = 1      
C    _ 
C    M (ij) = M(ij) + M(ie) * T(ej) = M(ij) + M(ei) * T(ej)
C    _
C    M (ij) = M(ij) - M(ei) * T(ej) if the perturbation is antisymmetric
C
C Iside = 2
C    _                _                       _
C    L (ij) = M(ji) + M(ei) * L(je) = M(ij) + M(ei) * L(je)
C    _                  _
C    L (ij) = - M(ij) + M(ei) * L(ej) if the perturbation is antisymmetric
C           _                                                         _
C Note that L(ji) which is what actually being evaluated is stored in L(ij).
C           _
C Moreover, L(ij) has no use in the rest of the program. 

      DO 10 ISPIN=1,1+IUHF
       LENOO=IRPDPD(IRREPX,20+ISPIN)
       LENVO=IRPDPD(IRREPX,8+ISPIN)
       I000=1
       CALL GETLST(SCR(I000),1,1,1,ISPIN,LISTT1)
       IOFFT0=I000
       IOFFOO0=1+(ISPIN-1)*IRPDPD(IRREPX,21)
       IOFFVO0=1+(ISPIN-1)*IRPDPD(IRREPX,9)
       DO 30 IRREPJ=1,NIRREP
        IRREPI=DIRPRD(IRREPJ,IRREPX)
        IRREPE=IRREPJ
        NUMJ=POP(IRREPJ,ISPIN)
        NUMI=POP(IRREPI,ISPIN)
        NUME=VRT(IRREPE,ISPIN)
        NROW=NUMI
        NCOL=NUMJ
        NSUM=NUME
        IOFFT=IOFFT0+(ISYMOFF(IRREPJ,1,8+ISPIN)-1)
        IOFFOO=IOFFOO0+(ISYMOFF(IRREPJ,IRREPX,20+ISPIN)-1)
        IOFFVO=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
        CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE*SIGN1,MAI(IOFFVO),
     &              NSUM,SCR(IOFFT),NSUM,ONE*SIGN2,MIJBAR(IOFFOO),NROW)
30     CONTINUE
10    CONTINUE
C
C DO VIRTUAL-VIRTUAL PART
C
C Iside = 1      
C    _
C    M (ab) = M(ba) - M(ma) * T(bm) = M(ab) - M(am) * T(bm) 
C    _
C    M (ab) = - M(ab) - (-) M(am) * T(bm) if the perturbation is antisymmetric
C 
C This is little tricky when the perturbation is not symmetric. The actual
C                           _                          _
C contraction evaluates the M(ba), but it is stored in M(ab) which is consistent
C with rest of the ACESII ordering of (virt,virt) quantities. But bare M is 
C ordered as M(ab) and note the negative sign in front of M(ab) when the
C perturbation is not symmetric.
C
C Iside = 2
C    _                _
C    L (ab) = M(ab) - M(am) * L(mb) 
C    _                _
C    L (ab) = M(ab) - M(am) * L(mb) if the perturbation is antisymmetric

      DO 110 ISPIN=1,1+IUHF
       LENVV=IRPDPD(IRREPX,18+ISPIN)
       LENVO=IRPDPD(IRREPX,8+ISPIN)
       I000=1
       I010=I000+NT(ISPIN)
       CALL GETLST(SCR(I000),1,1,1,ISPIN,LISTT1)
       IOFFT0=I000
       IOFFVV0=1+(ISPIN-1)*IRPDPD(IRREPX,19)
       IOFFVO0=1+(ISPIN-1)*IRPDPD(IRREPX,9)
       DO 130 IRREPB=1,NIRREP
        IRREPA=DIRPRD(IRREPB,IRREPX)
        IRREPM=IRREPB
        NUMA=VRT(IRREPA,ISPIN)
        NUMB=VRT(IRREPB,ISPIN)
        NUMM=POP(IRREPM,ISPIN)
        NROW=NUMA
        NCOL=NUMB
        NSUM=NUMM
        IOFFT=IOFFT0+(ISYMOFF(IRREPM,1,8+ISPIN)-1)
        IOFFVV=IOFFVV0+(ISYMOFF(IRREPB,IRREPX,18+ISPIN)-1)
        IOFFVO=IOFFVO0+(ISYMOFF(IRREPM,IRREPX,8+ISPIN)-1)
        CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM*SIGN1,MAI(IOFFVO),
     &             NROW,SCR(IOFFT),NCOL,ONE*SIGN1,MABBAR(IOFFVV),
     &             NROW)
130    CONTINUE
110   CONTINUE
C
C DO VIRTUAL-OCCUPIED PART
C
C Iside = 1
C    _                                       _             _
C    M (ai) = M(ai) + T(ae,im) M(me) + T(ei) M(ae) - T(am) M(mi)
C
C Iside = 2       
C    _                                       _             _
C    L (ia) = M(ia) + L(im,ae) M(em) + L(ie) M(ea) - L(ma) M(im)
C 
      DO 210 ISPIN=1,1+IUHF
       LENVV=IRPDPD(IRREPX,18+ISPIN)
       LENVO=IRPDPD(IRREPX,8+ISPIN)
       LENOO=IRPDPD(IRREPX,20+ISPIN)
       LENT1=NT(ISPIN)
       LISTT2AB=100*(ISIDE-1)+38-ISPIN
       LISTT2AA=100*(ISIDE-1)+33+ISPIN
       NT2DSZAB=IRPDPD(IRREPX,ISYTYP(1,LISTT2AB))
       NT2DISAB=IRPDPD(IRREPX,ISYTYP(2,LISTT2AB))
       NT2DSZAA=IRPDPD(IRREPX,ISYTYP(1,LISTT2AA))
       NT2DISAA=IRPDPD(IRREPX,ISYTYP(2,LISTT2AA))
       I000=1
       I0GV=I000+LENT1
       I0GO=I0GV+NFEA(ISPIN)
       I010=I0GO+NFMI(ISPIN)
       I020=I010+MAX(NT2DSZAA*NT2DISAA,NT2DSZAB*NT2DSZAB)
       IF(IUHF.EQ.0)I030=I020+NT2DSZAA*NT2DISAA
       CALL GETLST(SCR(I000),1,1,1,ISPIN,LISTT1)
       IF(ISIDE.EQ.2)THEN
        CALL GETLST(SCR(I0GO),1,1,1,ISPIN,191)
        CALL GETLST(SCR(I0GV),1,1,1,ISPIN,192)
       ENDIF
       IOFFT0=I000
       IOFFVV0=1+(ISPIN-1)*IRPDPD(IRREPX,19)
       IOFFOO0=1+(ISPIN-1)*IRPDPD(IRREPX,21)
       IOFFVO0=1+(ISPIN-1)*IRPDPD(IRREPX,9)
       DO 230 IRREPI=1,NIRREP
        IRREPA=DIRPRD(IRREPI,IRREPX)
        NUMA=VRT(IRREPA,ISPIN)
        NUMI=POP(IRREPI,ISPIN)
        NROW=NUMA
        NCOL=NUMI
        IOFFVO=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
C
C Iside = 1
C  
C    RHS(ai) = M(ae)  * T(ei) + M(ae) [Iside = 1]
C
C Iside = 2
C
C     LHS(ai) = M(ea) * L(ie) +  M(ia) = M(ea) * L(ie) + M(ai)
C     LHS(ai) = M(ea) * L(ie) - M(ai) if the perturbation is antisymmetric
C
C Refer to the earlier comments. 
C
        IRREPE=IRREPI
        NUME=VRT(IRREPE,ISPIN)
        NSUM=NUME
        IOFFT=IOFFT0+(ISYMOFF(IRREPI,1,8+ISPIN)-1)
        IOFFVV=IOFFVV0+(ISYMOFF(IRREPE,IRREPX,18+ISPIN)-1)
        CALL XGEMM('N','N',NROW,NCOL,NSUM,ONE,MAB(IOFFVV),
     &             NROW,SCR(IOFFT),NSUM,ONE*SIGN2,MAIBAR(IOFFVO),NROW)
C 
C Iside = 1
C                        _  
C    RHS(ai) = - T(am) * M(mi) + RHS(ai)
C          
C Iside = 2
C                        _                      
C    LHS(ai) = - L(ma) * M(im) + LHS(ai)
C
C  Refer to the earlier comments
C
        IRREPM=IRREPA
        NUMM=POP(IRREPM,ISPIN)
        NSUM=NUMM
        IOFFT=IOFFT0+(ISYMOFF(IRREPM,1,8+ISPIN)-1)
        IF(ISIDE.EQ.1)THEN
         IOFFOO=IOFFOO0+(ISYMOFF(IRREPI,IRREPX,20+ISPIN)-1)
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,SCR(IOFFT),
     &              NROW,MIJBAR(IOFFOO),NSUM,ONE,MAIBAR(IOFFVO),NROW)
        ELSE
         IOFFOO=IOFFOO0+(ISYMOFF(IRREPM,IRREPX,20+ISPIN)-1)
         CALL XGEMM('N','T',NROW,NCOL,NSUM,ONEM,SCR(IOFFT),
     &              NROW,MIJ(IOFFOO),NCOL,ONE,MAIBAR(IOFFVO),NROW)
        ENDIF
C
        IF(ISIDE.EQ.2)THEN
C    
C  LHS(ai) = G(ea) * M(ie) + LHS(ai) [Iside = 2]
C  LHS(ai) = M(ma) * G(im) + LHS(ai) [Iside = 2]
C
C Once again this is tricky for antisymetric perturbations. To do the
C matrix multiplication we need M(ei) or M(am) matrix elements. It looks
C like we need to negate them to take care of the antisymmetry, but in
C the program M(ie) and elements are ordered as M(ei) facilitating  direct
C multiplication.
C
         IRREPE=DIRPRD(IRREPI,IRREPX)
         NUME=VRT(IRREPE,ISPIN)
         IOFFG=I0GV+(ISYMOFF(IRREPA,1,18+ISPIN)-1)
         IOFFM=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
         IOFFTAR=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
         NSUM=NUME
         CALL XGEMM('T','N',NROW,NCOL,NSUM,ONE,SCR(IOFFG),
     &              NSUM,MIABAR(IOFFVO),NSUM,ONE,MAIBAR(IOFFTAR),NROW)
C
C  LHS(ai) = M(ma) * G(im) + LHS(ai) [Iside = 2]
C
         IRREPM=DIRPRD(IRREPA,IRREPX)
         NUMM=POP(IRREPM,ISPIN)
         IOFFG=I0GO+(ISYMOFF(IRREPI,1,20+ISPIN)-1)
         IOFFM=IOFFVO0+(ISYMOFF(IRREPM,IRREPX,8+ISPIN)-1)
         IOFFTAR=IOFFVO0+(ISYMOFF(IRREPI,IRREPX,8+ISPIN)-1)
         NSUM=NUMM
         CALL XGEMM('N','N',NROW,NCOL,NSUM,ONEM,MIABAR(IOFFM),NROW,
     &              SCR(IOFFG),NSUM,ONE,MAIBAR(IOFFTAR),NROW)
        ENDIF 
C
230    CONTINUE
C 
C [Iside = 1]
C
C  RHS(ai) = T(ai,em)*M(me) + RHS(ai) = T(ai,em)*M(em) + RHS(ai)
C
C  RHS(ai) = - T(ai,em)*M(em) + RHS(ai) If the perturbation is antisymmetric
C 
C [Iside = 2]
C                     _
C  LHS(ai) = L(ai,em)*M(em) + LHS(ai) 
C
C If ISIDE = 1 or (ISIDE = 2 and IEOMPROP not equal to 2) include the
C following contributions. 
C
       IF (.NOT. (ISIDE .EQ. 2 .AND. IEOMPROP .EQ. 2)) THEN
C
          IF(ISPIN.EQ.1)THEN
             IOFFVO1=1
             IOFFVO2=1+IUHF*IRPDPD(IRREPX,9)
          ELSE
             IOFFVO1=1+IUHF*IRPDPD(IRREPX,9)
             IOFFVO2=1
          ENDIF
C
          CALL GETLST(SCR(I010),1,NT2DISAB,1,IRREPX,LISTT2AB)
          IF(IUHF.EQ.0)THEN
             CALL GETLST(SCR(I020),1,NT2DISAA,1,IRREPX,
     &                   100*(ISIDE-1)+39)
             CALL SSCAL(NT2DSZAB*NT2DISAB,TWO,SCR(I010),1)
             CALL SAXPY(NT2DSZAB*NT2DISAB,ONEM,SCR(I020),1,SCR(I010),1)
          ENDIF
          CALL XGEMM('N','N',NT2DSZAB,1,NT2DISAB,ONE*SIGN1,SCR(I010),
     &                NT2DSZAB,MAI(IOFFVO2),NT2DISAB,ONE,
     &                MAIBAR(IOFFVO1),NT2DSZAB)
          IF(IUHF.NE.0)THEN
             CALL GETLST(SCR(I010),1,NT2DISAA,1,IRREPX,LISTT2AA)
             CALL XGEMM('N','N',NT2DSZAA,1,NT2DISAA,ONEM*SIGN1,
     &                   SCR(I010),NT2DSZAA,MAI(IOFFVO1),NT2DISAA,ONE,
     &                   MAIBAR(IOFFVO1),NT2DSZAA)
          ENDIF
       ENDIF
C
C For ISIDE = 2 (Left side) and IRREPX = 1 Add MBAR00*L(ai) 
C unlinked piece. We need to play with this so as to setup a another
C flag so that we can tell it when we want to include this contribution. 
C Added on 08/1994 by Ajith. IEOMPROP = 0 include these unlinked
C contributions and the results should be consistent with John's
C implementation. IEOMPROP = 1 do not include the unlinked terms and
C should be consistent with Hideo's modification.
C
       IF (ISIDE .EQ. 2 .AND. IRREPX .EQ. 1 .AND. IEOMPROP .EQ. 0) THEN
          CALL GETLST(SCR(I000),1,1,1,ISPIN,190)
          CALL SAXPY (NT(ISPIN),MBAR00,SCR(I000),1,MAIBAR(IOFFVO0),1)
       ENDIF
C
210   CONTINUE
C
#ifdef _DEBUG_LVL0
       Write(6,*)
       If (iside .eq.1)  call Checksum("oprbar1r-f:",
     +                                  MAIBAR(ioffvo0),Lenvo)

       If (iside .eq.2)  call Checksum("oprbar1l-f:",
     +                                  MAIBAR(ioffvo0),Lenvo)
#endif
      IF(PRINT)THEN
       WRITE(6,500)
500    FORMAT(T3,'@OPRBAR1R-I, Checksum of transformed matrix ',
     &           'elements: ')
       lenoo=irpdpd(irrepx,21)+iuhf*irpdpd(irrepx,22)
       lenvv=irpdpd(irrepx,19)+iuhf*irpdpd(irrepx,20)
       lenvo=irpdpd(irrepx,9)+iuhf*irpdpd(irrepx,10)
       write(6,*)' OO ',sdot(lenoo,mijbar,1,mijbar,1)
       write(6,*)' VV ',sdot(lenvv,mabbar,1,mabbar,1)
       write(6,*)' VO ',sdot(lenvo,maibar,1,maibar,1)
      ENDIF
C
      RETURN
      END
