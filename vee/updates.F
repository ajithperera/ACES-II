      SUBROUTINE UPDATES(IRREPX,S,LIST0,IOFFSIN,LIST1,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(*)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      I0SA=1
      IF(IUHF.NE.0)THEN
       I0SB=I0SA+IRPDPD(IRREPX,9)
       I0SAB=I0SB+IRPDPD(IRREPX,10)
       I0SBB=I0SAB+IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
       I0SAA=I0SBB+IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
      ELSE
       I0SAB=I0SA+IRPDPD(IRREPX,9)
      ENDIF
C
      CALL PUTLST(S(I0SA),1,1,1,1+IOFFSIN,LIST1)
      CALL PUTALL(S(I0SAB),1,IRREPX,LIST0+2)
      IF(IUHF.NE.0)THEN
       CALL PUTLST(S(I0SB),1,1,1,2+IOFFSIN,LIST1)
       CALL PUTALL(S(I0SAA),1,IRREPX,LIST0)
       CALL PUTALL(S(I0SBB),1,IRREPX,LIST0+1)
      ENDIF

#ifdef _DEBUG_LVL0
      Write(6,*) "Checksums @update"
      If (list0 .Eq. 461) Then 
      l1a= irpdpd(irrepx,9)
      l2ab=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
      call getlst(S(i0sa),1,1,1,1+ioffsin,490)
      call checksum("List-490",S(i0sa),l1a,B)
      call getall(s(i0sab),l2ab,irrepx,463)
      call checksum("List-463",S(i0sab),l2ab,B)
      call checksum("List-472",S,l2ab+l1a,B)
      If (Iuhf.Ne.0) Then
         l1b=irpdpd(irrepx,10)
         l2aa=IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
         l2bb=IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         l2ab=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         call getlst(S(i0sb),1,1,1,2+ioffsin,490)
         call checksum("List-490",S(i0sb),l1b,B)
         call getall(s(i0saa),l2aa,irrepx,461)
         call checksum("List-461",S(i0saa),l2aa,B)
         call getall(s(i0sbb),l2bb,irrepx,462)
         call checksum("List-462",S(i0sbb),l2bb,B)
         call getall(s(i0sab),l2ab,irrepx,463)
         call checksum("List-463",S(i0sab),l2ab,B)
         call checksum("List-472",S,l1a+l1b+l2aa+l2bb+l2ab,B)
      Endif 
      Elseif(list0 .eq.444) Then
      l1a= irpdpd(irrepx,9)
      l2ab=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
      call getlst(S(i0sa),1,1,1,1+ioffsin,490)
      call checksum("List-490",S(i0sa),l1a,B)
      call getall(s(i0sab),l2ab,irrepx,446)
      call checksum("List-446",S(i0sab),l2ab,B)
      call checksum("List-472",S,l2ab+l1a,B)
      If (Iuhf.Ne.0) Then
         l1b=irpdpd(irrepx,10)
         l2aa=IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
         l2bb=IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         l2ab=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         call getlst(S(i0sb),1,1,1,2+ioffsin,490)
         call checksum("List-490",S(i0sb),l1b,B)
         call getall(s(i0saa),l2aa,irrepx,444)
         call checksum("List-444",S(i0saa),l2aa,B)
         call getall(s(i0sbb),l2bb,irrepx,445)
         call checksum("List-445",S(i0sbb),l2bb,B)
         call getall(s(i0sab),l2ab,irrepx,446)
         call checksum("List-446",S(i0sab),l2ab,B)
         call checksum("List-472",S,l1a+l1b+l2aa+l2bb+l2ab,B)
      Endif
      Endif 
#endif 
C
      RETURN
      END
