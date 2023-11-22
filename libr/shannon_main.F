      Subroutine Shannon_main(Work, Maxcor, Iuhf)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION Work(MAXCOR)
C
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)

      If (Iuhf.NE.0) Then
C
        It1off = 1
        Call Getlst(Work(It1off), 1, 1, 1, 1, 90)
        Nt1absiz = IRPDPD(1,9)
#ifdef _DEBUG_LVLM1
         Write(6,*) Nt1absiz
         Write(6,"(a)") "T1 AA vector"
         Write(6,"(6(1x,F10.6))") (Work(It1off-1+i), i=1, Nt1absiz)
#endif
        It1off = It1off + Nt1absiz
        Call Getlst(Work(It1off), 1, 1, 1, 2, 90)
#ifdef _DEBUG_LVLM1
         Write(6,"(a)") "T1 BB vector"
         Write(6,"(6(1x,F10.6))") (Work(It1off-1+i), i=1,IRPDPD(1,10) )
#endif
        It1off = It1off + Nt1absiz
        Nt1siz = IRPDPD(1,9) + IRPDPD(1,10)

        It2off = 1 + Nt1siz 
        Nt2absiz = Isymsz(1,3)
        Call Getall(Work(it2off), Nt2absiz, 1, 44)

#ifdef _DEBUG_LVLM1
         Write(6,*)
         Write(6,"(a)") "T2 AA vector"
         Write(6,"(6(1x,F10.6))") (Work(It2off-1+i), i=1, Nt2absiz)
#endif 
        It2off = It2off + Nt2absiz
        Nt2absiz = Isymsz(2,4) 
        Call Getall(Work(it2off), Nt2absiz, 1, 45)
#ifdef _DEBUG_LVLM1
         Write(6,*)
         Write(6,"(a)") "T2 BB vector"
         Write(6,"(6(1x,F10.6))") (Work(It2off-1+i), i=1, Nt2absiz)
#endif 
        It2off = It2off + Nt2absiz
        Nt2absiz =  Isymsz(13,14)
        Call Getall(Work(it2off), Nt2absiz, 1, 46)
 
        Nt2siz = Isymsz(1,3) + Isymsz(2,4) + Isymsz(13,14)

        Shannon = Shannon12(Nt1siz, Nt2siz, Work(It1off),
     &                      Work(It2off))

        Write(6,*)
        Write(6,"(a,F13.7)") " Shannon Index =", Shannon
        Write(6,*)
         
      Else

         It1off = 1
         Nt1siz = IRPDPD(1,9)
         Call Getlst(Work(It1off), 1, 1, 1, 1, 90)
#ifdef _DEBUG_LVLM1
         Write(6,"(a)") "T1 vector"
         Write(6,"(F10.6)") (Work(It1off-1+i), i=1, Nt1siz)
#endif

         It2off = It1off + Nt1siz

         Nt2siz = Isymsz(13,14)
         Call Getall(Work(it2off), Nt2siz, 1, 46)
#ifdef _DEBUG_LVLM1
         Write(6,*)
         Write(6,"(a)") "T2 vector"
         Write(6,"(F10.6)") (Work(It2off-1+i), i=1, Nt2siz)
#endif 

         Shannon= Shannon12(Nt1siz, Nt2siz, Work(It1off),
     &                     Work(It2off))
         Write(6,*)
         Write(6,"(a,F13.7)") " Shannon Index =", Shannon
         Write(6,*) 

      Endif

      RETURN
      END

    
