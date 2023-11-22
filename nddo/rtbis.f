      subroutine RTBIS(out,X1,X2,XACC,L,condon,D)
      implicit double precision (a-h,o-z)
      PARAMETER (JMAX=40)
c
c     this routine was taken from numerical recipes and modified to determine
c     additive terms needed for integral calculations.  
c
c      FMID=FUNC(X2)
      call func(x2,fmid,L,condon,D)
c      F=FUNC(X1)
      call func(x1,f,L,condon,D)
      IF(F*FMID.GE.0.) PAUSE 'Root must be bracketed for bisection.'
      IF(F.LT.0.)THEN
c        RTBIS=X1
         out=x1
       DX=X2-X1
      ELSE
c        RTBIS=X2
         out=x2
       DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
c        XMID=RTBIS+DX
        XMID=out+DX
c        FMID=FUNC(XMID)
        call func(xmid,fmid,L,condon,D)
c        IF(FMID.LE.0.)RTBIS=XMID
        IF(FMID.LE.0.)out=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      PAUSE 'too many bisections'
      END
