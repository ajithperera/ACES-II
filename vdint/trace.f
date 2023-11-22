      SUBROUTINE TRACE
C
C Written 4-Dec-1983 hjaaj
C Revised 12-Jan-1984 hjaaj (use LIB$SIGNAL for VAX,
C                            traceback without exit)
C
      DATA A/1.0/,B/0.0/
      C = A/B
      RETURN
      END
