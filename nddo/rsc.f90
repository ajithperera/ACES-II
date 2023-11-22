      FUNCTION RSC (K,NA,EA,NB,EB,NC,EC,ND,ED)
use tables
!    *
!     CALCULATE THE RADIAL PART OF ONE-CENTER TWO-ELECTRON INTEGRALS
!     (SLATER-CONDON PARAMETER).
!     K     - TYPE OF INTEGRAL, CAN BE EQUAL TO 0,1,2,3,4 IN SPD-BASIS
!     NA,NB - PRINCIPLE QUANTUM NUMBER OF AO, ELECTRON 1
!     EA,EB - EXPONENTS OF AO, ELECTRON 1
!     NC,ND - PRINCIPLE QUANTUM NUMBER OF AO, ELECTRON 2
!     EC,ED - EXPONENTS OF AO, ELECTRON 2
!     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer,intent(in)::k,na,nb,nc,nd
      double precision,intent(in)::ea,eb,ec,ed
!      COMMON
!     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
!     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
!     ./PSC   / F(30),B(30,30)
      DATA A0    /  0.529167D0/
      DATA AFACT / 57.29577951308232D0/
      DATA EV    / 27.21D0/
      DATA EVCAL / 23.061D0/
      DATA PI    /  3.141592653589793D0/
      DATA W1    / 14.399D0/
      DATA W2    /  7.1995D0/
      DATA BIGEXP/ 50.0D0/
            DATA ZERO  /  0.0D0/
      DATA ONE   /  1.0D0/
      DATA TWO   /  2.0D0/
      DATA THREE /  3.0D0/
      DATA FOUR  /  4.0D0/
      DATA PT5   /  0.5D0/
      DATA PT25  /  0.25D0/
      AEA    = DLOG(EA)
      AEB    = DLOG(EB)
      AEC    = DLOG(EC)
      AED    = DLOG(ED)
      NAB    = NA+NB
      NCD    = NC+ND
      ECD    = EC+ED
      EAB    = EA+EB
      E      = ECD+EAB
      N      = NAB+NCD
      AE     = DLOG(E)
      A2     = DLOG(TWO)
      ACD    = DLOG(ECD)
      AAB    = DLOG(EAB)
      C      = DEXP(F(N)+NA*AEA+NB*AEB+NC*AEC+ND*AED &
                  +PT5*(AEA+AEB+AEC+AED)+A2*(N+2) &
                  -PT5*(F(2*NA+1)+F(2*NB+1) &
                  +F(2*NC+1)+F(2*ND+1))-AE*N)
      C      = C*EV
      S0     = ONE/E
      S1     = ZERO
      S2     = ZERO
      M      = NCD-K
      DO 10 I=1,M
      S0     = S0*E/ECD
   10 S1     = S1+S0*(B(NCD-K,I)-B(NCD+K+1,I))/B(N,I)
      M1     = M+1
      M2     = NCD+K+1
      DO 20 I=M1,M2
      S0     = S0*E/ECD
   20 S2     = S2+S0*B(M2,I)/B(N,I)
      S3     = DEXP(AE*N-ACD*M2-AAB*(NAB-K))/B(N,M2)
      RSC    = C*(S1-S2+S3)
      RETURN
      END function rsc
