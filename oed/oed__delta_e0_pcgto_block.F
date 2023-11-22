C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
         SUBROUTINE  OED__DELTA_E0_PCGTO_BLOCK
     +
     +                    ( LMNVAL,NBATCH,NINT1D,
     +                      ATOMIC,
     +                      MIJ,NIJ,NIJBEG,NIJEND,
     +                      NPGTOA,NPGTOB,
     +                      NXYZA,NXYZB,NXYZET,
     +                      SHELLA,SHELLB,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      XN,YN,ZN,
     +                      ABX,ABY,ABZ,
     +                      ALPHAA,ALPHAB,
     +                      PRIMA,PRIMB,
     +                      NORMA,NORMB,
     +                      RHOAB,
     +                      PAX,PAY,PAZ,PINVHF,SCALE,
     +                      INT1DX,INT1DY,INT1DZ,
     +                                BATCH )
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_E0_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL_1D_INTEGRALS
C                OED__OVL_INT1D_TO_E0
C  DESCRIPTION : This operation calculates a batch of unnormed delta
C                function integrals between primitive cartesian gaussians
C                for the shells A and B
C                and the block of ij exponent pairs. The total number
C                of delta integrasl  generated here is thus given by
C                the total number of cartesian monomials NXYZET times
C                the total number of exponent pairs MIJ in the present
C                block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,nxyzet)
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    integral batch
C                    NINT1D       =  space needed for each of the 1D
C                                    X,Y,Z integral arrays
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    NXYZx        =  # of cartesian monomials for csh x
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A and P=A+B
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B
C                    ABm          =  the m=x,y,z-coordinate differences
C                                    between centers A and B
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
C                    PRIMx        =  i,j labels of primitives for the
C                                    respective contraction shells
C                                    x = A,B
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B
C                    RHOAB        =  the complete set of NIJ exponential
C                                    prefactors between contraction
C                                    shells A and B
C                    PAx          =  will hold current MIJ coordinate
C                                    x=X,Y,Z differences P-A between
C                                    centers P and A
C                    PINVHF       =  will hold current MIJ values of
C                                    1/(2*P), where P are the exponent
C                                    sums for contraction shells A
C                                    and B
C                    SCALE        =  will hold current MIJ values of
C                                    scaling factors
C                    INT1Dx       =  will hold all current 1D integrals
C                                    for each cartesian component
C                                    (x = X,Y,Z)
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian overlap [E|0] integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC

         INTEGER     I,J,M
         INTEGER     IJ
         INTEGER     MIJ
         INTEGER     NBATCH,NINT1D
         INTEGER     NIJ,NIJBEG,NIJEND
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NXYZET
         INTEGER     NXYZA,NXYZB
         INTEGER     SHELLA,SHELLB

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  EXPA,EXPB
         DOUBLE PRECISION  PINV,PVAL
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  ZERO,HALF,ONE,ONEP5

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)

         DOUBLE PRECISION  RHOAB   (1:NIJ)

         DOUBLE PRECISION  BATCH   (1:NIJ,1:NXYZET)

         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  SCALE   (1:MIJ)

         DOUBLE PRECISION  INT1DX  (1:NINT1D)
         DOUBLE PRECISION  INT1DY  (1:NINT1D)
         DOUBLE PRECISION  INT1DZ  (1:NINT1D)

         INTEGER NUCLEI, LMNVAL(3,84)
         INTEGER LA, MA, NA, LB, MB, NB, NAB, NXYZ, COMPA, COMPB

         double precision  Xn, Yn, Zn
         DOUBLE PRECISION ANX, ANY, ANZ, BNX, BNY, BNZ, DIST_NA
         DOUBLE PRECISION DIST_NB, PREFACA, PREFACB, EXPONTA, EXPONTB
         DOUBLE PRECISION DEXP, PI3_2, FACT
C
         PARAMETER  (ZERO  = 0.D0)
         PARAMETER  (HALF  = 0.5D0)
         PARAMETER  (ONE   = 1.D0)
         PARAMETER  (ONEP5 = 1.5D0)
         PARAMETER  (PI3_2 = 5.568327996831707D0)
C
C------------------------------------------------------------------------
C
C             ...calculate the quantities needed to establish the
C                1D overlap integrals.
C
C     
#ifdef _DEBUG_LVL0
      Write(6,*) "OED__DELTA_E0_PCGTO_BLOCK"
      Write(6,*) "The shell and nuclear coordiantes"
      Write(6,*)
      Write(6,"(3(1x,F12.6))") xA, yA, zA
      Write(6,"(3(1x,F12.6))") xB, yB, zB
      Write(6,"(3(1x,F12.6))") xn, yn, zn
      Write(6,*) "The shell angular momentums"
      Write(6,"(2(1x,I3))") SHELLA, SHELLB
      Write(6,*) "The number of primitives per shell"
      Write(6,"(2(1x,I3))") NPGTOA, NPGTOB
      Write(6,*) "The number of monomials per shell and total ints."
      Write(6,"(3(1x,I3))") NXYZA, NXYZB, NXYZET
      Write(6,*) "The exponents A and B"
      Write(6,"(4(1x,F15.7))") (ALPHAA(i),i=1,NPGTOA)
      Write(6,"(4(1x,F15.7))") (ALPHAB(i),i=1,NPGTOB)
      Write(6,*)
#endif

        IF (SHELLA  .GT. 7 .OR. SHELLB .GT. 7) ThEN
           Write(*,"(a,a)")"There is a internal limit of angular",
     &                   " momentum 7 for the Delta function integral!"
           Write(*,*) "OED__OVL_E0_PCGTO_BLOCK"
           STOP
        Endif
          
        ANX = (XN - XA)
        ANY = (YN - YA)
        ANZ = (ZN - ZA)

        BNX = (XN - XA)
        BNY = (YN - YA)
        BNZ = (ZN - ZA)

        NXYZ = 0
        NAB  = 0
     
C The special cases, atomic and/or shell angular momentum is zero

        IF (SHELLA .EQ. 0 .AND. SHELLB .EQ. 0) THEN

           IF (ATOMIC) THEN

                  NXYZ = 1

                  DO NA = 1, NPGTOA
                     DO NB = 1, NPGTOB

                        NAB = NAB + 1

                        DIST_NA = ANX**2 + ANY**2 + ANZ**2
                        DIST_NB = BNX**2 + BNY**2 + BNZ**2

                        EXPONTA = DEXP(-ONE*EXPA*DIST_NA)
                        EXPONTB = DEXP(-ONE*EXPB*DIST_NB)

                        FACT = NORMA(NA)*NORMB(NB)/PI3_2

                        BATCH(NAB, NXYZ) = FACT*EXPONTA*EXPONTB

                      ENDDO
                  ENDDO
           ELSE

               DO NA = 1, NPGTOA
                  DO NB = 1, NPGTOB
  
                     NAB = NAB + 1

                     FACT = NORMA(NA)*NORMB(NB)/PI3_2
                     EXPA = ALPHAA(NA)
                     EXPB = ALPHAB(NB)

                     DIST_NA = ANX**2 + ANY**2 + ANZ**2
                     DIST_NB = BNX**2 + BNY**2 + BNZ**2

                     EXPONTA = DEXP(-ONE*EXPA*DIST_NA)
                     EXPONTB = DEXP(-ONE*EXPB*DIST_NB)

                     FACT   = NORMA(NA)*NORMB(NB)/PI3_2

                     BATCH(NAB, NXYZ) = EXPONTA*EXPONTB*FACT
                  ENDDO
               ENDDO
           ENDIF

        ELSE IF (SHELLA .EQ. 0) THEN 

            IF (ATOMIC) THEN

               DO COMPB = 1, NXYZB

                  NXYZ = NXYZ + 1

                  DO NA = 1, NPGTOA
                     DO NB = 1, NPGTOB

                        NAB = NAB + 1

                        FACT = NORMA(NA)*NORMB(NB)/PI3_2

                        BATCH(NAB, NXYZ) = FACT
                      
                      ENDDO
                  ENDDO
               ENDDO
            ELSE
               DO COMPB = 1, NXYZB

                  NXYZ = NXYZ + 1

                  DO NA = 1, NPGTOA
                     DO NB = 1, NPGTOB

                        NAB = NAB + 1
                     
                         FACT = NORMA(NA)*NORMB(NB)/PI3_2
                         EXPA = ALPHAA(NA)
                         EXPB = ALPHAB(NB)

                         DIST_NA = ANX**2 + ANY**2 + ANZ**2
                         DIST_NB = BNX**2 + BNY**2 + BNZ**2

                         EXPONTA = DEXP(-ONE*EXPA*DIST_NA)
                         EXPONTB = DEXP(-ONE*EXPB*DIST_NB)

                         BATCH(NAB, NXYZ) = EXPONTA*EXPONTB*FACT
                      
                      ENDDO
                  ENDDO
               ENDDO
            ENDIF

      ELSE IF (SHELLB .EQ. 0) THEN
     
            IF (ATOMIC) THEN

               DO COMPA = 1, NXYZB

                  NXYZ = NXYZ + 1

                  DO NA = 1, NPGTOA
                     DO NB = 1, NPGTOB

                        NAB = NAB + 1

                        FACT = NORMA(NA)*NORMB(NB)/PI3_2

                        BATCH(NAB, NXYZ) = FACT

                      ENDDO
                  ENDDO
               ENDDO

            ELSE

               DO COMPB = 1, NXYZB

                  NXYZ = NXYZ + 1

                  LB = LMNVAL(1, COMPB+1)
                  MB = LMNVAL(2, COMPB+1)
                  NB = LMNVAL(3, COMPB+1)

                  DO NA = 1, NPGTOA
                     DO NB = 1, NPGTOB

                        NAB = NAB + 1
 
                         FACT = NORMA(NA)*NORMB(NB)/PI3_2

                         EXPA = ALPHAA(NA)
                         EXPB = ALPHAB(NB)

                         DIST_NA = ANX**2 + ANY**2 + ANZ**2
                         DIST_NB = BNX**2 + BNY**2 + BNZ**2

                         PREFACA = ANX**LA * ANY**MA * ANZ**NA

                         EXPONTA = DEXP(-ONE*EXPA*DIST_NA)
                         EXPONTB = DEXP(-ONE*EXPB*DIST_NB)

                         BATCH(NAB, NXYZ) = PREFACA*
     &                                      EXPONTA*EXPONTB*FACT
                     ENDDO
                  ENDDO
               ENDDO
           ENDIF
           RETURN
        ENDIF 

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)")" The uncontracted <mu|Delta(r-Rn)|nu> shell and",
     &                 " primitive pair."
      Write(6,"(a,5(1x,I3))")" NAB, NIJ, NXYZ{A,B}, NXYZET :", NAB,NXYZ,
     &                         NXYZA, NXYZB, NXYZET
      Do i=1, NAB
         Write(6,"(1x,a,i4)") "The prim. pair: ", i
         write(6,"(6(1x,F20.10))") (BATCH(i,j), j=1, NXYZ)
      Enddo
#endif
       
C General Cases

        DO COMPA = 1, NXYZA
           
           LA = LMNVAL(1, COMPA)
           MA = LMNVAL(2, COMPA)
           NA = LMNVAL(3, COMPA)

           DO COMPB = 1, NXYZB
              
              LB = LMNVAL(1, COMPB)
              MB = LMNVAL(2, COMPB)
              NB = LMNVAL(3, COMPB)

              NXYZ = NXYZ + 1 
              
              DO NA = 1, NPGTOA
                 DO NB = 1, NPGTOB
               
                     NAB = NAB + 1

                     FACT = NORMA(NA)*NORMB(NB)/PI3_2
                     EXPA = ALPHAA(NA)
                     EXPB = ALPHAB(NB)

                     DIST_NA = ANX**2 + ANY**2 + ANZ**2  
                     DIST_NB = BNX**2 + BNY**2 + BNZ**2
 
                     PREFACA = ANX**LA * ANY**MA * ANZ**NA
                     PREFACB = ANX**LB * ANY**MB * ANZ**NB
             
                     EXPONTA = DEXP(-ONE*EXPA*DIST_NA)
                     EXPONTB = DEXP(-ONE*EXPB*DIST_NB)

                     BATCH(NAB, NXYZ) = PREFACA*PREFACB*
     &                                  EXPONTA*EXPONTB*FACT
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "Intermediates to accumulat <mu|Delta(r-Rn)|nu>"
      Write(6,"(a,a)") "EXPA, EXPB, PREFACA, PREFACB, EXPONTA, EXPONTB"
      Write(6,"(6(1x,F12.6))") EXPA, EXPB, DIST_NA, DIST_NB, 
     &                         EXPONTA, EXPONTB
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

#ifdef _DEBUG_LVL0
      Write(6,"(a,a)")" The uncontracted <mu|Delta(r-Rn)|nu> shell and",
     &                 " primitive pair."
      Write(6,"(a,4(1x,I3))")" NAB, NXYZ, NIJ, NXYZET :", NAB, NXYZ,
     &                         NIJ, NXYZET
      Do i=1, NAB
         Write(6,"(1x,a,i4)") "The prim. pair: ", i
         write(6,"(6(1x,F20.10))") (BATCH(i,j), j=1, NXYZ)
      Enddo
#endif
         RETURN
         END
