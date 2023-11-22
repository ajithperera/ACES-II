      Program test  
      Implicit Double Precision (A-H, O-Z)

      Parameter (Maxang = 7)
      
      Integer Rx, Sy, Tz, Rstval
      Dimension Ideg(0:Maxang), Istart(0:Maxang), Iend(0:Maxang),
     &          Imom(0:Maxang), Lmnval(3,84),
     &          Q(Maxang+2, Maxang+2),Fac2(-1:4*Maxang),
     &          Fact(0:2*Maxang),Faco(0:2*Maxang),
     &          Qnl(0:2*Maxang, 0:2*Maxang), Rstval(3,84),
     &          Fijk(0:4*Maxang,0:4*Maxang,0:4*Maxang),
     &          ACS_xyz(0:Maxang,3), BCS_xyz(0:Maxang,3),
     &          Sum1_xyz(0:Maxang,3),Bcoefs(0:2*Maxang,0:Maxang),
     &          Ang(0:2*Maxang,0:2*Maxang)
     
      common/qstore/Q,alpha,xk,t
      common/factorials/Fact, Fac2,Faco
      common /Bncoefs/Bcoefs
      common /Fints/Fijk

      La1 = 0
      Lb1 = 0
      
c Set lmnval array:
     
      II = 0
      DO LVAL = 0,MAXANG-1
         DO L = LVAL,0,-1
            LEFT = LVAL - L
            DO M = LEFT,0,-1
               N = LEFT - M
               II = II + 1
               LMNVAL(1,II) = L
               LMNVAL(2,II) = M
               LMNVAL(3,II) = N
             ENDDO
          ENDDO
      ENDDO
C
      Write(6,"(a)") "The lmn values"
      Do i=1, 3
         Write(6,"(10I3)") (lmnval(i, i1),i1=1,84)
         Write(6,*)
      Enddo
c
       
      Ideg(0)= 1
      II = 1
      Do I = 1, Maxang
         II = II + 1
         Ideg(I) = II*(II+1)/2
      Enddo
    
      Write(6,*) "I*(I+1)/2 values"
      Write(6,"(8I3)") (Ideg(I), I= 0, Maxang)
      Write(6,*)
      
      Istart(0) = 1
      Istart(1) = 2
      Itotal = 3 
      Do I = 2, Maxang
         Istart(I) = Istart(I-1)+ Ideg(I-1) 
         itotal = Itotal + Ideg(I)
      Enddo

      Write(6,*) "IStart values"
      Write(6,"(8I3)") (Istart(I), I= 0, Maxang)
      Write(6,*) "The size of the Ixyz array ", Itotal 
      Write(6,*)

      Iend(0) = 1
      Do I = 1, Maxang
         Iend(I) = Iend(I-1) + Ideg(I)
      Enddo

      Write(6,*) "Iend values"
      Write(6,"(8I4)") (IEnd(I), I= 0, Maxang)
      Write(6,*)

      Do I=1, Maxang
         Imom(I) = I-1
      Enddo 

      Write(6,*) "Angular momentum values"
      Write(6,"(8I3)") (Imom(I), I= 1, Maxang)
      Write(6,*)
C
      Indx = 0
      NX =1
      Ny =0
      NZ =0
      Nxx = Nx + 1
      Nyy = Ny + 1
      Nzz = Nz + 1
      Nxyz = Nxx + Nyy + Nzz + 1
      Do N = 1, Nxyz
         Do Rx = 1, Nxx
            Do Sy = 1, Nyy
               Do Tz = 1, Nzz
                  Write(6,*) Rx+Sy+Tz, N+2

                  If (Rx+Sy+Tz .EQ. (N+2)) Then
                     Indx = Indx + 1
                     Write(6,*) Rx, Sy, TZ
                     Rstval(1, Indx) = Rx - 1
                     Rstval(2, Indx) = Sy - 1
                     Rstval(3, Indx) = TZ - 1
                  Endif

                Enddo
             Enddo
         Enddo
      Enddo
C 
      Write(6,*) "Checking the setting of rst values"
      Do i=1, 3
         Write(6,"(10I3)") (Rstval(i, i1),i1=1,84)
         Write(6,*)
      Enddo

C Get No. of cartesian components of La1 and Lb1

      Ia_Ncart = Ideg(La1)
      Ib_Ncart = Ideg(Lb1)
      Ia_start = Istart(La1)
      Ia_End   = IEnd(La1)
      Ib_start = Istart(Lb1)
      Ib_End   = IEnd(Lb1)

      ltot1 = La1 + Lb1 
      ltot1 = 0
    
      call Factorial
      call Dfactorial
      call Factorialo
      call Binomial_coefs
C
      Write(6,"(a)") "Double factorial array"
      Write(6,"(6(1x,F10.1))") (fac2(j1), j1=-1, 2*maxang)
      Write(6,*)

      Write(6,"(a)") "Factorial array"
      Write(6,"(6(1x,F10.1))") (fact(j1), j1=0, 11)
      Write(6,*) 
      
      Write(6,"(a)") "Factorial_O array"
      Write(6,"(6(1x,F10.1))") (faco(j1), j1=0, 11)
      Write(6,*) 

      Write(6,"(a)") "Binomial expansion Coefs. array"
      do i1 = 0, 2*Maxang
          Write(6,"(6(1x,F10.2))") (Bcoefs(i1, j1), j1=0,i1)
          Write(6,*)
      enddo
      Write(6,*)

      PI  = Dacos(-1.0D0)
      Sqrt_Fpi = Dsqrt(4.0D0*PI)
           Fpi = 4.0D0*Pi
      Fijk(0,0,0) = Fpi
 
      Do I = 0, 4*Maxang, 2
         Do J = 0, 4*Maxang, 2
            Do K = 0, 4*Maxang, 2

                  If (I .NE. 0) Then

                  Fijk(I,J,K) = Fijk(I-2,J,K)*(I-1)/
     &                              (I+J+K+1.0D0)  
                  Elseif (J .NE. 0) Then

                      Fijk(I,J,K) = Fijk(I,J-2,K)*(J-1)/
     &                                   (I+J+K+1.0D0)
                  Elseif (K. NE. 0) Then

                      Fijk(I,J,K) = Fijk(I,J,K-2)*(K-1)/
     &                               (I+J+K+1.0D0)
                  Endif

            Enddo 
         Enddo
      Enddo
      
C      Write(6,*) "Checking Angular integral Table"
C      Do i=0, 4*Maxang, 2
C        Do j=0, 4*Maxang, 2
C           Write(6,"(4F15.7)") (Fijk(i,j,k),k=0, 4*Maxang, 2)
C        Enddo
C      Enddo
C      Write(6,*)

      Call cortab
C
C Radial integral can be done here. Pot(ij)= Sum{sum{j)dijr^(nj-2)exp(zetar^2)j
C     
C      Do I = 1, Ideg(La1)
           exp1 = 2.225000
C      Enddo

C      Do I = 1, Ideg(La1)
           exp2 = 2.225000
C      Enddo

C Coordiantes of A, B and C, C is where the ECP is

      xa = 0.0
      ya = 0.0
      za =  1.417294

      xb = 0.0
      yb = 0.0
      zb =  1.417294
 
      xc = 0.0
      yc = 0.0
      zc =  1.417294

      xca = xc - xa
      yca = yc - ya
      zca = zc - ya

      xcb = xc - xb
      ycb = yc - yb 
      zcb = zc - zb

C Common center of A and B (P) and the exponent and the factor
        
      exp12 = exp1 + exp2
      disab = (xa-xb)**2 + (ya-yb)*2 + (za-zb)**2
      prefa = exp1*exp2*disab/exp12

      xp = (exp1*xa + exp2*xb)/exp12
      yp = (exp1*ya + exp2*yb)/exp12
      zp = (exp1*za + exp2*zb)/exp12
C
      xpc = xp - xc
      ypc = yp - yc
      zpc = zp - zc

      dpc2 = xpc**2 + ypc**2 + zpc**2
      If (dpc2 .Gt. 0.0D0) Then
         dpc = dsqrt(dpc2)
         xph = xpc/dpc
         yph = ypc/dpc
         zph = zpc/dpc
      Else
         dpc = 0.0D0
         xph = 0.0D0
         yph = 0.0D0
         zph = 0.0D0
      Endif

      Write(6,"(a,3(1x,F10.5))") "X,Y,Z hats ", xph, yph, zph

      dexp12 = 2.0D0*exp12

C Radial integral can be computed here.  N_sum = 1 
      DO I = 1, Maxang
         DO J= 1, Maxang
            Q(I,J) = 0.0D0
         ENDDO
      ENddo

      Do Ipot = 1, N_sum
  
C Combined the gaussian with C (potential is a Gaussian)

         exp12c = exp12 + Zeta
         prefac = exp12*Zeta*pc2/exp12c
         Prefac = Exp(-prefac)
         

         Rk = dexp12*pc
         Ck = dexp12**2*Dpc2/exp12c
C
C Xval=Beta^2/4*Alpha
C
         Alpha = exp12c
         Beta  = Dexp12*DPc
         Xval  = ((exp12**2)*(Dpc2))/Alpha  
C
         t =   1.0
         xk   = 3.0
         Nlp  = 3
         Xval = t
         Ltot0 =6
         Ltot1 =7
         Beta = Xk
         
         call Radint_recur(Nlp, Ltot0, Alpha, Beta, Xval, Qnl)
         call recur1(Nlp,Ltot1)

         Write(6,*) "The standard"
         Do i=1, Maxang+2
            Write(6,"(4(1x,F10.7))") (Q(j,I), j=1, Maxang+2)
         Enddo

         Write(6,*)
         Write(6,*) "The test"
         Do i=0, 2*Maxang
            Write(6,"(4(1x,F10.7))") (Qnl(j,I), j=0, 2*Maxang)
         Enddo 
      
      Enddo
C outer most loops are over on La_Ncart and Lb_Ncart (s=1,p=3,d=6...)

      Do Ia_cart = Ia_start,  Ia_end

C Here extrct Nax, Nay and Naz for Ia_cart and Ib_cart; for example
C when La_Ncart=3, pfunctions, Ia_cart=px, Ia_cart=py, Ia_cart=3, pz.

         Nax =   Lmnval(1, Ia_cart)
         Nay =   Lmnval(2, Ia_cart)
         Naz =   Lmnval(3, Ia_cart)

         call Make_xyz(Nax, XCA, ACS_xyz(0,1))
         call Make_xyz(Nay, YCA, ACS_xyz(0,2))
         call Make_xyz(Naz, ZCA, ACS_xyz(0,3))
C
         Write(6,*) "The ACS_xyz array" 
         do j=1, 3
            Write(6,*) "The J value :", J
            if (j .eq. 1) k=Nax
            if (j .eq. 2) k=Nay
            if (j .eq. 3) k=Naz
            Write(6,"(6(1x,F10.5))") (Acs_xyz(i,j), i=0,k)
         enddo
      
C Here we can set up the r,s,t values for this angular momentum
C 
CSSS         Call set_rst(Nax, Nay, Naz, Lstval(1,1))

         Do Ib_cart = Ib_start, Ib_end

            Nbx =   Lmnval(1, Ib_cart)
            Nby =   Lmnval(2, Ib_cart)
            Nbz =   Lmnval(3, Ib_cart)

            call Make_xyz(Nbx, XCB, BCS_xyz(0,1))
            call Make_xyz(Nby, YCB, BCS_xyz(0,2))
            call Make_xyz(Nbz, ZCB, BCS_xyz(0,3))
C  
            Write(6,*) "The BCS_xyz array" 
            do j=1, 3
               Write(6,*) "The J value :", J
               if (j .eq. 1) k=Nbx
               if (j .eq. 2) k=Nby
               if (j .eq. 3) k=Nbz
               Write(6,"(6(1x,F10.5))") (Bcs_xyz(i,j), i=0,k)
            enddo
C
            call make_bc_cacb_xyz(Nax, Nbx, ACS_xyz(0,1),
     &                            BCS_xyz(0,1), Sum1_xyz(0,1))
            call make_bc_cacb_xyz(Nay, Nby, ACS_xyz(0,2),
     &                            BCS_xyz(0,2), Sum1_xyz(0,2))
            call make_bc_cacb_xyz(Naz, Nbz, ACS_xyz(0,3),
     &                            BCS_xyz(0,3), Sum1_xyz(0,3))

            Write(6,*) "The SUm1_xyz array" 
            do j=1, 3
               Write(6,*) "The J value :", J
               if (j .eq. 1) k=Nax
               if (j .eq. 2) k=Nay
               if (j .eq. 3) k=Naz
               Write(6,"(6(1x,F10.5))") (Sum1_xyz(i,j), i=0,k)
            enddo

C Lambda summation limits:(itype=0 for type1 and 1 for type 2

            Lab_hi = La1 + Lb1
   
            Write(6,"(a,6(1x,I3))") "Naxyz and Nbxyz", Nax, Nay, Naz,
     &                               Nbx, Nby, Nbz      
            Write(6,*)
C 
c Now we have loops over Nax, Nay, Naz, Nbx, Nby, Nbz
c
                     Do iabx = 0, Nax+Nbx
                        Do iaby = 0, Nay+Nby
                           Do iabz = 0, Naz+Nbz
C                             
                              Iabxyz = iabx+iaby+iabz
C 
                              Do Lam = 0, lab_hi

c Loops over mu (-la to +la and -lb to +lb)

                                 Ang_N_lam = 0.0D0
                                 Ang(Iabxyz, Lam) = 0.0D0
c
                                 Do mu = -Lam, Lam
C 
C The non-zero angular integral occurs for Lam < iabxyz and Lam+iabxyz
C even                                  
                                    If (Lam .LE. Iabxyz .and.  
     &                                  Mod(Lam+Iabxyz,2) .NE. 1)
     &                              Then

                                    Call Make_res(Iabx, Iaby, Iabz,
     &                                            Lam, Mu, Xph, Yph, 
     &                                            Zph, Res)
                                    Call Ang_int(0, 0, Iabx, Iaby,
     &                                           Iabz, Lam, Mu, 
     &                                           Ang_N_lam_mu) 

      Write(6,"(a,3(1x,F15.8))") "Before Sum:", Ang_N_lam,
     &                            ANg_N_lam_mu, Res

                              Ang_N_lam = Ang_N_lam + Ang_N_lam_mu*
     &                                    Res*Sqrt_fpi

C
      Write(6,*) "The angular part of the ints in Main"
      Write(6,"(5(1x,i2),(2(1x,F15.8)))"),Lam,mu,Iabx,Iaby,Iabz,
     &                                 Ang_N_Lam
c
                                    Endif
C
                                 Enddo

                      Ang(Iabxyz, Lam) = Ang(Iabxyz, Lam) +
     &                                   Ang_N_lam*
     &                                   Sum1_xyz(Iabxyz,1)*
     &                                   Sum1_xyz(Iabxyz,2)*
     &                                   Sum1_xyz(Iabxyz,3)

                      Write(6,"(a,F15.8)")"The angular int:",
     &                                     Ang(Iabxyz, Lam)
                      Write(6,*)
                      Enddo
            
                Enddo
            Enddo
          Enddo
C
        Enddo
      Enddo
C
      DPC = 1.0D0

      If (DPC .EQ. 0.0) Then

c PC = means A, B and C coincide, hence M_lam(0) = 1 (M is Bessel functions)

           Lamu = 1 
      Else
           Lamu = ltot1
      Endif

      Write(6,"(a,1x,I2)") "The value of Ltot1", Ltot1
      Do lam = 1, Lamu

         Nhi = Ltot1 - Mod(Ltot1 - Lam, 2)
C         Write(6,"(a,1x,i2)") "Mod of (Ltot1-lam,2)", 
C     &         Mod(Ltot1 - Lam, 2)
        
         Write(6,*) lam
         Do N = lam, nhi, 2

           Write(6, "(a,I2,2x,I2)") "Value of lambda and N ", lam, N

         Enddo
      Enddo
      Stop
      End
   
      subroutine recur1(npi, ltot1)
c     ----- controls computation of q(n,l) by recurrence -----
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
      common/qstore/q(9,9),alpha,xk,t

      data a0, a1s2, a1, a2, a3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/
      np1l=mod(npi,2)+1
      np1u=npi+1


      Write(6,*) alpha, xk, t,npi, ltot1

      if(xk.ne.a0) go to 24
c     ----- xk=0 -----
      f1sta=a1/alpha
      Write(6,*) "np1l and np1u", np1l, np1u

      if(np1l.eq.2) go to 12
      Write(6,*) pi, a1s2*dsqrt(pi/alpha)
      q(1,1)=a1s2*dsqrt(pi/alpha)
      fista=a1s2*f1sta
      go to 14

   12 q(1,1)=a1s2*f1sta
      Write(6,*) "Printing q(1,1)", q(1,1)
      fista=f1sta

      Write(6,*) "np1l and np1u", np1l, np1u
   14 if(np1l.eq.np1u) go to 18

      np1l=np1l+2
      do 16 np1=np1l,np1u,2
      q(1,1)=fista*q(1,1)
   16 fista=fista+f1sta

   18 if(ltot1.le.2) return
      do 20 lp1=3,ltot1,2
      Write(6,*) "ltot1", ltot1, lp1, fista, q(lp1-2,1)
      q(lp1,1)=fista*q(lp1-2,1)
      
   20 fista=fista+f1sta
      return

      Write(6,*) (ltot1 -2)

   24 if(ltot1-2) 28,26,30
c     ----- ltot1 1 or 2 -----
   26 q(2,2)=qcomp(npi+1,1)
   28 q(1,1)=qcomp(npi,0)
      return

   30 talph=alpha+alpha
      Write(6,*) "I am here1", npi
      Write(6,*) "I am here2", n2, np1l
      
      if(npi.gt.0) go to 32
      n2=3
      Write(6,*) "np1l when npi 0", np1l
      go to 34

   32 n2=1


   34 if(np1l.eq.2) go to 48

c     ----- npi even -----
      q(n2,1)=qcomp(2,0)
      Write(6,*) "I am taking this value", q(n2,1), n2

      lu=ltot1-n2
       
       Write(6,*) "I am here3", lu

      if(lu.eq.0) go to 38
c     ----- use (38a) -----

      do 36 l=1,lu
      Write(6,*) "I am here n2 and l+1 value", n2+l, l+1
      q(n2+l,l+1)=(xk/talph)*q(n2+l-1,l)
      Write(6,*) q(n2+l,l+1), q(n2+l-1,l)
   36 continue
      Write(6,*) "I am here 4", (npi-2), q(4,2) 
      if(npi-2) 38,68,60

c     ----- npi=0 -----
   38 if(t.gt.a3) go to 42

c This block compute ltot <= 4 (all the contrbutions for t< 3)
c     ----- recur down using (38c) etc. -----
      q(ltot1,ltot1)=qcomp(ltot1-1,ltot1-1)
      f2lp1=(ltot1+ltot1-3)
      q(ltot1-1,ltot1-1)=xk*(q(ltot1,ltot1-2)-q(ltot1,ltot1))/f2lp1
      do 40 lr=3,ltot1
      l=ltot1-lr
      f2lp1=f2lp1-a2
   40 q(l+1,l+1)=(talph*q(l+3,l+1)-xk*q(l+2,l+2))/f2lp1

      if(ltot1-4) 74,74,68

c     ----- recur up using (38b) etc. -----
c This block compute ltot <= 4 (all the contrbutions for t> 3)

   42 q(1,1)=qcomp(0,0)
      f2lm1=a1
      do 44 lp2=3,ltot1
      q(lp2-1,lp2-1)=(talph*q(lp2,lp2-2)-f2lm1*q(lp2-2,lp2-2))/xk
   44 f2lm1=f2lm1+a2
      Write(6,*) "f2lml value", f2lm1, q(ltot1,ltot1-2), 
     &            q(ltot1-1,ltot1-1)
      q(ltot1,ltot1)=q(ltot1,ltot1-2)-f2lm1*q(ltot1-1,ltot1-1)/xk

      if(ltot1-4) 74,74,68

c     ----- npi odd -----

   48 if(t.gt.a3) go to 52
c     ----- recur down using (38f) -----
      q(ltot1,ltot1)=qcomp(ltot1,ltot1-1)
      q(ltot1-1,ltot1-1)=qcomp(ltot1-1,ltot1-2)
      Write(6,*) qcomP(7,6), Qcomp(6,5)
      Write(6,*) "Starting vals", q(ltot1-1,ltot1-1), q(ltot1,ltot1)
      f2lp2=(ltot1+ltot1-2)
      do 50 lr=3,ltot1
      l=ltot1-lr
      f2lp2=f2lp2-a2
      Write(6,*) "Going down using 38f", l+1, f2lp2, q(l+3,l+3), 
     &            q(l+2,l+2)
      q(l+1,l+1)=(talph*q(l+3,l+3)-(xk-(f2lp2+a1)*(talph/xk))*
     1   q(l+2,l+2))/f2lp2
   50 continue
      go to 56

c     ----- recur up using (38e) -----
   52 q(1,1)=qcomp(1,0)
      q(2,2)=qcomp(2,1)
      f2lm2=a0
      do 54 lp1=3,ltot1
      f2lm2=f2lm2+a2
      Write(6,*) "odd; x> 3", q(lp1-2,lp1-2), q(lp1-1,lp1-1), f2lm2
      q(lp1,lp1)=(f2lm2*q(lp1-2,lp1-2)+(xk-(f2lm2+a1)*(talph/xk))*
     1   q(lp1-1,lp1-1))/talph
   54 continue
   56 if(npi.eq.1) go to 68

c     ----- npi greater than 2.  use (38d) etc. -----
   60 np1l=mod(npi-3,2)+4
      fnm3=(np1l-6)
      do 64 np1=np1l,np1u,2
      Write(6,*) "I am here for np1 loop", np1l, np1u, fnm3
      fnm3=fnm3+a2
      fnplm1=fnm3
      do 62 lp2=2,ltot1
      fnplm1=fnplm1+a2
      Write(6,*) "fnplm1", fnplm1, lp2-1, q(lp2-1,lp2-1), q(lp2,lp2)
      q(lp2-1,lp2-1)=(fnplm1*q(lp2-1,lp2-1)+xk*q(lp2,lp2))/talph
   62 continue
      Write(6,*) "what is in here", ltot1 , fnm3+a1, q(ltot1,ltot1),
     &           q(ltot1-1,ltot1-1), qcomp(ltot1+3,ltot1)
      q(ltot1,ltot1)=((fnm3+a1)*q(ltot1,ltot1)+xk*
     1   q(ltot1-1,ltot1-1))/talph
   64 continue
c     ----- use (38d) -----
C ltot1 > 4 
   68 nbnd=(ltot1-n2)/2
      lp1mx=(ltot1-n2)+1
      fnm3=(npi+n2-4)
      do 70 ibnd=1,nbnd
      fnm3=fnm3+a2
      fnplm1=fnm3
      n2=n2+2
      lp1mx=lp1mx-2
     
      do 70 lp1=1,lp1mx
      fnplm1=fnplm1+a2
      Write(6,*) "I am Here 68", n2+lp1-1, lp1, fnplm1, 
     &            q(n2+lp1-3,lp1), q(n2+lp1-2,lp1+1)
      q(n2+lp1-1,lp1)=(fnplm1*q(n2+lp1-3,lp1)+xk*q(n2+lp1-2,lp1+1))/
     1   talph
   70 continue
   74 return
      end
      function qcomp(n,l)
c     ----- computes q(n,l)                                -----
c     ----- scaled by exp(-t) to prevent overflows         -----
c     ----- arguments are alpha, xk, and t=xk**2/(4*alpha) -----
c     ----- no restriction on the magnitude of t           -----
c     ----- increase dfac array to raise n, l restrictions -----
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
      common /dfac/ dfac(29)
      common/qstore/dum1(81),alpha,xk,t
      dimension tmin(9)
      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1   18.0d0,15.0d0/
      data am1,   a0,   accpow, accasy, a1,   a2,   a4
     1    /-1.0d0,0.0d0,1.0d-14,1.0d-10,1.0d0,2.0d0,4.0d0/
      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c     ----- use alternating series (n+l.le.22.and.l.le.10) -----

      if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(xk/(alpha+alpha))**l
      prefac=sqpi2*xkp*dfac(n+l+1)/
     1   ((alpha+alpha)**((n-l)/2)*dsqrt(alpha+alpha)*dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
 10   if(num.eq.0) go to 20
      fnum=num
      term=term*fnum*t/(xden*xc)
      xc=xc+am1
      sum=sum+term
      num=num+2
      xden=xden+a2
      go to 10
   20 qcomp=prefac*sum
      return
   30 if(t.lt.tmin(min0(n,8)+1)) go to 60
c     ----- use asymptotic series (arbitrary n,l) -----
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp*sqpi2/((alpha+alpha)*dsqrt(alpha+alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(dabs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
c     ----- use power series (n+l.le.22.and.l.le.10) -----
   60 if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(xk/(alpha+alpha))**l
      prefac=dexp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) prefac=prefac*sqpi2/dsqrt(alpha+alpha)
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end
      subroutine cortab
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    common/fact/fac(13),fprod(7,7)
      common/fact/fac(17),fprod(9,9)
c*    common /zlmtab/ zlm(130)
      common /zlmtab/ zlm(377)
      data a0, a1s2, a1, a2, a3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/
c     ----- this value of pi may have to be changed if the     -----
c     ----- compiler will not accept extra significant figures -----
      pi=3.141592653589793238462643383279502884197169399d0
      sqpi2=sqrt(a1s2*pi)
      rt3=sqrt(a3)
      rt6=sqrt(6.0d0)
      rt5s2=sqrt(2.5d0)
      rt3s2=sqrt(1.5d0)
      dfac(1)=a1
      dfac(2)=a1
      fi=a0
c*    do 10 i=1,21
      do 10 i=1,27
      fi=fi+a1
   10 dfac(i+2)=fi*dfac(i)
      fac(1)=a1
      fi=a0
c*    do 20 i=1,12
      do 20 i=1,16
      fi=fi+a1
   20 fac(i+1)=fi*fac(i)
c*    do 30 l1=1,7
      do 30 l1=1,9
      do 30 k1=1,l1
   30 fprod(k1,l1)=fac(l1+(k1-1))/(fac(k1)*fac(l1-(k1-1)))
      Return 
      End

      Subroutine Make_BC_CACB_xyz(Nkaxyz, Nkbxyz, CASxyz, CBSxyz,
     &                            TCABxyz)
C
c This routine compute sums; Sum_a,b{a,b=0,n_i,n_j}{fac(ni,a)CA_i^^(n_i-a)
C CB_j^^(n_j-b)
C where i=x,y,z
C
      Parameter (Maxang = 7)
      Implicit Double Precision (A-H, O-Z)

      Dimension TCABxyz(0:Maxang), CASxyz(0:Maxang),
     &          CBSxyz(0:Maxang), Bcoefs(0:2*Maxang,0:2*Maxang)

      common /Bncoefs/Bcoefs

      Do I = 0, Maxang
        Tcabxyz(i) = 0.0D0
      Enddo

      Do Ka = 0, Nkaxyz
         Pre_Faca = Bcoefs(Nkaxyz, Ka)
         Do Kb = 0, Nkbxyz

            Pre_Facb = Bcoefs(Nkbxyz, Kb)

            TCABxyz(Ka+Kb) = TCABxyz(Ka+Kb) + Pre_Faca*Pre_Facb*
     &                    CASxyz(Nkaxyz-Ka)*CBSxyz(Nkbxyz-Kb)

            Write(6,"(2(1x,i3),5(1x,F10.4))"), Ka, Kb, 
     &            pre_faca, Pre_Facb, CASxyz(Nkaxyz-Ka), 
     &            CBSxyz(Nkbxyz-Kb), TCABxyz(kA+KB)
         Enddo
      Enddo

      Return
      End
     
      Subroutine Make_xyz(Nxyz, Cxyz, CS_xyz)
C
c This routine compute arrays PI{n_i=n,,k_i=0...n_i}CA_i^^(n_i-k_i) i=x,y,z
C where i=x,y,z
C

      Implicit Double Precision (A-H, O-Z)

      Parameter(Maxang=7)

      Dimension CS_xyz(0:Maxang)

      CS_xyz(0) = 1.0D0
      CS_xyz(1) = Cxyz
      
      IF (NXyz .EQ. 0 .OR. Nxyz .EQ. 1) Return

      Do I = 2, Nxyz

            CS_xyz(I) = CS_xyz(I-1)*CS_xyz(1)

      Enddo

      Return
      End
C
      Subroutine Binomial_coefs

      Implicit Integer (A-Z)
      Parameter(Maxang = 7)

      Double Precision Bcoefs(0:2*Maxang,0:2*Maxang)
      common /Bncoefs/Bcoefs

      Bcoefs(0, 0) = 1.0D0

      Do I1 = 1, 2*Maxang

           J2 = 0
           Bcoefs(I1, J2) = 1.D0

           Do J1 = 2,  I1
              J2 = J2 + 1
              Bcoefs(I1, J2) = Bcoefs(I1-1,J2-1) + Bcoefs(I1-1,J2)
           Enddo

              J2 = J2 + 1
              Bcoefs(I1, J2) = 1.D0

      Enddo

      Return
      End
C
      Subroutine Factorial
C
C Compute factorial of n! (n(n-1)...21), and 0! = 1
c
      Implicit Double Precision(A-H, O-Z)
      Parameter(Maxang=7)

      Dimension Fact(0:2*Maxang), Fac2(-1:4*Maxang),Faco(0:2*Maxang)

      common /factorials/Fact,Fac2,Faco

      fact(0)  = 1.0D0
      Fact(1)  = 1.0D0
      Fact(2)  = 2.0D0

      FI = 2.0D0
      Do I = 2, 2*Maxang
         FI = FI + 1.0D0
         Fact(I+1) =  FI*Fact(I)
      EndDo

      Return
      End

      Subroutine dfactorial
C
C Compute Double factorial of n!! (n(n-2)...531) or (n(n-2)..6421)
C depending on n is odd or even. 0!!,-1!! = 0
c
      Parameter (Maxang=7)
      Implicit Double Precision(A-H, O-Z)

      Dimension Fact(0:2*Maxang),Fac2(-1:4*Maxang),Faco(0:2*Maxang)

      common /factorials/Fact,Fac2,Faco

      Fac2(-1) = 1.0D0
      Fac2(0)  = 1.0D0
      Fac2(1)  = 1.0D0
      Fac2(2)  = 2.0D0

      FI = 2.0D0

      Do I = 2, 4*Maxang
         FI = FI + 1.0D0
         Fac2(I+1) =  FI*Fac2(I-1)
      EndDo

      Return
      End
      Subroutine Factorialo
C
C Compute factorial of (2n-1)! (2n-1..3.1), and 0! = 1
c
      Implicit Double Precision(A-H, O-Z)

      Parameter(Maxang=7)
      Dimension Fact(0:2*Maxang),Fac2(-1:4*Maxang),Faco(0:2*Maxang)

      common /factorials/Fact,Fac2,Faco

      faco(0)  = 1.0D0
      Faco(1)  = 1.0D0

      Do I = 2, 2*Maxang
         Faco(I) =  (2*I-1)*Faco(I-1)
      EndDo

      Return
      End

 
      Subroutine Radint_recur(Nlp, Ltot, Alpha, Beta, Xval, Qnl)
  
      Implicit Double Precision (A-H, O-Z)
      Logical Eve, Odd

      Parameter (Maxang = 7)
 
      Dimension Qnl(0:2*Maxang, 0:2*Maxang)

      Data A0, Ahalf, A1, A2, A3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/
      common /factorials/Fact(0:2*Maxang), Fac2(-1:2*Maxang),
     &                   Faco(0:2*Maxang)
C 
      Odd = ((Mod(Nlp,2)+1) .EQ. 2)
      Eve = ((Mod(Nlp,2)+1) .EQ. 1)
C
      Talpha = 2.0D0*Alpha
      PI     = DAcos(-1.0D0)
C
C Compute the radial integral, Int {0,inf} exp(-alpha r^n) M_l(beta r)
C employing recursion relations given in Mcmurchie and Davidson, 
c Journal of Computational Physics, 44, 289 (1981).
C
      If (Beta .Eq. 0.0D0) then
         
C Using standard integral for Bessel functions, Page 302, Text book
C of Mathematical functions (Abramowitz and Stegun, Eds. 1964, US
C BUreau of Standards).

         If (Nlp .Eq. 0) Then

            If (Ltot .LE. 1) Then

                Qnl(0,0) = Ahalf*Dsqrt(Pi/Alpha)

             Else

                Qnl(0,0) = Ahalf*Dsqrt(Pi/Alpha)

                Ieve = 1
                Do N = 2, Ltot, 2

                   NN = N - Ieve
                   Qnl(N,0) = Faco(NN)/((A2**(NN+1))*(Alpha**NN))*
     &                        (Dsqrt(Pi/Alpha))
                   Ieve  = Ieve + 1

                Enddo
             Endif

         Else

             If (Ltot .LE. 1) Then

                If ((Mod(Nlp,2)+1) .EQ. 1) Then

                   NN = Nlp/2
                   Qnl(0,0) = Faco(NN)/((A2**(NN+1))*(Alpha**NN))*
     &                     (Dsqrt(Pi/Alpha))
                Else

                   NN = (Nlp - 1)/2
                   Qnl(0,0) = Fact(NN)/(A2*(Alpha**(NN+1)))

                Endif

             Else

                Do N = 0, ltot, 2

                   If ((Mod(Nlp,2)+1) .EQ. 1) Then 

                      NN = (Nlp+N)/2
                      Write(6,*) "NN----", NN
                      Qnl(N,0) = Faco(NN)/((A2**(NN+1))*(Alpha**NN))*
     &                            (Dsqrt(Pi/Alpha)) + Qnl(N,0)

                   Else

                      Write(6,*) "NN----", NN, N
                       NN = (Nlp+N-1)/2
                       Qnl(N,0) = Fact(NN)/(A2*(Alpha**(NN+1)))

                   Endif    

                 Enddo

             Endif
         Endif

      Else 

         If (Ltot .LE. 1) Then
            Do L = 0, Ltot
               Do N = 0, L
C
C Only (N-L) positive and even is needed.
C
                   If ((N-L) .GE. 0) Qnl(N,L) = Qcomp(N+Nlp, L)

               Enddo
            Enddo
         Endif

         If (Ltot .LE. 1) Return
C
         If (Nlp .GT. 0) Then

            If (Eve) Then

C Npl > 0 and even 
C
C N = 0+2
               Qnl(0,0) = Qcomp(2, 0) 

               Do LL = 0, Ltot-1
C N = LL+3-2
                  N = LL+1
                  L = LL+1
                  Qnl(N,L) = Beta*(Qnl(N-1,L-1))/Talpha

               Enddo

               If ((Nlp-2) .EQ. 0) Then
C
C Use recursion 38D of Davidson paper.
C        
                   Lmax = Ltot 
                   Do NN = 0, Ltot-2, 2

                      Lmax = Lmax - 2
                      Do L= 0, Lmax

C NLM = N+2-1; N=0+2, 0+4,...

                         N = L + 2 + NN
                         NLM = N + L + 1
                         Write(6,*) "N and L", N, L , NLM, Qnl(N-2,L),
     &                               Qnl(N-1,L+1)
                         Qnl(N,L) = ((NLM)*Qnl(N-2,L) +
     &                                Beta*Qnl(N-1,L+1))/Talpha

                      Enddo
                   Enddo

               Else
C         
C Use recursion 38D in MD paper (multiple times).
C
                  Call Qcomp_38D(Maxang, Qnl, NlP, Ltot, Talpha, Beta)
             
               Endif

            Else

               If (Xval .GT. 3.0D0) Then
C      
C Npl > 0, and odd  Xval > 3.0
C N = 0+1, 1+1...
                   Qnl(0,0) = Qcomp(1,0)
                   Qnl(1,1) = Qcomp(2,1)
C
C Use recursion 38E in MD paper.
C
                   Do LL = 3, Ltot+1
        
                      L = LL - 1
                      N = LL - 1
                      Write(6,*) "X>3;odd", N, L, Qnl(N-2,L-2),
     &                            Qnl(N-1,L-1),N+L-2
                      Qnl(N,L) =  ((N+L-2)*Qnl(N-2,L-2) +
     &                             (Beta-(2*L-1)*(Talpha/Beta))*
     &                              Qnl(N-1,L-1))/Talpha
                   Enddo
C
C Use recursion 38D in MD paper (multiple times).
C
                   Call Qcomp_38D_odd(Maxang, Qnl, Nlp, Ltot, 
     &                                Talpha, Beta)

               Else
C      
C Npl > 0, and odd  Xval < 3.0
C N = Ltot+1, Ltot+1+....

                  Qnl(Ltot-1,Ltot-1)   = Qcomp(Ltot, Ltot-1)
                  Qnl(Ltot,Ltot) = Qcomp(Ltot+1, Ltot)
C
C Using 38F from MD paper.
C
                  Do LL = (Ltot-1), 1, -1
                     L = LL - 1
                     N = LL - 1
                     Write(6,*) "X<3;odd", N, L, Qnl(N+2,L+2),
     &               Qnl(N+1,L+1),N+L+2

                     Qnl(N,L) = (Talpha*Qnl(N+2,L+2) - (Beta-(2*L+3)*
     &                           (Talpha/Beta))*Qnl(N+1,L+1))/(N+L+2)
                  Enddo
C
C Use recursion 38D in MD paper (multiple times).
C
                  Call Qcomp_38D_Odd(Maxang, Qnl, Nlp, Ltot, 
     &                               Talpha, Beta)
                 
               Endif

            Endif          

         Else
C
C Nlp=0, Ltot >=2 bblock
C             
            If (Ltot .EQ. 2) Then
C
C Using prescription A in Fig. 1 of MD paper
c 
               Qnl(2,0) = Qcomp(2,0)  

               If (Xval .Gt. 3.0D0) Then

                  Qnl(0,0) = Qcomp(0,0)
C
C Using 38B of MD paper
C
                  Qnl(1,1) = (Talpha*Qnl(2,0) - Qnl(0,0))/Beta
C
C Combining 38A and 38B of MD paper,note the cancellation of constants.
C
CSSSS Qnl(3,1) = Beta*Q(2,0)/Talpha
C
                  Qnl(2,2) = Qnl(2,0) - 3.0D0*Qnl(1,1)/Beta

               Else
C
C Using prescription C in Fig. 1 and combining with 38C and 38A.
c 
                  Qnl(2,2) = Qcomp(2,2)

CSSS Qnl(1,1) = Beta*Qnl(2,0)/Talpha 
CSSS Qnl(3,1) = Beta*Qnl(2,0)/Talpha

                  Qnl(1,1) = Beta*(Qnl(2,0) - Qnl(2,2))/3.0D0
                  Qnl(0,0) = (Talpha*Qnl(2,0) - Beta*Qnl(1,1))

               Endif

            Else
C
C Nlp=0, Ltot>2 block
C Using recursion 38A of MD paper.
C
               Qnl(2,0) = Qcomp(2,0)

               Do L = 0, Ltot - 3

                  Qnl(3+L,L+1) = Beta*(Qnl(L+2,L))/Talpha

               Enddo

               If (Xval .Gt. 3.0D0) Then

                   Qnl(0,0) = Qcomp(0,0)

                   Do L = 1, Ltot - 1
C
C Using 38B of MD paper.
C
                      Qnl(L,L) = (Talpha*Qnl(L+1,L-1) - (2*L-1)*
     &                            Qnl(L-1,L-1))/Beta
                   Enddo
C
C Combining 38A and 38B of MD paper,note the cancellation of constants.
C
                   Qnl(Ltot,Ltot) = Qnl(Ltot,Ltot-2) - (2*Ltot-1)*
     &                              Qnl(Ltot-1, Ltot-1)/Beta
                   
C Using 38D of MD paper (Note the error in MD paper)
C
                   Lmax = Ltot - 2
                   Do NN = 0, Ltot-4, 2
                      Lmax = Lmax - 2
                      Do L= 0, Lmax

                         N = L + 4 + NN
                         Write(6,*) "N and L", N, L , Qnl(N-2,L),
     &                               Qnl(N-1,L+1)
                         Qnl(N,L) = ((N+L-1)*Qnl(N-2,L) +
     &                                Beta*Qnl(N-1,L+1))/Talpha

                      Enddo
                   Enddo
    
               Else
C
C Using prescription C in Fig. 1 and combining with 38C and 38A.
c
                  Qnl(2,2) = Qcomp(2,2)

CSSS         Qnl(1,1) = Beta*Qnl(2,0)/Talpha     
CSSS         Qnl(3,1) = Beta*Qnl(2,0)/Talpha

                  Qnl(1,1) = Beta*(Qnl(2,0) - Qnl(2,2))/3.0D0
                              
                  Do L = Ltot-3, 0, -1

                     Qnl(L,L) = (Talpha*Qnl(L+2,L) - Beta*Qnl(L+1,L+1)/
     &                          (2*L+1))
                  Enddo
C
C Using 38D of MD paper (Note the error in MD paper),
C
                   Lmax = Ltot - 2
                   Do NN = 0, Ltot-4, 2
                      Lmax = Lmax - 2
                      Do L= 0, Lmax

                         N = L + 4 + NN
                         Write(6,*) "N and L", N, L , Qnl(N-2,L),
     &                               Qnl(N-1,L+1)
                         Qnl(N,L) = ((N+L-1)*Qnl(N-2,L) +
     &                                Beta*Qnl(N-1,L+1))/Talpha

                      Enddo
                   Enddo

C
               Endif
          
            Endif
C
         Endif
C
      Endif

      Return
      End

      Subroutine Qcomp_38D_odd(Maxang, Qnl, Nlp, Ltot, Talpha, Beta)
C
      Implicit Double Precision(A-H, O-Z)

      Dimension Qnl(0:2*Maxang, 0:2*Maxang)
C
C Using 38D of MD paper (Note the error in MD paper)
C
      Write(6,*)
      Write(6,*) "The Qnl at entry to Qcomp_38D"
      Do i=0, 2*Maxang
         Write(6,"(4(1x,F10.7))") (Qnl(j,I), j=0, 2*Maxang)
      Enddo

      If (Nlp .EQ. 1) Then

         Lmax = Ltot
         Do NN = 0, Ltot-2, 2
            Lmax = Lmax - 2

            Do L = 0, Lmax

                N = L + 2 + NN
                NLM = N + NLP + L - 1
                Write(6,*) "in Qcomp_38D-0", N, L, NLM, Qnl(N-2,L),
     &                      Qnl(N-1,L+1)
                Qnl(N,L) = ((NLM)*Qnl(N-2,L) + Beta*Qnl(N-1,L+1))/
     &                       Talpha

            Enddo
         Enddo

      Else

         NP_Lo = Mod(Nlp-3,2) + 4
         NP_Hi = Nlp + 1
         Ifac  = 2
         
         Write(6,*) "in Qcomp 38D NP_LO,HI,LTOT", NP_LO, NP_hi, ltot
         Do NP_D = NP_Lo, NP_hi, 2
            Do LL = 2, Ltot+1
C
C N = LL (Note that what we are doing here is (Qnl(0+2+2..) where 
C Qnl(0+2) is stored in Qnl(0,0)
C
               N = LL - 2
               L = LL - 2
               NLM = N + L + Ifac
               Write(6,*) "In Qcomp_38D-1", N, L, NLM, Qnl(N,L),
     &                     Qnl(N+1,L+1)
               Qnl(N,L) = (NLM*Qnl(N,L) + Beta*Qnl(N+1,L+1))/Talpha

            Enddo

C N = Ltot + 1
               N = Ltot 
               L = Ltot 
               NLM = N + L + 3
               Write(6,*) "In Qcomp_38D-2", N, L, NLM, Qnl(N,L),
     &                     Qnl(N-1,L-1)
               Qnl(N,L) = (2*Qnl(N,L) + Beta*Qnl(N-1,L-1))/Talpha
               Ifac = Ifac + 2
         Enddo

         Lmax = Ltot
         Do NN = 0, Ltot-2, 2 
            Lmax = Lmax - 2
               
            Do L = 0, Lmax

                N = L + 2 + NN
                NLM = N + NLP + L - 1
                Write(6,*) "in Qcomp_38D-3", N, L, NLM, Qnl(N-2,L), 
     &                      Qnl(N-1,L+1)
                Qnl(N,L) = ((NLM)*Qnl(N-2,L) + Beta*Qnl(N-1,L+1))/
     &                       Talpha
     
            Enddo 
         Enddo
         
C                  
      Endif

      Return
      End

      Subroutine Qcomp_38D(Maxang, Qnl, Nlp, Ltot, Talpha, Beta)
C
      Implicit Double Precision(A-H, O-Z)

      Dimension Qnl(0:2*Maxang, 0:2*Maxang)
C
C Using 38D of MD paper (Note the error in MD paper)
C
      Write(6,*)
      Write(6,*) "The Qnl at entry to Qcomp_38D"
      Do i=0, 2*Maxang
         Write(6,"(4(1x,F10.7))") (Qnl(j,I), j=0, 2*Maxang)
      Enddo

      If (Nlp .EQ. 1) Then

         Lmax = Ltot
         Do NN = 0, Ltot-2, 2
            Lmax = Lmax - 2

            Do L = 0, Lmax

                N = L + 2 + NN
                NLM = N + NLP + L - 1
                Write(6,*) "in Qcomp_38D-3", N, L, NLM, Qnl(N-2,L),
     &                      Qnl(N-1,L+1)
                Qnl(N,L) = ((NLM)*Qnl(N-2,L) + Beta*Qnl(N-1,L+1))/
     &                       Talpha

            Enddo
         Enddo

      Else

         NP_Lo = Mod(Nlp-3,2) + 4
         NP_Hi = Nlp + 1
         Ifac  = 3
         
         Write(6,*) "in Qcomp 38D NP_LO,HI,LTOT", NP_LO, NP_hi, ltot
         Do NP_D = NP_Lo, NP_hi, 2
            Do LL = 2, Ltot+1
C
C N = LL (Note that what we are doing here is (Qnl(0+2+2..) where 
C Qnl(0+2) is stored in Qnl(0,0)
C
               N = LL - 2
               L = LL - 2
               NLM = N + L + Ifac
               Write(6,*) "In Qcomp_38D-1", N, L, NLM, Qnl(N,L),
     &                     Qnl(N+1,L+1)
               Qnl(N,L) = (NLM*Qnl(N,L) + Beta*Qnl(N+1,L+1))/Talpha

            Enddo

C N = Ltot + 1
               N = Ltot 
               L = Ltot 
               NLM = N + L + 3
               Write(6,*) "In Qcomp_38D-2", N, L, NLM, Qnl(N,L),
     &                     Qnl(N-1,L-1)
               Qnl(N,L) = (2*Qnl(N,L) + Beta*Qnl(N-1,L-1))/Talpha
               Ifac = Ifac + 2
         Enddo

         Lmax = Ltot
         Do NN = 0, Ltot-2, 2 
            Lmax = Lmax - 2
               
            Do L = 0, Lmax

                N = L + 2 + NN
                NLM = N + NLP + L - 1
                Write(6,*) "in Qcomp_38D-3", N, L, NLM, Qnl(N-2,L), 
     &                      Qnl(N-1,L+1)
                Qnl(N,L) = ((NLM)*Qnl(N-2,L) + Beta*Qnl(N-1,L+1))/
     &                       Talpha
     
            Enddo 
         Enddo
         
C                  
      Endif

      Return
      End

      Subroutine Ang_int(L, M, Ii, Jj, Kk, Lam, Mu, Angi)
     
      Parameter (Maxang = 7)

      Implicit Double Precision (A-H, O-Z)
      Dimension ILMF(122), ILMX(581),ILMY(581),ILMZ(581),ZLM(130)
      common /ztabcm/ZLM
      common /Fints/Fijk(0:4*Maxang,0:4*Maxang,0:4*Maxang)

      DATA ILMF/1, 2,3,4, 5,7,8,10,11, 12,14,16,18,20,22,23,
     &25,28,30,34,36,39,41,43,45, 47,50,53,57,61,64,67,70,72,76,78,
     &81,85,88,94,98,104,107,111,114,117,121,125,128,
     &131,135,139,145,151,157,163,167,171,175,178,184,188,194,197,
     &201,206,210,218,224,233,239,247,251,256,260,264,270,276,282,
     &288,292, 296,301,306,314,322,331,340,348,356,361,366,371,375,383,
     &389,398,404,412,416, 421,427,432,442,450,462,471,483,491,501,506,
     &512,517,522,530,538,547,556,564,572,577, 582/


      DATA ILMX/0, 1,0,0, 2,0,1,0,0,0,1, 3,1,2,0,1,1,0,0,0,0,1,2,0,
     &4,2,0,3,1,2,0,0,2,1,1,0,0,0,0,0,1,1,2,0,3,1, 5,3,1,0,2,4,3,1,1,3,
     &2,0,2,0,1,1,1,0,0,0,0,0,0,1,1,2,2,0,0,3,1,4,2,0, 6,4,2,0,5,3,1,4,
     &2,0,4,2,0,3,1,1,3,2,0,2,0,2,0,1,1,1,0,0,0,0,0,0,0,1,1,1,2,0,2,0,3,
     &1,3,1,4,2,0,5,3,1, 7,3,5,1,6,4,0,2,5,3,1,5,3,1,4,2,0,4,2,0,3,1,3,
     &1,3,1,2,0,2,0,2,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,2,0,2,0,2,0,3,1,3,
     &1,4,2,0,4,2,0,5,3,1,0,4,2,6, 8,6,4,2,0,7,5,3,1,6,4,2,0,6,4,2,0,5,
     &3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,2,0,2,0,2,0,2,0,1,1,1,1,0,
     &0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,
     &5,3,1,6,4,2,0,7,5,3,1, 9,7,5,3,1,8,6,4,2,0,7,5,3,1,7,5,3,1,6,4,2,
     &0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,3,1,2,0,
     &2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,2,0,
     &3,1,3,1,3,1,4,2,0,4,2,0,4,2,0,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,
     &8,6,4,2,0, 10,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,8,6,4,2,0,7,5,3,1,7,5,
     &3,1,6,4,2,0,6,4,2,0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,4,
     &0,2,3,1,3,1,3,1,3,1,2,0,2,0,2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,
     &0,0,0,0,1,1,1,1,1,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,0,2,4,0,2,4,0,2,
     &4,5,3,1,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,7,5,3,1,8,6,4,2,0,9,7,
     &5,3,1/

      DATA ILMY/0, 0,0,1, 0,2,0,0,0,1,1, 0,2,0,2,0,0,0,0,1,1,1,1,3,
     &0,2,4,0,2,0,2,2,0,0,0,0,0,0,1,1,1,1,1,3,1,3, 0,2,4,4,2,0,0,2,2,0,
     &0,2,0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,1,3,1,3,5, 0,2,4,6,0,2,4,0,
     &2,4,0,2,4,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,1,3,1,
     &3,1,3,1,3,5,1,3,5, 0,4,2,6,0,2,6,4,0,2,4,0,2,4,0,2,4,0,2,4,0,2,0,
     &2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,
     &3,1,3,5,1,3,5,1,3,5,7,3,5,1, 0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,6,0,
     &2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,
     &0,0,0,0,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,
     &1,3,5,1,3,5,7,1,3,5,7, 0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,
     &6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,
     &0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,
     &1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,
     &1,3,5,7,9, 0,2,4,6,8,10,0,2,4,6,8,0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,
     &4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,
     &4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,1,
     &1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,3,1,5,3,1,5,3,
     &1,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,9,1,3,
     &5,7,9/

      DATA ILMZ/0, 0,1,0, 0,0,1,2,0,1,0, 0,0,1,1,2,0,3,1,2,0,1,0,0,
     &0,0,0,1,1,2,2,0,0,3,1,4,2,0,3,1,2,0,1,1,0,0, 0,0,0,1,1,1,2,2,0,0,
     &3,3,1,1,4,2,0,5,3,1,4,2,0,3,1,2,0,2,0,1,1,0,0,0, 0,0,0,0,1,1,1,2,
     &2,2,0,0,0,3,3,1,1,4,4,2,2,0,0,5,3,1,6,4,2,0,5,3,1,4,2,0,3,3,1,1,2,
     &2,0,0,1,1,1,0,0,0, 0,0,0,0,1,1,1,1,2,2,2,0,0,0,3,3,3,1,1,1,4,4,2,
     &2,0,0,5,5,3,3,1,1,6,4,2,0,7,5,3,1,6,4,2,0,5,3,1,4,4,2,2,0,0,3,3,1,
     &1,2,2,2,0,0,0,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0,3,
     &3,3,1,1,1,4,4,4,2,2,2,0,0,0,5,5,3,3,1,1,6,6,4,4,2,2,0,0,7,5,3,1,8,
     &6,4,2,0,7,5,3,1,6,4,2,0,5,5,3,3,1,1,4,4,2,2,0,0,3,3,3,1,1,1,2,2,2,
     &0,0,0,1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,
     &3,1,1,1,1,4,4,4,2,2,2,0,0,0,5,5,5,3,3,3,1,1,1,6,6,4,4,2,2,0,0,7,7,
     &5,5,3,3,1,1,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,7,5,3,1,6,6,4,4,2,2,0,0,
     &5,5,3,3,1,1,4,4,4,2,2,2,0,0,0,3,3,3,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,
     &0,0,0,0,0, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,3,3,3,3,1,1,
     &1,1,4,4,4,4,2,2,2,2,0,0,0,0,5,5,5,3,3,3,1,1,1,6,6,6,4,4,4,2,2,2,0,
     &0,0,7,7,5,5,3,3,1,1,8,8,6,6,4,4,2,2,0,0,9,7,5,3,1,10,8,6,4,2,0,9,
     &7,5,3,1,8,6,4,2,0,7,7,5,5,3,3,1,1,6,6,4,4,2,2,0,0,5,5,5,3,3,3,1,1,
     &1,4,4,4,2,2,2,0,0,0,3,3,3,3,1,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,1,0,0,
     &0,0,0/

      Call make_ztab
C
      Id = L*(L+1) - M + 1

      IMn = ILmf(Id)
      Imx = ILmf(Id+1) - 1

      JD = Lam*(Lam+1) - Mu + 1

      JMn = ILmf(Jd)
      Jmx = ILmf(Jd+1) - 1
C
      Angi = 0.0D0
CSS      Write(6,"(a,10(1x,I3))") "At entry to Ang_int ", II,JJ,KK,Imn,
CSS     &                          lam,mu,Imx,Jmn,Jmx
      Do I = Imn, Imx
         Angj = 0.D0
         Do J = Jmn, Jmx
            Jx = ILmx(I) + Ii + ILmx(J)
            Jy = ILmy(I) + Jj + ILmy(J)
            Jz = ILmz(I) + Kk + ILmz(J)  
        
            Angj = Angj + Zlm(J)*Fijk(Jx, Jy, Jz)

CSSS      Write(6,*) "The angular part of Ints"
CSSS      Write(6,"(3(1x,i2),3(1x,F15.8))"), jx,jy,jz,
CSSS     &          Zlm(J),Fijk(Jx, Jy, Jz),angj

         Enddo
         Angi = Angi + Angj*Zlm(I)
      Enddo

      Return
      End

      Subroutine Make_ztab
 
      Implicit Double Precision (a-h,o-z)
   
      Dimension Zlm(130)
      common /ztabcm/ZLM

c This Routine set up the real spherical harmonics in the
C form of liner cartesian products. The table is generated
C using the expression provided by Richard J. Mathar. 
C Also, see
c
      Pi = Dacos(-1.0D0)
      fpi = 4.0D0*pi
C l=0,m=0

      zlm(1) = sqrt(1.0d0/fpi)
c
C l=1,m=+1,0,-1
c
      zlm(2) = sqrt(3.0d0/fpi)
      zlm(3) = zlm(2)
      zlm(4) = zlm(2)
C
C l=2,m=-2,-1....1,2
C
      zlm(5) = sqrt(15.0d0/fpi)/2.0d0
      zlm(6) = -zlm(5)
      zlm(7) = 2.0d0*zlm(5)
      zlm(8) = 3.0d0*sqrt(5.0d0/fpi)/2.0d0
      zlm(9) = -zlm(8)/3.0d0
      zlm(10) = zlm(7)
      zlm(11) = zlm(7)
c l=3
      zlm(12) = sqrt(35.0d0/(8.0d0*fpi))
      zlm(13) = -3.0d0*zlm(12)
      zlm(14) = sqrt(105.0d0/(4.0d0*fpi))
      zlm(15) = -zlm(14)
      zlm(16) = 5.0d0*sqrt(21.0d0/(8.0d0*fpi))
      zlm(17) = -zlm(16)/5.0d0
      zlm(18) = 5.0d0*sqrt(7.0d0/fpi)/2.0d0
      zlm(19) = -3.0d0*zlm(18)/5.0d0
      zlm(20) = zlm(16)
      zlm(21) = zlm(17)
      zlm(22) = 2.0d0*zlm(14)
      zlm(23) = -zlm(13)
      zlm(24) = -zlm(12)
c l=4
      zlm(25) = sqrt(315.0d0/(64.0d0*fpi))
      zlm(26) = -6.0d0*zlm(25)
      zlm(27) = zlm(25)
      zlm(28) = sqrt(315.0d0/(8.0d0*fpi))
      zlm(29) = -3.0d0*zlm(28)
      temp    = sqrt(45.0d0/fpi)/4.0d0
      zlm(30) = 7.0d0*temp
      zlm(31) = -zlm(30)
      zlm(32) = temp
      zlm(33) = -temp
      temp    = sqrt(45.0d0/(8.0d0*fpi))
      zlm(34) = 7.0d0*temp
      zlm(35) = -3.0d0*temp
      temp    = sqrt(9.0d0/fpi)/8.0d0
      zlm(36) = 35.0d0*temp
      zlm(37) = -30.0d0*temp
      zlm(38) = 3.0d0*temp
      zlm(39) = zlm(34)
      zlm(40) = zlm(35)
      temp    = sqrt(45.0d0/(4.0d0*fpi))
      zlm(41) = 7.0d0*temp
      zlm(42) = -temp
      zlm(43) = -zlm(29)
      zlm(44) = -zlm(28)
      zlm(45) = sqrt(315.0d0/(4.0d0*fpi))
      zlm(46) = -zlm(45)
c l=5
      zlm(47) = sqrt(693.0d0/(128.0d0*fpi))
      zlm(48) = -10.0d0*zlm(47)
      zlm(49) = 5.0d0*zlm(47)
      zlm(50) = sqrt(3465.0d0/(64.0d0*fpi))
      zlm(51) = -6.0d0*zlm(50)
      zlm(52) = zlm(50)
      temp    = sqrt(385.0d0/(128.0d0*fpi))
      zlm(53) = 9.0d0*temp
      zlm(54) = -27.0d0*temp
      zlm(55) = 3.0d0*temp
      zlm(56) = -temp
      temp    = sqrt(1155.0d0/fpi)/4.0d0
      zlm(57) = 3.0d0*temp
      zlm(58) = -zlm(57)
      zlm(59) = -temp
      zlm(60) = +temp
      temp    = sqrt(165.0d0/fpi)/8.0d0
      zlm(61) = 21.0d0*temp
      zlm(62) = -14.0d0*temp
      zlm(63) = temp
      temp    = sqrt(11.0d0/fpi)/8.0d0
      zlm(64) = 63.0d0*temp
      zlm(65) = -70.0d0*temp
      zlm(66) = 15.0d0*temp
      zlm(67) = zlm(61)
      zlm(68) = zlm(62)
      zlm(69) = zlm(63)
      temp    = sqrt(1155.0d0/fpi)/2.0d0
      zlm(70) = 3.0d0*temp
      zlm(71) = -temp
      zlm(72) = -zlm(54)
      zlm(73) = -zlm(55)
      zlm(74) = -zlm(53)
      zlm(75) = -zlm(56)
      zlm(76) = sqrt(3465.0d0/fpi)/2.0d0
      zlm(77) = -zlm(76)
      zlm(78) = zlm(49)
      zlm(79) = zlm(48)
      zlm(80) = zlm(47)
c l=6
      temp    = sqrt(3003.0d0/(512.0d0*fpi))
      zlm(81) = 6.0d0*temp
      zlm(82) = -20.0d0*temp
      zlm(83) = zlm(81)
      zlm(84) = sqrt(9009.0d0/(128.0d0*fpi))
      zlm(85) = -10.0d0*zlm(84)
      zlm(86) = 5.0d0*zlm(84)
      temp    = sqrt(819.0d0/(256.0d0*fpi))
      zlm(87) = 11.0d0*temp
      zlm(88) = -66.0d0*temp
      zlm(89) = zlm(87)
      zlm(90) = -temp
      zlm(91) = 6.0d0*temp
      zlm(92) = -temp
      temp    = sqrt(1365.0d0/(128.0d0*fpi))
      zlm(93) = 11.0d0*temp
      zlm(94) = -33.0d0*temp
      zlm(95) = 9.0d0*temp
      zlm(96) = -3.0d0*temp
      temp    = sqrt(1365.0d0/(512.0d0*fpi))
      zlm(97) = 33.0d0*temp
      zlm(98) = -zlm(97)
      zlm(99) = -18.0d0*temp
      zlm(100) = +18.0d0*temp
      zlm(101) = temp
      zlm(102) = -temp
      temp     = sqrt(273.0d0/fpi)/8.0d0
      zlm(103) = 33.0d0*temp
      zlm(104) = -30.0d0*temp
      zlm(105) = 5.0d0*temp
      temp     = sqrt(13.0d0/fpi)/16.0d0
      zlm(106) = 231.0d0*temp
      zlm(107) = -315.0d0*temp
      zlm(108) = 105.0d0*temp
      zlm(109) = -5.0d0*temp
      zlm(110) = zlm(103)
      zlm(111) = zlm(104)
      zlm(112) = zlm(105)
      temp     = sqrt(1365.0d0/(128.0d0*fpi))
      zlm(113) = 33.0d0*temp
      zlm(114) = -18.0d0*temp
      zlm(115) = temp
      zlm(116) = -zlm(94)
      zlm(117) = -zlm(93)
      zlm(118) = -zlm(95)
      zlm(119) = -zlm(96)
      temp     = sqrt(819.0d0/fpi)/4.0d0
      zlm(120) = 11.0d0*temp
      zlm(121) = -zlm(120)
      zlm(122) = -temp
      zlm(123) = temp
      zlm(124) = zlm(86)
      zlm(125) = zlm(85)
      zlm(126) = zlm(84)
      zlm(127) = sqrt(3003.0d0/(512.0d0*fpi))
      zlm(128) = -15.0d0*zlm(127)
      zlm(129) = -zlm(128)
      zlm(130) = -zlm(127)

      Return
      End

      Subroutine Make_res(Ii, Jj, Kk, Lam, Mu, xph, yph, zph, Res)
     
      Parameter (Maxang = 7)

      Implicit Double Precision (A-H, O-Z)
      Dimension ILMF(122), ILMX(581),ILMY(581),ILMZ(581),ZLM(130)
      common /ztabcm/ZLM
      common /Fints/Fijk(0:4*Maxang,0:4*Maxang,0:4*Maxang)

      DATA ILMF/1, 2,3,4, 5,7,8,10,11, 12,14,16,18,20,22,23,
     &25,28,30,34,36,39,41,43,45, 47,50,53,57,61,64,67,70,72,76,78,
     &81,85,88,94,98,104,107,111,114,117,121,125,128,
     &131,135,139,145,151,157,163,167,171,175,178,184,188,194,197,
     &201,206,210,218,224,233,239,247,251,256,260,264,270,276,282,
     &288,292, 296,301,306,314,322,331,340,348,356,361,366,371,375,383,
     &389,398,404,412,416, 421,427,432,442,450,462,471,483,491,501,506,
     &512,517,522,530,538,547,556,564,572,577, 582/


      DATA ILMX/0, 1,0,0, 2,0,1,0,0,0,1, 3,1,2,0,1,1,0,0,0,0,1,2,0,
     &4,2,0,3,1,2,0,0,2,1,1,0,0,0,0,0,1,1,2,0,3,1, 5,3,1,0,2,4,3,1,1,3,
     &2,0,2,0,1,1,1,0,0,0,0,0,0,1,1,2,2,0,0,3,1,4,2,0, 6,4,2,0,5,3,1,4,
     &2,0,4,2,0,3,1,1,3,2,0,2,0,2,0,1,1,1,0,0,0,0,0,0,0,1,1,1,2,0,2,0,3,
     &1,3,1,4,2,0,5,3,1, 7,3,5,1,6,4,0,2,5,3,1,5,3,1,4,2,0,4,2,0,3,1,3,
     &1,3,1,2,0,2,0,2,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,2,0,2,0,2,0,3,1,3,
     &1,4,2,0,4,2,0,5,3,1,0,4,2,6, 8,6,4,2,0,7,5,3,1,6,4,2,0,6,4,2,0,5,
     &3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,2,0,2,0,2,0,2,0,1,1,1,1,0,
     &0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,
     &5,3,1,6,4,2,0,7,5,3,1, 9,7,5,3,1,8,6,4,2,0,7,5,3,1,7,5,3,1,6,4,2,
     &0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,3,1,2,0,
     &2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,2,0,
     &3,1,3,1,3,1,4,2,0,4,2,0,4,2,0,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,
     &8,6,4,2,0, 10,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,8,6,4,2,0,7,5,3,1,7,5,
     &3,1,6,4,2,0,6,4,2,0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,4,
     &0,2,3,1,3,1,3,1,3,1,2,0,2,0,2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,
     &0,0,0,0,1,1,1,1,1,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,0,2,4,0,2,4,0,2,
     &4,5,3,1,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,7,5,3,1,8,6,4,2,0,9,7,
     &5,3,1/

      DATA ILMY/0, 0,0,1, 0,2,0,0,0,1,1, 0,2,0,2,0,0,0,0,1,1,1,1,3,
     &0,2,4,0,2,0,2,2,0,0,0,0,0,0,1,1,1,1,1,3,1,3, 0,2,4,4,2,0,0,2,2,0,
     &0,2,0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,1,3,1,3,5, 0,2,4,6,0,2,4,0,
     &2,4,0,2,4,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,1,3,1,
     &3,1,3,1,3,5,1,3,5, 0,4,2,6,0,2,6,4,0,2,4,0,2,4,0,2,4,0,2,4,0,2,0,
     &2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,
     &3,1,3,5,1,3,5,1,3,5,7,3,5,1, 0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,6,0,
     &2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,
     &0,0,0,0,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,
     &1,3,5,1,3,5,7,1,3,5,7, 0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,
     &6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,
     &0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,
     &1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,
     &1,3,5,7,9, 0,2,4,6,8,10,0,2,4,6,8,0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,
     &4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,
     &4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,1,
     &1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,3,1,5,3,1,5,3,
     &1,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,9,1,3,
     &5,7,9/

      DATA ILMZ/0, 0,1,0, 0,0,1,2,0,1,0, 0,0,1,1,2,0,3,1,2,0,1,0,0,
     &0,0,0,1,1,2,2,0,0,3,1,4,2,0,3,1,2,0,1,1,0,0, 0,0,0,1,1,1,2,2,0,0,
     &3,3,1,1,4,2,0,5,3,1,4,2,0,3,1,2,0,2,0,1,1,0,0,0, 0,0,0,0,1,1,1,2,
     &2,2,0,0,0,3,3,1,1,4,4,2,2,0,0,5,3,1,6,4,2,0,5,3,1,4,2,0,3,3,1,1,2,
     &2,0,0,1,1,1,0,0,0, 0,0,0,0,1,1,1,1,2,2,2,0,0,0,3,3,3,1,1,1,4,4,2,
     &2,0,0,5,5,3,3,1,1,6,4,2,0,7,5,3,1,6,4,2,0,5,3,1,4,4,2,2,0,0,3,3,1,
     &1,2,2,2,0,0,0,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0,3,
     &3,3,1,1,1,4,4,4,2,2,2,0,0,0,5,5,3,3,1,1,6,6,4,4,2,2,0,0,7,5,3,1,8,
     &6,4,2,0,7,5,3,1,6,4,2,0,5,5,3,3,1,1,4,4,2,2,0,0,3,3,3,1,1,1,2,2,2,
     &0,0,0,1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,
     &3,1,1,1,1,4,4,4,2,2,2,0,0,0,5,5,5,3,3,3,1,1,1,6,6,4,4,2,2,0,0,7,7,
     &5,5,3,3,1,1,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,7,5,3,1,6,6,4,4,2,2,0,0,
     &5,5,3,3,1,1,4,4,4,2,2,2,0,0,0,3,3,3,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,
     &0,0,0,0,0, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,3,3,3,3,1,1,
     &1,1,4,4,4,4,2,2,2,2,0,0,0,0,5,5,5,3,3,3,1,1,1,6,6,6,4,4,4,2,2,2,0,
     &0,0,7,7,5,5,3,3,1,1,8,8,6,6,4,4,2,2,0,0,9,7,5,3,1,10,8,6,4,2,0,9,
     &7,5,3,1,8,6,4,2,0,7,7,5,5,3,3,1,1,6,6,4,4,2,2,0,0,5,5,5,3,3,3,1,1,
     &1,4,4,4,2,2,2,0,0,0,3,3,3,3,1,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,1,0,0,
     &0,0,0/

      Call make_ztab

      JD = Lam*(Lam+1) - Mu + 1

      JMn = ILmf(Jd)
      Jmx = ILmf(Jd+1) - 1
C
      Res = 0.0D0
      Write(6,"(a,2(1x,i2))"),"Jmin and Jmax",Jmn,jmx

      Do J = Jmn, Jmx

         Jx = ILmx(J)
         Jy = ILmy(J)
         Jz = ILmz(J)  
        
          Res = Res + Zlm(J)*(Xph**Jx)*(Yph**Jy)*(Zph**Jz)
  
      Write(6,"(a,3(1x,i2),2(1x,F15.8))"), "The real Sph ",jx,jy,jz,
     &          Zlm(J),Res

      Enddo

      Return
      End

     
