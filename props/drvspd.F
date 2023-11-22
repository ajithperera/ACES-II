      SUBROUTINE DRVSPD(SCR, NATOM, IATNUM, NORBS, NUMATOM)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION SPENS,M_PE
      DIMENSION SCR(4*NORBS*NORBS), NATOM(NUMATOM)
      DIMENSION IATNUM(NUMATOM)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS

C Constants required to transform a.u. to MHz.

      au2ev = 27.2113961d+00      ! Conversion from au to eV
      ev2Mhz= 2.41798836d+08      ! Conversion from eV to MHz

      g_e      =  2.002319304386d+00  ! Electron g-factor
      C_light  = 137.0359895d+00      ! Speed of light in au
      M_pe     = 1836.152701d+00      ! Proton-electron mass ratio
      Beta_e   = 1.0d0/(2.0d0*C_light)
      Beta_N   = Beta_e/M_pe          ! Beta_N =1/2*c*Mpe
      pi       = acos(-1.0d0)
      Eight_pi_by_3 = (8.0d0*pi/3.0d0)

      fac = Eight_pi_by_3 * g_e * Beta_e * Beta_N
      con  = fac * au2ev * ev2Mhz

#ifdef _DEBUG_LVLM
      print *
      write(*,*)'Aiso = 8*pi/3  * g_e * g_N * Beta_e * Beta_N *CalcTe'
      Write(*,'(A,1x,F25.10)') 'Conversion from au to eV : ', au2ev
      Write(*,'(A,1x,F25.10)') 'Conversion from eV to Mhz: ', ev2Mhz
      Write(*,'(A,1x,F25.10)') 'g_e                      : ', g_e
      Write(*,'(A,1x,F25.10)') 'Beta_e                   : ', Beta_e
      Write(*,'(A,1x,F25.10)')'Beta_N                   : ', Beta_N
      print *
#endif 

      IMXATM=50
      IONE=1
      I000=1
      I010=I000+NORBS*NORBS
      I020=I010+NORBS*NORBS
      I030=I020+NORBS*NORBS
      I040=I030+NORBS*NORBS
      IATOM=0
      IRWND=0
      CALL GETREC(20,'JOBARC','MAP2ZMAT',NUMATOM,NATOM)
      CALL GETREC(20,'JOBARC','ATOMCHRG',NUMATOM,IATNUM)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") " Atom ordering and their charges"
      Write(6,"(6(1x,I4))") (natom(i),i=1,Numatom)
      Write(6,"(6(1x,I4))") (iatnum(i),i=1,Numatom)
#endif 
30    CALL SEEKLB('   DEN  ',IERR,IRWND)
      IF(IERR.NE.0)RETURN
      IRWND=1
      IF(IATOM.EQ.0)WRITE(LUOUT,1000)
      IATOM=IATOM+1
      CALL COMPPR(SPENS, SCR(I030), SCR(I010), NORBS, .FALSE.)
      CALL ATOM_GFAC(DBLE(IATNUM(NATOM(IATOM))),GFACTOR,ISOTOPE,
     $               ABUND,NISO)
#ifdef _DEBUG_LVL0
      Write(6,"(a,1x,I3,1x,F15.9)") "The g-factor for nuclei",
     $                               IATNUM(NATOM(IATOM)),GFACTOR
      Write(6,"(a,1x,F15.9)") "The conversion factor for nuclei",
     &                         CON*GFACTOR 
#endif   
      WRITE(LUOUT,2000) NATOM(IATOM),SPENS,SPENS*GFACTOR*CON
      GOTO 30

1000  FORMAT(T3,' Spin densities at atomic centers ',/,
     &       5X,'Z-matrix',20X,'Spin',20X, 'HFCC(A_iso)'/,5X,'center',
     &       21X,'Density',20X,'(MHz)') 
2000  FORMAT(5X,I3,20X,F15.10,10X,F15.10)
      END 
