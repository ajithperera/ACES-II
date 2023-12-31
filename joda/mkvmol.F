
C Purpose: Produce an VMOL .MOL file based on the ZMAT file
C
C Arguments:
C     PtGrp   Symbol for the point group symmetry to use in (input)
C             making the INP file.
C     NAtms   Number of atoms (input)
C     AtNr    Atomic numbers for each center (input)
C     NUnique Number of unique internal coordinates (input)
C     geom    Cartesian coordinates (input)
c             will not be the same as Q in FINDIF
C     GenBy   To hold the Genby list (internal)
C     IStat   Error codes (output)
C             = 0  Successful
C             = 1  Fatal error from dependent
C             = 3  I/O Error on AbiFil
C             = 5  I/O Error on ZFil
C             = 7  No non-dummy centers!
C             = 8  ZFil ends after Z-matrix
C             = 10 ZFil ends after JODA control info
C     BasName Holds the basis set name given in the ZMAT file
C
C Common blocks used:
C     LUnits  Sets unit names & numbers
C     Coord
C
C Internal Variables:
C     JobTtl  Holds the job title given in the ZMAT file
C     Wrk     Scratch string to copy lines form ZMAT to ABI
C
C Dependents
C     PitPtG  Translates a schoenflies symbol into Pitzer's form
C     WrPBas  Lookup & write basis set in Pitzer style
C
C Files used:
C     ZFil   Holds the basis specification & the "rest" of INP
C     AbiFil Gets written out here (kills previous one)
C     BasFil The basis set library
C
C Limitations:
C     Only looks up bases for the geometrically unique atoms.
c
C     Extra blank lines at the end of a "truncated" ZFil will cause
C     it to act like it encountered a real error on reading the file.

#include "flags.h"

      SUBROUTINE MKVMOL(geom, PTGRP,NATMS,NUNIQUE,ZSYM,ATNR,GENBY,
     $     SCRATCH, ISTAT, BasName)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "mxatms.par"
#include "io_units.par"
#include "linelen.par"
#include "baslims.par"
#include "molstrct.com"
#include "coord.com"

      character*(*) BasName

      integer    lusew
      parameter (lusew  =  15)
      character*(*) SewFil
      parameter    (SewFil = 'MOLCAS.INP')

      Character*(*) PtGrp, ZSym(NAtms)
      Integer AtNr(NAtms), GenBy(NAtms), IStat, AtoI,IZLoc(MxAtms)
      Character*(linelen) JobTtl, Wrk
      Character*3 Ka(3)
      Character*6 Slask,NullSt
      Character*1 Lb
      Character*2 Indxx(150)
      Character*2 Symbol

cYAU - parse ZMAT like fetchz()
      integer izl(2,7)
      logical bTmp

csb 1/97 Hold long basis set names
#include "baslen.par"
      character*(baslen) BasNam(MxAtms), BasNam_hold(MxAtms), BlnkBN,
     &                   ScrBas,EcpNam(MxAtms),EcpNam_hold(MxAtms)

      dimension icrcor(mxatms)
      common /turbo / iturbo,matom,ioffsh

CJDW  11/2/94. We now use IFLAGS2(108) to set TLA.
C     Parameter (TLA = 1.D-9)
CSSS      INTEGER SHELLANG,SHELLLOC,SHELLSIZ,SHELLPRM
CSSS      INTEGER SHELLORB,SHOFFSET,PRIMORBT,PROFFSET
      double precision geom(3*mxatms)
C
      LOGICAL XYZIN,NWFINDIF
      Logical Opn,bad123, GENBAS
CKJW 5-24-00
      logical seward,bCpBasis,bExist,bOpened
      Integer OldLu
      Dimension Nord(2*MxAtms),Scratch(NAtms)

      integer  atomnumb
C
C     Main OPTIM control data
C     IPRNT   Print level - not used yet by most routines
C     INR     Step-taking algorithm to use
C     IVEC    Eigenvector to follow (TS search)
C     IDIE    Ignore negative eigenvalues
C     ICURVY  Hessian is in curviliniear coordinates
C     IMXSTP  Maximum step size in millibohr
C     ISTCRT  Controls scaling of step
C     IVIB    Controls vibrational analysis
C     ICONTL  Negative base 10 log of convergence criterion.
C     IRECAL  Tells whether Hessian is recalculated on each cyc
C     INTTYP  Tells which integral program is to be used
C              = 0 Pitzer
C              = 1 VMol
C     XYZTol  Tolerance for comparison of cartesian coordinates
C
CSSS      COMMON /MOLSTR1/ SHELLANG(MXTNSH),SHELLLOC(MXTNSH),
CSSS     &                 SHELLSIZ(MXTNSH),SHELLPRM(MXTNSH),
CSSS     &                 BASISEXP(10000),BASISCNT(10000)
CSSS      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IIATOM
CSSS      COMMON /MOLSTR3/ SHELLORB(MXTNSH),SHOFFSET(MXTNSH),
CSSS     &                 PRIMORBT(MXTNSH),PROFFSET(MXTNSH)

      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON /INPTYP/XYZIN,NWFINDIF
      COMMON /FLAGS/IFLAGS(100),IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C     Length of BasNam is fixed by the basis library
C
C
C     This length is fixed by ABINITIO's input section
C
      Data Indxx /'01','02','03','04','05','06','07','08','09','00',
     & '0A','0B','0C','0D','0E','0F','0G','0H','0I','0J','0K','0L',
     & '0M','0N','0O','0P','0Q','0R','0S','0T','0U','0V','0W','0X',
     & '0Y','0Z','0a','0b','0c','0d','0e','0f','0g','0h','0i','0j',
     & '0k','0l','0m','0n','0o','0p','0q','0r','0s','0t','0u','0v',
     & '0w','0x','0y','0z',"11","12","13","14","15","16","17","18",
     & "19","10","21","22","23","24","25","26","27","28","29","20",
     & "31","32","33","34","35","36","37","38","39","30","41","42",
     & "43","44","45","46","47","48","49","40","51","52","53","54",
     & "55","56","57","58","59","50","61","62","63","64","65","67",
     & "68","69","60","71","72","73","74","75","76","77","78","79",
     & "70","81","82","83","84","85","86","87","88","89","80","91",
     & "92","93","94","95","96","97","98","99"/

csb
csb      DATA BASIS /'PBS     ','STO-3G  ','DZ      ',
csb     &            'VDZ     ','3-21G   ','6-31G   ',
csb     &            '4-31G   ','TZ      ','PVTZ    ',
csb     &            '6-311G  ','DZP     ','PVDZ    ',
csb     &            'TZP     ','PVQZ    ','TZ2P    ',
csb     &            'WMR     ','6-31G*  ','6-31G** ',
csb     &            '6-311G* ','6-311G**','PV5Z    ',
csb     &            'sv      ','svp     ','dz      ',
csb     &            'dzp     ','tzp     ','tz2p    ',
csb     &            'tzpl    ','tz2pl   ','qzp     ',
csb     &            'qz2p    ','pz2p    '/

c ----------------------------------------------------------------------

      external atomnumb
C
C Initialize some things.
C
      do i=1,baslen
        BlnkBN(i:i)=' '
      end do
      ScrBas=BlnkBN

      bad123 = .false.
      Slask='INTGRL'
      NullSt='      '
      LB='#'
      J=0
      I=0
      IGENER=1
      IF10=0
      IF18=1
      ID3=0
      IONE=1
      IRSTRT=0
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,INEWFD)
      IECP=IFLAGS(71)
      IECPU=0
      IPRINT=0
      INUTMP=0
      ISTATS=0
      ISTAT=0
      IOFFPRIM=1
      IOFFSHEL=1
C
      ICOUNT=0
c
c iturbo = 1: create an interfacefile to TURBOMOLE
c
      iturbo=IFLAGS(h_IFLAGS_turbomole)
      if (iturbo.eq.1) then
        matom=0
        ioffsh=0
      endif
C
C SEE IF A GENBAS FILE IS PRESENT
C
      INQUIRE(FILE='GENBAS',EXIST=GENBAS)
      IF (.NOT.GENBAS) INQUIRE(FILE='ZMAT.BAS',EXIST=GENBAS)
C
C     Cautiously open the ZMAT file
C
      Inquire (FILE=ZFil, OPENED=Opn)
      If (.not.Opn) Open (LuZ, FILE=ZFil, STATUS='OLD')
      Rewind LuZ

C     skip the header
      bTmp = .true.
      do while (bTmp)
         read (LuZ, '(A)', ERR=8500, END=8510) JobTtl
         call parsez(JobTtl,izl)
         i = izl(1,1)
         bTmp = ((i.eq.0).or.(JobTtl(i:i).EQ.'%'))
      end do

      IF ((.NOT.XYZIN).OR.(NWFINDIF)) THEN
         i = 1
         do while (i.le.natms)
            read (LuZ,'(A)', ERR=8500, END=8510) Wrk
            call parsez(wrk,izl)
            if (i.eq.3) then
               izz = AtoI(wrk(izl(1,2):izl(2,2)))
               bad123 = ( izz.eq.1            .and.
     &                    zsym(1)(1:1).ne.'X' .and.
     &                    zsym(2)(1:1).ne.'X'       )
            end if
            i = i + 1
         end do
      ELSE
         DO I = 1, NATMS
            READ (LUZ,'(A)', ERR=8500, END=8510) WRK
         END DO
      END IF
      IF (IFLAGS(68) .EQ. 0) THEN
         Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
         Do i = 1, NUnique
            Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
         End Do
      END IF
      Read (LuZ, '(A)', ERR=8500, END=8800) Wrk
C
C I think this is the best way to avoid the problem of Cartesian  
C frequency finite difference runs when the first two atoms are 
C not the same. This problem must have existed from the very 
C begining. It get unnoticed because we do not do lots and lots 
C of Cartesian Frqs (and then of course the first two atoms are
C the same, then no one notice it). 01/2016, Ajith Perera.

      IF (XYZIN) BAD123 = .FALSE.

c YAU : old
c     Get the optimizer settings & another blank line after
c      Call GtFlgs(0,IErr,
c     &            IPrnt, INR,   IContl,IVec,  IDie,
c     &            ICurvy,IStCrt,IMxStp,IVib,  IRecal,
c     &            IntTyp,IDisFD,IGrdFD,ICNTYP,ISym,
c     &            IBasis,IDFGHI,
c     &            BasName)
c      If (IErr.eq.2) Goto 8800
c      If (IErr.eq.1) Goto 8500
cC      Read (LuZ, '(11I3)', ERR=8500, END=8800) IPrnt, INR, IContl,
cC     $ IVec, IDie, ICurvy, IStCrt, IMxStp, IVib, IRecal, IntTyp
c      Read (LuZ, '(A)', ERR=8500, END=8810) Wrk
c YAU : new
      IPrnt  = iFlags(h_IFLAGS_joda_print)
      INR    = iFlags(h_IFLAGS_opt_method)
      IContl = iFlags(h_IFLAGS_convergence)
      IVec   = iFlags(h_IFLAGS_eigenvector)
      IDie   = iFlags(h_IFLAGS_negeval)
      ICurvy = iFlags(h_IFLAGS_curvilinear)
      IStCrt = iFlags(h_IFLAGS_stp_siz_ctl)
      IMxStp = iFlags(h_IFLAGS_max_step)
      IVib   = iFlags(h_IFLAGS_vibration)
      IRecal = iFlags(h_IFLAGS_eval_hess)
      IntTyp = iFlags(h_IFLAGS_integrals)
      IDisFD = iFlags(h_IFLAGS_fd_stepsize)
      IGrdFD = iFlags(h_IFLAGS_points)
      ICNTYP = iFlags(h_IFLAGS_contraction)
      ISym   = iFlags(h_IFLAGS_symmetry)
      IBasis = iFlags(h_IFLAGS_basis)
      IDFGHI = iFlags(h_IFLAGS_spherical)
c   o if ([basis=special]) then [scroll ZMAT down to the basis set definitions]
c     WARNING: THIS IS A HACK!
c     Basis set information handling should be totally revamped.
      if (iBasis.eq.0) then
         bTmp = .true.
         do while (bTmp)
            read(LuZ,'(A)',ERR=8500,END=8800) wrk
            bTmp = (wrk(2:2).ne.':'.and.wrk(3:3).ne.':')
         end do
c      o go back for the next read
         backspace(luz)
      end if
c YAU : end
      If (IContl.eq.0) IContl = 4
      IF (IPrnt.ge.100) THEN
         Write (LuOut, '(4X,A,I3)') 'Print Level: ', IPrnt
         Write (LuOut, '(4X,A,I3)') 'Step-taking: ', INR
         Write (LuOut, '(4X,A,I3)') 'E-Vec to follow: ',IVec
         Write (LuOut, '(4X,A,I3)') 'Die: ',IDie
         Write (LuOut, '(4X,A,I3)') 'Use curvilinear xform: ',ICurvy
         Write (LuOut, '(4X,A,I3)') 'Max Step: ',IMxStp
         Write (LuOut, '(4X,A,I3)') 'Step Contol: ', IStCrt
         Write (LuOut, '(4X,A,I3)') 'Vibrational Analysis: ', IVib
         Write (LuOut, '(4X,A,I3)') 'Recalculation of Hessian: ',IRecal
         Write (LuOut, '(4X,A,I3)') 'Convergence cutoff: ',IContl
         Write (LuOut, '(4X,A,I3)') 'Integral type: ',IntTyp
      END IF

C     Now read the basis sets - each line corresponds to a (non-dummy)
C     center in the Z-matrix, in the order of appearance.  We put this
C     into an array length NAtms, with blanks for dummies.  This makes
C     it easier to match the basis with the right coordinates later.
C     NOTE: blank line after basis names.
C     
C Having to type the same basis set name for every atom is error prone
C and time consuming for when the molecule is large (> 10 atoms). In
C order to aviod that this new logic will ask the user to input the 
C number of distinct atoms in the molecule and the use specify the 
C basis for those distinct atoms and the rest will get generated 
C automatically. If the user wants to have different basis sets for
C distinct atoms then he/she has to specify the basis set for each
C line correspond to a center in the Z-matrix. The distinct atom 
C basis pair do not need need to be in particular order. The 
C same applies for ECP basis sets. Ajith Perera, 06/2013.
C
C     
      If (Ibasis .EQ. 0) Then
          backspace(luz)
          read(LuZ,'(i4)',err=8500,end=8510,iostat=IOS) Num_unique
      Else
          Num_unique = 0
      Endif 

      If (Num_unique .EQ. Natms .OR. Num_unique .EQ. 0) Then

          Do 250 i = 1, NAtms
             BasNam(i) = BlnkBN
             if (AtNr(i).ne.0) then
               if (iBasis.eq.0) then

c   o read special
                 read(LuZ,'(A)',err=8500,end=8510,iostat=IOS) BasNam(i)

               else
C
csb 1/97 Allow arbitrary basis set names

                  BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
     &                      BasName(1:linblnk(BasName))
C
csb               izz=index(Basis(IBasis),'*')
csb               ixx=index(Basis(IBasis),'**')
csb               ipop=0
csbCJDW 8/29/94. 1.5 ---> 4.5 (what about He ?)
csb               if(ixx.ne.0.and.atmass(i).gt.4.5)ipop=1
csb               if(izz.ne.ixx.and.atmass(i).lt.4.5)ipop=1
csb               BasNam(i)=ZSym(i)(1:linblnk(ZSym(i)))//':'//
csb     &                   Basis(IBasis)(1:linblnk(Basis(IBasis))-ipop)
               end if

c        end if ([not dummy atom])
             end if
 250      Continue
     
      Else

          Do Iatm = 1,  Num_unique
             read(LuZ,'(A)',err=8500,end=8510,iostat=IOS) 
     &                      BasNam_hold(Iatm)
          Enddo
C
          Katm  = 0
          Do Iatm = 1, NAtms
             If (AtNr(Iatm).EQ.0) Katm = Katm + 1

             Do Jatm = 1, Num_unique
                Ilen   = Index(Basnam_hold(Jatm), ":")
                Symbol = Basnam_hold(Jatm)(1:Ilen-1)
                JAtnr  = ATOMNUMB(Symbol)
                If (AtNr(Iatm) .EQ. JAtNr) Then
                    Katm = Katm + 1
                    Basnam(Katm) = Basnam_hold(Jatm)
                Endif
             Enddo
          Enddo

      Endif 

      If(IBasis.eq.0)Read (LuZ,'(A)', ERR=8500, END=8510) Wrk

C ADDITION (JFS,4/90)
C Now take care of situation which occurs when 2--1--3 Z-matrix
C  specification and nonstandard basis set input is used.  In this case,
C  the basis sets for atoms #1 and #2 must be switched.
C
       CALL GETREC(-1,"JOBARC","12SWITCH",IONE,IAMBAD)
       BAD123 = (IAMBAD .EQ. 1)

      IF (BAD123 .AND. IBASIS .EQ.0 .and. INEWFD .NE. 0) then
          ScrBas=BasNam(1)
          BasNam(1)=BasNam(2)
          BasNam(2)=ScrBas
      ENDIF 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The basis set definitions"
      do i=1,natms
         if ( Atnr(i).ne.0) Write(6,"(1x,i3,1x,A20)") 
     +                      Atnr(i), basnam(i)
      enddo
#endif
C
C DEAL WITh ECPs
C
      IF (IECP.EQ.1) THEN

csb        IF (IBASIS.EQ.0) THEN
C
C  NOW READ THE NAMES OF THE ECPS. JODA TRIES TO READ THEM FROM
C  THE ZMAT FILE. THE USE OF ECPS REQUIRES "BASIS = SPECIAL"
C
        Read (LuZ,'(A)', ERR=8500, END=8510) Wrk
        If (Wrk(2:2) .EQ. ":" .OR. Wrk(3:3) .EQ. ":") Then
           backspace(Luz)
           Num_unique = 0
        Else
           backspace(Luz)
           read(LuZ,'(i4)',err=8500,end=8510,iostat=IOS) Num_unique
        ENdif 

        If (Num_unique .EQ. Natms .OR. Num_unique .EQ. 0) Then
           DO 255 I = 1, NATMS
              ECPNAM(I) = BLNKBN
              IF (ATNR(i) .NE. 0 .AND. IBASIS .EQ. 0) THEN
                 READ (LUZ, '(A)', ERR=8500, END=8510, IOSTAT=IOS)
     &           ECPNAM(I)
csb
              ElseIf(AtNr(i) .ne. 0 .and. IBasis .ne. 0)then

                 EcpNam(i)=BasName

              ENDIF
 255       CONTINUE

        Else

          Do Iatm = 1,  Num_unique

              read(LuZ,'(A)',err=8500,end=8510,iostat=IOS)
     &                   Ecpnam_hold(Iatm)
          Enddo
C
          Katm = 0
          Do Iatm = 1, NAtms
             If (AtNr(Iatm).EQ.0) Katm = Katm + 1
             Do Jatm = 1, Num_unique
                Ilen   = Index(Basnam_hold(Jatm), ":")
                Symbol = Ecpnam_hold(Jatm)(1:Ilen-1)
                JAtnr  = ATOMNUMB(Symbol)
                If (AtNr(Iatm) .EQ. JAtNr) Then
                   Katm  = Katm + 1
                   Ecpnam(Katm) = Ecpnam_hold(Jatm)
                Endif
             Enddo
          Enddo

        Endif

        IF(IBASIS.EQ.0)READ (LUZ,'(A)', ERR=8500, END=8510) WRK

C ADDITION (JFS,4/90)
C Now take care of situation which occurs when 2--1--3 Z-matrix
C specification and nonstandard basis set input is used.  In this case,
C the basis sets for atoms #1 and #2 must be switched.
C
       If(BAD123 .AND. IBASIS .eq.0 .AND. INEWFD .NE. 0) then
          SCRBAS=ECPNAM(1)
          ECPNAM(1)=ECPNAM(2)
          ECPNAM(2)=SCRBAS
        ENDIF

#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The ECP definitions" 
      do i=1,natms
         if ( Atnr(i).ne.0) Write(6,"(1x,i3,1x,A20)") 
     +                      Atnr(i), Ecpnam(i)
      enddo
      Write(6,*)
#endif

C
C NOW CORRECT THE CHARGES OF THE ATOMS BY THE CORE-ELECTRONS
C DESCRIBED BY THE ECP
C
        CALL CORCHR(ECPNAM,ATNR,ICRCOR,NATMS)
        DO 94 I=1,NATMS
          IF (ATNR(I).NE.0) THEN
            ICRCOR(i)=ATNR(I)-ICRCOR(i)
          ENDIF
   94   CONTINUE

C The following block of code is needed for inputs that has 
C dummy atoms and basis=special specification along with 
C number of unique atoms. 

        NREAL = 0
        DO IATMS = 1, NATMS
           IF (ATNR(IATMS) .NE. 0) THEN
              NREAL = NREAL + 1
              ECPNAM_HOLD(NREAL) = ECPNAM(IATMS)
           ENDIF
        ENDDO 
        CALL PUTREC(10,'JOBARC','ECPNAM',80*NREAL,ECPNAM_HOLD)
        CALL PUTREC(10,'JOBARC','NATOLD',1,NREAL)
C
csb        ELSE
csb          WRITE(*,*) 'INCOMPATIBLE KEYWORDS IN MKVMOL!'
csb          CALL ERREX
csb        ENDIF
      ENDIF
C
CCH END OF READING ECP NAMES

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Let's find out which centers are equivalent
C
      Call SymEqv ( NAtms, GenBy)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Count number of centers having *different* atomic numbers in this grou
C   Top half of Nord contains the number of symmetry equivalent atoms
C   for a given atomic number.
C
      Do I = 1, NAtms
         If (GenBy(i).eq.0.and.AtMass(i).gt.1.D-2) Then
            J = J + 1
            Scratch(J) = AtMass(I)
            Nord(J) = I
         End If
      End Do
      If (IVIB.ne.2) Then
c         Call PikSr2(J,Scratch,Nord)
         NonTyp = 1
         If (J.Eq.1) Then
            Nord(MxAtms+1) = 1
         Else
            Do I = 2, J
               Dif = Dabs(Scratch(i)-Scratch(i-1))
               If (NonTyp.eq.1) IBot=0
c               If (Dif.lt.1.D-14) GoTo 314
               Nord(MxAtms+Nontyp) = I - (IBot+1)
               IBot = IBot + Nord(MxAtms+Nontyp)
               NonTyp = NonTyp + 1
 314           If (I.eq.J) Nord(MxAtms+Nontyp) = I - Ibot
            End Do
         End If
      Else If (ivib.eq.2) Then
C
C For finite difference calcs, reset this so that it writes out the
C  atoms to the MOL
C  file in exactly the same order as they appear in the ZMAT.  This
C  ordering is assumed in FINDIF!
C
         NonTyp = J
         Do I = 1, J
            Nord(MxAtms+I) = 1
         End Do
      End If

C Now get symmetry information written by SYMMETRY to file VMINF.
      Opn = .False.
      OldLu = 0
      IF (PTGRP.EQ.'C1 ') THEN
         INQUIRE(FILE='VMLSYM',EXIST=Opn)
         IF (OPN) THEN
            OPEN(UNIT=30,FILE='VMLSYM',FORM='UNFORMATTED',
     &           STATUS='OLD')
            CLOSE(UNIT=30,STATUS='DELETE')
         END IF
         GOTO 5400
      END IF
      Inquire (FILE='VMLSYM', OPENED=Opn, NUMBER=OldLu)
      If (Opn) Close (OldLu)
C
      Open (30, FILE='VMLSYM', STATUS='OLD',FORM='UNFORMATTED')
      Rewind(30)
      Read (30,Err=9400) NSymOp, (KA(i),i=1,3)
      Close(30, Status='Keep')
5400  CONTINUE
      If (PtGrp.EQ.'C1 ') Then
         Do I = 1, 3
            KA(I)='   '
         End Do
         NSymOp = 0
      End If
#ifdef _DEBUG_LVL0
       Write(6,"(a,3(1x,I3))") "The nord array: ", (Nord(I),i=1,Natms)
#endif 
C
C     Now we've got (pretty much) all the information we need, so
C     we can start writing out the MOL file.
C     Begin with a very cautious opening of the file.
C
        Opn = .False. 
        Inquire (file=molfil, OPENED=Opn)
        If ( .not. Opn ) Open (LuVMol, FILE=molFil, FORM='FORMATTED',
     $     STATUS='UNKNOWN')
        Rewind LuVMol
C
C Write lines of VMOL deck.  All of this mysterious SLASK stuff, etc.
C
        idosph=1-IDFGHI
cYAU        Write(LuVMol,'(A6,4X,10I5)')Slask,IGENER,If10,If18,IECP,IECPU,
cYAU     &    IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LuVMol,'(A6,4X,10I5)')
     &     'INTGRL',IGENER,If10,If18,IECP,IECPU,
     &              IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LuVMol,'(8x,a)')
     &     ' *** ACES2 Program System (Release V0.1) *** '
        Write(LuVMol,'(A)',ERR=8300) JobTtl(:linblnk(JobTtl))
CJDW  11/2/94
        TLA = DBLE(10.0D+00**(-IFLAGS2(108)))
        Write(LuVMol,'(2I5,3A3,1X,E10.2,10X,2I5)')NonTyp,NSymop,Ka,Tla,
     &    ID3,IRstrt
        Zink=9999.0
        Zink2=3.0
        Write(LuVMol,'(2F10.2)',ERR=8300)Zink,Zink2
      if(inttyp.eq.4) then
        seward = .true.
        opn = .false.
        inquire(file=sewfil,opened=opn)
        if (.not. opn)open(lusew,file=sewfil,form='FORMATTED',
     &    status='UNKNOWN')
        rewind lusew
        write(lusew,'(A)')' &SEWARD  &END'
        write(lusew,'(A)')'ACES'
        write(lusew,'(A)')'TITLE'
        write(lusew,'(A)') JobTtl(:linblnk(JobTtl))
        write(lusew,'(A)')'MOLECULE'
        write(lusew,'(A)')'THRESHOLD'
        write(lusew,'(E10.2)') DBLE(10.0D+00**(-IFLAGS2(108)))
c kjw 7-27-00 bug fix for cases with no symmetry
        if(ka(1).eq.'   '.and.ka(2).eq.'   '.and.ka(3).eq.'   ')then
        else
          write(lusew,'(A)')'SYMMETRY'
          write(lusew,'(3a3)') Ka
        endif
      endif

c   o open a basis set file for rdbas and create ZMAT.BAS if necessary
      inquire(unit=12,opened=bOpened)
      inquire(file='ZMAT.BAS',exist=bExist)
      if (bExist) then
         if (.not.bOpened) then
            open(unit=12,file='ZMAT.BAS',form='formatted',status='old')
         end if
         bCpBasis = .false.
         rewind(12)
      else
         if (.not.bOpened) then
            open(unit=12,file='GENBAS',form='formatted',status='old')
         else
            print *, '@MKVMOL: Assertion failed.'
            print *, '         ZMAT.BAS does not exist but unit 12 is ',
     &               'already open.'
            call errex
         end if
         open(unit=13,file='ZMAT.BAS',form='formatted',status='new')
         rewind(13)
         bCpBasis = .true.
      end if

C
C     Charge & multiplicity
C
c      Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
C      Read (LuZ, '(A)', ERR=8500, END=8510) Wrk
C
C Now write the "real" input.
C
      IINFIL=0
      SHOFFSET(1)=1
      PROFFSET(1)=1
      Do 318 I = 1, Nontyp
       NStart = 1
       Do J = 1, I-1
          NStart = NStart + Nord(MxAtms+J)
       End Do
c       If(IntTyp.Eq.1)Then
c        Write(LuVMol,7000,ERR=8300)
c     &  NullSt,Float(IAtnum(Nord(NStart))),Nord(MxAtms+I)
c       Elseif(IntTyp.Eq.2)Then
c        WRITE(LUVMOL,7001,ERR=8300)
c     &  Float(IAtnum(Nord(NStart))),Nord(MxAtms+I)
c       ENDIF
c7000   Format(A6,F14.8,I5,$)
c7001   Format(F10.1,I5,$)
       SCRBAS=BASNAM(NORD(NSTART))
       IATTYP=1
       IIATOM=NORD(NSTART)
       IATOM=IATNUM(IIATOM)
       IF(IATOM.GT.2)IATTYP=2
       IF(IATOM.GT.10)IATTYP=3
       IF(GENBAS)THEN
        CALL RDBAS(SCRBAS,0,IFLAGS(64+IATTYP),LUVMOL,
     &             IntTyp,NullSt,IAtnum(Nord(NStart)),
     &             Nord(MxAtms+I),seward,.false.,13)
       ELSE
        Call WrPBas (LuBasL,BasFil,SCRBAS,
     &  LuVMol,0,1,IStat,IntTyp,Nullst,IAtnum(Nord(Nstart)),
     &  Nord(MxAtms+I))
       ENDIF
       Do J = 1, Nord(MxAtms+J)
          IBotm = -2 + 3*Nord(NStart-1+J)
          IINFIL = IINFIL + 1
          IZLOC(IINFIL) = NORD(NSTART+J-1)
          Icount = Icount + 1
          Write(LuVMol,'(A2,A2,3F20.12)',ERR=8300)
     &                  ZSym(Nord(NStart)),Indxx(ICOUNT),
     &                  (geom(IJ),IJ=IBotm,Ibotm+2)
       End Do
       IF(GENBAS)THEN
        CALL RDBAS(SCRBAS,1,IFLAGS(64+IATTYP),LUVMOL,
     &             ijunk,ijunk,ijunk,ijunk,seward,bCpBasis,13)
        if(inttyp.eq.4)then
          if(iflags(62).eq.0)then
            write(lusew,'(A)')'Cartesian all'
          else
            write(lusew,'(A)')'Spherical all'
          endif      
          wrk=zsym(nord(nstart))
          isew=index(wrk,' ')
          write(lusew,'(a6,3f20.12)') wrk(1:isew-1)//indxx(icount),
     &      (geom(isew),isew=ibotm,ibotm+2)
          write(lusew,'(A)')'End of Basis Set'
        endif
        SHOFFSET(I+1)=IOFFSHEL
        PROFFSET(I+1)=IOFFPRIM
        SHELLORB(I)=SHOFFSET(I+1)-SHOFFSET(I)
        PRIMORBT(I)=PROFFSET(I+1)-PROFFSET(I)
       ELSE
        Call WrPBas (LuBasL,BasFil,SCRBAS,
     &  LuVMol,1,1,IStat,ijunk,ijunk,ijunk,ijunk)
       ENDIF
C
C           Handle errors in WrPBas
C
            If (Mod(IStat,2) .eq. 1) then
               IStat = 1
               Write (LuErr, 9910)
               Close (LuZ)
               Close (LuVMol)
               Return
            ElseIf (IStat .eq. 8) then
               IStat = 3
               Write (LuErr, 9935) AbiFil, LuVMol
               Close (LuZ)
               Close (LuVMol)
               Return
            EndIf
 318    Continue
        Write(LuVMol,'(A6)',ERR=8300)'FINISH'
        if(inttyp.eq.4)then
          write(lusew,'(A)')'End of input'
        endif

c   o if we created ZMAT.BAS, then close GENBAS and connect ZMAT.BAS to unit 12
      if (bCpBasis) then
         close(12,status='keep')
         close(13,status='keep')
         open(unit=12,file='ZMAT.BAS',form='formatted',status='old')
         bCpBasis = .false.
      end if

c
c create the interfacefile to TURBOMOLE
c
      if (iturbo.eq.1) then
cYAU        ispher=iflags(62)
cYAU        call mkturb(basnam,geom,zsym,atmass,iatnum,natms,nontyp,
cYAU     &              ispher)
         write(*,*) '@MKVMOL: MkTurb was temporarily removed.'
         write(*,*) '         This calculation cannot complete.'
         call errex
      end if
C
C     Dump vector to disk which relates the atoms listed in the
C     MOL file to their position in the user supplied Z-matrix.
C
      CALL PUTREC(20,'JOBARC','MPVMZMAT',NATOMS,IZLOC)
#ifdef _DEBUG_LVL0
      Write(6,"(2a)") " The MPVMZMAT record (symmetry unique only,",
     &                " others are zero)"
      Write(6,"(6(1x,i3))") (Izloc(i),i=1,Natoms)
#endif 
C
 9910 Format (' @MKDECK-F, Dependent terminated with FATAL error.')
 9935 Format (' @MKDECK-F, Dependent reports error accessing file ',
     $   A,' on unit ',I3,'.')
C
C This better not happen.
C
      If (Nontyp .eq. 0) then
         IStat = 7
         Write (LuErr, 9970)
         Close (LuZ)
         Close (LuVmol)
         Return
 9970    Format (' @MKDECK-F, No non-dummy centers in Z-matrix.')
      EndIf
C
C     The rest of this just copied from ZMAT
C
C
C     Make sure every body gets closed & stuff
C
 8000 If(IntTyp.ne.0)Close (LuVMol)
      if(inttyp.eq.4)close (lusew)
      Close (LuZ)
C
C WRITE BASIS SET INFORMATION TO JOBARC
C
      IONE=1
      NTOTPRIM=IOFFPRIM-1
      NTOTSHEL=IOFFSHEL-1
      IF(IECP.EQ.1) THEN
       CALL PUTREC(20,'JOBARC','ATOMCHRG',NATOMS,ICRCOR)
       CALL PUTREC(20,'JOBARC','ATCHRORG',NATOMS,IATNUM)
      ELSE
       CALL PUTREC(20,'JOBARC','ATOMCHRG',NATOMS,IATNUM)
      ENDIF
      IF(INEWFD.EQ.0)THEN
       CALL PUTREC(20,'JOBARC','NTOTSHEL',IONE,NTOTSHEL)
       CALL PUTREC(20,'JOBARC','NTOTPRIM',IONE,NTOTPRIM)
       CALL PUTREC(20,'JOBARC','BASISEXP',NTOTPRIM*IINTFP,BASISEXP)
       CALL PUTREC(20,'JOBARC','BASISCNT',NTOTPRIM*IINTFP,BASISCNT)
       CALL PUTREC(20,'JOBARC','SHELLSIZ',NTOTSHEL,SHELLSIZ)
       CALL PUTREC(20,'JOBARC','SHELLPRM',NTOTSHEL,SHELLPRM)
       CALL PUTREC(20,'JOBARC','SHELLANG',NTOTSHEL,SHELLANG)
       CALL PUTREC(20,'JOBARC','SHELLLOC',NTOTSHEL,SHELLLOC)
       CALL PUTREC(20,'JOBARC','SHOFFSET',NONTYP,SHOFFSET)
       CALL PUTREC(20,'JOBARC','SHELLORB',NONTYP,SHELLORB)
       CALL PUTREC(20,'JOBARC','PROFFSET',NONTYP,PROFFSET)
       CALL PUTREC(20,'JOBARC','PRIMORBT',NONTYP,PRIMORBT)
      ENDIF
C
CJDW  6/ 6/95. Add JFS call to SHELLINF.
C
      call shellinf(nontyp,natoms,ntotshel,shellorb,shoffset,shellang,
     &              shellsiz)
C
      Return
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     Get here via an I/O error on LuAbi or LuVMol
C
 8300 IStat = 3
      If(IntTyp.ne.0)then
       Write (LuErr, 9930) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9930) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 9930 Format (' @MKDECK-F, I/O error on file ',A,'.')
C
C     Get here via an I/O error on VMLSYM
C
 9400 IStat=3
      Write (LuErr, 9930)'VMLSYM'
      Close(LuVMol)
      Close(LuZ)
      Return
C
C     Get here via an I/O error on LuZ
C
 8500 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9950) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9950) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 8510 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9955) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9955) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
C
C     If the ZFil ends after the Z-matrix or the optimization control
C     info, it probably means that JODA is being used "stand-alone"
C     and we want to allow for that possibility.
C
 8800 IStat = 8
      If(IntTyp.ne.0)then
       Write (LuOut, 9800) ZFil, MolFil
       Close (LuVMol)
      else
       Write (LuOut, 9800) ZFil, AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 8810 IStat = 10
      If(IntTyp.ne.0)then
       Write (LuOut, 9810) ZFil, MolFil
       Close (LuVMol)
      else
       Write (LuOut, 9810) ZFil, AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return
 9950 Format (' @MKDECK-F, I/O error on file ',A,'.')
 9955 Format (' @MKDECK-F, Premature end-of-file on ',A,'.')
 9800 Format (/1X,A,' ends after Z-matrix - cannot ',
     $   'finish making ',A,'.')
 9810 Format (/1X,A,' ends after JODA control info - cannot ',
     $   'finish making ',A,'.')
      End
