
C     THIS IS A MODIFICATION OF THE MKVMOL.F
C     FILE. THE PURPOSE IS TO WRITE A .MOL-LIKE
C     FILE WHICH INCLUDES ALL NON-DUMMY
C     ATOMS AND NOT JUST THE SYMMETRY UNIQUE CENTERS.
C     THIS ROUTINE IS ONLY EXECUTED WHEN THE
C     USER HAS CHOSEN AN NDDO INITIAL GUESS
C     -CARLOS TAYLOR

#include "flags.h"

      SUBROUTINE MKERDOED(geom, PTGRP,NATMS,NUNIQUE,ZSYM,ATNR,GENBY,
     $     SCRATCH, ISTAT, BasName)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

#include "baslims.par"
#include "mxatms.par"
#include "io_units.par"
#include "linelen.par"
#include "molstrct.com"
#include "coord.com"

C SET UNIT NUMBER CARLOS TAYLOR
      parameter(LUERDOED=40,LUERDOED_INPUT=50,ATOB=.529177d0)
      character*10 NddoFil
      character*4 CREF(3),CREFJ
      data CREF/' RHF',' UHF','ROHF'/,NddoFil/"NDDO.INPUT"/
      character*(*) BasName

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
      character*(baslen) BasNam(MxAtms),BlnkBN,ScrBas,EcpNam(MxAtms),
     &                   EcpNam_hold(MxAtms), BasNam_hold(MxAtms)

      dimension icrcor(mxatms)
      common /turbo / iturbo,matom,ioffsh

CJDW  11/2/94. We now use IFLAGS2(108) to set TLA.
C     Parameter (TLA = 1.D-9)
CSSS      INTEGER SHELLANG,SHELLLOC,SHELLSIZ,SHELLPRM
CSSS      INTEGER SHELLORB,SHOFFSET,PRIMORBT,PROFFSET
      double precision geom(3*mxatms)
C
      LOGICAL XYZIN,NWFINDIF
      Logical Opn,bad123,GENBAS 
CKJW 5-24-00
      logical seward,bCpBasis,bExist,bOpened
      data seward /.false./
      Integer OldLu
      Dimension Nord(2*MxAtms),Scratch(NAtms)
      Integer atomnumb 
C
CSSS      COMMON /MOLSTR1/ SHELLANG(100),SHELLLOC(100),SHELLSIZ(100),
CSSS     &                 SHELLPRM(100),BASISEXP(10000),BASISCNT(10000)
CSSS      COMMON /MOLSTR2/ NTOTPRIM,NTOTSHEL,IOFFPRIM,IOFFSHEL,IIATOM
CSSS      COMMON /MOLSTR3/ SHELLORB(100),SHOFFSET(100),PRIMORBT(100),
CSSS     &                 PROFFSET(100)

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
C
C     Cautiously open the ZMAT file
C
      Inquire (FILE=ZFil, OPENED=Opn)
      If (.not.Opn) Open (LuZ, FILE=ZFil, STATUS='OLD')
      Rewind LuZ

C     skip the header
      bTmp = .true.
      do while (bTmp)
         read (LuZ, '(A)') JobTtl
         call parsez(JobTtl,izl)
         i = izl(1,1)
         bTmp = ((i.eq.0).or.(JobTtl(i:i).EQ.'%'))
      end do

      IF ((.NOT.XYZIN).OR.(NWFINDIF)) THEN
         i = 1
         do while (i.le.natms)
            read (LuZ,'(A)') Wrk
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
            READ (LUZ,'(A)') WRK
         END DO
      END IF
      Read (LuZ,'(A)') Wrk
      IF (IFLAGS(68) .EQ. 0) THEN
         Do i = 1, NUnique
            Read (LuZ, '(A)') Wrk
         End Do
      END IF
      Read (LuZ, '(A)') Wrk
C
C I think this is the best way to avoid the problem of Cartesian
C frequency finite difference runs when the first two atoms are
C not the same. This problem must have existed from the very
C begining. It get unnoticed because we do not do lots and lots
C of Cartesian Frqs (and then of course the first two atoms are
C the same, then no one notice it). 01/2016, Ajith Perera.

      IF (XYZIN) BAD123 = .FALSE.

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
            read(LuZ,'(A)') wrk
            bTmp = (wrk(2:2).ne.':'.and.wrk(3:3).ne.':')
         end do
c      o go back for the next read
         backspace(luz)
      end if
c YAU : end
      If (IContl.eq.0) IContl = 4

C     Now read the basis sets - each line corresponds to a (non-dummy)
C     center in the Z-matrix, in the order of appearance.  We put this
C     into an array length NAtms, with blanks for dummies.  This makes
C     it easier to match the basis with the right coordinates later.
C     NOTE: blank line after basis names.
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

       If(BAD123 .AND. IBASIS .EQ.0 .AND. INEWFD .NE. 0) then
          ScrBas=BasNam(1)
          BasNam(1)=BasNam(2)
          BasNam(2)=ScrBas
       EndIf

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The basis set definitions"
      do i=1,natms
         If (Atnr(i) .ne.0) Write(6,"(1x,i3,1x,A20)") 
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
C  specification and nonstandard basis set input is used.  In this case,
C  the basis sets for atoms #1 and #2 must be switched.
C
       If(BAD123 .AND. IBASIS .EQ. 0 .AND. INEWFD .NE. 0) then
          SCRBAS=ECPNAM(1)
          ECPNAM(1)=ECPNAM(2)
          ECPNAM(2)=SCRBAS
        ENDIF
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The ECP definitions"
      do i=1,natms
         If (Atnr(i) .ne.0) Write(6,"(1x,i3,1x,A20)") 
     +                      Atnr(i), basnam(i)
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

C
C Count number of centers having *different* atomic numbers in this grou
C   Top half of Nord contains the number of symmetry equivalent atoms
C   for a given atomic number.
C
      Do I = 1, NAtms
C Write MOL file for all the atoms not just symmetry unique ones
CSSD         If (GenBy(i).eq.0.and.AtMass(i).gt.1.D-2) Then
            J = J + 1
            Scratch(J) = AtMass(I)
            Nord(J) = I
CSSS         End If
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
      Close(30, Status='Delete')
5400  CONTINUE
      If (PtGrp.EQ.'C1 ') Then
         Do I = 1, 3
            KA(I)='   '
         End Do
         NSymOp = 0
      End If
C
C     Now we've got (pretty much) all the information we need, so
C     we can start writing out the MOL file.
C     Begin with a very cautious opening of the file.
        Opn = .False.
        Open (LUERDOED, FILE='MOL_ERDOED', FORM='FORMATTED',
     $     STATUS='UNKNOWN')
        Rewind LUERDOED
C
C Write lines of VMOL deck.  All of this mysterious SLASK stuff, etc.
C
        idosph=1-IDFGHI
cYAU        Write(LUERDOED,'(A6,4X,10I5)')Slask,IGENER,If10,If18,IECP,IECPU,
cYAU     &    IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LUERDOED,'(A6,4X,10I5)')
     &     'INTGRL',IGENER,If10,If18,IECP,IECPU,
     &              IPRINT,INUTMP,ISTATS,IDOSPH
        Write(LUERDOED,'(8x,a)')
     &     ' *** ACES2 Program System (Release V0.1) *** '
        Write(LUERDOED,'(A)') JobTtl(:linblnk(JobTtl))
CJDW  11/2/94
        TLA = DBLE(10.0D+00**(-IFLAGS2(108)))
        Write(LUERDOED,'(2I5,3A3,1X,E10.2,10X,2I5)')NonTyp,NSymop,Ka,
     &    TLA,ID3,IRstrt
        Zink=9999.0
        Zink2=3.0
        Write(LUERDOED,'(2F10.2)')Zink,Zink2
C
C Now write the "real" input.
C
      IINFIL=0
      SHOFFSET(1)=1
      PROFFSET(1)=1
      Do 318 I = 1, Nontyp
C
C     THE FOLLOWING IF-THEN SOLVES THE PROBLEM WITH DUMMY ATOMS - CARLOS TAYLOR
C
       if(atnr(i).eq.0)goto 318
       NStart = 1
       Do J = 1, I-1
          NStart = NStart + Nord(MxAtms+J)
       End Do
       SCRBAS=BASNAM(NORD(NSTART))
      IATTYP=1
       IIATOM=NORD(NSTART)
       IATOM=IATNUM(IIATOM)
       IF(IATOM.GT.2)IATTYP=2
       IF(IATOM.GT.10)IATTYP=3
C
C     I'M NOT SURE ABOUT THE FOLLOWING CHANGE LUERDOED FOR LUVMOL
C
        CALL RDBAS(SCRBAS,0,IFLAGS(64+IATTYP),LUERDOED,
     &             IntTyp,NullSt,IAtnum(Nord(NStart)),
     &             Nord(MxAtms+I),seward,.false.,6)
       Do J = 1, Nord(MxAtms+J)
          IBotm = -2 + 3*Nord(NStart-1+J)
          IINFIL = IINFIL + 1
          IZLOC(IINFIL) = NORD(NSTART+J-1)
          Icount = Icount + 1
          Write(LUERDOED,'(A2,A2,3F20.12)')
     &                  ZSym(Nord(NStart)),Indxx(ICOUNT),
     &                  (geom(IJ),IJ=IBotm,Ibotm+2)
       End Do
        CALL RDBAS(SCRBAS,1,IFLAGS(64+IATTYP),LUERDOED,
     &             ijunk,ijunk,ijunk,ijunk,seward,.false.,6)
        SHOFFSET(I+1)=IOFFSHEL
        PROFFSET(I+1)=IOFFPRIM
        SHELLORB(I)=SHOFFSET(I+1)-SHOFFSET(I)
        PRIMORBT(I)=PROFFSET(I+1)-PROFFSET(I)
 318    Continue
        Write(LUERDOED,'(A6)')'FINISH'
C
C Also open the "NDDO.INPUT" file for Carlos Taylor's NDDO program. Note
C that projected NDDO is used as the initial guess for the SCF program when
C it is requested by the user (GUESS=NDDO). A. Perera, 10/2004
C
C     Determine if RHF or UHF (IREF=1 FOR RHF AND IREF=2 FOR UHF)
      IREF=IFLAGS(11)+1
      IMULT=IFLAGS(29)
      ICHG=IFLAGS(28)
      OPEN(UNIT=LUERDOED_INPUT, FILE=NddoFil,Form="Formatted",
     &     STATUS="UNKNOWN")
      WRITE(LUERDOED_INPUT,'(A40)') "'NDDO INPUT GENERATED VIA ACES2'"
      IF(IREF.EQ.1)WRITE(LUERDOED_INPUT,151)ICHG
      IF(IREF.EQ.2)WRITE(LUERDOED_INPUT,150)IMULT,ICHG
      CALL GETREC(20, 'JOBARC', 'ATOMCHRG', NATOMS, ICRCOR)
      DO IATOMS = 1, NATOMS
         ISTART = 3*(IATOMS - 1)
         WRITE(LUERDOED_INPUT, 99) ICRCOR(IATOMS), GEOM(ISTART+ 1)*
     &                           ATOB, GEOM(ISTART+ 2)*ATOB,
     &                           GEOM(ISTART+3)*ATOB
 99      FORMAT(I2,F10.5,F10.5,F10.5)
 150     FORMAT( "'AM1 XYZ UHF MULT=",I1," CHARGE=",I1,"'")
 151     FORMAT( "'AM1 XYZ RHF CHARGE=",I1,"'")
      END DO
      CLOSE(LUERDOED)
      CLOSE(LUERDOED_INPUT)
      Close (LuZ)
      Return
 9930 Format (' @MKDECK-F, I/O error on file ',A,'.')
 9950 Format (' @MKDECK-F, I/O error on file ',A,'.')
C
 9400 IStat=3
      Write (LuErr, 9930)'VMLSYM'
      Close(LuVMol)
      Close(LuZ)
      Return

C Get here via an I/O error on LuZ

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

C Get here via an I/O error on LuZ

 8510 IStat = 5
      If(IntTyp.ne.0)then
       Write (LuErr, 9950) MolFil
       Close (LuVMol)
      else
       Write (LuErr, 9950) AbiFil
       Close (LuAbi)
      EndIf
      Close (LuZ)
      Return

      End
