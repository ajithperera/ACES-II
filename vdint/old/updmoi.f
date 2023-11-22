

      SUBROUTINE UPDMOI(LSTDIS,DISSIZ,LSTSPN,LSTNUM,IENTER,IOFF)
C
C UPDATES THE MOIO AND MOIOWD VECTORS WHEN YOU WRITE A NEW LIST.
C
C     LSTDIS= NUMBER OF DISTRIBUTIONS IN THE LIST.
C     DISSIZ= THE SIZE OF THE INDIVIDUAL DISTRIBUTIONS (IN FP WORDS).
C     LSTSPN= "SPIN CASE" FOR LIST [I OF MOIO(I,J)].
C     LSTNUM= NUMBER FOR LIST [J OF MOIO(I,J)].
C     IENTER= 0 UNLESS THE MOINTS FILE IS NEW.  1 IS USED THEN.
C              AT THE END OF EXECUTION OF A LINK, CALL THIS
C              WITH A '2', AND IT WILL WRITE OUT TOTREC TO JOBARC.
C              AT BEGINNING, A '3' WILL INITIALIZE THE VALUE.
C     IOFF  = -1 BEGINS LIST ON A PHYSICAL RECORD BOUNDARY.
C                ANY OTHER VALUE BEGINS AT FIRST AVAILABLE WORD.
C
CEND
      IMPLICIT INTEGER (A-Z)
      CHARACTER*8 NAMES(5)
      CHARACTER*80 FNAME
      LOGICAL YESNO,YESNO1,YESNO2
      DIMENSION TOTREC(5),TOTWRD(5),ITOPRC(5),JUNK(8192)
      COMMON / / ICORE(1)
      COMMON /LISTS/ MOIO(10,500),MOIOWD(10,500),MOIOSZ(10,500),
     &               MOIODS(10,500),MOIOFL(10,500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /MACHSP2/MASK1,MASK2,ISHFSZ
      COMMON /FILSPC/ ILNBUF,IPRCLN,IPRCWD
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /IOPOS/ ICRSIZ,ICHCSZ,IOFFX(2),LENREC
      COMMON /CACHEINF/ CACHSIZ,CACHSZP1,CACHDIR(100),CACHPOS(100),
     &                  CACHFILE(100),CACHMOD(100),OLDEST
      DATA NAMES /'MOINTS  ','GAMLAM  ','MOABCD  ','DERINT  ',
     &            'DERGAM  '/
      SAVE TOTREC,TOTWRD
C
      IPACK(I,J)=OR(J,ISHFT(I-49,ISHFSZ))
      UPACKR(I) =AND(I,MASK1)
      UPACKF(I) =AND(ISHFT(I,-ISHFSZ),MASK2)+49
      CALL IZERO(JUNK,8192)
      DO 1001 I=1,5
       ITOPRC(I)=TOTREC(I)
1001  CONTINUE
      IFIVE=5
      IF(IENTER.EQ.2)THEN
       CALL PUTREC(20,'JOBARC','TOTRECMO',IFIVE,TOTREC)
       CALL PUTREC(20,'JOBARC','TOTWRDMO',IFIVE,TOTWRD)
       RETURN
      ENDIF
      IF(IENTER.EQ.3)THEN
       CALL GETREC(20,'JOBARC','TOTRECMO',IFIVE,TOTREC)
       CALL GETREC(20,'JOBARC','TOTWRDMO',IFIVE,TOTWRD)
       RETURN
      ENDIF
      ITYPE=1+(LSTNUM-1)/100
      IF(IENTER.EQ.1)THEN
       TOTREC(ITYPE)=1
       TOTWRD(ITYPE)=0
       ITOPRC(ITYPE)=0
      ENDIF
c
c  updmoi should not be used to change size of an exisiting list
c
      IF ( MOIO(LSTSPN,LSTNUM) .NE. 0 ) then
         if ( MOIOSZ(LSTSPN,LSTNUM).eq.DISSIZ .and.
     $        MOIODS(LSTSPN,LSTNUM) .eq. LSTDIS ) then
            RETURN
         else
            write (6,*) 'tried to change size of exisiting list: ',
     $         lstspn, ',', lstnum
            write (6,*) 'current: distributions -- ',
     $         moiods(lstspn,lstnum), ' size -- ', moiosz(lstspn,lstnum)
            write (6,*) 'requested: distributions -- ',
     $         lstdis, ' size -- ', dissiz
            call errex
         endif
      endif

C
C ASSIGN MOIO(LSTSPN,LSTNUM) AND MOIOWD(LSTSPN,LSTNUM) AND THE FILE NUMBER
C
      IFILE=ITYPE+49
      CALL GFNAME(NAMES(ITYPE),FNAME,ILENGTH)
      INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO1,OPENED=YESNO2)
      IF(.NOT.YESNO1.OR..NOT.YESNO2)THEN
       CALL OPNFIL(ITYPE)
      ENDIF
      MOIOSZ(LSTSPN,LSTNUM)=DISSIZ
      MOIODS(LSTSPN,LSTNUM)=LSTDIS
      MOIOFL(LSTSPN,LSTNUM)=IFILE
      IF(IOFF.EQ.-1)THEN
C
C LISTS STARTS ON A PHYSICAL RECORD BOUNDARY.
C
       MOIOWD(LSTSPN,LSTNUM)=1
       TOTREC(ITYPE)=TOTREC(ITYPE)+1
       MOIO(LSTSPN,LSTNUM)=TOTREC(ITYPE)
       TOTWRD(ITYPE)=0
      ELSE
       MOIOWD(LSTSPN,LSTNUM)=TOTWRD(ITYPE)+1
       IF(TOTWRD(ITYPE).EQ.IPRCWD)THEN
        MOIOWD(LSTSPN,LSTNUM)=1
        TOTREC(ITYPE)=TOTREC(ITYPE)+1
       ENDIF
       MOIO(LSTSPN,LSTNUM)=TOTREC(ITYPE)
      ENDIF
C
C COMPUTE TOTAL NUMBER OF WORDS IN LIST.
C
      NWORDS=LSTDIS*DISSIZ*IINTFP
C
C COMPUTE THE NUMBER OF FULL RECORDS THIS WILL TAKE UP.  USUALLY THIS
C   WILL BE ZERO.  ALSO COMPUTE HOW MANY WORDS OF A PARTIAL RECORD ARE
C   REQUIRED (OFTEN JUST NWORDS).
C
      NFULL=NWORDS/IPRCWD
      NPART=MOD(NWORDS,IPRCWD)
C
C INCREMENT TOTREC BY THE NUMBER OF FULL RECORDS THAT THIS WILL TAKE UP.
C
      TOTREC(ITYPE)=TOTREC(ITYPE)+NFULL
C
C COMPUTE POSITION OF LAST WORD WRITTEN IN PARTIAL RECORD.  IF THE LOGIC
C  RECORD GOES ACROSS A PHYSICAL RECORD BOUNDARY, INCREMENT TOTREC BY ON
C
      NLEFT=IPRCWD-TOTWRD(ITYPE) 
      IF(NLEFT.GE.NPART)THEN
       TOTWRD(ITYPE)=TOTWRD(ITYPE)+NPART
      ELSE
       TOTWRD(ITYPE)=NPART-NLEFT
       TOTREC(ITYPE)=TOTREC(ITYPE)+1
      ENDIF
C
C WRITE AN EMPTY BUFFER TO RECORD RECORD IF IT IS BEYOND EOF.
C
      IF(TOTREC(ITYPE).GT.ITOPRC(ITYPE))THEN
C
C FLUSH CURRENT CONTENTS OF CACHE #1
C
       IFILE0=CACHFILE(1)
       IREC  =UPACKR(CACHDIR(1))
       IF(IREC.NE.0.AND.IFILE0.NE.0)THEN
        CALL ACES_IO_WRITE(IFILE0,IREC,ICORE(CACHPOS(1)),ICHCSZ)
       ENDIF
       CALL IZERO(ICORE(CACHPOS(1)),ICHCSZ)
       CALL ACES_IO_WRITE(IFILE,TOTREC(ITYPE),ICORE(CACHPOS(1)),ICHCSZ)
       CACHDIR(1)=0
       CACHFILE(1)=0
      ENDIF
      RETURN
      END