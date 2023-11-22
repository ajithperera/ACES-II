      subroutine mkarhf1(amat,scratch,icore,maxcor,BrEDUNDAnt)
C
C   A(AI,BJ) = <AB//IJ> + <AJ//IB> + <Ab//Ij> + <Aj//Ib>
C
C            = <AB//IJ> - <JA//IB> + 2 <Ab//Ij>
      
      implicit double precision(a-h,o-z)
      integer dirprd,dissyw,pop,vrt,REFLIST
      LOGICAL bRedundant
      dimension amat(1),scratch(1),icore(1)
      common/info/nocco(2),nvrto(2)
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
c below common block is not consitent with rest of code
c need to be change I don't know what impact it might have on 
c rest of the code but i am doing it
       common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
       common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
c      COMMON/SYMPOP2/IRPDPD(8,22)
c      COMMON/SYMPOP/IRP_DM(8,22),ISYTYP(2,500),NTOT(18)
c      COMMON/SYM2/POP(8,2),VRT(8,2),NJUNK(6)
      common /dropgeo/ ndrgeo

c      COMMON /SHIFT/ ISHIFT 

c above ishift doesn't make any sence         
c i am removing it 

      DATA ONE,ONEM,TWO /1.0D+0,-1.0D+0,2.0D+0/

C   READ IN FIRST THE INTEGRALS <AB//IJ> (ORDERING AI,BJ) 
C   THEY ARE STORED ON LIST 19

c  below listw info has modifed commented one is the original
c one since we don't know the ishift info or we don't ishifht info
c thus it is remoced

c      LISTW=19  + ISHIFT 
c      end
      listw=19
      
      numsyw=irpdpd(1,isytyp(2,listw))
      dissyw=irpdpd(1,isytyp(1,listw))
   
C
c      IF(bRedundant) THEN
        CALL GETLST(AMAT,1,NUMSYW,2,1,LISTW)
c      ELSE
c        CALL GETLST_NR(AMAT,ICORE(1),MAXCOR,LISTW,1)
c      ENDIF
      
C
C   NOW ADD THE INTEGRALS <AJ//BI> (ORDERING AI,BJ)
C   THEY ARE STORED ON LIST 23, HOWEVER TAKE CARE OF
C   THE SIGN, SINCE THEY ARE ACTUALLY <JA//BI>
C

c      LISTW=23 + ISHIFT 
      listw=23

      NUMSYW=IRPDPD(1,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(1,ISYTYP(1,LISTW))

      CALL GETLST(SCRATCH,1,NUMSYW,2,1,LISTW)

      CALL SAXPY(NUMSYW*DISSYW,ONEM,SCRATCH,1,AMAT,1)

C   NOW ADD THE INTEGRALS <Ab//Ij> (ORDERING AI,bj)
C   WITH A FACTOR 2 TO A
C   THEY ARE STORED ON LIST 18
C
c      LISTW=18 + ISHIFT 
      listw=18     
      NUMSYW=IRPDPD(1,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(1,ISYTYP(1,LISTW))
C
c      IF(bRedundant) THEN
       CALL GETLST(SCRATCH,1,NUMSYW,2,1,LISTW)
c      ELSE
c       CALL GETLST_NR(SCRATCH,ICORE(1),MAXCOR,LISTW,1)
c      ENDIF

C
      CALL SAXPY(NUMSYW*DISSYW,TWO,SCRATCH,1,AMAT,1)    
C
C  ALL DONE
C
      RETURN
      END      
