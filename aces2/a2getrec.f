
c The job archive subsystem must be started and stopped before and after
c every record transaction; otherwise, changes to and from the archive
c file will not get transferred between xaces2 and the member executables.

      subroutine a2getrec(iFlag,szArchive,szRecord,iLength,iDest)
      implicit none

      integer iFlag, iLength, iDest(*)
      character*(*) szArchive, szRecord

      call aces_ja_init
      call getrec(iFlag,szArchive,szRecord,iLength,iDest)
      call aces_ja_fin

      return
      end

