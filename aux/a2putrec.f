
c The job archive subsystem must be started and stopped before and after
c every record transaction; otherwise, changes to and from the archive
c file will not get transferred between xaces2 and the member executables.

      subroutine a2putrec(iFlag,szArchive,szRecord,iLength,iSrc)
      implicit none

      integer iFlag, iLength, iSrc(*)
      character*(*) szArchive, szRecord

      call aces_ja_init
      call putrec(iFlag,szArchive,szRecord,iLength,iSrc)
      call aces_ja_fin

      return
      end

