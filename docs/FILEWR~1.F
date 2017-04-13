C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      SUBROUTINE FILE_WRAPPER (FILNAM,LUNIT,LSTAT,LFORM,IOERR)
C
C     This routine is a short wrapper around the file handler CHOPEN.
C     If an error occurs while a file is opened, the user may be able
C     to make corrections. In the Windows version there should be a
C     dialog box here...
C
      CHARACTER*255 FILNAM
      INTEGER*4 LUNIT,LSTAT,LFORM,IOERR
C
      IOERR = 0
      CALL CHOPEN(FILNAM,LUNIT,LSTAT,LFORM,IOERR)
      IF (IOERR .NE. 0) THEN
         WRITE (IOLOG, 1000) IOERR
         IOERR = ABS(IOERR)
      ENDIF
      RETURN
C
 1000 FORMAT(' Error while opening I/O channel',I3,': error status',I3)
      END
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
