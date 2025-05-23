C   THIS CODE MUST BE RUN BEFORE THE CODE Fill_SW_gaps.for  !
c
c   This code goes over the merged yearly 5-min avg solar wind data file (original 5-min OMNI, downloaded from CDAWEB)
c    and fills the gaps in the IMF coverage, using linear interpolation. The gaps are filled on condition
c    that they are not too wide (controlled by the upper limit LIMGAP on the number of missing 5-min records).
c
c    Author:  N. A. Tsyganenko

      SUBROUTINE FILLIMFGAPS(OMNIYEAR, length, DIRECTORY2)
c
      DIMENSION VXGSE(105430),VYGSE(105430),VZGSE(105430),TEMP(105430),
     +DEN(105430),BXGSM(105430),BYGSM(105430),BZGSM(105430),SYMH(105430)
     +,IDAY(105430),IHOUR(105430),MIN(105430), IMFFLAG(105430),
     +VGSE(105430)

      INTEGER*4 OMNIYEAR,length
      CHARACTER(length + 80) :: NAMEIN,NAMEOUT
      CHARACTER(length) :: DIRECTORY,DIRECTORY2
      character(len=256) :: dir
      CHARACTER(LEN=256) :: cmd, found_path
      INTEGER :: ierr, last_backslash

      
      ! Detect Windows by checking the environment variable "OS"
      CHARACTER(LEN=32) :: os_name
      CALL GET_ENVIRONMENT_VARIABLE("OS", os_name)

      IF (INDEX(os_name, "Windows") > 0) THEN
          
      ! Correctly formatted Windows search command
      cmd = 'where /R C:\ omni_5min_*.lst ' // '> found_pyd.txt'
          
      CALL EXECUTE_COMMAND_LINE(cmd, EXITSTAT=ierr)

      IF (ierr == 0) THEN
      OPEN(UNIT=10, FILE='found_pyd.txt', STATUS='OLD')
      READ(10, '(A)') found_path
      CLOSE(10, STATUS='DELETE')
      ! Initialize position of the last backslash
      last_backslash = 0

      ! Loop through the string from the end to find the first backslash
      DO i = LEN_TRIM(found_path), 1, -1
          IF (found_path(i:i) == '\') THEN
              last_backslash = i
              EXIT
          END IF
      END DO

      IF (last_backslash > 0) THEN
      found_path = found_path(1:last_backslash-1)  ! Keep only the directory part
      DIRECTORY = found_path
      DIRECTORY2 = DIRECTORY
      ELSE
      END IF
      ELSE
      END IF
      ELSE
      END IF

      DIRECTORY2 = DIRECTORY
c
c     DO 111 IYEAR=1995,2013
      DO 111 IYEAR=OMNIYEAR,OMNIYEAR

      DO I=1,105430
      IMFFLAG(I)=+1   !  INITIALLY SET AT +1, WHICH MEANS "IMF DATA AVAILABLE" FOR ALL 5-MIN INTERVALS
      END DO
      WRITE(NAMEIN, '(A, "omni_5min_", I0, ".lst")') TRIM(DIRECTORY) 
     & // "/", IYEAR

      WRITE(NAMEOUT, '(A, I0, "_IMF_gaps_le_3hrs_filled.txt")') 
     & TRIM(DIRECTORY)// "/", IYEAR

      LIMGAP=36           !  GAPS LARGER THAN 36*5=180 MIN = 3 HOURS WILL REMAIN UNFILLED
C
      OPEN (UNIT=1,FILE=NAMEIN,STATUS='OLD')

C     READ ALL THE DATA INTO RAM:

      I=0

 1    READ (1,502,END=2) IYEA,IDA,IHOU,MI,BX,BY,BZ,V,VX,VY,VZ,DE,T,ISYMH
      MINYEAR=(IDA-1)*1440+IHOU*60+MI
      IND=MINYEAR/5+1

      if (ind.lt.0) then
       print *, ' I=',I
       print *, ' IYEA=',IYEA
       print *, ' IDA=',IDA 
       print *, ' IHOU=',IHOU
       print *, ' MI=',MI
       print *, ' IND=',IND
       print *, ' MINYEAR=',MINYEAR
      endif

      IF(BX.EQ.9999.99.OR.BY.EQ.9999.99.OR.BZ.EQ.9999.99)IMFFLAG(IND)=-1
      IDAY (IND)=IDA
      IHOUR(IND)=IHOU
      MIN  (IND)=MI
      BXGSM(IND)=BX
      BYGSM(IND)=BY
      BZGSM(IND)=BZ
      VGSE(IND)=V
      VXGSE(IND)=VX
      VYGSE(IND)=VY
      VZGSE(IND)=VZ
CCCC      VTH  (IND)=0.15745*SQRT(T)  !  SQRT(3kT/M_p) in km/s    !  DISABLED HERE
      TEMP (IND)=T
      DEN  (IND)=DE
      SYMH (IND)=ISYMH
      I=I+1
      GOTO 1

 2    CLOSE(1)
      NREC=I

      !PRINT *, '  READING IN RAM FINISHED; NREC=',NREC    !  NUMBER OF RECORDS

 502  FORMAT(2I4,2I3,53X,F8.2,16X,2F8.2,16X,4F8.1,F7.2,F9.0,98X,I6)
 503  FORMAT(2I4,2I3,53X,F8.2,16X,2F8.2,16X,4F8.1,F7.2,F9.0,F6.0,3X,I2)

C       NOW MAKE THE INTERPOLATION FOR THE INTERNAL GAPS

       IFLAG=0    !  INDICATES THAT WE STILL DID NOT GET TO THE FIRST RECORD WITH VALID IMF DATA

       DO 21 I=1,NREC

          IF (IMFFLAG(I).EQ.-1.AND.IFLAG.EQ.0) GOTO 21
          IFLAG=1   !  GOT TO THE FIRST VALID RECORD

        IF (IMFFLAG(I).NE.-1.AND.IMFFLAG(I+1).EQ.-1) THEN  ! FIRST INTERNAL NO-IMF-DATA INTERVAL ENCOUNTERED
        IBEG=I
        BXBEG=BXGSM(IBEG)
        BYBEG=BYGSM(IBEG)
        BZBEG=BZGSM(IBEG)
        ENDIF

        IF (IMFFLAG(I).EQ.-1.AND.IMFFLAG(I+1).NE.-1) THEN    !  FIRST VALID RECORD AFTER THE IMF GAP ENCOUNTERED
        IEND=I+1

        IF (IEND-IBEG.GT.LIMGAP) GOTO 21

        BXEND=BXGSM(IEND)
        BYEND=BYGSM(IEND)
        BZEND=BZGSM(IEND)

        DO 22 K=IBEG+1,IEND-1
        A=(BXEND-BXBEG)/(IEND-IBEG)              !   LINEAR INTERPOLATION
        B=(BXBEG*IEND-BXEND*IBEG)/(IEND-IBEG)
        BXGSM(K)=A*K+B
        A=(BYEND-BYBEG)/(IEND-IBEG)
        B=(BYBEG*IEND-BYEND*IBEG)/(IEND-IBEG)
        BYGSM(K)=A*K+B
        A=(BZEND-BZBEG)/(IEND-IBEG)
        B=(BZBEG*IEND-BZEND*IBEG)/(IEND-IBEG)
        BZGSM(K)=A*K+B

        IMFFLAG(K)=2   
C   WILL BE WRITTEN INTO THE OUTPUT FILE RECORD, INDICATING THAT THE VALUE WAS
C                            OBTAINED USING THE LINEAR INTERPOLATION
 22     CONTINUE
        ENDIF

 21     CONTINUE

       OPEN (UNIT=1,FILE=NAMEOUT)

       DO 7 I=1,NREC
         WRITE (1,503) IYEA,IDAY(I),IHOUR(I),MIN(I),BXGSM(I),BYGSM(I),
     *   BZGSM(I),VGSE(I),VXGSE(I),VYGSE(I),VZGSE(I),DEN(I),TEMP(I),
     *   SYMH(I),IMFFLAG(I)
  7    CONTINUE

       CLOSE(1)

111    CONTINUE
      END SUBROUTINE FILLIMFGAPS


C   THIS CODE MUST BE RUN BEFORE THE CODE Fill_SW_gaps.for  !
c
c   This code goes over the merged yearly 5-min avg solar wind data file (original 5-min OMNI, downloaded from CDAWEB)
c    and fills the gaps in the IMF coverage, using linear interpolation. The gaps are filled on condition
c    that they are not too wide (controlled by the upper limit LIMGAP on the number of missing 5-min records).
c
c    Author:  N. A. Tsyganenko

      SUBROUTINE FILLIMFGAPSLINUX(OMNIYEAR, length, DIRECTORY)
c
      DIMENSION VXGSE(105430),VYGSE(105430),VZGSE(105430),TEMP(105430),
     +DEN(105430),BXGSM(105430),BYGSM(105430),BZGSM(105430),SYMH(105430)
     +,IDAY(105430),IHOUR(105430),MIN(105430), IMFFLAG(105430),
     +VGSE(105430)

      INTEGER*4 OMNIYEAR,length
      CHARACTER(length + 80) :: NAMEIN,NAMEOUT
      CHARACTER(length) :: DIRECTORY

c
c     DO 111 IYEAR=1995,2013
      DO 111 IYEAR=OMNIYEAR,OMNIYEAR

      DO I=1,105430
      IMFFLAG(I)=+1   !  INITIALLY SET AT +1, WHICH MEANS "IMF DATA AVAILABLE" FOR ALL 5-MIN INTERVALS
      END DO
      WRITE(NAMEIN, '(A, "omni_5min_", I0, ".lst")') TRIM(DIRECTORY) 
     & // "/", IYEAR

      WRITE(NAMEOUT, '(A, I0, "_IMF_gaps_le_3hrs_filled.txt")') 
     & TRIM(DIRECTORY)// "/", IYEAR

      LIMGAP=36           !  GAPS LARGER THAN 36*5=180 MIN = 3 HOURS WILL REMAIN UNFILLED
C
      OPEN (UNIT=1,FILE=NAMEIN,STATUS='OLD')

C     READ ALL THE DATA INTO RAM:

      I=0

 1    READ (1,502,END=2) IYEA,IDA,IHOU,MI,BX,BY,BZ,V,VX,VY,VZ,DE,T,ISYMH
      MINYEAR=(IDA-1)*1440+IHOU*60+MI
      IND=MINYEAR/5+1

      if (ind.lt.0) then
       print *, ' I=',I
       print *, ' IYEA=',IYEA
       print *, ' IDA=',IDA 
       print *, ' IHOU=',IHOU
       print *, ' MI=',MI
       print *, ' IND=',IND
       print *, ' MINYEAR=',MINYEAR
      endif

      IF(BX.EQ.9999.99.OR.BY.EQ.9999.99.OR.BZ.EQ.9999.99)IMFFLAG(IND)=-1
      IDAY (IND)=IDA
      IHOUR(IND)=IHOU
      MIN  (IND)=MI
      BXGSM(IND)=BX
      BYGSM(IND)=BY
      BZGSM(IND)=BZ
      VGSE(IND)=V
      VXGSE(IND)=VX
      VYGSE(IND)=VY
      VZGSE(IND)=VZ
CCCC      VTH  (IND)=0.15745*SQRT(T)  !  SQRT(3kT/M_p) in km/s    !  DISABLED HERE
      TEMP (IND)=T
      DEN  (IND)=DE
      SYMH (IND)=ISYMH
      I=I+1
      GOTO 1

 2    CLOSE(1)
      NREC=I

      !PRINT *, '  READING IN RAM FINISHED; NREC=',NREC    !  NUMBER OF RECORDS

 502  FORMAT(2I4,2I3,53X,F8.2,16X,2F8.2,16X,4F8.1,F7.2,F9.0,98X,I6)
 503  FORMAT(2I4,2I3,53X,F8.2,16X,2F8.2,16X,4F8.1,F7.2,F9.0,F6.0,3X,I2)

C       NOW MAKE THE INTERPOLATION FOR THE INTERNAL GAPS

       IFLAG=0    !  INDICATES THAT WE STILL DID NOT GET TO THE FIRST RECORD WITH VALID IMF DATA

       DO 21 I=1,NREC

          IF (IMFFLAG(I).EQ.-1.AND.IFLAG.EQ.0) GOTO 21
          IFLAG=1   !  GOT TO THE FIRST VALID RECORD

        IF (IMFFLAG(I).NE.-1.AND.IMFFLAG(I+1).EQ.-1) THEN  ! FIRST INTERNAL NO-IMF-DATA INTERVAL ENCOUNTERED
        IBEG=I
        BXBEG=BXGSM(IBEG)
        BYBEG=BYGSM(IBEG)
        BZBEG=BZGSM(IBEG)
        ENDIF

        IF (IMFFLAG(I).EQ.-1.AND.IMFFLAG(I+1).NE.-1) THEN    !  FIRST VALID RECORD AFTER THE IMF GAP ENCOUNTERED
        IEND=I+1

        IF (IEND-IBEG.GT.LIMGAP) GOTO 21

        BXEND=BXGSM(IEND)
        BYEND=BYGSM(IEND)
        BZEND=BZGSM(IEND)

        DO 22 K=IBEG+1,IEND-1
        A=(BXEND-BXBEG)/(IEND-IBEG)              !   LINEAR INTERPOLATION
        B=(BXBEG*IEND-BXEND*IBEG)/(IEND-IBEG)
        BXGSM(K)=A*K+B
        A=(BYEND-BYBEG)/(IEND-IBEG)
        B=(BYBEG*IEND-BYEND*IBEG)/(IEND-IBEG)
        BYGSM(K)=A*K+B
        A=(BZEND-BZBEG)/(IEND-IBEG)
        B=(BZBEG*IEND-BZEND*IBEG)/(IEND-IBEG)
        BZGSM(K)=A*K+B

        IMFFLAG(K)=2   
C   WILL BE WRITTEN INTO THE OUTPUT FILE RECORD, INDICATING THAT THE VALUE WAS
C                            OBTAINED USING THE LINEAR INTERPOLATION
 22     CONTINUE
        ENDIF

 21     CONTINUE

       OPEN (UNIT=1,FILE=NAMEOUT)

       DO 7 I=1,NREC
         WRITE (1,503) IYEA,IDAY(I),IHOUR(I),MIN(I),BXGSM(I),BYGSM(I),
     *   BZGSM(I),VGSE(I),VXGSE(I),VYGSE(I),VZGSE(I),DEN(I),TEMP(I),
     *   SYMH(I),IMFFLAG(I)
  7    CONTINUE

       CLOSE(1)

111    CONTINUE
      END SUBROUTINE FILLIMFGAPSLINUX