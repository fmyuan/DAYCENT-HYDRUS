
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C 3456789012345678901234567890123456789012345678901234567890123456789012 4567890
      CHARACTER*20 FUNCTION READSTRING (UNIT,buffr,lch,ststrt,stend)
C
C=======================================================================
C
C  Description
C    Parses character data separated by white spaces or hard delimiters.
C    see comments below.  Null record results in zero parsed field
C
C  Arguments
C
C  - Input
C    BUFFR, character*(*), input string
C    LCH   , integer,       is the index of the first character to parse
C                           updated by the routine as parsing occurs
C    UNIT,   integer,       IO unit to read input lines from
C
C  - Output
C    ststrt, integer,       index of string start
C    stend,  integer,       index of string end
C
C  History
C    Modified : 03/31/95 : Kendrick Killian
C              Forced a soft terminator at the buffer end.  Corrects the
C              bug which causes right justified strings to be dropped
C    Modified : 01/25/95 : Kendrick Killian
C              Modified routine to locate and return a single field
C    Written : 3/15/91 : Kendrick Killian 
C
C  Error Conditions
C    None
C
C  External References
C       subroutine READCLIN   reads an input line
C
C  Additional comments
C    1) a data field is established/bounded by a pair of upstream and
C       downstream delimiters.
C    2) commas are considered as hard delimiters:
C       - hard delimiters delimits field with no exception.
C       - consequtive hard delimiters constitute a null field.
C    3) White spaces (blanks and tabs) are soft delimeters:
C       - soft delimiters that preceed or follow an established data 
C         field are ignored.
C       - consecutive soft delimeters do not constitute a null field.
C    4) a record is assumed to be bounded by soft delimiters.
C
C      Copyright 1992 - 95    Kendrick Killian     All Rights Reserved
C========================== BEGIN EXECUTABLE CODE ======================
C
      implicit none
      CHARACTER    buffr*(*), CHAR*1
      INTEGER      UNIT, ststrt, stend, lch

      INTEGER I2, I, LENLIN
      LOGICAL ISWHITE, ISHARD
C
C  Logical statement functions to defines delimeters
C
      ISWHITE(CHAR) = ((CHAR.EQ.' ') .OR. (CHAR.EQ.'	'))
      ISHARD(CHAR) = (CHAR.EQ.',')
C
C Set Initial Values
C
      READSTRING = " "
      I = lch - 1
C
C Set status (I2) flag to a hard delimiter condition
C    I2 condition codes are:
C        -1  Field terminated with a white space
C         0  Last Field terminated with a hard delimiter
C        >0  tracking through a token; location of last hard delimeter
C
      I2 = 0
C
C determine record specified length
C
      LENLIN = LEN(BUFFR)
C
C read a record if required
C

20    I= I+1
      IF (I.GT.LENLIN) THEN
        IF (I2.GT.0) THEN
          ststrt = I2
          stend = LENLIN
          GOTO 90
        ELSE
          CALL READCLIN(UNIT,BUFFR,lch)
          I = lch
        ENDIF
      ENDIF
C
C Check to see if we have encountered a hard delimeter
C
        IF (ISHARD(BUFFR(I:I))) THEN
C
C Check if the previous field was terminated with a hard terminator
C  If so this is a null field
C
          IF (I2 .EQ. 0 ) THEN
            ststrt = lch
            stend  = lch
C
C Check to see if we are stepping over a data field
C If so we need to terminate the field
C
          ELSE IF (I2 .GT. 0 ) THEN
            ststrt = I2
            stend = I-1
          ENDIF
C
C  Found a field
          GOTO 90
C
C Have we encountered a white character with an intervening field
C
        ELSE IF ((I2 .GT. 0) .AND. ISWHITE(BUFFR(I:I))) THEN
          ststrt = I2
          stend = I-1
          GOTO 90
C
C If the character is not white and we are starting a new field and
C we need to store its location  
C (NOTE: don't worry about hard terminators since they won't get here
C
        ELSE IF ((I2 .LE. 0) .AND. .NOT.ISWHITE(BUFFR(I:I))) THEN
          I2 = I
        END IF
C
C If we are here we are either stepping over white spaces after 
C terminating a field or we are stepping over a valid field
C  In either case    KEEP GOING
C
      GOTO 20
C
C Clean up and end processing
C
   90 CONTINUE
      lch = I+1
      READSTRING = BUFFR(ststrt:stend)
      RETURN
      END

C2345678901234567890123456789012345678901234567890123456789012345678901234567890
      REAL FUNCTION READREAL(UNIT, BUFFR, LCH, ISTAT)
C
C=======================================================================
C
C  This is a general purpose routine for parsing a numeric character
C  string and converting it to a number. 
C
C  Arguments
C  - Input
C    UNIT , INTEGER,       unit to read string from if the pointer is past
C                          the buffr end  LCH > len(BUFFR).  a new BUFFR
C                          will be read and the conversion will start with
C                          the new line.  NOTE: if the buffer is empty and 
C                          LCH < the length of buffer no string is read.
C    BUFFR, character*(*), character string to be converted.
C    LCH  , INTEGER,       is the index of the first character to parse
C                          this is updated as parsing occurs
C    ISTAT, integer,       NOTE: if READREAL is called with ISTAT =-2 or -3,
C                                the previous real is returned.
C                                This allows access to the real decoded 
C                                in a failed integer read.
C  - Output
C    LCH  , INTEGER,    the index of the next token
C    ISTAT, integer,    conversion status
C                   < -1, Conversion error
C                      0, no numerical conversion (alpha field)
C                      1, valid integer (and real)
C                      2, valid floating point and integer approximation
C                      3, valid floating point and integer overflow
C
C  External References
C       subroutine READCLIN   reads an input line
C
C  History
C    Modified : 02/12/95 : Kendrick Killian
C              Substantial recoding to reduce the possibility of misreporting
C              unusual strings (strings starting with .,-,+,e)
C    Modified : 01/25/95 : Kendrick Killian
C              Modified arguments
C    Modified : 01/23/95 : Kendrick Killian
C              Removed FORTRAN77 non-compliant internal unformatted READ
C    Written  : 03/27/91 : Kendrick Killian
C
C      Copyright 1992 - 95    Kendrick Killian     All Rights Reserved
C========================== BEGIN EXECUTABLE CODE ======================
C
      implicit none
      CHARACTER BUFFR*(*),CHR*1
      INTEGER LCH, ISTAT, readint, UNIT
C
C  local variables
      LOGICAL NOSPAC, RDEC, INTEG, NEG, NEGE
      INTEGER I, J, LOCE, LENLIN, Cip, ND
      DOUBLE PRECISION Cfp
      SAVE Cfp
      INTEGER*2 Cex, Cdec
C
C  Overflow/underflow value for a 32 bit integer
      INTEGER ovrflw, Ndig
      PARAMETER (ovrflw = 2147483647, Ndig = 9)
C
C  Statement functions
      LOGICAL  ISUPC,ISLOWC,ISDIGIT,ISWHITE
C
      ISUPC(CHR) = CHR.GE.'A' .AND. CHR.LE.'Z'
      ISLOWC(CHR) = CHR.GE.'a' .AND. CHR.LE.'z'
      ISDIGIT(CHR) = CHR.GE.'0' .AND. CHR.LE.'9'
      ISWHITE(CHR)  = CHR.EQ.' ' .OR. CHR.EQ.'	'
C
      IF ((ISTAT .EQ. -2) .OR. (ISTAT .EQ. -3)) THEN
        READREAL = Cfp
        ISTAT=0
        RETURN
      ENDIF
      readreal  = 0.
      RDEC = .TRUE.
      GOTO 10

      ENTRY readint(UNIT, BUFFR, LCH, ISTAT)
      RDEC = .FALSE.
      readint = 0
C
 10   LENLIN = LEN(BUFFR)
      IF (LCH.GT.LENLIN) CALL READCLIN(UNIT,BUFFR,LCH)
C
      I     = LCH-1
      J     = 0
      ND    = 0
      Cex   = 0
      Cip   = 0
      Cfp   = 0.
      Cdec  = 0
      LOCE  = 0
      ISTAT = 0
      NEG    = .FALSE.
      NEGE   = .FALSE.
      INTEG  = .TRUE.
      NOSPAC = .TRUE.
C
C check for numeric entries.  Increment I until a non numeric character
C
20    I = I+1
      IF (I .LE. LENLIN) THEN
c      write (*,'(a,i3,3a,l2,f18.0,i10,3i4,3L3,2i4)')
c     &' character i=',I ,' /',buffr(I:I),'/',integ,
c     &   Cfp, Cip, ND, Cdec, Cex, NEG, NEGE, NOSPAC, J, loce
C   Accept a digit
        IF (ISDIGIT(BUFFR(I:I)) .AND. NOSPAC) THEN
          J = ICHAR(BUFFR(I:I)) - ICHAR('0')
          IF (LOCE.GT.0) THEN
C   record an exponent digit
            Cex = Cex * 10 + J
          ELSEIF (ND .GE. Ndig) THEN
C   Suppress excess precision  (record decimal point shift on integer)
            IF (INTEG) Cdec = Cdec+1
          ELSEIF (ND+J .EQ. 0)  THEN
C   Suppress leading zero's (record decimal shift if is a decimal
            ISTAT = 1
            IF (.NOT.INTEG) Cdec = Cdec-1
          ELSE
C   record a significant digit
            ISTAT = 1
            Cip = Cip * 10 + J
            ND = ND + 1
            IF (.NOT.INTEG) Cdec = Cdec-1
          ENDIF
          J=I
          NOSPAC = .TRUE.
C  Accept the first decimal point
        ELSE IF ((BUFFR(I:I).EQ.'.') .AND. INTEG .AND. NOSPAC) THEN
          INTEG = .FALSE.
          IF (J.EQ.0) J = I
C  Accept a plus sign
        ELSE IF ((BUFFR(I:I).EQ.'+') .AND. (J.LE.LOCE)) THEN
          IF (J.EQ.0) J = I
C  Accept a minus sign
        ELSE IF ((BUFFR(I:I).EQ.'-') .AND. (J.LE.LOCE)) THEN
          IF (J.EQ.0) J = I
          IF (LOCE.EQ.0) THEN
            NEG  = .TRUE.
          ELSE
            NEGE = .TRUE.
          ENDIF
C   Step over spaces
        ELSE IF (ISWHITE(BUFFR(I:I))) THEN
          NOSPAC = ((LOCE.ge.J).or.(J.eq.0))
C   Check for an E in an exponent field
        ELSE IF (LOCE.EQ.0 .AND. ISTAT.NE.0 .AND.
     &    ( (BUFFR(I:I).EQ.'E') .OR. (BUFFR(I:I).EQ.'e') .OR.
     &      (BUFFR(I:I).EQ.'D') .OR. (BUFFR(I:I).EQ.'d') ) ) THEN
          INTEG = .FALSE.
          NOSPAC = .TRUE.
          LOCE = I
        ELSE
c      WRITE (*,*) ' decoding',LCH,i,J,ND,LOCE,NOSPAC,ISTAT,BUFFR(I:I)
C         IF (I.GT.LENLIN .AND. J.EQ.0) GOTO 10
          GOTO 30
        ENDIF
C step to the next character
      GOTO 20
      ENDIF
C
C was anything accepted
30        LCH = I
          IF (ISTAT .EQ. 0) then
            IF (I.GT.LENLIN) GOTO 10
            IF(J.NE.0) LCH = J
            RETURN
          ENDIF
C  update the buffer pointer for special cases
          IF (LOCE .GT. J) LCH = LOCE
          IF (BUFFR(I:I).EQ.',') LCH = I + 1
C
CC  if you have reached here start the conversion to a real number
C         start by combining the mantissa and the exponents
c      WRITE (*,*) ' parsed',Cfp,Cip,Cdec,Cex,ND,NEG,NEGE,ISTAT
      IF (NEGE) Cex = -Cex
      Cex = Cex+Cdec
C
      Cfp =  DBLE(Cip) *10.d0**Cex
C
      IF (INTEG) THEN
        ISTAT = 1
      ELSE IF (Cfp.le.ovrflw) THEN
        Cip = Cfp
        ISTAT = 2
      ELSE
        Cip = ovrflw
        ISTAT = 3
      ENDIF
      IF (NEG) THEN
        Cfp = -Cfp
        Cip = -Cip
      ENDIF
C
      IF (RDEC) THEN
        readreal = Cfp
c        WRITE (*,*) ' --- Real Return --',readreal, LCH, ISTAT
      ELSE
        readint  = Cip
c        WRITE (*,*) ' --- Integer Return --',readint, LCH, ISTAT
      ENDIF
      RETURN
      END

      SUBROUTINE READCLIN(UNIT,BUFFR,LCH)
C
C=======================================================================
C
C  Description
C    Reads a character buffer from the specified unit.
C    The routine skips comment lines and upper cases alpha characters
C
C  - Input
C    UNIT,   integer,       IO unit to read input lines from
C
C  - Output
C    BUFFR, character*(*), input buffer
C    LCH   , integer,       is the index of the first character to parse
C
C  History
C    Modified : 03/14/95 : Kendrick Killian
C              added a 5 character over read to check for input overflow
C    Rewritten : 1/25/95 : Kendrick Killian 
C                          Recoded to remove special field processing
C
C  Error Conditions
C    None
C
C  External References
C    NONE
C
C  Additional comments
C    1) input commands are converted to upper case 
C    2) UNIX like comments can be inserted in the input stream
C       - # is the comment character
C       - comments extend to the end of the line
C
C      Copyright 1995   Colorado State University    All Rights Reserved
C=======================================================================
C
      implicit none
      CHARACTER BUFFR*(*), PAD*5, Comment*1 
      INTEGER LCH, I, UNIT
      PARAMETER (COMMENT='#')
C
10      READ(UNIT,'(a,a5)') BUFFR,PAD

C  Check for comments and buffer overflow
        I = INDEX(BUFFR,Comment)

C  error if non comment data extends past the buffer dimension
        IF (I.EQ.0) THEN
          IF (PAD.NE.'     ') THEN
            WRITE(*,'(a/4a)') ' Fatal Error: Input line to long:',
     &                      BUFFR,"/",PAD,"/"
            STOP 1
          ENDIF
        ELSE

C  Remove comments
          IF(I .EQ. 1) GOTO 10
          BUFFR = BUFFR(:I-1)
        ENDIF

C  check for a blank line
        LCH = 0
 15     LCH = LCH +1
C  NOTE: the second string is a TAB  This does not comply with FORTRAN 77
C        standards.  If this is a problem: ASCII tab is ICHAR(9)
        IF(LCH .GT. LEN(BUFFR)) GOTO 10
        IF(BUFFR(LCH:LCH).EQ.' ' .OR. BUFFR(LCH:LCH).EQ.'  ') GOTO 15

C  upper case the input line
      DO 50 I=LCH,LEN(BUFFR)
        IF (BUFFR(I:I).GE.'a' .AND. BUFFR(I:I).LE.'z') BUFFR(I:I) = 
     &               CHAR(ICHAR(BUFFR(I:I)) + (ICHAR('A')-ICHAR('a')))
 50   CONTINUE
      RETURN
      END
