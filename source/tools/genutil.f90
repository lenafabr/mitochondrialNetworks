MODULE GENUTIL
  ! generally useful utilities
  USE MT19937 ! mersenne random number generator

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793D0
  
CONTAINS  

  INTEGER FUNCTION STRING2NUM(STRINGIN,APPENDNUM)
    ! convert the string to a unique number based on ascii characters
    ! the characters SPACE, {,},(,),[,],",`,<,> and all nonprintable characters are ignored
    ! at most the last five characters (ignoring the unacceptable characters above) at the end of the string are used
    ! any leading "!" do not affect the final number (these map to 0)
    ! if APPENDNUM is specificied, only use the last 4 characters of the string as well as the additional number modulo 84

    IMPLICIT NONE
    CHARACTER(LEN=*) :: STRINGIN
    CHARACTER(LEN=5) :: STRING
    INTEGER, OPTIONAL :: APPENDNUM
    INTEGER :: DIGARRAY(5)
    INTEGER :: ALLOWED(84)
    INTEGER :: N, I, D, COUNT
    CHARACTER*84 :: ALLOWEDSTR

    ! set the allowed characters
    ALLOWED(1:6) = (/33,35,36,37,38,39/)
    ALLOWED(7:24) = (/(I,I=42,59)/)
    ALLOWED(25:27) = (/61,63,64/)
    ALLOWED(28:53) = (/(I, I=65,90)/)
    ALLOWED(54:56) = (/92,94,95/)
    ALLOWED(57:82) = (/(I, I=97,122)/)
    ALLOWED(83:84) = (/124,126/)

    N = LEN(STRINGIN)
    IF (PRESENT(APPENDNUM)) THEN
       STRING(1:4) = STRINGIN(N-3:N)
       STRING(5:5) = ACHAR(ALLOWED(MOD(APPENDNUM,84)+1))
    ELSE
       STRING = STRINGIN(N-4:N)
    ENDIF
    N =  5


    DO I = 1,84
       ALLOWEDSTR(I:I) = ACHAR(ALLOWED(I))
    ENDDO

    DIGARRAY = 0
    COUNT = 0
    DO I = 0,N-1
       D = INDEX(ALLOWEDSTR,STRING(N-I:N-I),.FALSE.)
       IF (D.EQ.0) THEN
          print*, 'Ignoring character:', D
          CYCLE
       ENDIF

       DIGARRAY(5-COUNT) = D-1
       COUNT = COUNT + 1
       IF (COUNT.GE.5) EXIT
    ENDDO

    STRING2NUM = BASE2DEC(DIGARRAY,5,84)
  END FUNCTION STRING2NUM

  INTEGER FUNCTION BASE2DEC(DIGARRAY,N,BASE)
  ! given a number in some integer base (specified as a list of digits)
  ! convert that number to a decimal integer
  ! N is the size of the list
  ! if resulting number is too large, wrap around to negative numbers
  ! starting from the right, only use as many of the digits as 
  ! will fit into the resulting integer between -HUGE and HUGE  
  ! if any digit is greater than base-1, print error and stop

  IMPLICIT NONE
  INTEGER, DIMENSION(N) :: DIGARRAY
  INTEGER, INTENT(IN) :: N, BASE
  INTEGER :: MAXDIG, I, D

  MAXDIG = INT(LOG(2*DBLE(HUGE(BASE))+2)/LOG(DBLE(BASE)))

  BASE2DEC = 0
  DO I = 0, MIN(N-1,MAXDIG-1)
     D = DIGARRAY(N-I)
     IF (D.EQ.0) CYCLE
     IF (D.GT.BASE-1) THEN
        PRINT*, 'ERROR in BASE2DEC: digit is bigger than base.', I, D, BASE
        STOP 1
     ENDIF
     
     BASE2DEC = BASE2DEC + D*BASE**I
  ENDDO
  
  END FUNCTION BASE2DEC

  SUBROUTINE INTERPARRAY(ARRAY,NA,COL,VAL,IND,INTERP)
    ! for an 2D array with dimensions NA
    ! use the values in column COL to interpolate for the value VAL
    ! return the index IND such that ARRAY(IND,COL)<VAL<ARRAY(IND+1,COL)
    ! and the interpolation of all other columns in INTERP
    ! If val is out of bounds returns IND=0 or IND=NL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NA(2), COL
    DOUBLE PRECISION, INTENT(IN) :: ARRAY(NA(1),NA(2))
    DOUBLE PRECISION, INTENT(IN) :: VAL
    INTEGER, INTENT(OUT) :: IND
    DOUBLE PRECISION, INTENT(OUT) :: INTERP(NA(2))
    DOUBLE PRECISION :: FRAC

    CALL INTERP1(ARRAY(:,COL),NA(1),VAL,IND)

    IF (IND.EQ.0.OR.IND.GE.NA(1)) RETURN

    FRAC = (VAL-ARRAY(IND,COL))/(ARRAY(IND+1,COL)-ARRAY(IND,COL))
    INTERP = (1-FRAC)*ARRAY(IND,:)+FRAC*ARRAY(IND+1,:)
       
  END SUBROUTINE INTERPARRAY

  SUBROUTINE INTERP1(LIST,NL,VAL,IND)
    ! for a monotonically increasing, double precision list
    ! find the index where LIST(IND) <VAL<LIST(IND+1)
    INTEGER, INTENT(IN) :: NL
    DOUBLE PRECISION, INTENT(IN) :: LIST(NL)
    DOUBLE PRECISION, INTENT(IN) :: VAL
    INTEGER, INTENT(OUT) :: IND
    DOUBLE PRECISION :: MINL, MAXL,PREVMAXL
    INTEGER :: MINI, MAXI,PREVMAXI
    LOGICAL :: VERBOSE = .FALSE.

    MINI = 1; MAXI = NL; MINL = LIST(1); MAXL = LIST(NL)
    IF (VAL.LT.MINL) THEN
       IND = 0; RETURN
    ELSEIF (VAL.EQ.MINL) THEN
       IND = 1; RETURN       
    ELSEIF (VAL.GT.MAXL) THEN
       IND = NL; RETURN
    ELSEIF (VAL.EQ.MAXL) THEN
       IND = NL-1; RETURN
    ENDIF

    DO WHILE (MAXI-MINI.GT.1.OR.MAXL.LT.VAL)
       IF (MAXL.GT.VAL) THEN
          PREVMAXL = MAXL; PREVMAXI = MAXI
          MAXI = MINI + (MAXI-MINI)/2
          MAXL = LIST(MAXI)
       ELSE
          MINI = MAXI; MAXI = PREVMAXI
          MINL = MAXL; MAXL = PREVMAXL
       ENDIF
       IF (VERBOSE) PRINT*, 'MINI, MAXI, MINL, MAXL', MINI, MAXI, MINL, MAXL,VAL
       if (maxi.eq.mini) then
          print*, 'something weird in interp1:', list(1), list(nl), val          
          stop 1
       endif
    ENDDO

    IF (.NOT.(MAXI.EQ.MINI+1.AND.LIST(MINI).LE.VAL.AND.LIST(MAXI).GE.VAL)) THEN
       PRINT*, 'SOMETHING IS WEIRD IN INTERP1', val, mini, maxi, list(mini), list(maxi)
       STOP 1
    ENDIF

    IND = MINI
  END SUBROUTINE INTERP1

  SUBROUTINE REPLACESUBSTR(INSTRING,C,REPL)
    ! replace the last instance of the substring C in INSTRING with REPL
    IMPLICIT NONE
    CHARACTER*100 :: INSTRING
    CHARACTER(LEN=*) :: C
    CHARACTER(LEN=*) :: REPL
    INTEGER :: LENC, IND

    INSTRING = ADJUSTL(INSTRING)

    LENC = LEN_TRIM(C)

    IND = INDEX(INSTRING,C,.TRUE.)
    IF (IND.GT.0) THEN! if * was found in the string

       INSTRING = INSTRING(1:IND-1) // TRIM(ADJUSTL(REPL)) // INSTRING(IND+LENC:100)   
    END IF
  END SUBROUTINE REPLACESUBSTR

  SUBROUTINE NORMALIZE(X)
    ! normalize a 3 dimensional vector

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT) :: X(3)
    DOUBLE PRECISION :: DX

    DX = SQRT(DOT_PRODUCT(X,X))
    X(:) = X(:)/DX

    RETURN
  END SUBROUTINE NORMALIZE

  DOUBLE PRECISION FUNCTION NORM(X)
    ! norm of 3D vector X

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: X(3)

    NORM = sqrt(DOT_PRODUCT(X,X))

  END FUNCTION NORM

  SUBROUTINE CROSS_PRODUCT(A, B, C)
    ! take the cross product of 3D vectors A and B; return result in C

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: A(3), B(3)
    DOUBLE PRECISION, INTENT(OUT) :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1)-A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)

    RETURN
  END SUBROUTINE CROSS_PRODUCT

  SUBROUTINE RANDOMAXIS(REFAX,CTRANGE,RANDAX)
    ! generate a random axis, within a certain range in cos(theta) 
    ! relative to the reference axis    
    DOUBLE PRECISION, INTENT(IN) :: REFAX(3), CTRANGE
    DOUBLE PRECISION, INTENT(OUT) :: RANDAX(3)
    DOUBLE PRECISION :: THETA, pHI, E1(3), E2(3), E3(3), X, Y, Z

    THETA = 1.0D0
    DO WHILE (THETA.EQ.1.0D0)
       ! get random number; MT19937 algorithm uses closed interval [0,1], 
       ! so ignore when R is exactly 1
       THETA = GRND() !get a random number
    ENDDO
    !CALL RANDOM_NUMBER(THETA)
    THETA = acos(1D0 - THETA*MAX(CTRANGE,2D0))

    PHI = 1.0D0
    DO WHILE (PHI.EQ.1.0D0)
       PHI = GRND() !get a random number
    ENDDO
    !CALL RANDOM_NUMBER(PHI)
    PHI = PHI*2*PI

    ! axis system relative to which angles are defined
    E3 = REFAX
    IF (E3(2) == 0 .AND. E3(3) == 0) THEN
       E2 = (/0D0,1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(E3, (/1D0,0D0,0D0/),E2)
    END IF
    CALL CROSS_PRODUCT(E2, E3, E1)

    CALL NORMALIZE(E1); CALL NORMALIZE(E2); CALL NORMALIZE(E3)

    ! generate the axis around which to rotate
    X = sin(THETA)*cos(PHI)
    Y = sin(THETA)*sin(PHI)
    Z = cos(THETA)

    RANDAX = X*E1 + Y*E2 + Z*E3

  END SUBROUTINE RANDOMAXIS

  SUBROUTINE GETPERP(V1,V2)
    ! get some unit vector perpendicular to V1 and store it in V2
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: V1(3)
    DOUBLE PRECISION, INTENT(OUT) :: V2(3)

    IF (V1(2).EQ.0.AND.V1(3).EQ.0) THEN
       V2 = (/0D0,1D0,0D0/)
    ELSE
       CALL CROSS_PRODUCT(V1,(/0D0,1D0,0D0/),V2)
       CALL NORMALIZE(V2)
    ENDIF
  END SUBROUTINE GETPERP

  SUBROUTINE LINELINEINTERSECT(P1,V1IN,P2,V2IN, T1, T2, PA1, PA2, INTPT)
    ! find the intersection point of a line defined by P1 + T1*V1
    ! (P1 is starting point, V1 is direction vector) and P2 + T2*V2
    ! T1, T2 = parameter coordinates along the 2 lines that mark the points
    ! of closest approach
    ! PA1, PA2 = points of closest approach
    ! INTPT = point halfway between the points of closest approach
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: P1, V1IN, P2, V2IN
    DOUBLE PRECISION, INTENT(OUT) :: T1, T2, INTPT(3), PA1(3), PA2(3)
    DOUBLE PRECISION :: V1(3), V2(3), DP(3), PV1, PV2, VV
    
    ! normalize the direction vectors
    V1 = V1IN/SQRT(SUM(V1IN**2))
    V2 = V2IN/SQRT(SUM(V2IN**2))
    
    DP = P2 - P1
    PV1 = DOT_PRODUCT(DP,V1)
    PV2 = DOT_PRODUCT(DP,V2)
    VV = DOT_PRODUCT(V1,V2)

    IF (ABS(ABS(VV)-1D0).LT.EPSILON(1D0)*2) THEN
       ! lines are parallel
       ! in this case, just take the point half way between the Ps     
       INTPT = (P1+P2)/2
       IF (VV.GT.0) THEN
          T1 = DOT_PRODUCT(INTPT-P1,V1)
          T2 = DOT_PRODUCT(INTPT-P2,V2)
       ELSE
          T1 = -DOT_PRODUCT(INTPT-P1,V1)
          T2 = -DOT_PRODUCT(INTPT-P2,V2)
       ENDIF
    ELSE
       ! find point of closest intersection

       T2 = (VV*PV1 - PV2)/(1-VV**2)
       T1 = PV1 + VV*T2
    ENDIF
    ! points of closest approach
    PA1 = P1 + V1*T1
    PA2 = P2 + V2*T2
    INTPT = (PA1+PA2)/2       
  END SUBROUTINE LINELINEINTERSECT


  SUBROUTINE SEGSEGINTERSECT(DIM,P1,P2,P3,P4, MUA, MUB, PA, PB, INTPT)
   ! find the intersection point of a segment defined by P1*(1-MUA) + P2*MUA
   ! (P1 is starting point, P2 is the ending point of first segment) and P3*(1-MUB) + P4*MUB
   ! MUA, MUB = parameter coordinates along the 2 lines that mark the points
   ! of closest approach
   ! PA, PB = points of closest approach
   ! INTPT = point halfway between the points of closest approach
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: DIM
   DOUBLE PRECISION, INTENT(IN) :: P1(DIM),P2(DIM),P3(DIM),P4(DIM)
   DOUBLE PRECISION, INTENT(OUT) :: MUA, MUB, PA(DIM), PB(DIM), INTPT(DIM)
   DOUBLE PRECISION :: D1343, D4321, D1321, D4343, D2121, DENOM
   DOUBLE PRECISION, DIMENSION(DIM) :: R21,R13,R43
   LOGICAL :: AOUT,BOUT

   AOUT=.FALSE.; BOUT=.FALSE.

   R21 = P2-P1
   R13 = P1-P3
   R43 = P4-P3
   D4321 = DOT_PRODUCT(R43,R21)
   D4343 = DOT_PRODUCT(R43,R43)
   D2121 = DOT_PRODUCT(R21,R21)
   D1343 = DOT_PRODUCT(R13,R43)
   D1321 = DOT_PRODUCT(R13,R21)
   DENOM = D2121*D4343 - D4321*D4321

   IF(DENOM.NE.0) THEN    
      MUA = (D1343*D4321 - D1321*D4343)/DENOM
      IF(MUA.GT.1D0) THEN
      AOUT = .TRUE.; MUA = 1D0
      ELSEIF(MUA.LT.0D0) THEN
      AOUT = .TRUE.; MUA = 0D0
      ENDIF
      MUB = (D1343 + MUA*D4321)/D4343
   
      IF(MUB.GT.1D0) THEN
         BOUT = .TRUE.; MUB = 1D0
      ELSEIF(MUB.LT.0D0) THEN
         BOUT = .TRUE.; MUB = 0D0
      ENDIF

      IF(BOUT) THEN
         MUA = (MUB*D4321 - D1321)/D2121
         IF(MUA.GT.1D0) THEN
         AOUT = .TRUE.; MUA = 1D0
         ELSEIF(MUA.LT.0D0) THEN
         AOUT = .TRUE.; MUA = 0D0
         ELSE
         AOUT = .FALSE.
         ENDIF
      ENDIF

   ELSE ! parallel lines
      MUA = 0D0
      MUB = D1343/D4343
      IF(MUB.GT.1D0) THEN
         BOUT = .TRUE.; MUB = 1D0
      ELSEIF(MUB.LT.0D0) THEN
         BOUT = .TRUE.; MUB = 0D0
      ENDIF

      IF(BOUT) THEN
         MUA = (MUB*D4321-D1321)/D2121
         IF(MUA.GT.1D0) THEN
         MUA = 1D0
         ELSEIF(MUA.LT.0D0) THEN
         MUA = 0D0
         ENDIF
      ENDIF
   ENDIF

   PA = (1-MUA)*P1 + MUA*P2
   PB = (1-MUB)*P3 + MUB*P4
   INTPT = 0.5*(PA + PB)
  END SUBROUTINE SEGSEGINTERSECT


  SUBROUTINE RAYRAYINTERSECT(P1,V1IN,P2,V2IN,RSTER, T1, T2, PA1, PA2, INTPT)
    ! find the intersection point of a ray defined by P1 + T1*V1 (t1>=0)
    ! (P1 is starting point, V1 is direction vector) and P2 + T2*V2 (t2>=0)
    ! T1, T2 = parameter coordinates along the 2 lines that mark the points
    ! of closest approach
    ! PA1, PA2 = points of closest approach
    ! INTPT = point halfway between the points of closest approach
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: P1, V1IN, P2, V2IN
    DOUBLE PRECISION, INTENT(IN) :: RSTER
    DOUBLE PRECISION, INTENT(OUT) :: T1, T2, INTPT(3), PA1(3), PA2(3)
    DOUBLE PRECISION :: V1(3), V2(3), DP(3), PV1, PV2, VV
    
    ! normalize the direction vectors
    V1 = V1IN/SQRT(SUM(V1IN**2))
    V2 = V2IN/SQRT(SUM(V2IN**2))
    
    DP = P2 - P1
    PV1 = DOT_PRODUCT(DP,V1)
    PV2 = DOT_PRODUCT(DP,V2)
    VV = DOT_PRODUCT(V1,V2)

    IF (ABS(ABS(VV)-1D0).LT.EPSILON(1D0)*2) THEN
       ! lines are parallel
       ! in this case, just take the point half way between the Ps     
       INTPT = (P1+P2)/2
       IF (VV.GT.0) THEN
          T1 = DOT_PRODUCT(INTPT-P1,V1)
          T2 = DOT_PRODUCT(INTPT-P2,V2)
       ELSE
          T1 = -DOT_PRODUCT(INTPT-P1,V1)
          T2 = -DOT_PRODUCT(INTPT-P2,V2)
       ENDIF

    ELSE
       ! find point of closest intersection       
       T2 = (VV*PV1 - PV2)/(1-VV**2)
       T1 = PV1 + VV*T2
    ENDIF

    IF (T1.LT.0) THEN
       T1 = 0D0
       ! closest approach of point p1 to line2
       T2 = MAX(-PV2,0D0)
    ELSE IF (T2.LT.0) THEN
       T2 = 0D0
       ! closest approach of point p2 to line1
       T1 = MAX(-PV1,0D0)
    ENDIF

    !T1 = MIN(T1,2*RSTER)
    !T2 = MIN(T2,2*RSTER)
    
    ! points of closest approach
    PA1 = P1 + V1*T1
    PA2 = P2 + V2*T2

    ! point in between
    INTPT = (PA1+PA2)/2
  END SUBROUTINE RAYRAYINTERSECT

  SUBROUTINE RANDOMSPHEREPOINT(R, SHELL, POINT)
     ! generate a random point in a sphere
     ! If SPHERE is true, choose the point on the surface, else distribute in the volume
     ! R = sphere radius
     ! POINT = (X, Y, Z) position vector
     IMPLICIT NONE
 
     DOUBLE PRECISION, INTENT(IN) :: R
     LOGICAL, INTENT(IN) :: SHELL
     DOUBLE PRECISION, INTENT(OUT) :: POINT(3)
     DOUBLE PRECISION :: PHI, COSTHETA, RAD, SINTHETA

     PHI = GRND()*2*PI
     COSTHETA = 2*GRND()-1
     IF(SHELL) THEN
        RAD = R
     ELSE
        RAD = R*(GRND()**(1D0/3))
     ENDIF
     SINTHETA = SQRT(1-COSTHETA**2)

     POINT = RAD*(/ SINTHETA*COS(PHI), SINTHETA*SIN(PHI), COSTHETA /)
  END SUBROUTINE RANDOMSPHEREPOINT

END MODULE GENUTIL
