module fermi

    implicit none
	
	contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Integrale de Fermi d'ordre -1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


real(8) function FDM0P5(XVALUE)


!   DESCRIPTION:
!
!      This function computes the Fermi-Dirac function of
!      order -1/2, defined as
!
!                     Int{0 to inf} t**(-1/2) / (1+exp(t-x)) dt
!         FDM0P5(x) = -----------------------------------------
!                                 Gamma(1/2)
!
!      The function uses Chebyshev expansions which are given to
!      16 decimal places for x <= 2, but only 10 decimal places 
!      for x > 2.
!
!
!   ERROR RETURNS:
!    
!      None.
!
!
!   MACHINE-DEPENDENT CONSTANTS:
!
!      NTERMS1 - INTEGER - The number of terms used from the array
!                          ARRFD1. The recommended value is such that
!                               ABS(ARRFD1(NTERMS1)) < EPS/10
!                          subject to 1 <= NTERMS1 <= 14.
!
!      NTERMS2 - INTEGER - The number of terms used from the array
!                          ARRFD2. The recommended value is such that
!                               ABS(ARRFD2(NTERMS2)) < EPS/10
!                          subject to 1 <= NTERMS1 <= 23.
!
!      NTERMS3 - INTEGER - The number of terms used from the array
!                          ARRFD3. The recommended value is such that
!                               ABS(ARRFD3(NTERMS3)) < EPS/10
!                          subject to 1 <= NTERMS3 <= 28.
!
!      XMIN1 - REAL - The value of x below which
!                         FDM0P5(x) = exp(x)
!                     to machine precision. The recommended value
!                     is    LN ( SQRT(2) * EPSNEG )
!
!      XMIN2 - REAL - The value of x below which
!                         FDM0P5(x) = 0.0 
!                     to machine precision. The recommended value
!                     is    LN ( XMIN )
!
!      XHIGH - REAL - The value of x above which
!                         FDM0P5(x) = 2 sqrt (x/pi) 
!                     to machine precision. The recommended value
!                     is    1 / sqrt( 2 * EPSNEG )
!
!      For values of EPS, EPSNEG, and XMIN the user should refer to the
!      paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
!   
!      This code is provided with single and double precision values
!      of the machine-dependent parameters, suitable for machines
!      which satisfy the IEEE floating-point standard.
!
!
!   AUTHOR:
!          DR. ALLAN MACLEOD,
!          DEPT. OF MATHEMATICS AND STATISTICS,
!          UNIVERSITY OF PAISLEY,
!          HIGH ST.,
!          PAISLEY,
!          SCOTLAND
!          PA1 2BE
!
!          (e-mail: macl-ms0@paisley.ac.uk )
!
!
!   LATEST UPDATE:
!                 20 NOVEMBER, 1996
!

      INTEGER NTERM1,NTERM2,NTERM3
      
	  REAL(8)  ARRFD1(0:14),ARRFD2(0:23),ARRFD3(0:58), &
			CHV,EXPX,FIFTY,FORTY2,GAM1P5,ONE,T,THREE,TWO,TWOE, &
			X,XHIGH,XMIN1,XMIN2,XSQ,XVALUE,ZERO


      DATA ARRFD1/	+1.7863596385102264d0,  &
					-0.999372007632333d-1,  &
					+0.64144652216054d-2,   &
					-0.4356415371345d-3,    &
					+0.305216700310d-4,     &
					-0.21810648110d-5,      &
					+0.1580050781d-6,       &
					-0.115620570d-7,        &
					+0.8525860d-9,          &
					-0.632529d-10,          &
					+0.47159d-11,           &
					-0.3530d-12,            &
					+0.265d-13,             &
					-0.20d-14,              &
					+0.2d-15                /  


      DATA ARRFD2( 0) /  1.6877111526052352d0   /
      DATA ARRFD2( 1) /  0.5978360226336983d0   / 
      DATA ARRFD2( 2) /  0.357226004541669d-1   /
      DATA ARRFD2( 3) / -0.132144786506426d-1   /
      DATA ARRFD2( 4) / -0.4040134207447d-3     /
      DATA ARRFD2( 5) /  0.5330011846887d-3     /
      DATA ARRFD2( 6) / -0.148923504863d-4      /
      DATA ARRFD2( 7) / -0.218863822916d-4      /
      DATA ARRFD2( 8) /  0.19652084277d-5       / 
      DATA ARRFD2( 9) /  0.8565830466d-6        /
      DATA ARRFD2(10) / -0.1407723133d-6        /
      DATA ARRFD2(11) / -0.305175803d-7         /
      DATA ARRFD2(12) /  0.83524532d-8          /
      DATA ARRFD2(13) /  0.9025750d-9           /
      DATA ARRFD2(14) / -0.4455471d-9           /
      DATA ARRFD2(15) / -0.148342d-10           /
      DATA ARRFD2(16) /  0.219266d-10           /
      DATA ARRFD2(17) / -0.6579d-12             /
      DATA ARRFD2(18) / -0.10009d-11            /
      DATA ARRFD2(19) /  0.936d-13              /
      DATA ARRFD2(20) /  0.420d-13              /
      DATA ARRFD2(21) / -0.71d-14               /
      DATA ARRFD2(22) / -0.16d-14               /
      DATA ARRFD2(23) /  0.4d-15                /


      DATA ARRFD3(0)  /   0.8707195029590563d0	/
      DATA ARRFD3(1)  /   0.59833110231733d-2	/
      DATA ARRFD3(2)  /  -0.432670470895746d-1	/
      DATA ARRFD3(3)  /  -0.393083681608590d-1	/
      DATA ARRFD3(4)  /  -0.191482688045932d-1	/
      DATA ARRFD3(5)  /  -0.65582880980158d-2	/
      DATA ARRFD3(6)  /  -0.22276691516312d-2	/
      DATA ARRFD3(7)  /  -0.8466786936178d-3	/
      DATA ARRFD3(8)  /  -0.2807459489219d-3	/
      DATA ARRFD3(9)  /  -0.955575024348d-4		/
      DATA ARRFD3(10) /  -0.362367662803d-4		/
      DATA ARRFD3(11) /  -0.109158468869d-4		/
      DATA ARRFD3(12) /  -0.39356701000d-5		/
      DATA ARRFD3(13) /  -0.13108192725d-5		/
      DATA ARRFD3(14) /  -0.2468816388d-6		/
      DATA ARRFD3(15) /  -0.1048380311d-6		/
      DATA ARRFD3(16) /   0.236181487d-7		/
      DATA ARRFD3(17) /   0.227145359d-7		/
      DATA ARRFD3(18) /   0.145775174d-7		/
      DATA ARRFD3(19) /   0.153926767d-7		/
      DATA ARRFD3(20) /   0.56924772d-8			/
      DATA ARRFD3(21) /   0.50623068d-8			/
      DATA ARRFD3(22) /   0.23426075d-8			/	
      DATA ARRFD3(23) /   0.12652275d-8			/
      DATA ARRFD3(24) /   0.8927773d-9			/
      DATA ARRFD3(25) /   0.2994501d-9			/
      DATA ARRFD3(26) /   0.2822785d-9			/
      DATA ARRFD3(27) /   0.910685d-10			/
      DATA ARRFD3(28) /   0.696285d-10			/
      DATA ARRFD3(29) /   0.366225d-10			/
      DATA ARRFD3(30) /   0.124351d-10			/
      DATA ARRFD3(31) /   0.145019d-10			/
      DATA ARRFD3(32) /   0.16645d-11			/
      DATA ARRFD3(33) /   0.45856d-11			/
      DATA ARRFD3(34) /   0.6092d-12			/
      DATA ARRFD3(35) /   0.9331d-12			/
      DATA ARRFD3(36) /   0.5238d-12			/
      DATA ARRFD3(37) /  -0.56d-14				/
      DATA ARRFD3(38) /   0.3170d-12			/
      DATA ARRFD3(39) /  -0.926d-13				/
      DATA ARRFD3(40) /   0.1265d-12			/
      DATA ARRFD3(41) /  -0.327d-13				/
      DATA ARRFD3(42) /   0.276d-13				/
      DATA ARRFD3(43) /   0.33d-14				/
      DATA ARRFD3(44) /  -0.42d-14				/
      DATA ARRFD3(45) /   0.101d-13				/
      DATA ARRFD3(46) /  -0.73d-14				/
      DATA ARRFD3(47) /   0.64d-14				/
      DATA ARRFD3(48) /  -0.37d-14				/
      DATA ARRFD3(49) /   0.23d-14				/
      DATA ARRFD3(50) /  -0.9d-15				/
      DATA ARRFD3(51) /   0.2d-15				/
      DATA ARRFD3(52) /   0.2d-15				/
      DATA ARRFD3(53) /  -0.3d-15				/
      DATA ARRFD3(54) /   0.4d-15				/
      DATA ARRFD3(55) /  -0.3d-15				/
      DATA ARRFD3(56) /   0.2d-15				/
      DATA ARRFD3(57) /  -0.1d-15				/
      DATA ARRFD3(58) /   0.1d-15				/


      DATA ZERO,ONE,TWO/ 0.0d0, 1.0d0, 2.0d0/
      DATA THREE,FORTY2,FIFTY/ 3.0D0, 42.0d0, 50.0d0/
      DATA GAM1P5/ 0.8862269254527580d0/
      DATA TWOE/ 5.4365636569180905d0/
!
!   Machine-dependent constants
!
      DATA NTERM1,NTERM2,NTERM3/14,23,58/
      DATA XMIN1,XMIN2,XHIGH/-36.39023d0,-708.39641d0,67108864.0d0/
!
!   Start calculation
!
      X=XVALUE
!
!   Code for x < -1
!
      IF ( X .LT. -ONE ) THEN
         IF ( X .GT. XMIN1 ) THEN
            EXPX = EXP(X)
            T = TWOE * EXPX - ONE
            FDM0P5 = EXPX * CHEVAL(NTERM1, ARRFD1, T)
         ELSE
            IF ( X .LT. XMIN2 ) THEN
               FDM0P5 = ZERO
            ELSE
               FDM0P5 = EXP(X)
            ENDIF
         ENDIF
      ELSE
!
!   Code for -1 <= x <= 2
!
         IF ( X .LE. TWO ) THEN
            T = ( TWO * X - ONE ) / THREE
            FDM0P5 = CHEVAL(NTERM2, ARRFD2, T)
         ELSE
!
!   Code for x > 2
!
            FDM0P5 = SQRT(X) / GAM1P5
            IF ( X .LE. XHIGH ) THEN 
               XSQ = X * X
               T = ( FIFTY - XSQ ) / ( FORTY2 + XSQ )
               CHV = CHEVAL(NTERM3, ARRFD3, T)
               FDM0P5 = FDM0P5 * ( ONE - CHV / XSQ )
            ENDIF
         ENDIF
      ENDIF

	return

end function FDM0P5


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Integrale de Fermi d'ordre 1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function FDP0P5(XVALUE)

!
!   DESCRIPTION:
!
!      This function computes the Fermi-Dirac function of
!      order 1/2, defined as
!
!                     Int{0 to inf} t**(1/2) / (1+exp(t-x)) dt
!         FDP0P5(x) = -----------------------------------------
!                                 Gamma(3/2)
!
!      The function uses Chebyshev expansions which are given to
!      16 decimal places for x <= 2, but only 10 decimal places 
!      for x > 2.
!
!
!   ERROR RETURNS:
!    
!      If XVALUE too large and positive, the function value
!      will overflow. An error message is printed and the function
!      returns the value 0.0.
!
!
!   MACHINE-DEPENDENT CONSTANTS:
!
!      NTERMS1 - INTEGER - The number of terms used from the array
!                          ARRFD1. The recommended value is such that
!                               ABS(ARRFD1(NTERMS1)) < EPS/10
!                          subject to 1 <= NTERMS1 <= 13.
!
!      NTERMS2 - INTEGER - The number of terms used from the array
!                          ARRFD2. The recommended value is such that
!                               ABS(ARRFD2(NTERMS2)) < EPS/10
!                          subject to 1 <= NTERMS1 <= 23.
!
!      NTERMS3 - INTEGER - The number of terms used from the array
!                          ARRFD3. The recommended value is such that
!                               ABS(ARRFD3(NTERMS3)) < EPS/10
!                          subject to 1 <= NTERMS3 <= 32.
!
!      XMIN1 - REAL - The value of x below which
!                         FDP0P5(x) = exp(x)
!                     to machine precision. The recommended value
!                     is   1.5*LN(2) + LN(EPSNEG)
!
!      XMIN2 - REAL - The value of x below which
!                         FDP0P5(x) = 0.0 
!                     to machine precision. The recommended value
!                     is    LN ( XMIN )
!
!      XHIGH1 - REAL - The value of x above which
!                         FDP0P5(x) = x**(3/2)/GAMMA(5/2)
!                     to machine precision. The recommended value
!                     is   pi / SQRT(8*EPS)
!
!      XHIGH2 - REAL - The value of x above which FDP0P5 would 
!                      overflow. The reommended value is
!                              (1.329*XMAX)**(2/3)
!
!      For values of EPS, EPSNEG, and XMIN the user should refer to the
!      paper by Cody in ACM. Trans. Math. Soft. Vol. 14 (1988) p303-311.
!   
!      This code is provided with single and double precision values
!      of the machine-dependent parameters, suitable for machines
!      which satisfy the IEEE floating-point standard.
!
!
!   AUTHOR:
!          DR. ALLAN MACLEOD,
!          DEPT. OF MATHEMATICS AND STATISTICS,
!          UNIVERSITY OF PAISLEY,
!          HIGH ST.,
!          PAISLEY,
!          SCOTLAND
!          PA1 2BE
!
!          (e-mail: macl-ms0@paisley.ac.uk )
!
!
!   LATEST UPDATE:
!                 20 NOVEMBER, 1996
!


      INTEGER NTERM1,NTERM2,NTERM3
      REAL(8)	ARRFD1(0:13),ARRFD2(0:23),ARRFD3(0:53), &
						CHV,EXPX,FIFTY,FORTY2,GAM2P5,ONE,T, &
						THREE,TWO,TWOE,X,XHIGH1,XHIGH2,XMIN1,XMIN2, &
						XSQ,XVALUE,ZERO


      DATA ARRFD1/1.8862968392734597d0,		&
                 -0.543580817644053d-1,		&
                  0.23644975439720d-2,		&
                 -0.1216929365880d-3,		&	
                  0.68695130622d-5,			&
                 -0.4112076172d-6,			&
                  0.256351628d-7,			&
                 -0.16465008d-8,			& 
                  0.1081948d-9,				&
                 -0.72392d-11,				&
                  0.4915d-12,				&
                 -0.338d-13,				&
                  0.23d-14,					&
                 -0.2d-15					/


      DATA ARRFD2( 0)  /  2.6982492788170612d0	/
      DATA ARRFD2( 1)  /  1.2389914141133012d0	/
      DATA ARRFD2( 2)  /  0.2291439379816278d0	/
      DATA ARRFD2( 3)  /  0.90316534687279d-2	/
      DATA ARRFD2( 4)  / -0.25776524691246d-2	/
      DATA ARRFD2( 5)  / -0.583681605388d-4		/
      DATA ARRFD2( 6)  /  0.693609458725d-4		/
      DATA ARRFD2( 7)  / -0.18061670265d-5		/
      DATA ARRFD2( 8)  / -0.21321530005d-5		/
      DATA ARRFD2( 9)  /  0.1754983951d-6		/
      DATA ARRFD2(10)  /  0.665325470d-7		/
      DATA ARRFD2(11)  / -0.101675977d-7		/
      DATA ARRFD2(12)  / -0.19637597d-8			/
      DATA ARRFD2(13)  /  0.5075769d-9			/
      DATA ARRFD2(14)  /  0.491469d-10			/
      DATA ARRFD2(15)  / -0.233737d-10			/
      DATA ARRFD2(16)  / -0.6645d-12			/
      DATA ARRFD2(17)  /  0.10115d-11			/
      DATA ARRFD2(18)  / -0.313d-13				/
      DATA ARRFD2(19)  / -0.412d-13				/
      DATA ARRFD2(20)  /  0.38d-14				/
      DATA ARRFD2(21)  /  0.16d-14				/
      DATA ARRFD2(22)  / -0.3d-15				/
      DATA ARRFD2(23)  / -0.1d-15				/


      DATA ARRFD3(0)   /  2.5484384198009122d0	/
      DATA ARRFD3(1)   /  0.510439408960652d-1	/
      DATA ARRFD3(2)   /  0.77493527628294d-2	/
      DATA ARRFD3(3)   / -0.75041656584953d-2	/
      DATA ARRFD3(4)   / -0.77540826320296d-2	/
      DATA ARRFD3(5)   / -0.45810844539977d-2	/
      DATA ARRFD3(6)   / -0.23431641587363d-2	/
      DATA ARRFD3(7)   / -0.11788049513591d-2	/
      DATA ARRFD3(8)   / -0.5802739359702d-3	/
      DATA ARRFD3(9)   / -0.2825350700537d-3	/
      DATA ARRFD3(10)  / -0.1388136651799d-3	/
      DATA ARRFD3(11)  / -0.680695084875d-4		/
      DATA ARRFD3(12)  / -0.335356350608d-4		/
      DATA ARRFD3(13)  / -0.166533018734d-4		/
      DATA ARRFD3(14)  / -0.82714908266d-5		/
      DATA ARRFD3(15)  / -0.41425714409d-5		/
      DATA ARRFD3(16)  / -0.20805255294d-5		/
      DATA ARRFD3(17)  / -0.10479767478d-5		/
      DATA ARRFD3(18)  / -0.5315273802d-6		/
      DATA ARRFD3(19)  / -0.2694061178d-6		/
      DATA ARRFD3(20)  / -0.1374878749d-6		/
      DATA ARRFD3(21)  / -0.702308887d-7		/
      DATA ARRFD3(22)  / -0.359543942d-7		/
      DATA ARRFD3(23)  / -0.185106126d-7		/
      DATA ARRFD3(24)  / -0.95023937d-8			/
      DATA ARRFD3(25)  / -0.49184811d-8			/
      DATA ARRFD3(26)  / -0.25371950d-8			/
      DATA ARRFD3(27)  / -0.13151532d-8			/
      DATA ARRFD3(28)  / -0.6835168d-9			/
      DATA ARRFD3(29)  / -0.3538244d-9			/
      DATA ARRFD3(30)  / -0.1853182d-9			/
      DATA ARRFD3(31)  / -0.958983d-10			/
      DATA ARRFD3(32)  / -0.504083d-10			/
      DATA ARRFD3(33)  / -0.262238d-10			/
      DATA ARRFD3(34)  / -0.137255d-10			/
      DATA ARRFD3(35)  / -0.72340d-11			/
      DATA ARRFD3(36)  / -0.37429d-11			/
      DATA ARRFD3(37)  / -0.20059d-11			/
      DATA ARRFD3(38)  / -0.10269d-11			/
      DATA ARRFD3(39)  / -0.5551d-12			/
      DATA ARRFD3(40)  / -0.2857d-12			/
      DATA ARRFD3(41)  / -0.1520d-12			/
      DATA ARRFD3(42)  / -0.811d-13				/
      DATA ARRFD3(43)  / -0.410d-13				/
      DATA ARRFD3(44)  / -0.234d-13				/
      DATA ARRFD3(45)  / -0.110d-13				/
      DATA ARRFD3(46)  / -0.67d-14				/
      DATA ARRFD3(47)  / -0.30d-14				/
      DATA ARRFD3(48)  / -0.19d-14				/
      DATA ARRFD3(49)  / -0.9d-15				/
      DATA ARRFD3(50)  / -0.5d-15				/
      DATA ARRFD3(51)  / -0.3d-15				/
      DATA ARRFD3(52)  / -0.1d-15				/
      DATA ARRFD3(53)  / -0.1d-15				/	


      DATA ZERO,ONE,TWO/ 0.0d0, 1.0d0 , 2.0d0/
      DATA THREE,FORTY2,FIFTY/ 3.0d0, 42.0d0, 50.0d0/
      DATA GAM2P5/0.1329340388179137d1/
      DATA TWOE/5.4365636569180905d0/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
      DATA NTERM1,NTERM2,NTERM3/13,23,53/
      DATA XMIN1,XMIN2/-35.7d0,-708.394d0/
      DATA XHIGH1,XHIGH2/7.45467d7,3.8392996d205/
!
!   Start calculation
!
      X=XVALUE
!
!   Test for error condition
!
      IF ( X .GT. XHIGH2 ) THEN
         PRINT*,'** ERROR ** - X TOO LARGE FOR FDP0P5',X
         FDP0P5 = ZERO
         RETURN
      ENDIF    
!
!   Code for x < -1
!
      IF ( X .LT. -ONE ) THEN
         IF ( X .GT. XMIN1 ) THEN
            EXPX = EXP(X)
            T = TWOE * EXPX - ONE
            FDP0P5 = EXPX * CHEVAL ( NTERM1 , ARRFD1 , T )
         ELSE
            IF ( X .LT. XMIN2 ) THEN
               FDP0P5 = ZERO
            ELSE
               FDP0P5 = EXP(X)
            ENDIF
         ENDIF
      ELSE
!
!   Code for -1 <= x <= 2
!
         IF ( X .LE. TWO ) THEN
            T = ( TWO * X - ONE ) / THREE
            FDP0P5 = CHEVAL ( NTERM2 , ARRFD2 , T )
         ELSE
!
!   Code for x > 2
!
            FDP0P5 = X * SQRT(X) / GAM2P5
            IF ( X .LE. XHIGH1 ) THEN 
               XSQ = X * X
               T = ( FIFTY - XSQ ) / ( FORTY2 + XSQ )
               CHV = CHEVAL ( NTERM3 , ARRFD3 , T )
               FDP0P5 = FDP0P5 * ( ONE + CHV / XSQ )
            ENDIF
         ENDIF
      ENDIF


	return

end function FDP0P5


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8)  function CHEVAL(N,A,T)

!
!   DESCRIPTION:
!
!      This function evaluates a Chebyshev series, using the
!      Clenshaw method with Reinsch modification, as analysed
!      in the paper by Oliver.
!
!
!   INPUT PARAMETERS
!
!       N - INTEGER - The no. of terms in the sequence
!
!       A - REAL ARRAY, dimension 0 to N - The coefficients of
!           the Chebyshev series
!
!       T - REAL - The value at which the series is to be
!           evaluated
!
!
!   REFERENCES
!
!        "An error analysis of the modified Clenshaw method for
!         evaluating Chebyshev and Fourier series" J. Oliver,
!         J.I.M.A., vol. 20, 1977, pp379-391
!
!
!   MACHINE-DEPENDENT CONSTANTS: NONE
!
!
!   INTRINSIC FUNCTIONS USED;
!
!      ABS
!
!
!    AUTHOR:  Dr. Allan J. MacLeod,
!             Dept. of Mathematics and Statistics,
!             University of Paisley ,
!             High St.,
!             PAISLEY,
!             SCOTLAND
!             ( e-mail:  macl-ms0@paisley.ac.uk )
!
!
!   LATEST MODIFICATION:
!                       21 September , 1995
!
!

      INTEGER I,N
      REAL(8) A(0:N),D1,D2,HALF,T,TEST,TT,TWO,U0,U1,U2,ZERO
      DATA ZERO,HALF/ 0.0d0 , 0.5d0 /
      DATA TEST,TWO/ 0.6d0 , 2.0d0 /
      U1 = ZERO
!
!   If ABS ( T )  < 0.6 use the standard Clenshaw method
!
      IF ( ABS( T ) .LT. TEST ) THEN
         U0 = ZERO
         TT = T + T
         DO 100 I = N , 0 , -1
            U2 = U1
            U1 = U0
            U0 = TT * U1 + A( I ) - U2
 100     CONTINUE
         CHEVAL =  ( U0 - U2 ) / TWO
      ELSE
	  
!
!   If ABS ( T )  > =  0.6 use the Reinsch modification
!
         D1 = ZERO
!
!   T > =  0.6 code
!
         IF ( T .GT. ZERO ) THEN
            TT =  ( T - HALF ) - HALF
            TT = TT + TT
            DO 200 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) + D2
               U1 = D1 + U2
 200        CONTINUE
            CHEVAL =  ( D1 + D2 ) / TWO
         ELSE
!
!   T < =  -0.6 code
!
            TT =  ( T + HALF ) + HALF
            TT = TT + TT
            DO 300 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) - D2
               U1 = D1 - U2
 300        CONTINUE
            CHEVAL =  ( D1 - D2 ) / TWO
         ENDIF
      ENDIF


	return

end function CHEVAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INVERSE FERMI  INTEGRAL OF ORDER 1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(8) function IFDP0P5(f)

	implicit none

!..this routine applies a rational function expansion to get the inverse
!..fermi-dirac integral of order 1/2 when it is equal to f.
!..maximum error is 4.19d-9.   reference: antia apjs 84,101 1993
!..
!..declare
 
	integer          i,m1,k1,m2,k2
    REAL(8) f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff

!..load the coefficients of the expansion
    data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
    data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3, &
						 6.610132843877d2,   3.818838129486d1, 1.0d0/
	
	
    data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
                          9.130355392717d1,  -1.670718177489d0/
    
	data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
                         -4.262314235106d-1,  4.997559426872d-1, &
                         -1.285579118012d0,  -3.930805454272d-1, &
                          1.0d0/
						  
    data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
                        -3.299466243260d-1,  4.077841975923d-1, &
                        -1.145531476975d0,  -6.067091689181d-2/


    if (f .lt. 4.0d0) then
	
		rn = f + a1(m1)
	
		do i=m1-1,1,-1
			rn = rn*f + a1(i)
		enddo
		
		den = b1(k1+1)
       
		do i=k1,1,-1
			den = den*f + b1(i)
		enddo
       
		IFDP0P5 = log(f * rn/den)

    else
	  
		ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
		
		rn = ff + a2(m2)
       
		do i=m2-1,1,-1
			rn = rn*ff + a2(i)
		enddo
       
		den = b2(k2+1)
       
		do i=k2,1,-1
			den = den*ff + b2(i)
		enddo
       
		IFDP0P5 = rn/(den*ff)
    
	end if
    
	return
    
	end function IFDP0P5


end module fermi
