MODULE my_module 
    USE iso_fortran_env, ONLY: li=>int64, dp=>real64 !iso_fortran_end is part of Fortran 2008 onwards. My aliases: li=Long Integer, dp=Double Precision
  IMPLICIT NONE
    TYPE :: func_result
        REAL(KIND=dp)                 :: min_x
        REAL(KIND=dp)                 :: max_x
        REAL(KIND=dp)                 :: root
        REAL(KIND=dp)                 :: mfr_typ
        CHARACTER(LEN=:), ALLOCATABLE :: latex
        CHARACTER(LEN=:), ALLOCATABLE :: descript
        REAL(KIND=dp)                 :: fx ! =f(x) evaluation
    END TYPE func_result

    REAL (KIND=dp), PARAMETER :: onePi = ACOS(-1.0_dp)
    REAL (KIND=dp), PARAMETER :: dbl_epsilon = EPSILON(1.0_dp) !Interval machine epsilon for double-precision: the difference between 1 and the next larger floating-point number
    INTEGER (KIND=li) :: funcN !Function number
    TYPE(func_result) :: fRes
    CONTAINS

    SUBROUTINE Algo1_Bisection(minX, maxX, xTol, dimEval, nEval, xret, fret)
!       Conventional bisection algorithm based on TOMS algoritm nr. 4 "Bisec", which was written in an early form of ALGOL
!        (e.g. ALGOL60) on a Royal Precision LGP-30 computer. The code here is Fortran 90 translation of "Bisec", but with
!        more modern termination criteria: 
!        (i) Stricter requirement for zero function value, i.e. fc=0, not abs(fc) <= user-defined epsilon
!       (ii) The "Birect" routine may never terminate if 1<<|xc| and input user-defined epsilon is less than machine precision
!            Here xTol functions as a *relative* tolerance. When the root is |xc|=0, this will never occur, so prevent infinite iteratin, we stop on nEval=dimEval
!      .Passed arguments:
        REAL (KIND=dp), INTENT(IN) :: minX !Left bound
        REAL (KIND=dp), INTENT(IN) :: maxX !Right bound
        REAL (KIND=dp), INTENT(IN) :: xTol !Tolerance. xTol must not be less than machine epsilon (dbl_epsilon), otherwise too many iterations
        INTEGER (KIND=li), INTENT(IN) :: dimEval !Maximum allowed number of function evaluations
        INTEGER (KIND=li), INTENT(OUT) :: nEval !actual number of function evaluations
        REAL (KIND=dp), INTENT(OUT) :: xRet !Returned x-value, where f(x)=0
        REAL (KIND=dp), INTENT(OUT) :: fRet !f(x), should be zero, or at least near-zero
!      .Internal variables:
        REAL (KIND=dp):: xa,xb,xc,fa,fb,fc

        xa = minX
        xb = maxX        
        fRes = TestFunction(funcN, xa); fa=fRes%fx
        fRes = TestFunction(funcN, xb); fb=fRes%fx
        nEval = 2

        xc = xa !wild guess in case fa and fb same sign
        fc = fa
        IF (fa * fb < 0.0_dp) THEN
            DO WHILE (fc /= 0.0_dp .AND. nEval < dimEval) !Termination criterion
                xc = xa * 0.5_dp + xb * 0.5_dp
                IF (xb - xa < xtol * ABS(xc)) EXIT
                fRes = TestFunction(funcN, xc); fc=fRes%fx
                nEval = nEval + 1
                !IF (0.0_dp < fa * fc) THEN !This is prone to underflow when the root is near zero
                IF (sign(1.0_dp,fa)==sign(1.0_dp,fc)) THEN !fix
                    xa = xc
                    fa = fc
                ELSE
                    xb = xc
                    fb = fc
                END IF
            END DO
        END IF
        xret = xc
        fret = fc
    END SUBROUTINE

    SUBROUTINE Algo2_BandBisection(minX, maxX, xtol, nEval, xRet, fRet)
!   This is the most accurate form of root finding possible. Features:
!   - Calculates root to the limit of floating-point precision, with no limit on the number of evaluations!
!   - Takes into account when f(x)=0 for a band of x-values. Brackets the zero-band from left & right sides, and calculates the median.
!   - Argument xtol is used to flag error when the zero-band is suspiciously wide. recommended value is 1E-14. Set xtol>=1 to avoid checking
!   - Suitable for calculating the "ground truth" root for evaluating other root-finding algorithms

!      .Passed arguments:
        REAL (KIND=dp), INTENT(IN) :: minX !Left bound
        REAL (KIND=dp), INTENT(IN) :: maxX !Right bound
        REAL (KIND=dp), INTENT(IN) :: xTol !Tolerance
        INTEGER (KIND=li), INTENT(OUT) :: nEval !actual number of function evaluations        
        REAL (KIND=dp), INTENT(OUT) :: xRet !Returned x-value, where f(x)≈0
        REAL (KIND=dp), INTENT(OUT) :: fRet !f(x), should be zero, or at least near-zero
!      .Internal variables:
        REAL (KIND=dp):: xa1,xb1,xc1,fa1,fb1,fc1,xa2,xb2,xc2,fc2,Xleft,Xright,maxWidth            
    
        xa1 = minX
        xb1 = maxX
        fRes = TestFunction(funcN, xa1); fa1=fRes%fx
        fRes = TestFunction(funcN, xb1); fb1=fRes%fx
        nEval = 2
    
        xc1 = xa1
        fc1 = fa1
        IF (SIGN(1.0_dp,fa1) /= SIGN(1.0_dp,fb1)) THEN
!           (1) Find zero
            DO
                IF (SIGN(1.0_dp,fa1) == SIGN(1.0_dp,fc1)) THEN !safer than 0<fa1*fc1 due to risk of underflow or overflow
                    xa1 = xc1
                    fa1 = fc1
                ELSE
                    xb1 = xc1
                END IF
                xc1 = xa1 + (xb1 - xa1) * 0.5_dp !Bisection, minimizes risk of overflow and underflow, and little risk of loss of precision when xba and xa1 are same magnitude                
                fRes = TestFunction(funcN, xc1); fc1=fRes%fx
                nEval = nEval + 1
                IF (xc1 == xa1 .OR. xc1 == xb1 .OR. fc1 == 0.0_dp) EXIT !End at limit of finite precision (xc= xa or xb) or is mid-point lands in a zero-band (fc=0)
            END DO
        
            IF (fc1 == 0.0_dp) THEN !The zero has finite width, so find the left and right edges of the zero band
!               (2a) Bracket left edge of zero-band
                xa2 = xa1
                xc2 = xc1
                fc2 = fc1
                DO
                    IF (fc2 == 0.0_dp) THEN; xb2 = xc2; ELSE; xa2 = xc2; END IF
                    xc2 = xa2 + (xb2 - xa2) * 0.5_dp
                    IF (xc2 == xa2 .OR. xc2 == xb2) EXIT !Reached limit of finite precision
                    fRes = TestFunction(funcN, xc2); fc2=fRes%fx
                    nEval = nEval + 1
                END DO
                Xleft = xc2

!               (2b) Bracket right edge of zero-band
                xb2 = xb1
                xc2 = xc1
                fc2 = fc1
                DO
                    IF (fc2 == 0.0_dp) THEN; xa2 = xc2; ELSE; xb2 = xc2; END IF
                    xc2 = xb2 + (xa2 - xb2) * 0.5_dp
                    IF (xc2 == xa2 .OR. xc2 == xb2) EXIT !Reached limit of finite precision
                    fRes = TestFunction(funcN, xc2); fc2=fRes%fx
                    nEval = nEval + 1
                END DO
                Xright = xc2

!               (2c) best-guess of root is median value
                xc1 = Xleft + (Xright - Xleft) * 0.5_dp
                fc1 = fc2 !no need to recalculate fc2 as it is zero everywhere within the band anyway

!               Optional additional check for band-width; it should be small. If it is large, there may be multiple zeros
                IF (1.0_dp < ABS(xc1)) THEN; maxWidth = xtol * ABS(xc1); ELSE; maxWidth = xtol; END IF !1E-14 recommended
                IF (maxWidth < Xright - Xleft) THEN
                    STOP !Suspiciously wide zero!
                ENDIF
            END IF
        END IF
        xret = xc1
        fret = fc1
    END SUBROUTINE

    SUBROUTINE ConfirmMachineEpsilon()
!       Calculates the machine's interval-epsilon = ULP(1.0).
!        Here, we use the industry-wide informal consensus on machine epsilon: the difference between 1 and the next larger floating-point number, https://en.wikipedia.org/wiki/Machine_epsilon
!       NOTE that the "relative rounding error" is half of the interval when unit roundoff applies the 'round-to-nearest' method, because real numbers in the lower half of the interval 1 to 1+epsilon are rounded down to 1, and numbers in the upper half of the interval are rounded up to 1+epsilon
!       IEEE 754 binary64: ULP(1.0) = 2**(-52) = 2.22E-16 (https://en.wikipedia.org/wiki/Double-precision_floating-point_format)
!       This harmonizes with programming languages C (DBL_EPSILON constant from float.h), Python (sys.float_info.epsilon), Matlab (eps)
!       Ref1: D. Goldberg (1991) "What every computer scientist should know about floating-point arithmetic", ACM Computing Surveys. Vol.23, No.1, pp 5–48, https://doi.org/10.1145/103162.103163
!       Ref2: https://en.wikipedia.org/wiki/Machine_epsilon
!        Below are shown 3 ways of checking the value of dbl_epsilon on this machine.
        REAL (KIND=dp):: eps, a, b 

        PRINT *, '2^(-52) =', 2.0_dp**(-52) !Quick way of calculating interval machine epsilon
        PRINT *, 'EPSILON(dp)=', EPSILON(1.0_dp) !Fortran90 function, returns smallest number eps of the same kind, such that 1 < 1+eps
        !Below are 3 ways of calculating/checking the value of dbl_epsilon on this machine:
        !Method 1:
        eps = 1.0_dp
        DO WHILE (1.0_dp < 1.0_dp+eps/2.0_dp)
            eps = eps / 2.0_dp
        END DO
        PRINT *, '1) Interval machine epsilon =', eps
        !Method 2:
        a = 1.0_dp + 1.1E-16
        b = 1.0_dp + 1.12E-16
        PRINT *, '2) ULP(1) =', b - a !a gets rounded down; b bets rounded up
        !Method 3:
        PRINT *, '3) Both must be T:', (1.0_dp == 1.0_dp+eps/2.0_dp), (1.0_dp < 1.0_dp+eps) !the first is rounded down; the second is rounded up
    END SUBROUTINE
    
        FUNCTION TestFunction(f_N, x) RESULT(res)
        !This is a carefully curated library of one-dimensional functions for testing root-finding algorithms.
        !ARGUMENT:  DESCRIPTION:
        !f_N      Function number. +ve to calculate f(x), -ve to get bounds and
        !x          If f_N<=0 then x is a data index [0,3]: x=0 output number of variants, x=1 output lower x-bound, x=2 output upper x-bound, x=3 output hardcoded root value of subtype Decimal
        !x          If 0<f_N then x is the value to calculate f(x). The input value of x should be within the bounds
        !Note 1: A function can be plotted by using as user-defined function in cells in a spreadsheet, e.g. =TestFunction(1, 1, 1.234)
        !Note 2: Many of the harcoded root values are coded as the sum of two numbers, each of 18 digits. This is to gain more accuracy, as numbers in Visual Basic are stored in higher precision but can only be iput as 15 digits
        !Note 3: Some functions that are tangential to the x-axis, giving a wide range of x-values with f(x)=0. This is solved here by approximating with a linear function near the zero, giving a precise zero where the line crosses the x-axis
        !Note 4: Hardcoded roots are only approximate, and are provided for reference/vetting only, deviate up to ±5E-15 from exact roots. Roots for comparing algorithms should be calculated using bisection with no limit to number of evaluations.
        !Author: P.G.Schild, 2023

!      .Passed arguments:
        INTEGER (KIND=li), INTENT(IN):: f_N
        REAL (KIND=dp), INTENT(IN):: x
        TYPE(func_result):: res
!       .Internal variables:
        REAL (KIND=dp):: a_

!-------------------------------------
!       Functions with a "broad" zero
!-------------------------------------        
        
        SELECT CASE (ABS(f_N))
        CASE (1)
            !REF: type III in [Bus75]. Used in f17 & f18 in [Suhadolnik12] with different bounds
            !Computation complexity = 6
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = 10.0_dp
                res%root     = 0.0_dp
                res%mfr_typ  = 3.3_dp
                res%descript = "Monotonic (strictly), 34%-constant-gradient, Smooth, Non-linear near root"
                res%latex    = "f(x) = x^{25}"
            ELSE !calculate f(x)
                res%fx = x ** 25
            END IF
    
        CASE (2)
            !REF: #12 in [Alefeld95], f48-f50 in [FernandezDiaz22]
            !NOTE: The analytical function is smooth, but finite precision gives it a fine-grained staircase pattern. The linear approximtion near the root below does not fully alleviate the issue, just ensures non-constant where it crosses the x-axis
            !Computation complexity = 19
            IF (f_N < 0) THEN !output properties
                res%min_x    = 1.0_dp
                res%max_x    = 100.0_dp
                res%root     = 33.0_dp
                res%mfr_typ  = 3.3_dp
                res%descript = "Monotonic (strictly), Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = x ** (1.0_dp / 33.0_dp) - 33.0_dp ** (1.0_dp / 33.0_dp)
            END IF
        
        CASE (3)
            !REF: f8 in [Thukral13] with bounds here by P.G.Schild 2023
            !Computation complexity = 46
            IF (f_N < 0) THEN !output properties
                res%min_x    = -3.0_dp
                res%max_x    = 3.0_dp
                res%root     = -1.57248441639688204179771917429_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Monotonic (weakly, with plateaus), 82%-constant-value, Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = (6.0_dp * x ** 5 + 5.0_dp * x ** 4 - 4.0_dp * x ** 3 + 3.0_dp * x ** 2 - 2.0_dp * x + 1.0_dp) ** 55
            END IF
        
        CASE (4)
            !REF: [Rice69], A15-A17 in [Nerinckx76], P7 in [Wu05]
            !Computation complexity = 9
            IF (f_N < 0) THEN !output properties
                res%min_x    = -onePi / 5.0_dp
                res%max_x    = 0.4_dp * onePi
                res%root     = 0.0_dp
                res%mfr_typ  = 3.3_dp
                res%descript = "Monotonic (strictly with 1 inflection), Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = TAN(x) - x
            END IF
        
        CASE (5)
            !REF: f4 in [Thukral12] and f4 in [Thukral13] with n_=3 & 95 respectively, with bounds here by P.G.Schild 2023
            !Computation complexity = 14
            IF (f_N < 0) THEN !output properties
                res%min_x    = -3.0_dp
                res%max_x    = 4.0_dp
                res%root     = 2.84243895378444706781658594015_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Monotonic (strictly), 42%-constant-value, 54%-constant-gradient, Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = (EXP(x) + x - 20.0_dp) ** 95
            END IF

        CASE (6)
            !REF: Example in [Brent71], Type IV in [Bus75], #13 in [Alefeld95], "Function 7" in [Chandrupatla97]; f3 in [BuzziFerraris14]
            !Computation complexity = 18
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = 4.0_dp
                res%root     = 0.0_dp
                res%mfr_typ  = -3.3_dp
                res%descript = "Monotonic (weakly, with plateaus with 1 inflection), Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = 0.0
                IF (x /= 0.0) THEN !avoid divide-by-zero
                    IF (dbl_epsilon < ABS(x)) res%fx = x * EXP(-1.0_dp / (x * x)) !avoid overflow
                END IF
            END IF
        
        CASE (7)
            !REF: f7 in [Thukral12] and f7 in [Thukral13] with n_=6 & 30 respectively, with bounds here [-1, 1.5] by P.G.Schild 2023
            !NOTE: Thukral used even-valued exponent n_, corrected here to odd-values n_=1, 3 & 29 to ensure that f(x) crosses x-axis.
            !Computation complexity = 30
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = 1.5_dp
                res%root     = 1.23941705309120820380838768087_dp
                res%mfr_typ  = -3.3_dp
                res%descript = "Non-monotonic (concave), 91%-constant-gradient, Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = (EXP(-x * x) - EXP(x * x) - x ** 8 + 10.0_dp) ** 29
            END IF
        
        CASE (8)
            !REF: f10 in [Thukral12] and f10 in [Thukral13] with n_=3 & 31 respectively with bounds by P.G.Schild 2023
            !Computation complexity = 24
            IF (f_N < 0) THEN !output properties
                res%min_x    = -3.0_dp
                res%max_x    = 6.0_dp
                res%root     = 5.46901233591014209815739616587_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Monotonic (strictly), 37%-constant-gradient, Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = (LOG(x * x + 3.0_dp * x + 5.0_dp) - 2.0_dp * x + 7.0_dp) ** 31
            END IF
        
        CASE (9)
            !REF: "Circles" in [Oliveira20]; "f13" in [Ganchovski23]
            !SHAPE: 
            !Computation complexity = 23
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = 1.0_dp
                res%root     = -1.0_dp/3.0_dp
                res%mfr_typ  = 4.3_dp
                res%descript = "Monotonic (strictly with 1 inflection), Discontinuous convexity (curvature-jolt), Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                a_ = 3 * x + 1
                res%fx = SIGN(1.0_dp,a_) * (1.0_dp - SQRT(1.0_dp - a_ * a_ / 81.0_dp))
            END IF
        
        CASE (10)
            !REF: Originated here by P.G.Schild, 2023
            !Computation complexity = 13
            IF (f_N < 0) THEN !output properties
                res%min_x    = -0.99_dp * onePi
                res%max_x    = 0.5_dp * onePi
                res%root     = 0.0_dp
                res%mfr_typ  = 3.3_dp
                res%descript = "Monotonic (strictly with 2 inflections), Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = SIGN(1.0_dp,x) * (1.0_dp - COS(x))
            END IF
        
        CASE (11)
            !REF: Originated here by P.G.Schild, 2023
            !Computation complexity = 19
            IF (f_N < 0) THEN !output properties
                res%min_x    = 0.0_dp
                res%max_x    = 1.0_dp + SQRT(0.5_dp)
                res%root     = 1.0_dp
                res%mfr_typ  = 5.3_dp
                res%descript = "Monotonic (strictly with 1 inflection), Cusp (gradient-jump) at left bound, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = (SQRT(2.0_dp * x - x * x) - 1.0_dp) * SIGN(1.0_dp,x - 1.0_dp)
            END IF
        
        CASE (12)
            !REF: Originated here by P.G.Schild, 2023
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = 0.0_dp
                res%max_x    = 10.0_dp
                res%root     = 3.6_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Monotonic (weakly, with plateaus), 68%-constant-value, Smooth, Non-linear near root"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = LOG(1.0_dp + EXP(-10.0_dp * x)) - LOG(1.0_dp + EXP(-36.0_dp))
            END IF

!-------------------------------------
!       Step-functions (discontinuous)
!-------------------------------------

        CASE (13)
            !REF: Originated here by P.G.Schild, 2023
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = onePi
                res%root     = 0.0_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Step function"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                IF (x < 0.0_dp) THEN; res%fx = -1.0_dp; ELSE; res%fx = 1.0_dp; END IF
            END IF

        CASE (14)
            !REF: Originated here by P.G.Schild, 2023
            !SHAPE: 
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = onePi
                res%root     = 0.0_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Step function "
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = 0
                IF (x < 0.0_dp) THEN; res%fx = -1.0_dp; ELSE If (0.0_dp < x) THEN; res%fx = 1.0_dp; END IF
            END IF

        CASE (15)
            !REF: Originated here by P.G.Schild, 2023
            !SHAPE: Step function 
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = onePi
                res%root     = 0.0_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Step function"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                IF (x <= 0.0_dp) THEN; res%fx = -1.0_dp; ELSE; res%fx = 1.0_dp; END IF
            END IF

        CASE (16)
            !REF: Originated here by P.G.Schild, 2023
            !SHAPE: Step function 
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = onePi
                res%root     = 1.0_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Step function"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                IF (x < 1.0_dp) THEN; res%fx = -1.0_dp; ELSE; res%fx = 1.0_dp; END IF
            END IF

        CASE (17)
            !REF: Originated here by P.G.Schild, 2023
            !SHAPE: Step function 
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = onePi
                res%root     = 1.0_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Step function"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = 0
                IF (x < 1.0_dp) THEN; res%fx = -1.0_dp; ELSE IF (1.0_dp < x) THEN; res%fx = 1.0_dp; END IF
            END IF

        CASE (18)
            !REF: Originated here by P.G.Schild, 2023
            !SHAPE: Step function 
            !Computation complexity = 26
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = onePi
                res%root     = 1.0_dp
                res%mfr_typ  = -1.3_dp
                res%descript = "Step function"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                IF (x <= 1.0_dp) THEN; res%fx = -1.0_dp; ELSE; res%fx = 1.0_dp; END IF
                END IF

!--------------------------------------
!       Reference function (super easy)
!--------------------------------------

        CASE (19)
            !REF: Originated here by P.G.Schild, 2023
            !Computation complexity = 1
            IF (f_N < 0) THEN !output properties
                res%min_x    = -1.0_dp
                res%max_x    = 3.0_dp*onePi
                res%root     = 1.01_dp
                res%mfr_typ  = 1.1_dp
                res%descript = "Linear function"
                res%latex    = "f(x) = "
            ELSE !calculate f(x)
                res%fx = x - 1.01_dp
            END IF
            
        END SELECT    
    END FUNCTION

    SUBROUTINE ShellSort(arr, Nused)
        !Shell Sort for in-place sorting of an array. Shell Sort is almost as fast as QuickSort for small arrays, with O(n^1.25)
        !Input arr(1 to Nused).
        !Advantages of ShellSort: Small code size, no use of call stack, reasonably fast, useful where memory is at a premium such as embedded and older mainframe applications.
        !Best case speed O(NlogN) and worst case O(Nlog²N) cannot be achieved together. With best case speed O(NlogN), best worst case is N^(4/3).

!      .Passed arguments:
        REAL (KIND=dp), INTENT(INOUT):: arr(*) !Left bound
        INTEGER (KIND=li), INTENT(IN):: Nused !Dimension of arr(1:Nused)
!      .Internal variables:
        REAL (KIND=dp):: tmp
        INTEGER (KIND=li):: gapi,gapn,k1,k2
        INTEGER (KIND=li), DIMENSION(9):: gaps
        gaps=(/1750, 701, 301, 132, 57, 23, 10, 4, 1/)
    
        DO gapi = 1,9 !Number of gaps in Marcin Ciura's Shell sort: 1750, 701, 301, 132, 57, 23, 10, 4, 1
            gapn=gaps(gapi)
            DO k1 = gapn + 1, Nused !This inner loop is same as insertion sort, which is the quickest sort method for small arrays (N < 30)
                tmp = arr(k1)
                k2 = k1
                DO WHILE (tmp < arr(k2 - gapn)) !+<+ for ascending sort; -<- for descending sort
                    arr(k2) = arr(k2 - gapn)
                    k2 = k2 - gapn
                    IF (k2 <= gapn) EXIT
                END DO
                arr(k2) = tmp
            END DO
        END DO
    END SUBROUTINE
        
END MODULE my_module


PROGRAM run_me
!   Test program for two bisection algorithms:
!   (1) Conventional bisection, TOMS Algorithm 4 (https://calgo.acm.org/), published in Comm ACM 3,3 (Mar 1960) p.174
!   (2) Band-bisection, i.e. that bounds the band of zero values, then assumes that the true root is in the middle
!   Author: Peter G. Schild, Oslo Metropilitan University, 2023
!
    USE my_module
    
!   Constants:
    INTEGER (KIND=li), PARAMETER :: dimEval = 10000_li
    INTEGER (KIND=li), PARAMETER :: dimFns = 19_li !number of test functions in module TestFunction_module
    INTEGER (KIND=li), PARAMETER :: dimTests = 100_li !Number of repeat tests for each function (with different bounds)
    INTEGER (KIND=li), PARAMETER :: rowRepeat = dimFns + 1_li
    CHARACTER (LEN=1), PARAMETER :: tab = CHAR(9)
!   Integer variables:
    INTEGER (KIND=li) :: iTest !Test/repeat counter
    INTEGER (KIND=li) :: iAlg !Algorithm counter 
    INTEGER (KIND=li) :: nEval !actual number of function evaluations
    INTEGER (KIND=li) :: sumEval !sum of nEval for all tests
    INTEGER (KIND=li) :: maxEval !max nEval for all tests
!   Real variables:
    REAL (KIND=dp) :: xa1,xb1,xc1,fa1,fb1,fc1
    REAL (KIND=dp) :: xa2,xb2,xc2,fa2,fb2,fc2
    REAL (KIND=dp) :: Xleft,Xright
    REAL (KIND=dp) :: X0,X1
    REAL (KIND=dp) :: ff
    REAL (KIND=dp) :: xret,fret !returned root estimate and its function f(x)
    REAL (KIND=dp) :: thisErr
    REAL (KIND=dp) :: ErrL1,ErrL2,ErrLinf
    REAL (KIND=dp) :: minXret,maxXret
    REAL (KIND=dp) :: mfr_type !function type code
    REAL (KIND=dp) :: mo_type !+ve= strictly nonotonic, -ve=non-monotonic or weakly monotonic
    REAL (KIND=dp) :: fn_type !'±1= Linear, ±2= Quadratic, ±3= Other smooth, ±4= Jolts, i.e. discontinuous f"(x), ±5=Kinks, i.e. discontinuous f'(x), ±6=Jumps, i.e. discontinuous f(x)
    REAL (KIND=dp) :: rt_type !.1= linear at root, .2= inflection at root, .3= if discontinuous (neither linear nor inflection) ar root
    TYPE(func_result) :: fDat

    OPEN(UNIT=1,FILE="TOMS4_output.xls",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
    WRITE(1,'("Func.",A,"Xroot",A,"Xleft-Xroot",A,"Xright-Xroot",A,"Alg.",A,"meanEvals",A,"maxEvals",A,"L1-norm",A,"L2-norm",A,"Linf-norm",A,"Xmin-Xroot",A,"Xmax-Xroot",A,"(Xmin-Xleft)/(Xright-Xleft)"A,"(Xmax-Xleft)/(Xright-Xleft)")') &
            tab,tab,tab,tab,tab,tab,tab,tab,tab,tab,tab,tab,tab
    CALL ConfirmMachineEpsilon() !Just to confirm that REAL64 in this PC/compiler is IEEE 754 precision; should be = 2**(-52) = 2.22E-16
    
    DO funcN = 1, dimFns

!---------------------------------
!       1. Read precise root
!---------------------------------
        fDat = TestFunction(-funcN, 1.0_dp)
        mo_type = SIGN(1.0_dp,fDat%mfr_typ)
        fn_type = INT(ABS(fDat%mfr_typ))
        rt_type = ABS(fDat%mfr_typ) - INT(ABS(fDat%mfr_typ))

!-------------------------
!   2: Bracket zero (may have finite width due to finite precision)
!-------------------------

        xa1 = fDat%min_x
        xb1 = fDat%max_x
        fRes = TestFunction(funcN, xa1); fa1=fRes%fx
        fRes = TestFunction(funcN, xb1); fb1=fRes%fx
        PRINT *,"funcN=",funcN,xa1,xb1,fa1,fb1
        IF (fa1==0.0) THEN
            STOP !The left bound is a zero. Inapproprate test function
        ELSE IF (fb1==0) THEN
            STOP !The right bound is a zero. Inapproprate test function
        ELSE IF (SIGN(1.0_dp,fa1) == SIGN(1.0_dp,fb1)) THEN
            print *, "funcn=", funcn
            print *, fDat%min_x, fDat%max_x
            print *, xa1, xb1
            print *, fa1, fb1
            print *, SIGN(1.0_dp,fa1), SIGN(1.0_dp,fb1)
            STOP !f(a) & f(b) have same sign. Inapproprate test function
        END IF
        
        !Find zero
        nEval = 2
        xc1 = xa1
        fc1 = fa1
        DO
            nEval = nEval + 1
            IF (SIGN(1.0_dp,fa1) == SIGN(1.0_dp,fc1)) THEN
                xa1 = xc1
                fa1 = fc1
            ELSE
                xb1 = xc1 !Note: no need to update fb1                
            END IF
            xc1 = xa1 + (xb1 - xa1) * 0.5_dp !Bisection, minimizes risk of overflow and underflow, and little risk of loss of precision when xba and xa1 are same magnitude
            fRes = TestFunction(funcN, xc1); fc1=fRes%fx
            IF (xc1 == xa1 .OR. xc1 == xb1 .OR. fc1 == 0) EXIT
        END DO

        IF (fc1 == 0.0_dp) THEN !zero has finite band width, so find left and right edges of zero-band
            IF (fa1 == 0.0_dp) THEN
                STOP
            ELSE IF (fb1 == 0.0_dp) THEN
                STOP
            ELSE IF (0.0_dp < fa1 * fb1) THEN
                STOP
            END IF

            !Bracket left edge of zero
            xa2 = xa1
            xc2 = xc1
            fc2 = fc1
            DO
                nEval = nEval + 1
                IF (fc2 == 0.0_dp) THEN; xb2 = xc2; ELSE; xa2 = xc2; END IF
                xc2 = xa2 + (xb2 - xa2) * 0.5_dp
                IF (xc2 == xa2 .OR. xc2 == xb2) EXIT
                fRes = TestFunction(funcN, xc2); fc2=fRes%fx
            END DO
            Xleft = xc2

            !Bracket right edge of zero
            xb2 = xb1
            xc2 = xc1
            fc2 = fc1
            DO
                nEval = nEval + 1
                IF (fc2 == 0.0_dp) THEN; xa2 = xc2; ELSE; xb2 = xc2; END IF
                xc2 = xb2 + (xa2 - xb2) * 0.5_dp
                IF (xc2 == xa2 .OR. xc2 == xb2) EXIT
                fRes = TestFunction(funcN, xc2); fc2=fRes%fx
            END DO
            Xright = xc2

            xc2 = Xleft + (Xright - Xleft) * 0.5_dp !median value
        ELSE !root has zero width
            Xleft = xa1
            xc2 = xc1
            Xright = xb1
        END IF

!----------------------------------
!   3. Test root-finding algorithms
!----------------------------------

        DO iAlg = 1, 2
            PRINT *, funcN, iAlg
            ErrL1 = 0.0_dp
            ErrL2 = 0.0_dp
            ErrLinf = 0.0_dp
            sumEval = 0_li
            maxEval = 0_li
            DO iTest = 1, dimTests
                ff = (iTest - 1.0_dp) / dimTests
                X0 = fDat%min_x + 0.5_dp * (fDat%root - fDat%min_x) * ff
                X1 = fDat%max_x - 0.5_dp * (fDat%max_x - fDat%root) * (1.0_dp - ff)
        
                SELECT CASE (iAlg)
                CASE (1); Call Algo1_Bisection(X0, X1, dbl_epsilon, dimEval, nEval, xret, fret)
                CASE (2); Call Algo2_BandBisection(X0, X1, 0.1_dp, nEval, xret, fret)
                CASE DEFAULT; STOP
                END SELECT
                
                thisErr = ABS(xret - fDat%root)
                ErrL1 = ErrL1 + (thisErr - ErrL1) / iTest
                ErrL2 = ErrL2 + (thisErr * thisErr - ErrL2) / iTest
                IF (ErrLinf < thisErr) ErrLinf = thisErr
                sumEval = sumEval + nEval
                IF (maxEval < nEval) maxEval = nEval
                IF (iTest == 1) THEN
                    minXret = xret
                    maxXret = xret
                ELSE
                    IF (xret < minXret) THEN
                        minXret = xret 
                    ELSE IF (maxXret < xret) THEN 
                        maxXret = xret
                    END IF
                END IF
            END DO !iTest loop

            IF (Xright == Xleft) THEN                    
                WRITE(1,'(I0,3(A,ES25.15E3),A,I0,A,ES25.15E3,A,I0,5(A,ES25.15E3))') &            
                    funcN,tab,fDat%root,tab,Xleft - fDat%root,tab,Xright - fDat%root,tab,iAlg,tab,sumEval/REAL(dimTests),tab,maxEval,tab,ErrL1 / dimTests,tab,SQRT(ErrL2 / (dimTests)),tab,ErrLinf,tab,minXret - fDat%root,tab,maxXret - fDat%root
            ELSE
                WRITE(1,'(I0,3(A,ES25.15E3),A,I0,A,ES25.15E3,A,I0,7(A,ES25.15E3))') &
                    funcN,tab,fDat%root,tab,Xleft - fDat%root,tab,Xright - fDat%root,tab,iAlg,tab,sumEval/REAL(dimTests),tab,maxEval,tab,ErrL1 / dimTests,tab,SQRT(ErrL2 / (dimTests)),tab,ErrLinf,tab,minXret - fDat%root,tab,maxXret - fDat%root,tab,(minXret - Xleft) / (Xright - Xleft),tab,(maxXret - Xleft) / (Xright - Xleft)
            END IF
        END DO !iAlg loop
    END DO !funcN loop
    PRINT *,'finished'
    WRITE(*,*) 'finished'
    CLOSE(1)
END PROGRAM run_me
