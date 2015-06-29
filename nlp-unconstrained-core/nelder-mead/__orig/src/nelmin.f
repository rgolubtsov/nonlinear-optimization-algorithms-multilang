C =============================================================================
C NLP-UNCONSTRAINED-CORE/NELDER-MEAD/__ORIG/SRC/NELMIN.F
C =============================================================================
C NONLINEAR OPTIMIZATION ALGORITHMS MULTILANG. VERSION 0.1
C =============================================================================
C NONLINEAR PROGRAMMING ALGORITHMS AS THE (UN-)CONSTRAINED MINIMIZATION
C PROBLEMS WITH THE FOCUS ON THEIR NUMERICAL EXPRESSION USING VARIOUS
C PROGRAMMING LANGUAGES.
C
C THIS IS THE NELDER-MEAD NONLINEAR UNCONSTRAINED MINIMIZATION ALGORITHM.
C =============================================================================

C     === THE USER-SUPPLIED OBJECTIVE FUNCTION F(X).
      FUNCTION F(X)
          IMPLICIT NONE

C         CONSTANT. THE MAXIMUM NUMBER OF VARIABLES.
          INTEGER    NV
          PARAMETER (NV = 250)

C         RETURN VAL. THE OBJECTIVE FUNCTION VALUE.
          DOUBLE PRECISION F

C         ARG. THE POINT AT WHICH F(X) SHOULD BE EVALUATED.
          DOUBLE PRECISION X(NV)

#ifndef WOODS
C     ROSENBROCK'S CLASSIC PARABOLIC VALLEY ('BANANA') FUNCTION.
          DOUBLE PRECISION A
          DOUBLE PRECISION B

          A = X(2) - X(1) * X(1)
          B = 1    - X(1)

          F = 100 * (A * A) + (B * B)
#else
C     WOODS -- A LA MORE, GARBOW AND HILLSTROM (TOMS ALGORITHM 566).
          DOUBLE PRECISION S1
          DOUBLE PRECISION S2
          DOUBLE PRECISION S3
          DOUBLE PRECISION T1
          DOUBLE PRECISION T2
          DOUBLE PRECISION T3
          DOUBLE PRECISION T4
          DOUBLE PRECISION T5

          S1 = X(2) - X(1) * X(1)
          S2 = 1    - X(1)
          S3 = X(2) - 1

          T1 = X(4) - X(3) * X(3)
          T2 = 1    - X(3)
          T3 = X(4) - 1

          T4 = S3 + T3
          T5 = S3 - T3

          F = 100 * (S1 * S1) + (S2 * S2)
     *       + 90 * (T1 * T1) + (T2 * T2)
     *       + 10 * (T4 * T4) + (T5 * T5) / 10
#endif
          RETURN
      END

C     === MAIN OPTIMIZATION SUBROUTINE.
C     THE NELMIN SUBROUTINE ITSELF (NELDER-MEAD MINIMIZATION).
      SUBROUTINE NELMIN(F, N, START, XMIN, YNEWLO, REQMIN, STEP, KONVGE,
     *                  KCOUNT, ICOUNT, NUMRES, IFAULT)

          IMPLICIT NONE

C         CONSTANT. THE MAXIMUM NUMBER OF VARIABLES.
          INTEGER    NV
          PARAMETER (NV = 250)

C         ARG. THE OBJECTIVE FUNCTION F(X).
          DOUBLE PRECISION F
          EXTERNAL         F

C         ARG. THE NUMBER OF VARIABLES.
          INTEGER N

C         ARG. START
          DOUBLE PRECISION START(NV)

C         ARG. XMIN
          DOUBLE PRECISION XMIN(NV)

C         ARG. YNEWLO
          DOUBLE PRECISION YNEWLO

C         ARG. REQMIN
          DOUBLE PRECISION REQMIN

C         ARG. STEP
          DOUBLE PRECISION STEP(NV)

C         ARG. KONVGE
          INTEGER KONVGE

C         ARG. KCOUNT
          INTEGER KCOUNT

C         ARG. ICOUNT
          INTEGER ICOUNT

C         ARG. NUMRES
          INTEGER NUMRES

C         ARG. IFAULT
          INTEGER IFAULT

C         TODO: IMPLEMENT THE SUBROUTINE BODY.

          RETURN
      END

C     === MAIN PROGRAM.
      PROGRAM AMOEBA
          IMPLICIT NONE

C         CONSTANT. THE MAXIMUM NUMBER OF VARIABLES.
          INTEGER    NV
          PARAMETER (NV = 250)

          INTEGER N
          INTEGER KONVGE
          INTEGER KCOUNT
          INTEGER I
          INTEGER ICOUNT
          INTEGER NUMRES
          INTEGER IFAULT

          DOUBLE PRECISION START(NV)
          DOUBLE PRECISION REQMIN
          DOUBLE PRECISION STEP(NV)
          DOUBLE PRECISION YNEWLO
          DOUBLE PRECISION XMIN(NV)

C         PROTO REF. THE OBJECTIVE FUNCTION F(X).
          DOUBLE PRECISION F

#ifndef WOODS
C     STARTING GUESS FOR ROSENBROCK'S TEST FUNCTION.
          PRINT 10
10        FORMAT (/,   'TEST01',
     *            /2x, 'Apply NELMIN to ROSENBROCK function.')

          N        =  2
          START(1) = -1.2
          START(2) =  1.0
#else
C     STARTING GUESS TEST PROBLEM 'WOODS'.
          PRINT 50
50        FORMAT (/,   'TEST05',
     *            /2x, 'Apply NELMIN to WOODS function.')

          N        =  4
          START(1) = -3.0
          START(2) = -1.0
          START(3) = -3.0
          START(4) = -1.0
#endif

          REQMIN  = 1.0D-08
          STEP(1) = 1.0
          STEP(2) = 1.0
          KONVGE  = 10
          KCOUNT  = 500

          PRINT 60
60        FORMAT (/2x, 'Starting point X:', /)

          DO 70 I = 1, N
              PRINT 80, START(I)
80            FORMAT (2x, g14.6)
70        CONTINUE

          YNEWLO = F(START)

          PRINT 90, YNEWLO
90        FORMAT (/2x, 'F(X) =', 1x, g14.6)

          CALL NELMIN(F, N, START, XMIN, YNEWLO, REQMIN, STEP, KONVGE,
     *                KCOUNT, ICOUNT, NUMRES, IFAULT)

          PRINT 100, IFAULT
100       FORMAT (/2x, 'Return code IFAULT =', 1x, i8)

          PRINT 110
110       FORMAT (/2x, 'Estimate of minimizing value X*:', /)

          DO 120 I = 1, N
              PRINT 130, XMIN(I)
130           FORMAT (2x, g14.6)
120       CONTINUE

          PRINT 140, YNEWLO
140       FORMAT (/2x, 'F(X*) =', 1x, g14.6)

          PRINT 150, ICOUNT
150       FORMAT (/2x, 'Number of iterations =', 1x, i8)

          PRINT 160, NUMRES
160       FORMAT (2x, 'Number of restarts   =', 1x, i8)
      END

C =============================================================================
