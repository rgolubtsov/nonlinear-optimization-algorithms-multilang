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
          INTEGER    VARS
          PARAMETER (VARS = 20)

C         RETURN VAL. THE OBJECTIVE FUNCTION VALUE.
          DOUBLE PRECISION F

C         ARG. THE POINT AT WHICH F(X) SHOULD BE EVALUATED.
          DOUBLE PRECISION X(VARS)

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
          SUBROUTINE NELMIN(  F,      N,  START,   XMIN, YNEWLO, REQMIN,
     *                     STEP, KONVGE, KCOUNT, ICOUNT, NUMRES, IFAULT)

          IMPLICIT NONE

C         CONSTANT. THE MAXIMUM NUMBER OF VARIABLES.
          INTEGER    VARS
          PARAMETER (VARS = 20)

C         CONSTANT. THE REFLECTION COEFFICIENT.
          DOUBLE PRECISION RCOEFF
          PARAMETER       (RCOEFF = 1.0D+00)

C         CONSTANT. THE EXTENSION COEFFICIENT.
          DOUBLE PRECISION ECOEFF
          PARAMETER       (ECOEFF = 2.0D+00)

C         CONSTANT. THE CONTRACTION COEFFICIENT.
          DOUBLE PRECISION CCOEFF
          PARAMETER       (CCOEFF = 0.5D+00)

C         CONSTANT. THE OPTIMALITY FACTOR.
          DOUBLE PRECISION EPS
          PARAMETER       (EPS = 0.001D+00)

C         ARG. THE OBJECTIVE FUNCTION F(X).
          DOUBLE PRECISION F
          EXTERNAL         F

C         ARG. THE NUMBER OF VARIABLES.
          INTEGER N

C         ARG. THE STARTING POINT FOR THE ITERATION.
          DOUBLE PRECISION START(VARS)

C         ARG. THE COORDINATES OF THE POINT WHICH IS ESTIMATED
C              TO MINIMIZE THE FUNCTION.
          DOUBLE PRECISION XMIN(VARS)

C         ARG. THE MINIMUM VALUE OF THE FUNCTION.
          DOUBLE PRECISION YNEWLO

C         ARG. THE TERMINATING LIMIT FOR THE VARIANCE
C              OF FUNCTION VALUES.
          DOUBLE PRECISION REQMIN

C         ARG. THE SIZE AND SHAPE OF THE INITIAL SIMPLEX.
          DOUBLE PRECISION STEP(VARS)

C         ARG. THE CONVERGENCE CHECK.
          INTEGER KONVGE

C         ARG. THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS.
          INTEGER KCOUNT

C         ARG. THE NUMBER OF FUNCTION EVALUATIONS USED.
          INTEGER ICOUNT

C         ARG. THE NUMBER OF RESTARTS.
          INTEGER NUMRES

C         ARG. THE ERROR INDICATOR.
          INTEGER IFAULT

          INTEGER JCOUNT
          INTEGER NN
          INTEGER I
          INTEGER J
          INTEGER ILO
          INTEGER IHI
          INTEGER L

          DOUBLE PRECISION DN
          DOUBLE PRECISION DNN
          DOUBLE PRECISION DEL
          DOUBLE PRECISION RQ
          DOUBLE PRECISION P(VARS, VARS + 1)
          DOUBLE PRECISION Y(VARS + 1)
          DOUBLE PRECISION X
          DOUBLE PRECISION YLO
          DOUBLE PRECISION Z
          DOUBLE PRECISION PBAR(VARS)
          DOUBLE PRECISION PSTAR(VARS)
          DOUBLE PRECISION YSTAR
          DOUBLE PRECISION P2STAR(VARS)
          DOUBLE PRECISION Y2STAR

C         CHECK THE INPUT PARAMETERS.
          IF (REQMIN .LE. 0.0D+00) THEN
              IFAULT = 1

              RETURN
          END IF

          IF (N .LT. 1) THEN
              IFAULT = 1

              RETURN
          END IF

          IF (VARS .LT. N) THEN
              IFAULT = 1

              RETURN
          END IF

          IF (KONVGE .LT. 1) THEN
              IFAULT = 1

              RETURN
          END IF

          ICOUNT = 0
          NUMRES = 0

          JCOUNT = KONVGE
          DN     = DBLE(N)
          NN     = N + 1
          DNN    = DBLE(NN)
          DEL    = 1.0D+00
          RQ     = REQMIN * DN

C         CONSTRUCTION OF INITIAL SIMPLEX.
1000      CONTINUE

          DO 10 I = 1, N
              P(I,NN) = START(I)
10        CONTINUE

          Y(NN) = F(START)

          DO 20 J = 1, N
              X        = START(J)
              START(J) = START(J) + STEP(J) * DEL

              DO 30 I = 1, N
                  P(I, J) = START(I)
30            CONTINUE

              Y(J) = F(START)

              START(J) = X
20        CONTINUE

          ICOUNT = ICOUNT + NN

C         THE SIMPLEX CONSTRUCTION IS COMPLETE.

C         FIND HIGHEST AND LOWEST Y VALUES.
C         YNEWLO = Y(IHI) INDICATES THE VERTEX OF THE SIMPLEX
C         TO BE REPLACED.
          YLO = Y(1)
          ILO = 1

          DO 40 I = 2, NN
              IF (Y(I) .LT. YLO) THEN
                  YLO = Y(I)
                  ILO = I
              END IF
40        CONTINUE

2000      CONTINUE

          YNEWLO = Y(1)
          IHI    = 1

          DO 50 I = 2, NN
              IF (YNEWLO .LT. Y(I)) THEN
                  YNEWLO = Y(I)
                  IHI    = I
              END IF
50        CONTINUE

C         CALCULATE PBAR, THE CENTROID OF THE SIMPLEX VERTICES
C         EXCEPTING THE VERTEX WITH Y VALUE YNEWLO.
          DO 60 I = 1, N
              Z = 0.0D+00

              DO 70 J = 1, NN
                  Z = Z + P(I, J)
70            CONTINUE

              Z       = Z - P(I, IHI)
              PBAR(I) = Z / DN
60        CONTINUE

C         REFLECTION THROUGH THE CENTROID.
          DO 80 I = 1, N
              PSTAR(I) = PBAR(I) + RCOEFF * (PBAR(I) - P(I, IHI))
80        CONTINUE

          YSTAR = F(PSTAR)

          ICOUNT = ICOUNT + 1

C         SUCCESSFUL REFLECTION, SO EXTENSION.
          IF (YSTAR .LT. YLO) THEN
              DO 90 I = 1, N
                  P2STAR(I) = PBAR(I) + ECOEFF * (PSTAR(I) - PBAR(I))
90            CONTINUE

              Y2STAR = F(P2STAR)

              ICOUNT = ICOUNT + 1

C             CHECK EXTENSION.
              IF (YSTAR .LT. Y2STAR) THEN
                  DO 100 I = 1, N
                      P(I, IHI) = PSTAR(I)
100               CONTINUE

                  Y(IHI) = YSTAR
C             RETAIN EXTENSION OR CONTRACTION.
              ELSE
                  DO 110 I = 1, N
                      P(I, IHI) = P2STAR(I)
110               CONTINUE

                  Y(IHI) = Y2STAR
              END IF
C         NO EXTENSION.
          ELSE
              L = 0

              DO 120 I = 1, NN
                  IF (YSTAR .LT. Y(I)) THEN
                      L = L + 1
                  END IF
120           CONTINUE

              IF (1 .LT. L) THEN
                  DO 130 I = 1, N
                      P(I, IHI) = PSTAR(I)
130               CONTINUE

                  Y(IHI) = YSTAR
C             CONTRACTION ON THE Y(IHI) SIDE OF THE CENTROID.
              ELSE IF (L .EQ. 0) THEN
                  DO 140 I = 1, N
                      P2STAR(I) = PBAR(I)
     *                          + CCOEFF * (P(I, IHI) - PBAR(I))
140               CONTINUE

                  Y2STAR = F(P2STAR)

                  ICOUNT = ICOUNT + 1

C                 CONTRACT THE WHOLE SIMPLEX.
                  IF (Y(IHI) .LT. Y2STAR) THEN
                      DO 150 J = 1, NN
                          DO 160 I = 1, N
                              P(I, J) = (P(I, J) + P(I, ILO)) * 0.5D+00
                              XMIN(I) =  P(I, J)
160                       CONTINUE

                          Y(J) = F(XMIN)
150                   CONTINUE

                      ICOUNT = ICOUNT + NN

                      IF (KCOUNT .LT. ICOUNT) THEN
                          GO TO 3000
                      END IF

                      YLO = Y(1)
                      ILO = 1

                      DO 170 I = 2, NN
                          IF (Y(I) .LT. YLO) THEN
                              YLO = Y(I)
                              ILO = I
                          END IF
170                   CONTINUE

                      GO TO 2000
C                 RETAIN CONTRACTION.
                  ELSE
                      DO 180 I = 1, N
                          P(I, IHI) = P2STAR(I)
180                   CONTINUE

                      Y(IHI) = Y2STAR
                  END IF
C             CONTRACTION ON THE REFLECTION SIDE OF THE CENTROID.
              ELSE IF (L .EQ. 1) THEN
                  DO 190 I = 1, N
                      P2STAR(I) = PBAR(I)
     *                          + CCOEFF * (PSTAR(I) - PBAR(I))
190               CONTINUE

                  Y2STAR = F(P2STAR)

                  ICOUNT = ICOUNT + 1

C                 RETAIN REFLECTION?
                  IF (Y2STAR .LE. YSTAR) THEN
                      DO 200 I = 1, N
                          P(I, IHI) = P2STAR(I)
200                   CONTINUE

                      Y(IHI) = Y2STAR
                  ELSE
                      DO 210 I = 1, N
                          P(I, IHI) = PSTAR(I)
210                   CONTINUE

                      Y(IHI) = YSTAR
                  END IF
              END IF
          END IF

C         CHECK IF YLO IMPROVED.
          IF (Y(IHI) .LT. YLO) THEN
              YLO = Y(IHI)
              ILO = IHI
          END IF

          JCOUNT = JCOUNT - 1

          IF (JCOUNT .NE. 0) THEN
              GO TO 2000
          END IF

C         CHECK TO SEE IF MINIMUM REACHED.
          IF (ICOUNT .LE. KCOUNT) THEN
              JCOUNT = KONVGE
              Z      = 0.0D+00

              DO 220 I = 1, NN
                  Z = Z + Y(I)
220           CONTINUE

              X = Z / DNN
              Z = 0.0D+00

              DO 230 I = 1, NN
                  Z = Z + (Y(I) - X) ** 2
230           CONTINUE

              IF (RQ .LT. Z) THEN
                  GO TO 2000
              END IF
          END IF

C         FACTORIAL TESTS TO CHECK THAT YNEWLO IS A LOCAL MINIMUM.
3000      CONTINUE

          DO 240 I = 1, N
              XMIN(I) = P(I, ILO)
240       CONTINUE

          YNEWLO = Y(ILO)

          IF (KCOUNT .LT. ICOUNT) THEN
              IFAULT = 2

              RETURN
          END IF

          IFAULT = 0

          DO 250 I = 1, N
              DEL     = STEP(I) * EPS
              XMIN(I) = XMIN(I) + DEL

              Z = F(XMIN)

              ICOUNT = ICOUNT + 1

              IF (Z .LT. YNEWLO) THEN
                  IFAULT = 2

                  GO TO 4000
              END IF

              XMIN(I) = XMIN(I) - DEL - DEL

              Z = F(XMIN)

              ICOUNT = ICOUNT + 1

              IF (Z .LT. YNEWLO) THEN
                  IFAULT = 2

                  GO TO 4000
              END IF

              XMIN(I) = XMIN(I) + DEL
250       CONTINUE

4000      CONTINUE

          IF (IFAULT .EQ. 0) THEN
              RETURN
          END IF

C         RESTART THE PROCEDURE.
          DO 260 I = 1, N
              START(I) = XMIN(I)
260       CONTINUE

          DEL    = EPS
          NUMRES = NUMRES + 1

          GO TO 1000
      END

C     === MAIN PROGRAM.
      PROGRAM AMOEBA
          IMPLICIT NONE

C         CONSTANT. THE MAXIMUM NUMBER OF VARIABLES.
          INTEGER    VARS
          PARAMETER (VARS = 20)

          INTEGER N
          INTEGER KONVGE
          INTEGER KCOUNT
          INTEGER I
          INTEGER ICOUNT
          INTEGER NUMRES
          INTEGER IFAULT

          DOUBLE PRECISION START(VARS)
          DOUBLE PRECISION REQMIN
          DOUBLE PRECISION STEP(VARS)
          DOUBLE PRECISION YNEWLO
          DOUBLE PRECISION XMIN(VARS)

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

          CALL NELMIN(   F,      N,  START,   XMIN, YNEWLO, REQMIN,
     *                STEP, KONVGE, KCOUNT, ICOUNT, NUMRES, IFAULT)

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
