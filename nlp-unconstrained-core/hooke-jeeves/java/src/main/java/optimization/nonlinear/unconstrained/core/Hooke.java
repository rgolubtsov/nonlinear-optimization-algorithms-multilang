/*
 * nlp-unconstrained-core/hooke-jeeves/java/src/main/java/
 * optimization/nonlinear/unconstrained/core/Hooke.java
 * ============================================================================
 * Nonlinear Optimization Algorithms Multilang. Version 0.1
 * ============================================================================
 * Nonlinear programming algorithms as the (un-)constrained minimization
 * problems with the focus on their numerical expression using various
 * programming languages.
 *
 * This is the Hooke and Jeeves nonlinear unconstrained minimization algorithm.
 * ============================================================================
 */

package optimization.nonlinear.unconstrained.core;

public class Hooke {
    public static final int MAXIMUM_NUMBER_OF_VARIABLES = 250;

    public static final double ENDING_VALUE_OF_STEPSIZE = 1E-6;

    public static final int MAXIMUM_NUMBER_OF_ITERATIONS = 5000;

    public static final int INDEX_ZERO = 0;
    public static final int INDEX_ONE = 1;
    private static final double ZERO_POINT_FIVE = 0.5;

    private int numberOfFunctionEvaluations = 0;

    /**
     * Helper method.
     * <br />
     * <br />Given a point, look for a better one nearby, one coord at a time.
     *
     * @param delta                        The delta between <code>previousBestValuedCoordinate</code>
     *                                     and <code>point</code>.
     * @param point                        The coordinate from where to begin.
     * @param previousBestValuedCoordinate The previous best-valued coordinate.
     * @param numberOfVariables            The number of variables.
     * @param objFunCls                    The class in which the objective function is defined.
     * @return The objective function value at a nearby.
     */
    private double bestNearby(double[] delta,
                              double[] point,
                              double previousBestValuedCoordinate,
                              int numberOfVariables,
                              ObjectiveFunction objFunCls) {

        double minumuValueOfFunction;
        double[] currentSearchPoint = new double[MAXIMUM_NUMBER_OF_VARIABLES];
        double currentFunctionValue;

        int iterationNumber;

        minumuValueOfFunction = previousBestValuedCoordinate;

        for (iterationNumber = 0; iterationNumber < numberOfVariables; iterationNumber++) {
            currentSearchPoint[iterationNumber] = point[iterationNumber];
        }

        for (iterationNumber = 0; iterationNumber < numberOfVariables; iterationNumber++) {
            currentSearchPoint[iterationNumber] = point[iterationNumber] + delta[iterationNumber];
            currentFunctionValue = objFunCls.findValueForArguments(currentSearchPoint);
            numberOfFunctionEvaluations++;

            if (currentFunctionValue < minumuValueOfFunction) {
                minumuValueOfFunction = currentFunctionValue;
            } else {
                delta[iterationNumber] = 0.0 - delta[iterationNumber];
                currentSearchPoint[iterationNumber] = point[iterationNumber] + delta[iterationNumber];
                currentFunctionValue = objFunCls.findValueForArguments(currentSearchPoint);
                numberOfFunctionEvaluations++;

                if (currentFunctionValue < minumuValueOfFunction) {
                    minumuValueOfFunction = currentFunctionValue;
                } else {
                    currentSearchPoint[iterationNumber] = point[iterationNumber];
                }
            }
        }

        for (iterationNumber = 0; iterationNumber < numberOfVariables; iterationNumber++) {
            point[iterationNumber] = currentSearchPoint[iterationNumber];
        }

        return minumuValueOfFunction;
    }

    public int findMinimum(int numberOfVariables,
                           double[] startingPointCoordinates,
                           double[] endingPointCoordinates,
                           double rhoStepsizeGeometricShrink,
                           double epsilonEndingValueOfStepsize,
                           int maximumNumberOfIterations,
                           ObjectiveFunction objectiveFunction) {


        double[] newFunctionArguments = new double[MAXIMUM_NUMBER_OF_VARIABLES];
        double[] previousFunctionArguments = new double[MAXIMUM_NUMBER_OF_VARIABLES];
        double[] delta = new double[MAXIMUM_NUMBER_OF_VARIABLES];
        double previousFunctionValue;

        for (int currentCoordinateNumber = 0; currentCoordinateNumber < numberOfVariables; currentCoordinateNumber++) {
            previousFunctionArguments[currentCoordinateNumber] = startingPointCoordinates[currentCoordinateNumber];
            newFunctionArguments[currentCoordinateNumber] = previousFunctionArguments[currentCoordinateNumber];

            delta[currentCoordinateNumber] = Math.abs(startingPointCoordinates[currentCoordinateNumber] * rhoStepsizeGeometricShrink);

            if (delta[currentCoordinateNumber] == 0.0) {
                delta[currentCoordinateNumber] = rhoStepsizeGeometricShrink;
            }
        }



        int numberOfIterations = 0;

        previousFunctionValue = objectiveFunction.findValueForArguments(newFunctionArguments);
        numberOfFunctionEvaluations++;

        double newFunctionValue = previousFunctionValue;
        double stepLength = rhoStepsizeGeometricShrink;
        while ((numberOfIterations < maximumNumberOfIterations) && (stepLength > epsilonEndingValueOfStepsize)) {
            numberOfIterations++;

            System.out.printf(
                    "\nAfter %5d funevals, f(x) =  %.4e at\n", numberOfFunctionEvaluations, previousFunctionValue
            );

            for (int  currentNumberOfVariable = 0; currentNumberOfVariable < numberOfVariables; currentNumberOfVariable++) {
                System.out.printf("   x[%2d] = %.4e\n", currentNumberOfVariable, previousFunctionArguments[currentNumberOfVariable]);
            }

            findBestNewPoint(numberOfVariables, newFunctionArguments, previousFunctionArguments);

            newFunctionValue = bestNearby(delta, newFunctionArguments, previousFunctionValue, numberOfVariables, objectiveFunction);

            // If we made some improvements, pursue that direction.
            int keep = 1;

            NewDirectionPursuer newDirectionPursuer = new NewDirectionPursuer(numberOfVariables, objectiveFunction, newFunctionArguments, previousFunctionArguments, delta, previousFunctionValue, newFunctionValue, keep).invoke();
            newFunctionValue = newDirectionPursuer.getNewFunctionValue();
            previousFunctionValue = newDirectionPursuer.getPreviousFunctionValue();

            if ((stepLength >= epsilonEndingValueOfStepsize) && (newFunctionValue >= previousFunctionValue)) {
                stepLength = stepLength * rhoStepsizeGeometricShrink;

                for (int currentVariableIndex = 0; currentVariableIndex < numberOfVariables; currentVariableIndex++) {
                    delta[currentVariableIndex] *= rhoStepsizeGeometricShrink;
                }
            }
        }

        findBestNewPoint(numberOfVariables, endingPointCoordinates, previousFunctionArguments);

        return numberOfIterations;
    }

    private void findBestNewPoint(int numberOfVariables, double[] newFunctionArguments, double[] previousFunctionArguments) {
        for (int newPointCoordinates = 0; newPointCoordinates < numberOfVariables; newPointCoordinates++) {
            newFunctionArguments[newPointCoordinates] = previousFunctionArguments[newPointCoordinates];
        }
    }

    private class NewDirectionPursuer {
        private int numberOfVariables;
        private ObjectiveFunction objectiveFunction;
        private double[] newFunctionArguments;
        private double[] previousFunctionArguments;
        private double[] delta;
        private double previousFunctionValue;
        private double newFunctionValue;
        private int keep;

        public NewDirectionPursuer(int numberOfVariables, ObjectiveFunction objectiveFunction, double[] newFunctionArguments, double[] previousFunctionArguments, double[] delta, double previousFunctionValue, double newFunctionValue, int keep) {
            this.numberOfVariables = numberOfVariables;
            this.objectiveFunction = objectiveFunction;
            this.newFunctionArguments = newFunctionArguments;
            this.previousFunctionArguments = previousFunctionArguments;
            this.delta = delta;
            this.previousFunctionValue = previousFunctionValue;
            this.newFunctionValue = newFunctionValue;
            this.keep = keep;
        }

        public double getPreviousFunctionValue() {
            return previousFunctionValue;
        }

        public double getNewFunctionValue() {
            return newFunctionValue;
        }

        public NewDirectionPursuer invoke() {
            double tmp;
            while ((newFunctionValue < previousFunctionValue) && (keep == 1)) {

                for (int i = 0; i < numberOfVariables; i++) {
                    // Firstly, arrange the sign of delta[].
                    if (newFunctionArguments[i] <= previousFunctionArguments[i]) {
                        delta[i] = 0.0 - Math.abs(delta[i]);
                    } else {
                        delta[i] = Math.abs(delta[i]);
                    }

                    // Now, move further in this direction.
                    tmp = previousFunctionArguments[i];
                    previousFunctionArguments[i] = newFunctionArguments[i];
                    newFunctionArguments[i] = newFunctionArguments[i] + newFunctionArguments[i] - tmp;
                }

                previousFunctionValue = newFunctionValue;

                newFunctionValue = bestNearby(delta, newFunctionArguments, previousFunctionValue, numberOfVariables, objectiveFunction);

                // If the further (optimistic) move was bad....
                if (newFunctionValue >= previousFunctionValue) {
                    break;
                }

                /*
                 * Make sure that the differences between the new and the old
                 * points are due to actual displacements; beware of roundoff
                 * errors that might cause newFunctionValue < previousFunctionValue.
                 */
                keep = 0;

                for (int i = 0; i < numberOfVariables; i++) {
                    keep = 1;

                    if (Math.abs(newFunctionArguments[i] - previousFunctionArguments[i])
                            > (ZERO_POINT_FIVE * Math.abs(delta[i]))) {

                        break;
                    } else {
                        keep = 0;
                    }
                }
            }
            return this;
        }
    }
}
