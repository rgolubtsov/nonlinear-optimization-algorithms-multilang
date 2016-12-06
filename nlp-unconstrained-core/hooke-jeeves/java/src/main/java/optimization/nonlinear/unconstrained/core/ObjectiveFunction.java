package optimization.nonlinear.unconstrained.core;

public interface ObjectiveFunction {

    default public double findValueForArguments(double[] parameters){
        return 0.0;
    }
}
