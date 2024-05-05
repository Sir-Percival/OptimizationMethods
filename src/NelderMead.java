import java.util.Arrays;
import java.util.Comparator;
import java.util.Scanner;
import java.util.stream.IntStream;

public class NelderMead
{
    private static final double ALPHA = 1.0; // reflection coefficient
    private static final double BETA = 0.5; // contraction coefficient
    private static final double GAMMA = 2.0; // expansion coefficient
    private static final double DELTA = 0.5; // reduction coefficient
    private static final double EPSILON = 0.00000001; // tolerance

    private interface Function
    {
        double calculate(double[] x);
    }

    public static void main(String[] args)
    {
        int n = 11;
        Function function2D = new Function() {
            @Override
            public double calculate(double[] x) {
                return Math.pow(2-x[0], 2) + (n+2)*Math.pow(x[1]-Math.pow(x[0],2), 2);
            }
        };

        Function function3D = new Function() {
            @Override
            public double calculate(double[] x) {
                return Math.pow(1-x[0], 2) + n*Math.pow(x[1]-Math.pow(x[0],2), 2) + Math.pow(1-x[1], 2) + n*Math.pow(x[2]-Math.pow(x[1],2), 2);
            }
        };

        Scanner scanner = new Scanner(System.in);
        int dimensions = 0;
        Function function = null;

        while (true)
        {
            System.out.print("Number of dimensions: ");
            dimensions = scanner.nextInt();
            if(dimensions == 2)
            {
                function = function2D;
                break;
            }
            else if (dimensions == 3)
            {
                function = function3D;
                break;
            }
            else
            {
                System.out.println("Enter 2 or 3!");
            }
        }

        minimize(function, dimensions);
    }

    private static double[][] sortPoints(double[] values, double[][] points)
    {
        Comparator<Integer> comparator = Comparator.comparingDouble(i -> values[i]);
        int[] indices = IntStream.range(0, values.length).boxed().sorted(comparator).mapToInt(i -> i).toArray();

        Arrays.sort(values);

        double[][] sortedPoints = new double[points.length][];
        for (int i = 0; i < points.length; i++) {
            sortedPoints[i] = points[indices[i]];
        }

        return sortedPoints;
    }

    private static void minimize(Function function, int dimensions)
    {
        // ========== 0 - Build initial simplex =========

        double[][] simplex = buildSimplex(dimensions);
        double[] fValues = new double[dimensions + 1];

        while (true)
        {
            // ========= 1 - Numbering of points =========
            for(int i = 0; i <= dimensions; i++)
            {
                fValues[i] = function.calculate(simplex[i]);
            }

            simplex = sortPoints(fValues, simplex);

            // ========= Calculate centroid =========

            double[] centroid = new double[dimensions];
            for(int i=0; i<dimensions; i++)
            {
                centroid[i] = 0.0;
                for(int j=0; j<dimensions; j++)
                {
                    centroid[i] += simplex[j][i];
                }
                centroid[i] /= dimensions;
            }

            // ========= 2 - Find reflection point =========

            double[] reflected = new double[dimensions];
            for(int i=0; i<dimensions; i++)
            {
                reflected[i] = (1 + ALPHA) * centroid[i] - ALPHA * simplex[simplex.length - 1][i];
            }
            double fReflected = function.calculate(reflected);

            // ========= 3 - If f(xr) < f(x1) find xe - expansion =========

            if(fReflected < fValues[0])
            {
                double[] expanded = new double[dimensions];
                for(int i=0; i<dimensions; i++)
                {
                    expanded[i] = GAMMA * reflected[i] + (1 - GAMMA) * centroid[i];
                }
                double fExpanded = function.calculate(expanded);

                // ========= 3A - If f(xe) < f(x1) - xn+1 = xe
                if(fExpanded < fValues[0])
                {
                    simplex[simplex.length - 1] = expanded;
                    fValues[fValues.length - 1] = fExpanded;
                }
                else    // ========= 3B - If f(xe) >= f(x1) - xn+1 = xr
                {
                    simplex[simplex.length - 1] = reflected;
                    fValues[fValues.length - 1] = fReflected;
                }
            }
            // ========= 4 - If f(x1) <= f(xr) <= f(xn) - xn+1 = xr =========
            else if (fValues[0] <= fReflected && fReflected <= fValues[fValues.length-2])
            {
                simplex[simplex.length - 1] = reflected;
                fValues[fValues.length - 1] = fReflected;
            }
            // ========= 5 - If f(xn) < f(xr) < f(xn+1) - xn+1 = xr + Contraction =========
            else if(fValues[fValues.length-2] < fReflected && fReflected < fValues[fValues.length-1])
            {
                simplex[simplex.length - 1] = reflected;
                fValues[fValues.length - 1] = fReflected;

                double[] contraction = new double[dimensions];
                for(int i=0; i<dimensions; i++)
                {
                    contraction[i] = BETA * simplex[simplex.length-1][i] + (1-BETA)*centroid[i];
                }
                double fContraction = function.calculate(contraction);

                // ========== 5A - If f(xc) < f(xn+1) - xn+1 = xc ==========
                if(fContraction < fValues[fValues.length-1])
                {
                    simplex[simplex.length - 1] = contraction;
                    fValues[fValues.length - 1] = fContraction;
                }
                // ========== 5B - Else - Reduction ==========
                else
                {
                    for(int i=0; i <= dimensions; i++)
                    {
                        for(int j=0; j < dimensions; j++)
                        {
                            simplex[i][j] = simplex[i][j] + simplex[0][j] * DELTA;
                        }
                        fValues[i] = function.calculate(simplex[i]);
                    }
                }
            }
            // ========== 6 - If f(xr) >= f(xn+1) - Contraction ==========
            else
            {
                double[] contraction = new double[dimensions];
                for(int i=0; i<dimensions; i++)
                {
                    contraction[i] = BETA * simplex[simplex.length-1][i] + (1-BETA)*centroid[i];
                }
                double fContraction = function.calculate(contraction);

                // ========== 6A - If f(xc) < f(xn+1) - xn+1 = xc ==========
                if(fContraction < fValues[fValues.length-1])
                {
                    simplex[simplex.length - 1] = contraction;
                    fValues[fValues.length - 1] = fContraction;
                }
                // ========== 6B - Else - Reduction ==========
                else
                {
                    for(int i=0; i <= dimensions; i++)
                    {
                        for(int j=0; j < dimensions; j++)
                        {
                            simplex[i][j] = simplex[i][j] + simplex[0][j] * DELTA;
                        }
                        fValues[i] = function.calculate(simplex[i]);
                    }
                }
            }

            // ========== 7 - Check stop criterion ==========
            double value = 0.0;
            double fCentroid = function.calculate(centroid);

            for (int i=0; i <= dimensions; i++)
            {
                value += Math.pow(fValues[i] - fCentroid, 2);
            }

            value = Math.sqrt(value / (dimensions+1));

            if(value <= EPSILON)
            {
                break;
            }
        }

        double[] minimum = simplex[0];
        System.out.println("Minimum point: " + Arrays.toString(minimum));
        System.out.printf("Minimum value: %.15f%n", fValues[0]);
    }

    private static double[][] buildSimplex(int dimensions)
    {
        System.out.println("===== CREATE INITIAL SIMPLEX =====");
        Scanner scanner = new Scanner(System.in);
        double[] firstPoint = new double[dimensions];
        for(int i=0; i<dimensions; i++)
        {
            System.out.print("Enter value #" + i + ": ");
            firstPoint[i] = scanner.nextDouble();
        }

        double[][] simplex = new double[dimensions + 1][dimensions];
        for (int i = 0; i <= dimensions; i++) {
            simplex[i] = firstPoint.clone();
            if (i > 0) {
                simplex[i][i - 1] += 1.0;
            }
        }

//        System.out.println("========= INITIAL SIMPLEX ========");
//        for(var a : simplex)
//        {
//            System.out.println(Arrays.toString(a));
//        }
//        System.out.println("==================================");

        return simplex;
    }
}
