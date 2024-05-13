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
    private static final double EPSILON = 1e-8; // tolerance

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

        int iteration = 1;

        while (true)
        {
            // ========= 1 - Numbering of points =========
            for(int i = 0; i <= dimensions; i++)
            {
                fValues[i] = function.calculate(simplex[i]);
            }

            System.out.println("========== ITERATION #" + iteration + " ==========");
            System.out.println("Simplex: " + Arrays.deepToString(simplex));
            System.out.println("fValues: " + Arrays.toString(fValues));

            simplex = sortPoints(fValues, simplex);


            System.out.println("Sorted simplex:");
            System.out.println(Arrays.deepToString(simplex));
            System.out.println("Sorted fValues: " + Arrays.toString(fValues));

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

            System.out.println("===== CENTROID =====");
            System.out.println("x0 = " + Arrays.toString(centroid));

            // ========= 2 - Find reflection point =========

            double[] reflected = new double[dimensions];
            for(int i=0; i<dimensions; i++)
            {
                reflected[i] = (1 + ALPHA) * centroid[i] - ALPHA * simplex[simplex.length - 1][i];
            }
            double fReflected = function.calculate(reflected);

            System.out.println("===== REFLECTION =====");
            System.out.println("Reflection point: " + Arrays.toString(reflected) + " = " + fReflected);

            // ========= 3 - If f(xr) < f(x1) find xe - expansion =========

            if(fReflected < fValues[0])
            {
                System.out.println("===== (3) f(x0) < f(x1) -> EXPANSION =====");

                double[] expanded = new double[dimensions];
                for(int i=0; i<dimensions; i++)
                {
                    expanded[i] = GAMMA * reflected[i] + (1 - GAMMA) * centroid[i];
                }
                double fExpanded = function.calculate(expanded);

                System.out.println("Expansion point: " + Arrays.toString(expanded) + " = " + fExpanded);

                // ========= 3A - If f(xe) < f(x1) - xn+1 = xe
                if(fExpanded < fValues[0])
                {
                    simplex[simplex.length - 1] = expanded;
                    fValues[fValues.length - 1] = fExpanded;
                    System.out.println("===== (3A) f(xe) < f(x1) -> xn+1 = xe =====");
                    System.out.println("New simplex:");
                    System.out.println(Arrays.deepToString(simplex));
                    System.out.println("New fValues: " + Arrays.toString(fValues));
                }
                else    // ========= 3B - If f(xe) >= f(x1) - xn+1 = xr
                {
                    simplex[simplex.length - 1] = reflected;
                    fValues[fValues.length - 1] = fReflected;
                    System.out.println("===== (3B) f(xe) >= f(x1) -> xn+1 = xr =====");
                    System.out.println("New simplex:");
                    System.out.println(Arrays.deepToString(simplex));
                    System.out.println("New fValues: " + Arrays.toString(fValues));
                }
            }
            // ========= 4 - If f(x1) <= f(xr) <= f(xn) - xn+1 = xr =========
            else if (fValues[0] <= fReflected && fReflected <= fValues[fValues.length-2])
            {
                simplex[simplex.length - 1] = reflected;
                fValues[fValues.length - 1] = fReflected;
                System.out.println("===== (4) f(x1) <= f(xr) <= f(xn) -> xn+1 = xr =====");
                System.out.println("New simplex:");
                System.out.println(Arrays.deepToString(simplex));
                System.out.println("New fValues: " + Arrays.toString(fValues));
            }
            // ========= 5 - If f(xn) < f(xr) < f(xn+1) - xn+1 = xr + Contraction =========
            else if(fValues[fValues.length-2] < fReflected && fReflected < fValues[fValues.length-1])
            {
                simplex[simplex.length - 1] = reflected;
                fValues[fValues.length - 1] = fReflected;

                System.out.println("===== (5) f(xn) < f(xr) < f(xn+1) -> xn+1 = xr =====");
                System.out.println("New simplex:");
                System.out.println(Arrays.deepToString(simplex));
                System.out.println("New fValues: " + Arrays.toString(fValues));

                double[] contraction = new double[dimensions];
                for(int i=0; i<dimensions; i++)
                {
                    contraction[i] = BETA * simplex[simplex.length-1][i] + (1-BETA)*centroid[i];
                }
                double fContraction = function.calculate(contraction);

                System.out.println("===== CONTRACTION =====");
                System.out.println("Contraction point: " + Arrays.toString(contraction) + " = " + fContraction);

                // ========== 5A - If f(xc) < f(xn+1) - xn+1 = xc ==========
                if(fContraction < fValues[fValues.length-1])
                {
                    System.out.println("===== (5A) f(xc) < f(xn+1) -> xn+1 = xc =====");
                    simplex[simplex.length - 1] = contraction;
                    fValues[fValues.length - 1] = fContraction;
                    System.out.println("New simplex:");
                    System.out.println(Arrays.deepToString(simplex));
                    System.out.println("New fValues: " + Arrays.toString(fValues));
                }
                // ========== 5B - Else - Reduction ==========
                else
                {
                    System.out.println("===== (5B) REDUCTION =====");
                    for(int i=0; i <= dimensions; i++)
                    {
                        for(int j=0; j < dimensions; j++)
                        {
                            simplex[i][j] = simplex[i][j] + simplex[0][j] * DELTA;
                        }
                        fValues[i] = function.calculate(simplex[i]);
                    }

                    System.out.println("New simplex:");
                    System.out.println(Arrays.deepToString(simplex));
                    System.out.println("New fValues: " + Arrays.toString(fValues));
                }
            }
            // ========== 6 - If f(xr) >= f(xn+1) - Contraction ==========
            else
            {
                System.out.println("===== (6) f(xr) >= f(xn+1) -> CONTRACTION =====");
                double[] contraction = new double[dimensions];
                for(int i=0; i<dimensions; i++)
                {
                    contraction[i] = BETA * simplex[simplex.length-1][i] + (1-BETA)*centroid[i];
                }
                double fContraction = function.calculate(contraction);

                System.out.println("Contraction point: " + Arrays.toString(contraction) + " = " + fContraction);

                // ========== 6A - If f(xc) < f(xn+1) - xn+1 = xc ==========
                System.out.println("===== (6A) f(xc) < f(xn+1) -> xn+1 = xc =====");
                if(fContraction < fValues[fValues.length-1])
                {
                    simplex[simplex.length - 1] = contraction;
                    fValues[fValues.length - 1] = fContraction;
                    System.out.println("New simplex:");
                    System.out.println(Arrays.deepToString(simplex));
                    System.out.println("New fValues: " + Arrays.toString(fValues));
                }
                // ========== 6B - Else - Reduction ==========
                else
                {
                    System.out.println("===== (6B) REDUCTION =====");
                    for(int i=0; i <= dimensions; i++)
                    {
                        for(int j=0; j < dimensions; j++)
                        {
                            simplex[i][j] = simplex[i][j] + simplex[0][j] * DELTA;
                        }
                        fValues[i] = function.calculate(simplex[i]);
                    }

                    System.out.println("New simplex:");
                    System.out.println(Arrays.deepToString(simplex));
                    System.out.println("New fValues: " + Arrays.toString(fValues));
                }
            }

            // ========== 7 - Check stop criterion ==========
            double value = 0.0;
            double fCentroid = function.calculate(centroid);

            for (int i=0; i <= dimensions; i++)
            {
                value += Math.pow(fValues[i] - fCentroid, 2);
            }

            System.out.println("===== STOP CRITERION =====");

            value = Math.sqrt(value / (dimensions+1));

            System.out.println("Best point: " + Arrays.toString(simplex[0]) + " = " + fValues[0]);
            System.out.println("Value = " + value + " (< " + EPSILON + ")");
            System.out.println("========== END OF ITERATION ==========\n");
            iteration++;

            if(value <= EPSILON)
            {
                break;
            }
        }

        double[] minimum = simplex[0];
        System.out.println("Final simplex: " + Arrays.deepToString(simplex));
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
