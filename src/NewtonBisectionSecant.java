import java.util.function.Function;

public class NewtonBisectionSecant
{
    public static void main(String[] args)
    {
        Function<Double, Double> f1 = (Double x) -> Math.pow(x, 3) - (2 * Math.pow(x, 2)) + 7 * x - 9;
        Function<Double, Double> f1Prime = (Double x) -> Math.pow(x, 2) * 3 - 4 * x + 7;

//        double zero = Newton(1, f1, f1Prime, 0.000001, 1000);
//        zero = Bisection(1, 5, f1, 0.000001, 1000);
//        zero = Secant(f1, 1, 5, 0.000001, 1000);

        Function<Double, Double> f2 = (Double x) -> 2 * Math.sin(x) - Math.log(x);
        Function<Double, Double> f2Prime = (Double x) -> 2 * Math.cos(x) - (1/x);

//        double zero = Newton(2, f2, f2Prime, 0.000001, 1000);
//        zero = Bisection(2, 3, f2, 0.000001, 1000);
//        zero = Secant(f2, 2, 3, 0.000001, 1000);

        Function<Double, Double> f3 = (Double x) -> (1 / (1 + Math.pow(Math.atan(x), 2))) - Math.pow(Math.E, Math.sqrt(x));
        Function<Double, Double> f3Prime = (Double x) -> ((-2*Math.atan(x))/((Math.pow(x, 2)+1)*Math.pow((1 + Math.pow(Math.atan(x), 2)), 2)))
                - (Math.pow(Math.E, Math.sqrt(x))/(2*Math.sqrt(x)));

//        double zero = Bisection(0, 1, f3, 0.000001, 1000);
//        zero = Newton(0, f3, f3Prime, 0.000001, 1000);
//        zero = Secant(f3, 0, 1, 0.000001, 1000);

        Function<Double, Double> f4 = (Double x) -> Math.atan(Math.log(Math.pow(x, 7)+5)-9);
        Function<Double, Double> f4Prime = (Double x) -> (7*Math.pow(x,6))/((Math.pow(x,7)+5) * ((Math.pow(Math.log(Math.pow(x,7)+5)-9,2))+1));

//        double zero = Bisection(3, 4, f4, 0.000001, 1000);
//        zero = Newton(3, f4, f4Prime, 0.000001, 1000);
//        zero = Secant(f4, 3, 4, 0.000001, 1000);

        Function<Double, Double> f5 = (Double x) -> Math.sin(1/x)-2*x;
        Function<Double, Double> f5Prime = (Double x) -> (-(Math.cos(1/x)/Math.pow(x, 2))) - 2;

        double zero = Bisection(-0.2, -0.12, f5, 0.000001, 1000);
        zero = Newton(-0.2, f5, f5Prime, 0.000001, 1000);
        zero = Secant(f5, -0.2, -0.12, 0.000001, 1000);

    }

    private static double Newton(double x0, Function<Double, Double> f, Function<Double, Double> fPrime, double epsilon, int max_iterations)
    {
        System.out.println("========== NEWTON ==========");
        for(int i=0; i<max_iterations; i++)
        {
            double y = f.apply(x0);
            double yPrime = fPrime.apply(x0);

            System.out.printf("Iter: %d, x: %f, f: %f%n", i+1, x0, y);

            if(Math.abs(yPrime) < 0)
            {
                System.out.println("Nie znaleziono");
                break;
            }

            double x1 = x0 - (y/yPrime);

            if(Math.abs(x1 - x0) <= epsilon)
            {
                System.out.println("Miejsce zerowe: " + x1);
                return x1;
            }

            x0 = x1;
        }

        System.out.println("Nie znaleziono");
        return 0.0;
    }

    private static double Bisection(double start, double stop, Function<Double, Double> f, double epsilon, int max_iterations)
    {
        System.out.println("========== BISECTION ==========");
        for(int i=0; i<max_iterations; i++)
        {
            double x1 = (start + stop) / 2;

            System.out.printf("Iter: %d, x: %f, f: %f%n", i+1, x1, f.apply(x1));

            if(Math.abs(f.apply(x1)) <= epsilon)
            {
                System.out.println("Miejsce zerowe: " + x1);
                return x1;
            }
            else if(f.apply(start) * f.apply(x1) < 0)
            {
                stop = x1;
            }
            else
            {
                start = x1;
            }
        }
        System.out.println("Nie znaleziono");
        return 0.0;
    }

    private static double Secant(Function<Double, Double> f, double x0, double x1, double epsilon, int max_iterations)
    {
        System.out.println("========== SECANT ==========");
        for(int i=0; i<max_iterations; i++)
        {
            double x2 = (f.apply(x1)*x0 - f.apply(x0)*x1) / (f.apply(x1) - f.apply(x0));
            double f2 = f.apply(x2);

            System.out.printf("Iter: %d, x: %f, f: %f%n", i+1, x2, f2);

            if(Math.abs(f2) < epsilon)
            {
                System.out.println("Miejsce zerowe: " + x2);
                return x2;
            }

            x0 = x1;
            x1 = x2;
        }

        System.out.println("Nie znaleziono");
        return 0.0;
    }
}