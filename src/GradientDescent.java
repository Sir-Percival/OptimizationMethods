public class GradientDescent
{
    public static void main(String[] args)
    {
//        simpleGradientDescent(1.4, 1.0, 0.1, 0.000001, 1000);   // #1
//        simpleGradientDescent(0.0, 0.0, 0.1, 0.000001, 1000);   // #2
        simpleGradientDescent(1.0, 1.0, 0.01, 0.000001, 1000);   // #3
    }

    private static void simpleGradientDescent(double x0, double y0, double learning_rate, double epsilon, int max_iterations)
    {
        double gradientX;
        double gradientY;

        for(int i=1; i<=max_iterations; i++)
        {
            gradientX = df_dx(x0, y0);
            gradientY = df_dy(x0, y0);

            double nextX = x0 + learning_rate * -gradientX;
            double nextY = y0 + learning_rate * -gradientY;

            System.out.printf("%d\t%f\t%f\t%f\t%f%n", i, x0, y0, f(x0, y0), learning_rate);

            double gradientNorm = Math.sqrt(gradientX * gradientX + gradientY * gradientY);
            double positionDifferenceNorm = Math.sqrt(Math.pow(nextX - x0, 2) + Math.pow(nextY - y0, 2));

            if(gradientNorm <= epsilon || positionDifferenceNorm <= epsilon)
            {
                break;
            }

            if(f(nextX, nextY) >= f(x0, y0))
            {
                learning_rate *= 0.5;
            }
            else
            {
                x0 = nextX;
                y0 = nextY;
            }
        }

        System.out.println("Znalezione minimum:");
        System.out.printf("f(%f,%f) = %f", x0, y0, f(x0, y0));
    }

    private static double f(double x, double y)
    {
//        return Math.pow(x-2*y, 2) + Math.pow(x-2, 4);   // #1
//        return Math.sin(x) * Math.cos(y);   // #2
        return - Math.pow(Math.tan(x*y), 2);
    }

    private static double df_dx(double x, double y)
    {
//        return 2 * (x-2*y) + 4 * Math.pow(x-2, 3);  // #1
//        return Math.cos(x) * Math.cos(y);   // #2
        return (-2 * y * Math.tan(x*y)) / Math.pow(Math.cos(x*y),2);
    }

    private static double df_dy(double x, double y)
    {
//        return -4 * (x-2*y);    // #1
//        return - Math.sin(x)*Math.sin(y);   // #2
        return (-2 * x * Math.tan(x*y)) / Math.pow(Math.cos(x*y),2);
    }
}
