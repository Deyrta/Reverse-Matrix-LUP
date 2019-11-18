using System;
using static System.Console;
using static System.Math;

namespace ConsoleApp6
{
    class Program
    {
        static void Main(string[] args)
        {
            SLAR_Gauss m = new SLAR_Gauss();
            // m.SetParametes();
            m.DoTriangle();
            m.Revers();
            m.OutputNewMatrix();
            SLAR_Gauss.InvertMatrix(); // send matrix here...
            ReadKey();
        }
    }
    class SLAR_Gauss
    {
        static int n = 3;//= int.Parse(ReadLine());
        double[,] Matrix = new double[3, 4] { { 10, -7, 0, 7 }, { -3, 2, 6, 4 }, { 5, -1, 5, 6 } };
        double[] Temporary = new double[n + 1];
        double[] Result = new double[n];

        double MaxElement = 0;
        double koef;

        int I = 1;
        int J = 0;
        int MaxI;
        int Count = 1;

        public static double[] LUPSolve(double[][] LU, int[] pi, double[] b)
        {
            int n = LU.Length - 1;
            double[] x = new double[n + 1];
            double[] y = new double[n + 1];
            double suml = 0;
            double sumu = 0;
            double lij = 0;

            /*
            * Solve for y using formward substitution
            * */
            for (int i = 0; i <= n; i++)
            {
                suml = 0;
                for (int j = 0; j <= i - 1; j++)
                {
                    /*
                    * Since we've taken L and U as a singular matrix as an input
                    * the value for L at index i and j will be 1 when i equals j, not LU[i][j], since
                    * the diagonal values are all 1 for L.
                    * */
                    if (i == j)
                    {
                        lij = 1;
                    }
                    else
                    {
                        lij = LU[i][j];
                    }
                    suml = suml + (lij * y[j]);
                }
                y[i] = b[pi[i]] - suml;
            }
            //Solve for x by using back substitution
            for (int i = n; i >= 0; i--)
            {
                sumu = 0;
                for (int j = i + 1; j <= n; j++)
                {
                    sumu = sumu + (LU[i][j] * x[j]);
                }
                x[i] = (y[i] - sumu) / LU[i][i];
            }
            return x;
        }

        public static Tuple<double[][], int[]> LUPDecomposition(double[][] A)
        {
            int n = A.Length - 1;
            /*
            * pi represents the permutation matrix.  We implement it as an array
            * whose value indicates which column the 1 would appear.  We use it to avoid 
            * dividing by zero or small numbers.
            * */
            int[] pi = new int[n + 1];
            double p = 0;
            int kp = 0;
            int pik = 0;
            int pikp = 0;
            double aki = 0;
            double akpi = 0;

            //Initialize the permutation matrix, will be the identity matrix
            for (int j = 0; j <= n; j++)
            {
                pi[j] = j;
            }

            for (int k = 0; k <= n; k++)
            {
                /*
                * In finding the permutation matrix p that avoids dividing by zero
                * we take a slightly different approach.  For numerical stability
                * We find the element with the largest 
                * absolute value of those in the current first column (column k).  If all elements in
                * the current first column are zero then the matrix is singluar and throw an
                * error.
                * */
                p = 0;
                for (int i = k; i <= n; i++)
                {
                    if (Math.Abs(A[i][k]) > p)
                    {
                        p = Math.Abs(A[i][k]);
                        kp = i;
                    }
                }
                if (p == 0)
                {
                    throw new Exception("singular matrix");
                }
                /*
                * These lines update the pivot array (which represents the pivot matrix)
                * by exchanging pi[k] and pi[kp].
                * */
                pik = pi[k];
                pikp = pi[kp];
                pi[k] = pikp;
                pi[kp] = pik;

                /*
                * Exchange rows k and kpi as determined by the pivot
                * */
                for (int i = 0; i <= n; i++)
                {
                    aki = A[k][i];
                    akpi = A[kp][i];
                    A[k][i] = akpi;
                    A[kp][i] = aki;
                }

                /*
                    * Compute the Schur complement
                    * */
                for (int i = k + 1; i <= n; i++)
                {
                    A[i][k] = A[i][k] / A[k][k];
                    for (int j = k + 1; j <= n; j++)
                    {
                        A[i][j] = A[i][j] - (A[i][k] * A[k][j]);
                    }
                }
            }
            return Tuple.Create(A, pi);
        }

        public static double[][] InvertMatrix(double[][] A)
        {
            int n = A.Length;
            //e will represent each column in the identity matrix
            double[] e;
            //x will hold the inverse matrix to be returned
            double[][] x = new double[n][];
            for (int i = 0; i < n; i++)
            {
                x[i] = new double[A[i].Length];
            }
            /*
            * solve will contain the vector solution for the LUP decomposition as we solve
            * for each vector of x.  We will combine the solutions into the double[][] array x.
            * */
            double[] solve;

            //Get the LU matrix and P matrix (as an array)
            Tuple<double[][], int[]> results = LUPDecomposition(A);

            double[][] LU = results.Item1;
            int[] P = results.Item2;

            /*
            * Solve AX = e for each column ei of the identity matrix using LUP decomposition
            * */
            for (int i = 0; i < n; i++)
            {
                e = new double[A[i].Length];
                e[i] = 1;
                solve = LUPSolve(LU, P, e);
                for (int j = 0; j < solve.Length; j++)
                {
                    x[j][i] = solve[j];
                }
            }
            return x;
        }


        public void SetParametes()
        {
            WriteLine("Enter the parameters");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    WriteLine("e[{0}][{1}]", i, j);
                    Matrix[i, j] = double.Parse(ReadLine());
                }
            }
        }

        public void SearchMaxElement()
        {
            MaxElement = Matrix[I - 1, J];
            for (int i = I - 1; i < n; i++)
            {
                if (Abs(MaxElement) < Abs(Matrix[i, J]))
                {
                    MaxI = i;
                    MaxElement = Matrix[i, J];
                }
            }
        }

        public void Change()
        {
            for (int j = J; j < n + 1; j++)
            {
                Temporary[j] = Matrix[I - 1, j];
                Matrix[I - 1, j] = Matrix[MaxI, j];
                Matrix[MaxI, j] = Temporary[j];
            }
        }
        public void DoTriangle()
        {
            while (Count < n + 1)
            {
                SearchMaxElement();
                Change();
                for (int i = I; i < n; i++)
                {
                    koef = Matrix[I, J] / Matrix[J, J];
                    for (int j = J; j < n + 1; j++)
                    {
                        Matrix[i, j] = Matrix[i, j] - Matrix[J, j] * koef;
                    }
                    I++;
                }
                Count++;
                J++;
                I = J + 1;
            }
        }
        public void OutputNewMatrix()
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    Write(Matrix[i, j]);
                    Write("\t");
                }
                WriteLine("\n");
            }
            for (int j = 0; j < n; j++)
            {
                Write("x{0}: ", j);
                WriteLine(Round(Result[j], 6));
            }

        }

        public void Revers()
        {
            for (int i = n - 1; i >= 0; i--)
            {
                double summ = 0;
                for (int j = i + 1; j < n; j++)
                {
                    summ += Matrix[i, j] * Result[j];
                }
                summ = Matrix[i, n] - summ;
                Result[i] = summ / Matrix[i, i];
            }
        }
    }
}
