using System.Diagnostics;
using System.Collections.Concurrent;

namespace RelaxationMethod
{
    class Program
    {
        private static Matrix A;
        private static Vector b;
        private static int n = 1000;
        private static double epsilon = 0.01;

        class Matrix
        {
            public double[,] matrix;
            private int lines;
            private int columns;

            public Matrix(int lines, int columns)
            {
                if (lines < 1 || columns < 1)
                {
                    throw new Exception("Неверный размер");
                }
                this.lines = lines;
                this.columns = columns;
                matrix = new double[lines, columns];
            }

            public Matrix(Matrix init) : this(init.lines, init.columns)
            {
                for (int i = 0; i < lines; i++)
                {
                    for (int j = 0; j < columns; j++)
                    {
                        this.matrix[i, j] = init.matrix[i, j];
                    }
                }
            }

            public int getLines()
            {
                return this.lines;
            }
            public int getColumns()
            {
                return this.columns;
            }

            public void swap(int fi, int fj, int si, int sj)
            {
                double tmp = matrix[fi, fj];
                matrix[fi, fj] = matrix[si, sj];
                matrix[si, sj] = tmp;
            }

            public void generateMatrix(int n)
            {
                double[,] matrix = new double[n, n];
                Random random = new Random();
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        int randInt = random.Next(1, 5);
                        matrix[i, j] = randInt * random.NextDouble();
                    }
                }

                this.matrix = matrix;
                this.columns = n;
                this.lines = n;
            }

            public Vector mul(Vector vector)
            {
                if (columns != vector.getLength())
                {
                    throw new Exception("Неверная матрица или вектор");
                }
                Vector result = new Vector(vector.getLength());
                for (int i = 0; i < lines; i++)
                {
                    result.vector[i] = 0;
                    for (int j = 0; j < columns; j++)
                    {
                        result.vector[i] += matrix[i,j] * vector.vector[j];
                    }
                }
                return result;
            }

            public Matrix mul(Matrix mtr)
            {
                if (columns != mtr.getLines())
                {
                    throw new Exception("Неверная матрица");
                }
                Matrix result = new Matrix(lines, mtr.getColumns());
                for (int i = 0; i < result.getLines(); i++)
                {
                    for (int j = 0; j < result.getColumns(); j++)
                    {
                        result.matrix[i,j] = 0;
                        for (int k = 0; k < columns; k++)
                        {
                            result.matrix[i,j] += this.matrix[i,k] * mtr.matrix[k,j];
                        }
                    }
                }
                return result;
            }

            public Matrix transpose()
            {
                if (lines != columns)
                {
                    throw new Exception("Неверная матрица");
                }
                Matrix result = new Matrix(this);
                for (int i = 0; i < lines; i++)
                {
                    for (int j = i + 1; j < columns; j++)
                    {
                        result.swap(i, j, j, i);
                    }
                }
                return result;
            }
        }


        class Vector
        {
            public double[] vector;
            private int length;

            public Vector(int length)
            {
                if (length < 1)
                {
                    throw new Exception("Неверный размер");
                }
                this.length = length;
                vector = new double[length];
            }

            public Vector(Vector init) : this(init.getLength())
            {
                for (int i = 0; i < length; i++)
                {
                    this.vector[i] = init.vector[i];
                }
            }

            public int getLength()
            {
                return length;
            }

            public void print(bool exponent)
            {
                foreach (double v in vector)
                {
                    if (exponent)
                    {
                        Console.WriteLine("%e\n", v);
                    }
                    else
                    {
                        Console.WriteLine(v);
                    }
                }
            }

            public void generateVector(int n)
            {
                double[] vector = new double[n];
                Random random = new Random();
                for (int i = 0; i < n; i++)
                {
                    int randInt = random.Next(1, 5);
                    vector[i] = randInt * random.NextDouble();
                }
                this.vector = vector;
                this.length = n;
            }

            public Vector subtract(Vector sub)
            {
                if (length != sub.getLength())
                {
                    throw new Exception("Неверный вектор");
                }
                Vector result = new Vector(length);
                for (int i = 0; i < length; i++)
                {
                    result.vector[i] = this.vector[i] - sub.vector[i];
                }
                return result;
            }

            // параллельная реализация
            public Vector subtract(Vector sub, int threadNum)
            {
                ParallelOptions Options = new ParallelOptions();
                if (length != sub.getLength())
                {
                    throw new Exception("Неверный вектор");
                }
                Vector result = new Vector(length);
               
                Options.MaxDegreeOfParallelism = threadNum;
                Parallel.ForEach(Partitioner.Create(0, length), Options, (range, loop) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        result.vector[i] = this.vector[i] - sub.vector[i];
                    }
                });

                return result;
            }

            public double normI()
            {
                double max = Math.Abs(vector[0]);
                for (int i = 1; i < length; i++)
                {
                    if (Math.Abs(vector[i]) > max)
                    {
                        max = Math.Abs(vector[i]);
                    }
                }
                return max;
            }

            // параллельная реализация
            public double normI(int threadNum)
            {
                ParallelOptions Options = new ParallelOptions();
                double max = Math.Abs(vector[0]);
                Options.MaxDegreeOfParallelism = threadNum;
                Parallel.For(1, length, Options, i =>
                {
                    if (Math.Abs(vector[i]) > max)
                    {
                        max = Math.Abs(vector[i]);
                    }
                });


                return max;
            }
        }

        // последовательная реализация
        private static Vector upperRelaxation()
        {
            Vector prevX;
            Vector nextX = new Vector(n);
            double sum;
            double omega = 1.1;
            int counter = 0;
            for (int i = 0; i < n; i++)
            {
                nextX.vector[i] = b.vector[i] / A.matrix[i,i];
            }
            do
            {
                prevX = new Vector(nextX);
                for (int i = 0; i < n; i++)
                {
                    sum = 0;
                    nextX.vector[i] = (1 - omega) * prevX.vector[i];
                    for (int j = 0; j < i; j++)
                    {
                        sum += (A.matrix[i,j] / A.matrix[i,i]) * nextX.vector[j];
                    }
                    nextX.vector[i] -= omega * sum;
                    sum = 0;
                    for (int j = i + 1; j < n; j++)
                    {
                        sum += (A.matrix[i,j] / A.matrix[i,i]) * prevX.vector[j];
                    }
                    nextX.vector[i] -= omega * sum;
                    nextX.vector[i] += omega * (b.vector[i] / A.matrix[i,i]);
                }
                counter++;
            } while (nextX.subtract(prevX).normI() > omega * epsilon);
            Console.WriteLine("Количество итераций: " + counter);
            return nextX;
        }

        public static double Add(ref double location1, double value)
        {
            double newCurrentValue = location1;
            while (true)
            {
                double currentValue = newCurrentValue;
                double newValue = currentValue + value;
                newCurrentValue = Interlocked.CompareExchange(ref location1, newValue, currentValue);
                if (newCurrentValue == currentValue)
                    return newValue;
            }
        }

        // параллельная реализация
        private static Vector upperRelaxationParallel(int threadNum)
        {
            Vector prevX;
            Vector nextX = new Vector(n);
            double omega = 1.1;
            int counter = 0;
            ParallelOptions Options = new ParallelOptions();
            Options.MaxDegreeOfParallelism = threadNum;
            Parallel.For(0, n, Options, i =>
            {
                nextX.vector[i] = b.vector[i] / A.matrix[i, i];
            });

            do
            {
                prevX = new Vector(nextX);
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    nextX.vector[i] = (1 - omega) * prevX.vector[i];
                    Parallel.For<double>(0, i, Options, () => 0, (j, loop, subtotal) =>
                    {
                        subtotal += (A.matrix[i, j] / A.matrix[i, i]) * nextX.vector[j];
                        return subtotal;
                    },
                    (subtotal) => Add(ref sum, subtotal));
                    nextX.vector[i] -= omega * sum;
                    sum = 0;
                    Parallel.For<double>(i + 1, n, Options, () => 0, (j, loop, subtotal) =>
                      {
                          subtotal += (A.matrix[i, j] / A.matrix[i, i]) * prevX.vector[j];
                          return subtotal;
                      },
                    (subtotal) => Add(ref sum, subtotal));
                    nextX.vector[i] -= omega * sum;
                    nextX.vector[i] += omega * (b.vector[i] / A.matrix[i, i]);

                }
                counter++;
            } while (nextX.subtract(prevX, threadNum).normI(threadNum) > omega * epsilon);
            Console.WriteLine("Количество итераций: " + counter);
            return nextX;
        }
       
        public static void Main(string[] args)
        {
            Vector x;
            Vector xPar;
            Stopwatch Watch = new Stopwatch();
            TimeSpan Time;

            Stopwatch Watch1 = new Stopwatch();
            TimeSpan Time1;

            try
            {
                A = new Matrix(n, n);
                b = new Vector(n);
                A.generateMatrix(n);
                b.generateVector(n);
                b = A.transpose().mul(b);
                A = A.transpose().mul(A);
                Watch.Start();
                x = upperRelaxation();
                Watch.Stop();
                Time = Watch.Elapsed;
                Console.WriteLine
                        (
                            "\nRunTime: " +
                            String.Format
                            (
                                "{0:00}:{1:00}:{2:00}.{3:00}",
                                Time.Hours,
                                Time.Minutes,
                                Time.Seconds,
                                Time.Milliseconds / 10
                            )
                        );
                Console.WriteLine("Ticks: " + Time.Ticks);
                Console.WriteLine("Вектор Х:");
                x.print(false);

                Console.WriteLine();
                Watch1.Start();
                xPar = upperRelaxationParallel(Environment.ProcessorCount);
                Console.WriteLine("Время параллельного выполнения");
                Watch1.Stop();
                Time1 = Watch1.Elapsed;
                Console.WriteLine
                        (
                            "\nRunTime: " +
                            String.Format
                            (
                                "{0:00}:{1:00}:{2:00}.{3:00}",
                                Time.Hours,
                                Time.Minutes,
                                Time.Seconds,
                                Time.Milliseconds / 10
                            )
                        );
                Console.WriteLine("Ticks: " + Time1.Ticks);
                Console.WriteLine("Вектор X:");
                xPar.print(false);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.ToString());
            }
        }
    }
}