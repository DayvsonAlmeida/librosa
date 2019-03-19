using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace FftPack{
    class Transform{
        // Implementação realizada apenas do tipo 2.
        // TODO: Implementar tipo 1, 3 e 4.
        /// <summary>
        /// Transformação Discreta de Cosseno 2-D
        /// </summary>
        /// <param name="x">Array 2D de entrada</param>
        /// <param name="axis">Eixo de normalização</param>
        /// <param name="type">Tipo da transformada {1,2,3,4}</param>
        /// <param name="norm">null ou ortho Default: NULL</param>
        /// <returns></returns>
        public static double[,] DCT2D(double[,] x, int axis = 0, int type = 2, string norm = null)
        {
            int rows = x.GetLength(0), columns = x.GetLength(1), N;
            double[,] y;
            double sum;
            y = new double[rows, columns];

            
            //Implementação do tipo 2
            //Transformar nas colunas
            if (axis == 0)
            {
                N = rows;
                for (int j = 0; j < columns; j++)
                {
                    //Itera sobre todas as casas da coluna j de Y
                    for (int k = 0; k < rows; k++)
                    {
                        sum = 0.0;
                        //Itera sobre todas as casas da coluna j de X
                        for (int n = 0; n < rows; n++)
                            sum += x[n, j] * Math.Cos(Math.PI * k * (2.0 * n + 1) / (2.0 * N));
                        sum *= 2;
                        y[k, j] = sum;
                        if (norm != null && norm.Equals("ortho"))
                        {
                            if (k == 0 )
                                y[k, j] *= Math.Sqrt(1.0 / (4.0 * N));
                            else
                                y[k, j] *= Math.Sqrt(1.0 / (2.0 * N));
                        }
                    }
                }
            }
            //Transformar nas linhas
            else
            {
                N = columns;
                for (int i = 0; i < rows; i++)
                {
                    //Itera sobre todas as casas da coluna j de Y
                    for (int k = 0; k < columns; k++)
                    {
                        sum = 0.0;
                        //Itera sobre todas as casas da coluna j de X
                        for (int n = 0; n < columns; n++)
                            sum += x[i, n] * Math.Cos(Math.PI * k * (2.0 * n + 1) / (2.0 * N));
                        sum *= 2;
                        y[i, k] = sum;
                        if (norm != null && norm.Equals("ortho"))
                        {
                            if (k == 0)
                                y[i, k] *= Math.Sqrt(1.0 / (4.0 * N));
                            else
                                y[i, k] *= Math.Sqrt(1.0 / (2.0 * N));
                        }
                    }
                }
            }
            return y;
        }

        // Implementação realizada apenas do tipo 2.
        // TODO: Implementar tipo 1, 3 e 4.
        /// <summary>
        /// Transformação Discreta de Cosseno 1-D
        /// </summary>
        /// <param name="x">Array de Entrada</param>
        /// <param name="type">Tipo da transformada {1,2,3,4}</param>
        /// <param name="norm">null ou ortho Default: NULL</param>
        /// <returns></returns>
        public static double[] DCT(double[] x, int type = 2, string norm = null)
        {
            double[] y;
            int N = x.Length;

            y = new double[N];
            double sum;
            for (int k = 0; k < N; k++)
            {
                sum = 0.0;
                for (int n = 0; n < N; n++)
                    sum += x[n] * Math.Cos(Math.PI * k * (2.0 * n + 1) / (2.0 * N));
                sum *= 2;
                y[k] = sum;
                if (norm != null && norm.Equals("ortho"))
                {
                    if (k == 0)
                        y[k] *= Math.Sqrt(1.0 / (4.0 * N));
                    else
                        y[k] *= Math.Sqrt(1.0 / (2.0 * N));
                }
            }

            return y;
        }
    }
}
