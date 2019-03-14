using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
namespace Librosa
{
    class Program
    {

        static void Main(string[] args)
        {
            double[][] m = Time_frequency.mel(sr:16000,n_fft:2048,n_mels:128);

            Console.Read();

        }


    }

    class Time_frequency
    {
        public static IEnumerable<double> Arange(double start, int count)
        {
            return Enumerable.Range((int)start, count).Select(v => (double)v);
        }

        /// <summary>
        ///     Retorna um vetor <double> com espaçamento linear entre o início e o fim especificados
        /// </summary>
        /// <param name="start"></param>
        /// <param name="stop"></param>
        /// <param name="num"></param>
        /// <param name="endpoint"></param>
        /// <returns></returns>
        public static IEnumerable<double> LinSpace(double start, double stop, int num, bool endpoint = true)
        {
            var result = new List<double>();
            if (num <= 0)
            {
                return result;
            }

            if (endpoint)
            {
                if (num == 1)
                {
                    return new List<double>() { start };
                }

                var step = (stop - start) / ((double)num - 1.0d);
                result = Arange(0, num).Select(v => (v * step) + start).ToList();
            }
            else
            {
                var step = (stop - start) / (double)num;
                result = Arange(0, num).Select(v => (v * step) + start).ToList();
            }

            return result;
        }

        /// <summary>
        ///     Vetor com espaçamento linear de frequencias
        /// </summary>
        /// <param name="sr"></param>
        /// <param name="n_fft"></param>
        /// <returns></returns>
        public static double[] fft_frequencies(double sr = 16000, int n_fft = 2048)
        {

            return LinSpace(0, sr / 2, (int)(1 + n_fft / 2), true).ToArray<double>();

        }

        /// <summary>
        ///     hertz para mels (double para double)
        /// </summary>
        /// <param name="frequencies"></param>
        /// <param name="htk"></param>
        /// <returns></returns>
        public static double hz_to_mel(double frequencies, bool htk = false)
        {
            if (htk)
            {
                return 2595 * Math.Log10(1 + frequencies / 700);
            }

            //numeros magicos pt1
            double f_min = 0;
            double f_sp = 200 / 3;
            double mels = (frequencies - f_min) / f_sp;

            double min_log_hz = 1000;
            double min_log_mel = (min_log_hz - f_min) / f_sp;

            double logstep = Math.Log(6.4) / 27;         

            if (frequencies >= min_log_hz)
            {
                mels = min_log_mel + Math.Log(frequencies / min_log_hz) / logstep;
            }
            return mels;
        }

        /// <summary>
        ///     mels (vetor de double) para hertz (double) 
        /// </summary>
        /// <param name="mels"></param>
        /// <param name="htk"></param>
        /// <returns></returns>
        public static double[] mel_to_hz(double[] mels, bool htk = false)
        {
            int i = 0;

            if (htk)
            {
                double[] m = new double[mels.Length];
                foreach (double d in mels)
                {
                    double v = 700 * (Math.Pow(10, d / 2595) - 1);
                    m[i] = v;
                    i++;
                }
                return m;
            }

            //numeros magicos pt2
            double f_min = 0;
            double f_sp = 200 / 3;
            double[] freqs = new double[mels.Length];
            foreach (double d in mels)
            {
                double v = f_min + f_sp * d;
                freqs[i] = v;
                i++;
            }


            double min_log_hz = 1000;
            double min_log_mel = (min_log_hz - f_min) / f_sp;

            double logstep = Math.Log(6.4) / 27; 

            if (mels.Length > 1)
            {
                i = 0;
                foreach (double d in mels)
                {
                    if (d >= min_log_mel)
                    {
                        freqs[i] = min_log_hz * Math.Exp(logstep * (mels[i] - min_log_mel));
                    }
                    i++;
                }
            }
            /*if (mels >= min_log_mel)
            {
                freqs = min_log_hz * Math.Exp(logstep * (mels - min_log_mel));
            }*/

            return freqs;
        }

        /// <summary>
        ///     Vetor de mels
        /// </summary>
        /// <param name="n_mels"></param>
        /// <param name="fmin"></param>
        /// <param name="fmax"></param>
        /// <param name="htk"></param>
        /// <returns></returns>
        public static double[] mel_frequencies(int n_mels = 128, double fmin = 0, double fmax = 11025, bool htk = false)
        {
            double min_mel = hz_to_mel(fmin, htk);
            double max_mel = hz_to_mel(fmax, htk);

            double[] mels = LinSpace(min_mel, max_mel, n_mels).ToArray<double>();

            return mel_to_hz(mels, htk);
        }

        /// <summary>
        ///     retorna um vetor de tamanho length-1, tal que: out[i]=a[i+1]-a[i]
        /// </summary>
        /// <param name="vet"></param>
        /// <returns></returns>
        public static double[] diff(double[] vet)
        {
            double[] v = new double[vet.Length - 1];
            for (int i = 0; i < vet.Length - 1; i++)
            {
                v[i] = vet[i + 1] - vet[i];
            }

            return v;
        }

        /// <summary>
        ///     numpy subtract outer: saida[i,j]
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static double[,] subtract_outer(double[] v1, double[] v2)
        {
            double[,] result = new double[v1.Length, v2.Length];
            for (int i = 0; i < v1.Length; i++)
            {
                for (int j = 0; j < v2.Length; j++)
                {
                    result[i, j] = v1[i] - v2[j];
                }
            }

            return result;
        }
        /// <summary>
        ///     overload de subtract_outer(double[], double[]) 
        /// </summary>
        /// <param name="vet"></param>
        /// <returns></returns>
        public static double[,] subtract_outer(double[] vet)
        {
            return subtract_outer(vet, vet);
        }

        public static double[] maximum(double[] a1, double[] a2)
        {
            double[] res = new double[a1.Length];
            if (a1.Length != a2.Length)
            {
                return null;
            }
            for (int i = 0; i < a1.Length; i++)
            {
                if (a1[i] > a2[i])
                {
                    res[i] = a1[i];
                }
                else
                {
                    res[i] = a2[i];
                }
            }

            return res;
        }
        public static double maximum(double v1, double v2)
        {
            if (v1 > v2)
            {
                return v1;
            }
            else
            {
                return v2;
            }
        }
        public static double[] maximum(double[] ar, double v)
        {
            double[] res = new double[ar.Length];
            for (int i = 0; i < ar.Length; i++)
            {
                if (ar[i] > v)
                {
                    res[i] = ar[i];
                }
                else
                {
                    res[i] = v;
                }
            }
            return res;
        }
        public static double[] maximum(double v, double[] ar)
        {
            double[] res = new double[ar.Length];
            for (int i = 0; i < ar.Length; i++)
            {
                if (ar[i] > v)
                {
                    res[i] = ar[i];
                }
                else
                {
                    res[i] = v;
                }
            }
            return res;
        }
        public static double[] minimum(double[] a1, double[] a2)
        {
            double[] res = new double[a1.Length];
            if (a1.Length != a2.Length)
            {
                return null;
            }
            for (int i = 0; i < a1.Length; i++)
            {
                if (a1[i] < a2[i])
                {
                    res[i] = a1[i];
                }
                else
                {
                    res[i] = a2[i];
                }
            }

            return res;
        }
        public static double minimum(double v1, double v2)
        {
            if (v1 < v2)
            {
                return v1;
            }
            else
            {
                return v2;
            }
        }
        public static double[] minimum(double[] ar, double v)
        {
            double[] res = new double[ar.Length];
            for (int i = 0; i < ar.Length; i++)
            {
                if (ar[i] < v)
                {
                    res[i] = ar[i];
                }
                else
                {
                    res[i] = v;
                }
            }
            return res;
        }
        public static double[] minimum(double v, double[] ar)
        {
            double[] res = new double[ar.Length];
            for (int i = 0; i < ar.Length; i++)
            {
                if (ar[i] < v)
                {
                    res[i] = ar[i];
                }
                else
                {
                    res[i] = v;
                }
            }
            return res;
        }

        /// <summary>
        ///     Checa se há um canal vazio em algum lugar.   
        /// </summary>
        /// <returns></returns>
        public static bool mel_f_check(double[] v1, double[] v2, int size)
        {
            bool[] vet = new bool[size];

            for (int i = 0; i < size; i++)
            {
                if (v1[i] == 0 || v2[i] > 0)
                {
                    vet[i] = true;
                }
                else
                {
                    vet[i] = false;
                }
            }
            for (int i = 0; i < size; i++)
            {
                if (vet[i] == false)
                {
                    return false;
                }
            }
            return true;

        }

        /// <summary>
        ///     retorna um vetor com os maiores valores do eixo y da matriz
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static double[] getMaxYOfMatrix(double[][] mat, int x, int y) {
            double[] vet = new double[y];

            for (int i = 0; i < x; i++)
            {
                vet[i] = -99999;
                for (int j = 0; j < y; j++)
                {
                    if (vet[i] < mat[i][j])
                    {
                        vet[i] = mat[i][j];
                    }
                }
            }
            return vet;

        }

        /// <summary>
        /// Escreve o array na tela(apenas para efeito de testes de funcionamento das funções)
        /// </summary>
        /// <param name="array"></param>
        /// <param name="size"></param>
        /// <param name="nome"></param>
        public static void write_array(double[] array, int size, string nome)
        {
            Console.WriteLine(nome + ": ");
            for (int i = 0; i < size; i++)
            {
                Console.Write(array[i] + " ");
            }
            Console.WriteLine();
            Console.WriteLine("tamanho: " + size);
        }

        public static void write_mat(double[,] mat, int x, int y, string nome){
            Console.WriteLine(nome + ": ");
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    Console.Write(mat[i,j] + " ");
                }
            }
            Console.WriteLine();
            Console.WriteLine("tamanho: " + x+" "+y);
        }

        /// <summary>
        ///    2/2 função melspectrogram
        /// </summary>
        /// <param name="sr"></param>
        /// <param name="n_fft"></param>
        /// <param name="n_mels"></param>
        /// <param name="fmin"></param>
        /// <param name="htk"></param>
        /// <param name="norm"></param>
        /// <returns></returns>
        //Parametros padrão sao os definidos como "melhores"
        public static double[][] mel(int sr = 16000, int n_fft = 2048, int n_mels = 128, double fmin = 0, bool htk = false, int norm = 1)
        {
            double fmax = (double)sr / 2;
            double[][] weights = new double[n_mels][];

            double[] fftfreqs = fft_frequencies().ToArray<double>();
            //write_array(fftfreqs, fftfreqs.Length, "fftfreqs");

            double[] mel_f = mel_frequencies(n_mels:n_mels+2,fmin:fmin,fmax:fmax);
            //write_array(mel_f, mel_f.Length, "mel_f");

            double[] fdiff = diff(mel_f);
            //write_array(fdiff, fdiff.Length, "fdiff");

            double[,] ramps = subtract_outer(mel_f, fftfreqs);
            //write_mat(ramps, ramps.GetLength(0), ramps.GetLength(1), "ramps");

            double[] lower = new double[ramps.GetLength(1)];
            double[] upper = new double[ramps.GetLength(1)];
            
            double zero = 0.0d;
            int tamanho = (int)(1 + (n_fft / 2));
            for (int i = 0; i < n_mels; i++)
            {
                weights[i] = new double[tamanho];
                for (int j = 0; j < ramps.GetLength(1); j++)
                {
                    lower[j] = -ramps[i,j]/fdiff[i];
                    upper[j] = ramps[i + 2, j] / fdiff[i];
                }
                weights[i] = maximum(zero, minimum(lower,upper));
                
            }


            if (norm == 1)
            {
                double[] enorm = new double[n_mels];
                for (int i = 2; i < n_mels+2; i++)
                {
                    enorm[i-2] = 2.0 / (mel_f[i] - mel_f[i - 2]);
                }
                for (int i = 0; i < n_mels; i++)
                {
                    for (int j = 0; j < tamanho; j++)
                    {
                               
                        weights[i][j] *= enorm[i];
                        
                    }
                    Console.WriteLine();
                }
            }
            //Console.WriteLine(weights[127][1021]);

            if ((mel_f_check(mel_f, getMaxYOfMatrix(weights, n_mels, tamanho), n_mels))!=true){
                Console.WriteLine("filtros vazios detectados na base mel-freq.");
                Console.WriteLine("Alguns canais podem produzir respostas vazias");
                Console.WriteLine("tente incrementar a taxa de exemplos (e fmax) ou reduzir n_mels");
            }

            return weights;
        }

    }
}
