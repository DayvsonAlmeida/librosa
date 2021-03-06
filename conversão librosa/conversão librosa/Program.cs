﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using Librosa;
using FftPack;
using MathNet.Numerics.IntegralTransforms;
using System.Numerics;

namespace LibrosaConv
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("INICIANDO");
            double[] y = Util.Load2("C:/Users/paulo/Documents/Autoleitura/balanceado16PCM_14-02-2019/0/0_F_17-01-2019_19_23_36.wav");
            int sr = 16000;

            var mels = Spectral.Melspectrogram(y: y, sr: sr);
            for (int i = 0; i < mels.GetLength(0); i++)
            {
                for (int j = 0; j < mels.GetLength(1); j++)
                {
                    Console.WriteLine(i + " " + j + " " + mels[i, j]);
                }
            }

            //var r = Spectrum.power_to_db(Spectral.Melspectrogram(y:y, sr:22050));

            //for (int i = 0; i < r.GetLength(0); i++)
            //{
            //    for (int j = 0; j < r.GetLength(1); j++)
            //    {
            //        Console.WriteLine(i+" "+j+" "+r[i,j]);
            //    }
            //}
            Console.Read();

        }
    }
    class Spectral
    {
        public static double[,] mfcc(double[] y, int sr = 16000, double[,] S = null, int n_mfcc = 20, int dct_type = 2, string norm = "ortho")
        {
            if (S == null)
            {
                S = Spectrum.power_to_db(Melspectrogram(y:y,sr:sr));
            }
            double[,] S_DCT2D = FftPack.Transform.DCT2D(S,norm:norm);
            double[,] out_put = new double[n_mfcc, S_DCT2D.GetLength(1)];
            for (int i = 0; i < n_mfcc; i++)
                for (int j = 0; j < S_DCT2D.GetLength(1); j++)
                    out_put[i, j] = S_DCT2D[i, j];
            return out_put;
        }

        public static double[,] Melspectrogram(double[] y, int sr = 16000, double[,] S = null, int n_fft = 2048, int hop_length = 512, int win_length = 0, bool center = true, double power = 2.0, int n_mels = 128)
        {
            S = Spectrum._spectrogram(y, ref n_fft, S, hop_length, power, win_length, center);
            var mel_basis = Librosa.Time_frequency.mel(sr, n_fft, n_mels);


            double[,] output = new double[mel_basis.Length, S.GetLength(1)];

            for (int i = 0; i < mel_basis.Length; i++)
            {
                for (int j = 0; j < S.GetLength(1); j++)
                {
                    double tmp = 0.0;
                    for (int k = 0; k < S.GetLength(0); k++)
                    {
                        tmp += mel_basis[i][k] * S[k, j];
                    }
                    output[i, j] = tmp;
                }
            }

            return output;
        }
    }
    class Spectrum
    {


        public static double[,] power_to_db(double[,] S, double reff = 1.0, double amin = 0.0000000001, double top_db = 80.0)
        {
            if (amin <= 0)
            {
                Console.WriteLine("parâmetro deve ser estritamente positivo");
            }

            double[,] magnitude = S;

            double ref_value = Math.Abs(reff);

            double[,] log_spec = utils.maximum(magnitude, amin);
            for (int i = 0; i < log_spec.GetLength(0); i++)
            {
                for (int j = 0; j < log_spec.GetLength(1); j++)
                {
                    log_spec[i, j] = 10 * Math.Log10(log_spec[i, j]) - 10 * Math.Log10(utils.maximum(amin, ref_value));
                }
            }

            if (top_db < 0)
            {
                Console.WriteLine("top_db precisa ser não negativo ;)");
            }
            log_spec = utils.maximum(log_spec, utils.maxValueOfMatrix(log_spec) - top_db);

            return log_spec;
        }


        /// <summary>
        /// Gera uma janela de Hann com o comprimento especificado
        /// </summary>
        /// <param name="len">Comprimento da Janela de Hann</param>
        /// <returns></returns>
    	public static double[] Get_window(int len)
        {
            double[] window;
            window = new double[len];
            for (int i = 0; i < len; i++)
                window[i] = point_hann(i, len);
            return window;
        }
        private static double point_hann(int n, int N)
        {
            return 0.5 * (1 - Math.Cos((2 * Math.PI * n) / (N - 1)));
        }

        /// <summary>
        /// Uma Simples implementação da Discrete Fourier Transform
        /// </summary>
        /// <param name="input"></param>
        /// <param name="partials"></param>
        /// <returns></returns>
        private static Tuple<double[], double[]> DFT(double[] input, int partials = 0)
        {
            int len = input.Length;
            double[] cosDFT = new double[len / 2 + 1];
            double[] sinDFT = new double[len / 2 + 1];

            if (partials == 0)
                partials = len / 2;

            for (int n = 0; n <= partials; n++)
            {
                double cos = 0.0;
                double sin = 0.0;

                for (int i = 0; i < len; i++)
                {
                    cos += input[i] * Math.Cos(2 * Math.PI * n / len * i);
                    sin += input[i] * Math.Sin(2 * Math.PI * n / len * i);
                }

                cosDFT[n] = cos;
                sinDFT[n] = sin;
            }

            return new Tuple<double[], double[]>(cosDFT, sinDFT);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="n_fft"></param>
        /// <param name="hop_length"></param>
        /// <param name="win_length"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        public static double[,,] STFT(double[] y, int n_fft = 2048, int hop_length = 0, int win_length = 0, bool center = true)
        {
            //Console.WriteLine("iniciando stft");
            if (win_length <= 0)
                win_length = n_fft;

            if (hop_length <= 0)
                hop_length = win_length / 4;

            var fft_window = Util.Pad_center(Spectrum.Get_window(win_length), n_fft);

            //var fft_window = Spectrum.Get_window(win_length);
            //fft_window = Util.Pad_center(fft_window, n_fft);
            var fft_window_column = Util.Line_to_Column(fft_window);

            double[] y_local = new double[y.Length];
            y.CopyTo(y_local, 0);
            if (center)
                y_local = Util.Pad_reflect(y_local, n_fft / 2);

            var y_frames = Util.Frame(y_local, n_fft, hop_length);

            double[,,] stft_matrix = new double[(1 + n_fft / 2), y_frames.GetLength(1), 2];
            const int MAX_MEM_BLOCK = 262144;
            int n_columns = (MAX_MEM_BLOCK / (stft_matrix.GetLength(0) * sizeof(double)));

            for (int bl_s = 0; bl_s < stft_matrix.GetLength(1); bl_s += n_columns)
            {
                var bl_t = Math.Min(bl_s + n_columns, stft_matrix.GetLength(1));
                var tmp = Util.GetColumns(y_frames, bl_s, bl_t);


                //Convolução da Janela de Hann com os Frames
                for (int pos = 0; pos < fft_window_column.GetLength(0); pos++)
                {
                    for (int col = 0; col < tmp.GetLength(1); col++)
                    {
                        tmp[pos, col] *= fft_window_column[pos, 0];
                    }
                }


                //Frames Transformados estão saindo nas linhas
                var tmp_rfft = Spectrum.FFT_frames(tmp); ;//Spectrum.RFFT(tmp);

                for (int i = 0; i < stft_matrix.GetLength(0); i++)
                {
                    for (int j = bl_s; j < bl_t; j++)
                    {
                        //Parte Real
                        stft_matrix[i, j, 0] = tmp_rfft[j - bl_s, i].Real;
                        //Parte Imaginária
                        stft_matrix[i, j, 1] = tmp_rfft[j - bl_s, i].Magnitude;
                    }
                }
            }
            return stft_matrix;
        }

        //Os Frames estão saindo nas linhas
        /*  Frames Antes da Transformada
            [f1][f2][f3]
            [f1][f2][f3]
            [f1][f2][f3]

            Frames Depois da Transformada no Python
            [f1][f2][f3]
            [f1][f2][f3]

            Frames Depois da Transformada no C#
            [f1][f1]
            [f2][f2]
            [f3][f3]              */
        private static Tuple<double[], double[]>[] RFFT(double[,] frames)
        {
            //Função deve realizar a transformada em cada coluna dos frames passados
            var current_frame = new double[frames.GetLength(0)];
            Tuple<double[], double[]>[] output = new Tuple<double[], double[]>[frames.GetLength(1)];

            for (int j = 0; j < frames.GetLength(1); j++)
            {
                for (int i = 0; i < frames.GetLength(0); i++)
                    current_frame[i] = frames[i, j];
                output[j] = Spectrum.DFT(current_frame);
            }

            return output;
        }
        private static Complex[,] FFT_frames(double[,] frames)
        {
            int count_frames = frames.GetLength(1);
            int length_frames = frames.GetLength(0);

            Complex[,] frame_mass_complex = new Complex[count_frames, length_frames]; //для хранения результатов FFT каждого фрейма в комплексном виде
            Complex[] FFT_frame = new Complex[length_frames];     //спектр одного фрейма
            for (int k = 0; k < count_frames; k++)
            {
                for (int i = 0; i < length_frames; i++) FFT_frame[i] = frames[i, k];
                Fourier.Forward(FFT_frame, FourierOptions.Matlab);
                int len_fft_frame = FFT_frame.Length;
                for (int i = 0; i < len_fft_frame; i++) frame_mass_complex[k, i] = FFT_frame[i];
            }
            return frame_mass_complex;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="y"></param>
        /// <param name="n_fft">Default = 2048</param>
        /// <param name="S"></param>
        /// <param name="hop_length"></param>
        /// <param name="power"></param>
        /// <param name="win_length"></param>
        /// <param name="center"></param>
        /// <returns></returns>
        public static double[,] _spectrogram(double[] y, ref int n_fft, double[,] S = null, int hop_length = 512, double power = 2.0, int win_length = 0, bool center = true)
        {
            //Console.WriteLine("iniciando spectrogram");

            double[,] S_out;
            if (S != null)
            {
                n_fft = 2 * (S.GetLength(0) - 1);
                S_out = S;
            }
            else
            {
                var tmp = Spectrum.STFT(y, n_fft, hop_length, win_length, center);
                S_out = new double[tmp.GetLength(0), tmp.GetLength(1)];
                for (int i = 0; i < tmp.GetLength(0); i++)
                {
                    for (int j = 0; j < tmp.GetLength(1); j++)
                    {
                        S_out[i, j] = Math.Pow(tmp[i, j, 0], 2) + Math.Pow(tmp[i, j, 1], 2);
                        S_out[i, j] = Math.Sqrt(S_out[i, j]);

                    }

                }
            }

            return S_out;
        }
    }
    class Util
    {
        /// <summary>
        /// Centraliza o array de origem em um array de tamanho size preenchido com zeros 
        /// </summary>
        /// <param name="data">Source</param>
        /// <param name="size">Tamanho do array de saída</param>
        /// <returns>Um array de tamanho size com o source centralizado com pad preenchido com zeros</returns>
        public static double[] Pad_center(double[] data, int size)
        {
            double[] new_window;
            new_window = new double[size];
            int start = (new_window.Length - data.Length) / 2;
            Array.Copy(data, 0, new_window, start, data.Length);
            return new_window;
        }


        /// <summary>
        /// Realiza o padding do array de origem com o reflexo do mesmo.
        /// </summary>
        /// <param name="data">Source</param>
        /// <param name="size">Tamanho do array de saída</param>
        /// <returns></returns>
        public static double[] Pad_reflect(double[] data, int size)
        {
            double[] new_window;

            //Copiando centralizada do array
            int total_len = data.Length + 2 * size;
            new_window = new double[total_len];
            Array.Copy(data, 0, new_window, size, data.Length);

            //Refletindo o array
            int left = size - 1;
            int right = size + data.Length;

            int i = 1, j = data.Length - 2; //i: primeiro elemento a ser refletido para à esquerda; j: primeiro elemento a ser refletido à direita
            int qtd = size;
            int i_direction = 1;
            int j_direction = -1;
            while (qtd > 0)
            {
                //Refletindo o item atual
                new_window[left] = data[i];
                new_window[right] = data[j];

                //Alternando direção de reflexo
                if (i + i_direction >= data.Length || i + i_direction < 0)
                    i_direction *= -1;
                if (j + j_direction >= data.Length || j + j_direction < 0)
                    j_direction *= -1;

                i += i_direction;
                j += j_direction;
                left--;
                right++;
                qtd--;
            }
            return new_window;
        }

        /// <summary>
        /// Cria uma perspectiva de janelas a partir de um conjunto de dados. 
        /// </summary>
        /// <param name="y">Vetor de origem dos dados</param>
        /// <param name="frame_length">Comprimento da janela</param>
        /// <param name="hop_length">Deslocamento entre janelas</param>
        /// <returns>Um array 2D com as janelas</returns>
        public static double[,] Frame(double[] y, int frame_length = 2048, int hop_length = 512)
        {
            if (y.Length < frame_length)
                return null;

            //Quantidade de Janelas que serão utilizadas [Valor Truncado]
            var n_frames = 1 + (int)((y.Length - frame_length) / hop_length);
            double[,] y_frames;
            y_frames = new double[frame_length, n_frames];
            //for (int i=0; i<frame_length; i++)
            //    y_frames[i] = new double [n_frames];

            //Extraindo n_frames Janelas do Time Series y
            for (int i = 0; i < frame_length; i++)
            {
                for (int j = 0; j < n_frames; j++)
                {
                    y_frames[i, j] = y[j * hop_length + i]; //Mapeamento vetor x janela
                }
            }
            return y_frames;
        }

        /// <summary>
        /// Converte um vetor linha para um vetor coluna
        /// </summary>
        /// <param name="data">Source</param>
        /// <returns>Retorna um vetor coluna com os elementos de data</returns>
        public static double[,] Line_to_Column(double[] data)
        {
            double[,] new_data;
            new_data = new double[data.Length, 1];
            for (int i = 0; i < data.Length; i++)
                new_data[i, 0] = data[i];
            return new_data;
        }

        /// <summary>
        /// Faz um slice do source
        /// </summary>
        /// <param name="source"></param>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        public static double[,] GetColumns(double[,] source, int start, int end)
        {
            int len = source.GetLength(0);
            int n_columns = end - start;

            double[,] output = new double[len, n_columns];

            for (int i = 0; i < len; i++)
            {
                for (int j = start; j < end; j++)
                {
                    output[i, j - start] = source[i, j];
                }
            }

            return output;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static double[] Load2(string path)
        {
            double[] y;
            using (var streamReader = new StreamReader(path))
            {
                var buffer = default(byte[]);
                using (var memstream = new MemoryStream())
                {
                    streamReader.BaseStream.CopyTo(memstream);
                    buffer = memstream.ToArray();
                    int read = buffer.Length;
                    short[] sampleBuffer = new short[read / 2];
                    y = new double[sampleBuffer.Length - 22];
                    Buffer.BlockCopy(buffer, 0, sampleBuffer, 0, read);


                    for (int i = 22; i < sampleBuffer.Length; i++)
                    {
                        y[i - 22] = (double)sampleBuffer[i] / 32768.0;

                    }

                }
            }
            return y;
        }
    }
}
