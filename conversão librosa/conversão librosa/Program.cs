using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace conversão_librosa{
    class Program{
        static void Main(string[] args){
            double []window;
            window = Spectrum.Get_window(10);
            for(int i=0; i<window.Length; i++)
                Console.Write(window[i].ToString()+"| ");
            Console.Write("\nResizing...\n");

            var saida = Spectrum.DFT(window);
            for (int i = 0; i < saida.Item1.Length; i++) {
                Console.WriteLine(saida.Item1.GetValue(i).ToString()+"+"+saida.Item2.GetValue(i).ToString()+"j");
            }
            /*double[] new_window;
            new_window = Util.Pad_reflect(window, 3);
            for (int i = 0; i < new_window.Length; i++) {
                Console.Write(new_window[i].ToString()+" | ");
            }*/
            Console.Read();
        }
    }

    class Spectrum{
        /// <summary>
        /// Gera uma janela de Hann com o comprimento especificado
        /// </summary>
        /// <param name="len">Comprimento da Janela de Hann</param>
        /// <returns></returns>
    	public static double[] Get_window(int len){
            double[] window;
            window = new double [len];
            for (int i = 0; i < len; i++)
                window[i] = point_hann(i, len);
            return window;
    	}
    	private static double point_hann(int n, int N){
            return 0.5 * (1 - Math.Cos((2 * Math.PI * n) / (N - 1)));
    	}

        /// <summary>
        /// Uma Simples implementação da Discrete Fourier Transform
        /// </summary>
        /// <param name="input"></param>
        /// <param name="partials"></param>
        /// <returns></returns>
        public static Tuple<double[], double[]> DFT(double[] input, int partials = 0)
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
    }
    class Util{
        /// <summary>
        /// Centraliza o array de origem em um array de tamanho size preenchido com zeros 
        /// </summary>
        /// <param name="data">Source</param>
        /// <param name="size">Tamanho do array de saída</param>
        /// <returns>Um array de tamanho size com o source centralizado com pad preenchido com zeros</returns>
        public static double[] Pad_center(double[] data, int size){
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
        public static double[] Pad_reflect(double[] data, int size){
            double[] new_window;
            
            //Copiando centralizada do array
            int total_len = data.Length+2*size;
            new_window = new double[total_len];
            Array.Copy(data, 0, new_window, size, data.Length);

            //Refletindo o array
            int left = size-1;
            int right = size+data.Length;
            
            int i = 1, j = data.Length - 2; //i: primeiro elemento a ser refletido para à esquerda; j: primeiro elemento a ser refletido à direita
            int qtd = size;
            int i_direction = 1;
            int j_direction = -1;
            while (qtd>0){
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
        public static double[,] Frame(double[] y, int frame_length=2048, int hop_length= 512){
            if (y.Length < frame_length)
                return null;

            //Quantidade de Janelas que serão utilizadas [Valor Truncado]
            var n_frames = 1 + ((y.Length - frame_length) / hop_length);

            double[,] y_frames;
            y_frames = new double[frame_length,n_frames];
            //for (int i=0; i<frame_length; i++)
            //    y_frames[i] = new double [n_frames];

            //Extraindo n_frames Janelas do Time Series y
            for (int i = 0; i < frame_length; i++){
                for (int j = 0; j < n_frames; j++){
                    y_frames[i,j] = y[j * hop_length + i]; //Mapeamento vetor x janela
                }
            }
            return y_frames;
        }
        

    }

}
