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
            for(int i=0; i<10; i++)
                Console.Write(window[i].ToString()+"| ");
            Console.Write("\nResizing...\n");

            int size = 12;
            double[] new_window;
            new_window = Util.Pad_center(window, size);
            var y_frames = Util.Frame(window, 3, 2);
            
            for (int i = 0; i < y_frames.GetLength(0); i++){
                for (int j = 0; j < y_frames.GetLength(1); j++) {
                    Console.Write(y_frames[i,j].ToString()+" | ");
                }
                Console.WriteLine();
            }
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
