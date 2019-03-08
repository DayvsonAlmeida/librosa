using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace conversão_librosa{
    class Program{
        static void Main(string[] args){
            double []window;
            window = Spectrum.get_window(10);
            for(int i=0; i<10; i++)
                Console.Write(window[i].ToString()+"| ");
            Console.Write("\nResizing...\n");

            int size = 12;
            
            double[] new_window;
            new_window = util.pad_center(window, size);
            for (int i = 0; i < size; i++)
                Console.Write(new_window[i].ToString() + "| ");

            Console.Read();
        }
    }

    class Spectrum{
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
    class util{
        public static double[] Pad_center(double[] data, int size){
            double[] new_window;
            new_window = new double[size];
            int start = (new_window.Length - data.Length) / 2;
            Array.Copy(data, 0, new_window, start, data.Length);
            return new_window;
        }

        public static double[][] Frame(){
            return null;
        }
    }
}
