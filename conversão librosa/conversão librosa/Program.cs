using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace conversão_librosa{
    class Program{
        static void Main(string[] args){
            double []window;
            window = Spectrum.get_window(2048);
            for(int i=0; i<10; i++)
                Console.Write(window[i].ToString()+", ");
            Console.Read();
        }
    }

    class Spectrum{
    	public static double[] get_window(int len){
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
}
