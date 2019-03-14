using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Numerics;
using System.Threading.Tasks;
using Accord.Audio;
using Accord.DirectSound;
using NAudio.Wave;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Filtering;

//svm
using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.ML;
using Emgu.CV.Structure;
using Emgu.CV.Util;

namespace SpeechRecognition
{
    class Program
    {
        static void Main(string[] args)
        {
            
            String path = @"C:\Users\PauloRenato\Documents\base_balanceada_14-02-2019";
            double[][] data = AudioFeatures.getDataFromBase(path, "");

            //double[][] data;
            //AudioFeatures.LoadData(@"C:\Projeto\autoleitura\paulorenato\AudioCSharp\Caracteristicas\4125_mfcc.txt", out data, 4125, 101); ;
            //AudioFeatures.scaleData(data);
            //double[][] trainData, testData;

            //AudioFeatures.MakeTrainTestByRef(data, 42,out trainData,out testData);
            //Train(trainData, @"C:\Projeto\autoleitura\paulorenato\AudioCSharp\Modelo\audioSVMTrain.xml");
            //Test(testData, @"C:\Projeto\autoleitura\paulorenato\AudioCSharp\Modelo\audioSVMTrain.xml");

            //Console.WriteLine(RecognizeAudio(@"C:\Users\PauloRenato\Documents\base_balanceada_14-02-2019\5\5_F_14-01-2019_21_22_23.wav",  @"C:\Projeto\autoleitura\paulorenato\AudioCSharp\Modelo\audioSVMTrain.xml"));


            //Console.ReadKey();
            

        }
        
        static void Train(double[][] trainData, string pathModel)
        {
            Console.WriteLine("Montando treino e target");
            SVMAudio mySVM = new SVMAudio();
            int qtdAudio = trainData.Length;
            int featuresSize = trainData[0].Length - 1;//menos a classe
            Image<Gray, float> features = new Image<Gray, float>(featuresSize, qtdAudio);
            Image<Gray, int> target = new Image<Gray, int>(1, qtdAudio);

            for(int audioIndex = 0; audioIndex < trainData.Length; audioIndex++)
            {
                for(int j = 0; j < featuresSize; j++)
                {
                    features.Data[audioIndex, j, 0] = (float)trainData[audioIndex][j];
                }
                target.Data[audioIndex, 0, 0] = (int)trainData[audioIndex][featuresSize]; // posição da classe
            }
            Console.WriteLine("Treinando");
            mySVM.Train(features, target);
            mySVM.SaveSVMToFile(pathModel);
        }

        static void Test(double[][] testData, string pathModel)
        {
            Console.WriteLine("Montando treino e target");
            SVMAudio mySVM = new SVMAudio();
            mySVM.LoadSVMFromFile(pathModel);
            int qtdAudio = testData.Length;
            int featuresSize = testData[0].Length - 1;//menos a classe
            Image<Gray, float> features = new Image<Gray, float>(featuresSize, qtdAudio);
            Image<Gray, int> target = new Image<Gray, int>(1, qtdAudio);
            Image<Gray, float> results = new Image<Gray, float>(1, qtdAudio);

            for (int audioIndex = 0; audioIndex < testData.Length; audioIndex++)
            {
                for (int j = 0; j < featuresSize; j++)
                {
                    features.Data[audioIndex, j, 0] = (float)testData[audioIndex][j];
                }
                target.Data[audioIndex, 0, 0] = (int)testData[audioIndex][featuresSize]; // posição da classe
            }
            Console.WriteLine("Testando");
            mySVM.Predict(features, results);
            int cont = 0;
            int[] qtd = new int[11];

            for (int i = 0; i < qtd.Length; i++) qtd[i] = 0;

            for (int i = 0; i < qtdAudio; i++)
            {
                Console.WriteLine("Correct: " + target.Data[i, 0, 0] + " Predicted: " + results.Data[i, 0, 0]);
                if (target.Data[i, 0, 0] == (int)results.Data[i, 0, 0])
                {
                    
                    cont++;
                }
                qtd[target.Data[i, 0, 0]]++;
            }
            Console.WriteLine("Acc: " + (float)cont / qtdAudio);

            for (int i = 0; i < qtd.Length; i++) Console.WriteLine("Quantidade "+i+": "+qtd[i]);
            Console.ReadKey();
        }
        static float RecognizeAudio(string path, string svmModel)
        {
            SVMAudio mySVM = new SVMAudio();
            mySVM.LoadSVMFromFile(svmModel);
            double[] y = AudioFeatures.load2(path);            
            double[,] mfcc_mass = MFCC_calculating.MFCC_20_calculation(y);
            double[] features_vector = FeaturesFromMFCC(mfcc_mass);
            int features_size = 100;
            Image<Gray, float> features = new Image<Gray, float>(features_size, 1);
            Image<Gray, float> results = new Image<Gray, float>(1, 1);

            AudioFeatures.NormalizeFromFile(@"C:\Projeto\autoleitura\paulorenato\AudioCSharp\Modelo\scaleSVMTrain.txt", features_vector);

            for (int i = 0; i < features_size; i++)
            {
                features.Data[0, i, 0] = (float)features_vector[i];
            }

            mySVM.Predict(features, results);
            return results.Data[0, 0, 0];
        }

        static double[] FeaturesFromMFCC(double[,] mfcc_mass)
        {            
            int count_frames = mfcc_mass.Length / 20;

            int features_size = 100;
            double[] features = new double[features_size];

            int feature_index = 0;
            for (int i = 0; i < 20; i++)
            {
                double media = 0.0;
                double max = Double.MinValue;
                double min = Double.MaxValue;
                double[] values = new double[count_frames];
                for (int k = 0; k < count_frames; k++)
                {
                    if (mfcc_mass[k, i] < min)
                    {
                        min = mfcc_mass[k, i];
                    }
                    if (mfcc_mass[k, i] > max)
                    {
                        max = mfcc_mass[k, i];
                    }
                    media += mfcc_mass[k, i];
                    values[k] = mfcc_mass[k, i];
                }
                media = media / count_frames;
                double variance = AudioFeatures.Variance(values, media);
                double mediana = AudioFeatures.GetMedian(values);

                features[feature_index++] = (float)mediana;
                features[feature_index++] = (float)media;
                features[feature_index++] = (float)variance;
                features[feature_index++] = (float)min;
                features[feature_index++] = (float)max;
            }
            return features;
        }
    }

    public class AudioFeatures
    {
        public static void LoadData(string path, out double[][] data, int qtdAudio, int features) {
            String input = File.ReadAllText(path);

            int i = 0, j = 0;
            data = new double[qtdAudio][];
            foreach (var row in input.Split('\n'))
            {
                j = 0;
                data[i] = new double[features];
                foreach (var col in row.Trim().Split(','))
                {
                    data[i][j] = double.Parse(col.Trim(), System.Globalization.CultureInfo.InvariantCulture);
                    j++;
                }
                i++;
            }

        }
        public static void NormalizeFromFile(string pathWeights, double[] features_vector){
            StreamReader pathStream = new StreamReader(pathWeights);
            string line;
            int index = 0;
            while (( line = pathStream.ReadLine()) != null)
            {
                
                //Console.WriteLine(features_vector[index]+" novo:    "+ features_vector[index] / Double.Parse(line));
                features_vector[index] = features_vector[index] / Double.Parse(line);
                index++;
            }

        }
        public static double[] load(string path)
        {
            double[] y;
            using (WaveFileReader reader = new WaveFileReader(path))
            {
                byte[] buffer = new byte[reader.Length];
                int read = reader.Read(buffer, 0, buffer.Length);
                short[] sampleBuffer = new short[read / 2];
                y = new double[sampleBuffer.Length];
                Buffer.BlockCopy(buffer, 0, sampleBuffer, 0, read);                

                int i = 0;
                foreach (short a in sampleBuffer)
                {
                    y[i] = (double)a / 32768.0;
                    i++;
                }
                
            }
            return y;
        }

        public static double[] load2(string path)
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
                    y = new double[sampleBuffer.Length-22];
                    Buffer.BlockCopy(buffer, 0, sampleBuffer, 0, read);

                    
                    for(int i = 22; i<sampleBuffer.Length; i++)
                    {
                        y[i-22] = (double)sampleBuffer[i] / 32768.0;
                        i++;
                    }

                }
            }
            return y;
        }

        /*public static double[][] mfcc(double[] y, int sr, int n_mfcc)
        {
            double[][] mfcc;
            return mfcc;
        }*/
        public static double GetMedian(double[] sourceNumbers)
        {
            //Framework 2.0 version of this method. there is an easier way in F4        
            if (sourceNumbers == null || sourceNumbers.Length == 0)
                throw new System.Exception("Median of empty array not defined.");

            //make sure the list is sorted, but use a new array
            double[] sortedPNumbers = (double[])sourceNumbers.Clone();
            Array.Sort(sortedPNumbers);

            //get the median
            int size = sortedPNumbers.Length;
            int mid = size / 2;
            double median = (size % 2 != 0) ? (double)sortedPNumbers[mid] : ((double)sortedPNumbers[mid] + (double)sortedPNumbers[mid - 1]) / 2;
            return median;
        }
        public static double Variance(double[] values, double mean, int start = 0)
        {
            double variance = 0;
            int end = values.Length;

            for (int i = start; i < end; i++)
            {
                variance += Math.Pow((values[i] - mean), 2);
            }

            int n = end - start;
            if (start > 0) n -= 1;

            return variance / (n);
        }

        public static double[][] getDataFromBase(String path, string fileName)
        {
            DirectoryInfo dir = new DirectoryInfo(path);
            DirectoryInfo[] subDirs = dir.GetDirectories();
            double[][] data = new double[4125][];

            int audioIndex = 0;

            //Salvar Características
            string file_full_path = Path.Combine(@"C:\Projeto\autoleitura\paulorenato\AudioCSharp\Caracteristicas", fileName);
            StreamWriter fileStream = new StreamWriter(file_full_path);
            foreach (DirectoryInfo sub in subDirs)
            {
                FileInfo[] files = sub.GetFiles("*.wav");
                Console.WriteLine("Current working directory: " + sub.Name);
                foreach (FileInfo file in files)
                {
                    data[audioIndex] = new double[101];
                    int featureIndex = 0;

                    double[] y = load2(file.FullName);

                    double[,] mfcc_mass = MFCC_calculating.MFCC_20_calculation(y);
                    int count_frames = mfcc_mass.Length / 20;

                    for (int i = 0; i < 20; i++)
                    {
                        double media = 0.0;
                        double max = Double.MinValue;
                        double min = Double.MaxValue;
                        double[] values = new double[count_frames];
                        for (int k = 0; k < count_frames; k++)
                        {
                            if (mfcc_mass[k, i] < min)
                            {
                                min = mfcc_mass[k, i];
                            }
                            if (mfcc_mass[k, i] > max)
                            {
                                max = mfcc_mass[k, i];
                            }
                            media += mfcc_mass[k, i];
                            values[k] = mfcc_mass[k, i];
                        }
                        media = media / count_frames;
                        double variance = Variance(values, media);
                        double mediana = GetMedian(values);

                        data[audioIndex][featureIndex++] = mediana;
                        data[audioIndex][featureIndex++] = media;
                        data[audioIndex][featureIndex++] = variance;
                        data[audioIndex][featureIndex++] = min;
                        data[audioIndex][featureIndex++] = max;
                        
                        fileStream.Write(mediana.ToString("0.000000000", System.Globalization.CultureInfo.InvariantCulture) + ",");
                        fileStream.Write(media.ToString("0.000000000", System.Globalization.CultureInfo.InvariantCulture) + ",");
                        fileStream.Write(variance.ToString("0.00000000000", System.Globalization.CultureInfo.InvariantCulture) + ",");
                        fileStream.Write(min.ToString("0.000000000", System.Globalization.CultureInfo.InvariantCulture) + ",");
                        fileStream.Write(max.ToString("0.000000000", System.Globalization.CultureInfo.InvariantCulture) + ",");
                        
                    }

                    data[audioIndex][featureIndex] = (double)int.Parse(sub.Name);
                    
                    fileStream.WriteLine(int.Parse(sub.Name));
                    fileStream.Flush();
                    
                    audioIndex++;
                }
            }
            fileStream.Flush();
            return data;
        }

        public static void scaleData(double[][] data, string path = @"C:\Users\PauloRenato\Documents\scaleSVMTrain.txt")
        {
            int featuresSize = data[0].Length - 1;
            StreamWriter fileStream = new StreamWriter(path);
            for (int j = 0; j < featuresSize; j++)
            {
                double max = double.MinValue;
                for (int i = 0; i < data.Length; i++)
                {
                    if (Math.Abs(data[i][j]) > max)
                    {
                        max = Math.Abs(data[i][j]);
                    }
                }
                for (int i = 0; i < data.Length; i++)
                {
                    data[i][j] = data[i][j] / max;
                }
                if (j != featuresSize - 1)
                {
                    fileStream.WriteLine(max);
                }
                else
                {
                    fileStream.Write(max);
                }

            }
            fileStream.Flush();
        }

        public static void MakeTrainTestByRef(double[][] allData, int seed, out double[][] trainData, out double[][] testData, double train_perc = 0.9)
        {
            Random rnd = new Random(seed);
            int totRows = allData.Length;

            int numTrainRows = (int)(totRows * train_perc);
            int numTestRows = totRows - numTrainRows;
            trainData = new double[numTrainRows][];
            testData = new double[numTestRows][];

            double[][] copy = new double[allData.Length][];
            for (int i = 0; i < copy.Length; ++i)
                copy[i] = allData[i];
            for (int i = 0; i < copy.Length; ++i)
            {
                int r = rnd.Next(i, copy.Length);
                double[] tmp = copy[r];
                copy[r] = copy[i];
                copy[i] = tmp;
            }

            for (int i = 0; i < numTrainRows; ++i)
                trainData[i] = copy[i];

            for (int i = 0; i < numTestRows; ++i)
                testData[i] = copy[i + numTrainRows];

        } // MakeTrainTestByRef
    }

    public class MFCC_calculating
    {
        public static double[] frame;        //один фрейм
        public static double[,] frame_mass;  //массив всех фреймов по 2048 отсчетов или 128 мс        
        public static Complex[,] frame_mass_FFT;     //массив результатов FFT для всех фреймов


        static int[] filter_points = {6,18,31,46,63,82,103,127,154,184,218,
                              257,299,348,402,463,531,608,695,792,901,1023};//массив опорных точек для фильтрации спекрта фрейма
        static double[,] H = new double[20, 1024];     //массив из 20-ти фильтров для каждого MFCC

        static double[] MFCC = new double[20];     //массив MFCC для данной речевой выборки   <<<<<<<<<<<<<<<<<<<<

        /// <summary>
        /// Функция для расчета MFCC для сигнала с частотой дискретизации 16кГц
        /// </summary>
        /// <param name="wav_PCM">Массив значений амплитуд аудиосигнала</param>
        /// <returns>Массив из 20-ти MFCC</returns>
        public static double[,] MFCC_20_calculation(double[] wav_PCM)
        {
            int count_frames = (wav_PCM.Length * 2 / 2048) + 1; //количество отрезков в сигнале
            //int count_frames = (wav_PCM.Length/ 256) +1;
            RMS_gate(wav_PCM);          //применение noise gate
            Normalize(wav_PCM);         //нормализация
            frame_mass = Set_Frames(wav_PCM);       //формирование массива фреймов
            Hamming_window(frame_mass, count_frames);        //окно Хэмминга для каждого отрезка
            frame_mass_FFT = FFT_frames(frame_mass, count_frames);       //FFT для каждого фрейма


            double[,] MFCC_mass = new double[count_frames, 20];         //массив наборов MFCC для каждого фрейма

            //***********   Расчет гребенчатых фильтров спектра:    *************
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 1024; j++)
                {
                    if (j < filter_points[i]) H[i, j] = 0;
                    if ((filter_points[i] <= j) & (j <= filter_points[i + 1]))
                        H[i, j] = ((double)(j - filter_points[i]) / (filter_points[i + 1] - filter_points[i]));
                    if ((filter_points[i + 1] <= j) & (j <= filter_points[i + 2]))
                        H[i, j] = ((double)(filter_points[i + 2] - j) / (filter_points[i + 2] - filter_points[i + 1]));
                    if (j > filter_points[i + 2]) H[i, j] = 0;
                }

            for (int k = 0; k < count_frames; k++)
            {
                //**********    Применение фильтров и логарифмирование энергии спектра для каждого фрейма   ***********
                double[] S = new double[20];
                for (int i = 0; i < 20; i++)
                {
                    for (int j = 0; j < 1024; j++)
                    {
                        S[i] += Math.Pow(frame_mass_FFT[k, j].Magnitude, 2) * H[i, j];
                    }
                    if (S[i] != 0) S[i] = Math.Log(S[i], Math.E);
                }

                //**********    DCT и массив MFCC для каждого фрейма на выходе     ***********
                for (int l = 0; l < 20; l++)
                    for (int i = 0; i < 20; i++) MFCC_mass[k, l] += S[i] * Math.Cos(Math.PI * l * ((i * 0.5) / 20));
            }

            //***********   Рассчет конечных MFCC для всей речевой выборки    ***********   

            return MFCC_mass;

            /*for (int i = 0; i < 20; i++)
            {
                for (int k = 0; k < count_frames; k++) MFCC[i] += MFCC_mass[k, i];
                MFCC[i] = MFCC[i] / count_frames;
            }

            return MFCC;*/
        }


        /// <summary>
        /// Функция для подавления шума по среднекравратичному уровню
        /// </summary>
        /// <param name="wav_PCM">Массив значений амплитуд аудиосигнала</param>
        private static void RMS_gate(double[] wav_PCM)
        {
            int k = 0;
            double[] buf_rms = new double[50];
            double RMS = 0;

            for (int j = 0; j < wav_PCM.Length; j++)
            {
                if (k < 100)
                {
                    RMS += Math.Pow((wav_PCM[j]), 2);
                    k++;
                }
                else
                {
                    if (Math.Sqrt(RMS / 100) < 0.005)
                        for (int i = j - 100; i <= j; i++) wav_PCM[i] = 0;
                    k = 0; RMS = 0;
                }
            }
        }

        /// <summary>
        /// Функция нормализации сигнала
        /// </summary>
        /// <param name="wav_PCM">Массив значений амплитуд аудиосигнала</param>
        private static void Normalize(double[] wav_PCM)
        {
            double[] abs_wav_buf = new double[wav_PCM.Length];
            for (int i = 0; i < wav_PCM.Length; i++)
                if (wav_PCM[i] < 0) abs_wav_buf[i] = -wav_PCM[i];   //приводим все значения амплитуд к абсолютной величине 
                else abs_wav_buf[i] = wav_PCM[i];                    //для определения максимального пика
            double max = abs_wav_buf.Max();
            double k = 1f / max;        //получаем коэффициент нормализации            

            for (int i = 0; i < wav_PCM.Length; i++)    //записываем нормализованные значения в исходный массив амплитуд
            {
                wav_PCM[i] = wav_PCM[i] * k;
            }
        }

        /// <summary>
        /// Функция для формирования двумерного массива отрезков сигнала длиной по 128мс.
        /// При этом начало каждого следующего отрезка делит предыдущий пополам
        /// </summary>
        /// <param name="wav_PCM">Массив значений амплитуд аудиосигнала</param>
        private static double[,] Set_Frames(double[] wav_PCM)
        {
            double[,] frame_mass_1;  //массив всех фреймов по 2048 отсчетов или 128 мс
            int count_frames = 0;
            int count_samp = 0;

            frame_mass_1 = new double[(wav_PCM.Length *2/2048) + 1, 2048];
            for (int j = 0; j < wav_PCM.Length; j++)
            {
                if (j >= 1024)      //запись фреймов в массив
                {
                    count_samp++;
                    if (count_samp >= 2049)
                    {
                        count_frames += 2;
                        count_samp = 1;
                    }
                    frame_mass_1[count_frames, count_samp - 1] = wav_PCM[j - 1024];
                    frame_mass_1[count_frames + 1, count_samp - 1] = wav_PCM[j];
                }
            }
            return frame_mass_1;
        }


        /// <summary>
        /// Оконная функция Хэмминга
        /// </summary>
        /// <param name="frames">Двумерный массив отрезвов аудиосигнала</param>
        /// <param name="wav_PCM">Массив значений амплитуд аудиосигнала</param>
        private static void Hamming_window(double[,] frames, int count_frames)
        {
            double omega = 2.0 * Math.PI / (2048f);
            for (int i = 0; i < count_frames; i++)
                for (int j = 0; j < 2048; j++)
                    frames[i, j] = (0.54 - 0.46 * Math.Cos(omega * (j))) * frames[i, j];
        }


        /// <summary>
        /// Быстрое преобразование фурье для набора отрезков
        /// </summary>
        /// <param name="frames">Двумерный массив отрезвов аудиосигнала</param>
        /// <param name="wav_PCM">Массив значений амплитуд аудиосигнала</param>
        private static Complex[,] FFT_frames(double[,] frames, int count_frames)
        {
            Complex[,] frame_mass_complex =
                new Complex[count_frames, 2048]; //для хранения результатов FFT каждого фрейма в комплексном виде
            Complex[] FFT_frame = new Complex[2048];     //спектр одного фрейма
            for (int k = 0; k < count_frames; k++)
            {
                for (int i = 0; i < 2048; i++) FFT_frame[i] = frames[k, i];
                Fourier.Forward(FFT_frame, FourierOptions.Matlab);
                for (int i = 0; i < 2048; i++) frame_mass_complex[k, i] = FFT_frame[i];
            }
            return frame_mass_complex;
        }
    }

    public class SVMAudio
    {
        public SVM audioSVM;
        public SVMAudio()
        {
            audioSVM = new SVM();
            audioSVM.SetKernel(SVM.SvmKernelType.Rbf);
            audioSVM.Type = SVM.SvmType.CSvc;
            audioSVM.C = 8;
            audioSVM.Gamma = 0.01;
            audioSVM.TermCriteria = new MCvTermCriteria(1000, 0.000001);

        }
        public void Train(Image<Gray, float> samples, Image<Gray, int> responses)
        {
            TrainData data = new TrainData(samples, Emgu.CV.ML.MlEnum.DataLayoutType.RowSample, responses);
            audioSVM.TrainAuto(data);
            //audioSVM.Train(samples, Emgu.CV.ML.MlEnum.DataLayoutType.RowSample, responses);
        }
        public void SaveSVMToFile(string path)
        {
            audioSVM.Save(path);
        }
        public void LoadSVMFromFile(string path)
        {
            FileStorage fs = new FileStorage(path, FileStorage.Mode.Read);
            audioSVM.Read(fs.GetFirstTopLevelNode());
            fs.Dispose();
        }
        public Image<Gray, float> Predict(Image<Gray, float> testData, Image<Gray, float> results)
        {            
            audioSVM.Predict(testData, results);
            return results;
        }
    }

}
