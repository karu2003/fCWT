//
//  main.cpp
//  fCWT
//
//  Created by Lukas Arts on 21/12/2020.
//  Copyright © 2021 Lukas Arts.
/*Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "main.h"
#include <fstream>
#include <algorithm>
#include <random>

using namespace std;

void saveRealDataToFile(const std::vector<float>& data, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Ошибка при открытии файла для записи: " << filename << std::endl;
        return;
    }
    int index = 0;
    for (const auto& val : data) {
        // add index to the file
        outfile << index++ << " " << val << std::endl;
    }
    outfile.close();
}

void saveComplexSignalToFile(const std::vector<std::complex<float>>& sig, const std::string& filename) {
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        for (size_t index = 0; index < sig.size(); ++index) {
            outFile << index << " " << sig[index].real() << " " << sig[index].imag() << std::endl;
        }
        outFile.close();
        std::cout << "Data saved to " << filename << " with index, real and imaginary parts." << std::endl;
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
}

void saveTFMToFile(const std::vector<std::complex<float>>& tfm, const std::string& filename, int n, float fn, float f0, float f1) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Unable to open file for writing." << std::endl;
        return;
    }

    outFile << n << " " << fn << " " << f0 << " " << f1 << std::endl;

    for (const auto& value : tfm) {
        outFile << value.real() << "," << value.imag() << std::endl;
    }
}

void generate_chirp_signal(int duration_points, int sample_rate, float start_freq, float stop_freq, float* output) {
    float t;
    float k = (stop_freq - start_freq) / duration_points;  // Коэффициент изменения частоты
    for (int i = 0; i < duration_points; i++) {
        t = (float)i / sample_rate;
        output[i] = sin(2 * M_PI * (start_freq * t + 0.5 * k * t * t));
    }
}

void generate_complex_chirp_signal(int duration_points, int sample_rate, float start_freq, float stop_freq, std::vector<std::complex<float>>* output) {
    output->resize(duration_points);
    float t1 = static_cast<float>(duration_points) / sample_rate;
    float k = (stop_freq - start_freq) / t1;

    for (int i = 0; i < duration_points; ++i) {
        float t = static_cast<float>(i) / sample_rate;
        float phase = 2.0f * M_PI * (start_freq * t + 0.5f * k * t * t);
        float real_part = std::cos(phase);
        float imag_part = std::sin(phase);
        (*output)[i] = std::complex<float>(real_part, imag_part);
    }
}

void add_noise_and_clamp(std::vector<float>& signal, float noise_amplitude) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1, 1);  // Генерация шума с амплитудой, которая поможет избежать выхода за пределы [-1, 1]

    // Наложение шума и ограничение результата диапазоном [-1, 1]
    std::transform(signal.begin(), signal.end(), signal.begin(), [noise_amplitude, &dis, &gen](float val) {
        float noisy_val = val + (dis(gen) * noise_amplitude);  // Добавление шума
        return std::max(-1.0f, std::min(1.0f, noisy_val));     // Ограничение значения диапазоном [-1, 1]
    });
}

void add_noise_and_clamp_c(std::vector<std::complex<float>>& vec, float noise_amplitude) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1, 1);  // Генерация шума в диапазоне [-0.5, 0.5]

    for (auto& val : vec) {
        // Добавление шума к действительной и мнимой частям
        float realPart = val.real() + dis(gen) * noise_amplitude;
        float imagPart = val.imag() + dis(gen) * noise_amplitude;

        // Ограничение действительной и мнимой частей диапазоном [-1, 1]
        realPart = std::max(-1.0f, std::min(1.0f, realPart));
        imagPart = std::max(-1.0f, std::min(1.0f, imagPart));

        // Обновление значения вектора
        val = std::complex<float>(realPart, imagPart);
    }
}

std::vector<float> generateSineSignal(int durationPoints, int sampleRate, float frequency, float amplitude = 1.0f, float phase = 0.0f) {
    std::vector<float> signal(durationPoints);
    float angularFrequency = 2.0f * M_PI * frequency;

    for (int i = 0; i < durationPoints; ++i) {
        float t = static_cast<float>(i) / sampleRate;
        signal[i] = amplitude * std::sin(angularFrequency * t + phase);
    }

    return signal;
}

std::vector<std::complex<float>> generateComplexSineSignal(int durationPoints, int sampleRate, float frequency, float amplitude = 1.0f, float phase = 0.0f) {
    std::vector<std::complex<float>> signal(durationPoints);
    float angularFrequency = 2.0f * M_PI * frequency;

    for (int i = 0; i < durationPoints; ++i) {
        float t = static_cast<float>(i) / sampleRate;
        float realPart = amplitude * std::cos(angularFrequency * t + phase);
        float imagPart = amplitude * std::sin(angularFrequency * t + phase);
        signal[i] = std::complex<float>(realPart, imagPart);
    }

    return signal;
}

int main(int argc, char* argv[]) {
    // int n = 1000; //signal length
    // const int fs = 1000; //sampling frequency
    // float twopi = 2.0*3.1415;

    // //3000 frequencies spread logartihmically between 1 and 32 Hz
    // const float f0 = 0.1;
    // const float f1 = 20;
    // const int fn = 200;

    int n = 1000;           // signal length 1000
    const int fs = 192000;  // sampling frequency 192000
    float noise = 1.0;  // noise amplitude
    const float wvl = 2.0f;  // wavelet sigma
    const float f0 = 3400;
    const float f1 = 34000;
    const int fn = 20;  // 200
    int chirp_n = 500;
    const float fstart = 7000;
    const float fend = 17000;

    float frequency = 17000.0f;  // Частота сигнала в Гц
    float amplitude = 1.0f;      // Амплитуда сигнала
    float phase = 0.0f;          // Начальная фаза сигнала в радианах

    // Define number of threads for multithreaded use
    const int nthreads = 1; // 8

    std::random_device rd;                 // Получаем случайное начальное число
    std::mt19937 gen(rd());                // Инициализируем генератор случайных чисел
    std::normal_distribution<> dis(0, 1);  // Нормальное распределение с mean = -1 и std dev = 1

    // input: n real numbers and fill with zeros
    std::vector<float> sig(n);
    std::fill(sig.begin(), sig.end(), 0.0f);

    // input: n complex numbers and fill with zeros
    std::vector<complex<float>> sigc(n);
    std::fill(std::begin(sigc), std::end(sigc), std::complex<float>(0.0, 0.0));

    // output: n x scales x 2 (complex numbers consist of two parts)
    std::vector<complex<float>> tfm(n * fn);

    // Generate chirp signal
    std::vector<float> chirp(chirp_n);
    std::vector<complex<float>> chirpc(chirp_n);

    // chirp = generateSineSignal(chirp_n, fs, frequency, amplitude, phase);
    // chirpc = generateComplexSineSignal(chirp_n, fs, frequency, amplitude, phase);

    generate_chirp_signal(chirp_n, fs, fstart, fend, &chirp[0]);
    generate_complex_chirp_signal(chirp_n, fs, fstart, fend, &chirpc);

    std::copy(chirp.begin(), chirp.end(), sig.begin());
    std::copy(chirpc.begin(), chirpc.end(), sigc.begin());

    add_noise_and_clamp(sig, noise);
    add_noise_and_clamp_c(sigc, noise);

    // initialize with 1 Hz cosine wave
    // for (auto& el : sig) {
    //     el = cos(twopi * ((float)(&el - &sig[0]) / (float)fs));
    // }

    // // initialize with 1 Hz cosine wave
    // for (auto& el : sigc) {
    //     el = complex<float>(cos(twopi * ((float)(&el - &sigc[0]) / (float)fs)), 0.0f);
    // }

    // Start timing
    auto start = chrono::high_resolution_clock::now();

    // Create a wavelet object
    Wavelet* wavelet;

    // Initialize a Morlet wavelet having sigma=1.0;
    Morlet morl(wvl);
    wavelet = &morl;

    // Other wavelets are also possible
    // DOG dog(int order);
    // Paul paul(int order);

    // Create the continuous wavelet transform object
    // constructor(wavelet, nthreads, optplan)
    //
    // Arguments
    // wavelet   - pointer to wavelet object
    // nthreads  - number of threads to use
    // optplan   - use FFTW optimization plans if true
    FCWT fcwt(wavelet, nthreads, true, false);

    // Generate frequencies
    // constructor(wavelet, dist, fs, f0, f1, fn)
    //
    // Arguments
    // wavelet   - pointer to wavelet object
    // dist      - FCWT_LOGSCALES | FCWT_LINSCALES for logarithmic or linear distribution of scales across frequency range
    // fs        - sample frequency
    // f0        - beginning of frequency range
    // f1        - end of frequency range
    // fn        - number of wavelets to generate across frequency range
    Scales scs(wavelet, FCWT_LINFREQS, fs, f0, f1, fn);

    // Инициализация вектора частот
    std::vector<float> frequencies(fn);

    // Получение вектора частот
    scs.getFrequencies(frequencies.data(), fn);

    // Вывод частот
    // for (int i = 0; i < fn; i++) {
    //     std::cout << "Frequency " << i << ": " << frequencies[i] << " Hz\n";
    // }

    // Perform a CWT
    // cwt(input, length, output, scales)
    //
    // Arguments:
    // input     - floating pointer to input array
    // length    - integer signal length
    // output    - floating pointer to output array
    // scales    - pointer to scales object
    fcwt.cwt(&sigc[0], n, &tfm[0], &scs);

    // End timing
    auto finish = chrono::high_resolution_clock::now();

    // Calculate total duration
    chrono::duration<double> elapsed = finish - start;

    saveRealDataToFile(sig, "sig.dat");
    saveComplexSignalToFile(sigc, "sigc.dat");
    saveTFMToFile(tfm, "tfm.cwt", n, fn, f0, f1); //frequencies[0])

    cout << "=== fCWT example ===" << endl;
    cout << "Calculate CWT of a 100k sample sinusodial signal using a [" << f0 << "-" << f1 << "] Hz linear frequency range and " << fn << " wavelets." << endl;
    cout << "====================" << endl;
    cout << "fCWT finished in " << elapsed.count() << "s" << endl;

    return 0;
}
