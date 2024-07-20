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

std::vector<float> generate_chirp_signal(int duration_points, int sample_rate, int start_freq, int stop_freq) {
    std::vector<float> signal(duration_points);
    double t;
    for (int i = 0; i < duration_points; ++i) {
        t = static_cast<double>(i) / sample_rate;
        double freq = start_freq + (stop_freq - start_freq) * (t * sample_rate / duration_points);
        signal[i] = std::sin(2 * M_PI * freq * t);
    }
    return signal;
}

std::vector<std::complex<float>> generate_complex_chirp_signal(int duration_points, int sample_rate, int start_freq, int stop_freq) {
    std::vector<std::complex<float>> signal(duration_points);
    double t;
    for (int i = 0; i < duration_points; ++i) {
        t = static_cast<double>(i) / sample_rate;
        double freq = start_freq + (stop_freq - start_freq) * (t * sample_rate / duration_points);
        signal[i] = std::complex<float>(std::cos(2 * M_PI * freq * t), std::sin(2 * M_PI * freq * t));
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
    float twopi = 2.0 * 3.1415;

    // 3000 frequencies spread logartihmically between 1 and 32 Hz
    const float f0 = 100;
    const float f1 = 80000;
    const int fn = 200;  // 200
    int chirp_n = 500;
    const float fstart = 7000;
    const float fend = 17000;

    // Define number of threads for multithreaded use
    const int nthreads = 8;

    // input: n real numbers
    std::vector<float> sig(n);
    std::fill(sig.begin(), sig.end(), 0.0f);

    std::random_device rd;                 // Получаем случайное начальное число
    std::mt19937 gen(rd());                // Инициализируем генератор случайных чисел
    std::normal_distribution<> dis(0, 1);  // Нормальное распределение с mean = 0 и std dev = 1

    std::generate(sig.begin(), sig.end(), [&]() { return dis(gen); });

    // input: n complex numbers
    std::vector<complex<float>> sigc(n);
    std::fill(std::begin(sigc), std::end(sigc), std::complex<float>(0.0, 0.0));

    std::generate(sigc.begin(), sigc.end(), [&]() {
        return std::complex<float>(dis(gen), dis(gen));
    });

    // output: n x scales x 2 (complex numbers consist of two parts)
    std::vector<complex<float>> tfm(n * fn);

    std::vector<float> chirp(chirp_n);
    std::vector<complex<float>> chirpc(chirp_n);
    chirp = generate_chirp_signal(chirp_n, fs, fstart, fend);
    chirpc = generate_complex_chirp_signal(chirp_n, fs, fstart, fend);

    // std::copy(chirp.begin(), chirp.end(), sig.begin());
    // std::copy(chirpc.begin(), chirpc.end(), sigc.begin());

    // Суммирование chirp с sig
    for (size_t i = 0; i < chirp.size(); ++i) {
        sig[i] += chirp[i];
    }

    // Суммирование chirpc с sigc
    for (size_t i = 0; i < chirpc.size(); ++i) {
        sigc[i] += chirpc[i];
    }

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
    Morlet morl(2.0f);
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
    saveTFMToFile(tfm, "tfm.cwt", n, fn, f0, f1);

    cout << "=== fCWT example ===" << endl;
    cout << "Calculate CWT of a 100k sample sinusodial signal using a [" << f0 << "-" << f1 << "] Hz linear frequency range and " << fn << " wavelets." << endl;
    cout << "====================" << endl;
    cout << "fCWT finished in " << elapsed.count() << "s" << endl;

    return 0;
}
