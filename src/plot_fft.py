import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.fft import fft, fftfreq

current_path = os.getcwd()
print(current_path)

# Находим все файлы .dat в каталоге ./build/ и его подкаталогах
dat_files = glob.glob("**/*.dat", recursive=True)
directory_path = "lib/fCWT/build/"

# Находим все файлы .dat в указанном каталоге
dat_files = glob.glob(f"{directory_path}/*.fft")

sample_rate = 192000
N = 1024

haft = int(1024 / 2)

xf = fftfreq(N, 1 / sample_rate)[: N // 2]

# xtick_labels = np.around(np.linspace(0, 192000/2, haft), 2)
# xtick_positions = np.linspace(0, haft, haft)

xtick_positions = np.linspace(0, haft - 1, haft // 10)
xtick_labels = np.around(np.linspace(0, sample_rate / 2, haft // 10), 2)

# Перебираем все найденные файлы
for file in dat_files:
    data = np.loadtxt(file)

    # Разделяем данные на x (номера отчетов) и y (значения)
    x = data[:, 0]  # Номера отчетов
    complex_data = data[:, 1] + 1j * data[:, 2]
    power_spectrum = np.abs(complex_data) ** 2
    power_spectrum = power_spectrum[: len(power_spectrum) // 2]

    # Находим максимум спектра и соответствующую частоту
    max_power_index = np.argmax(power_spectrum)
    max_power = power_spectrum[max_power_index]
    max_frequency = np.linspace(0, sample_rate / 2, len(power_spectrum))[
        max_power_index
    ]

    # power_spectrum = 10 * np.log10(power_spectrum)
    # Визуализация мощности спектра
    plt.figure(figsize=(10, 6))
    plt.plot(power_spectrum, label="Power Spectrum")
    # plt.scatter(
    #     max_frequency, max_power, color="red", label="Max Power"
    # )  # Отмечаем максимум
    plt.annotate(
        f"Max: {max_power:.2f} dB\nFreq: {max_frequency:.2f} Hz",
        (max_frequency, max_power),
        textcoords="offset points",
        xytext=(10, -15),
        ha="center",
    )

    plt.title(f"Power Spectrum of the Chirp Signal {file}")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power (dB)")
    plt.xticks(xtick_positions, xtick_labels, rotation="vertical", fontsize=8)
    plt.legend()

plt.show()
