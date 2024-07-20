import numpy as np
import matplotlib.pyplot as plt
import glob
import os

current_path = os.getcwd()
print(current_path)

# Находим все файлы .dat в каталоге ./build/ и его подкаталогах
dat_files = glob.glob("**/*.dat", recursive=True)
directory_path = "lib/fCWT/build/"

# Находим все файлы .dat в указанном каталоге
dat_files = glob.glob(f"{directory_path}/*.fft")

# Перебираем все найденные файлы
for file in dat_files:
    data = np.loadtxt(file)

    # Разделяем данные на x (номера отчетов) и y (значения)
    x = data[:, 0]  # Номера отчетов
    complex_data = data[:, 1] + 1j * data[:, 2]
    power_spectrum = np.abs(complex_data) ** 2
    # Визуализация мощности спектра
    plt.figure(figsize=(10, 6))
    plt.plot(power_spectrum, label='Мощность спектра')
    plt.title(f"График {file}")
    plt.xlabel('Номер отчета')
    plt.ylabel('Мощность')
    plt.legend()
plt.show()
