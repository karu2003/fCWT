import numpy as np
import matplotlib.pyplot as plt
import os

current_path = os.getcwd()
directory_path = "lib/fCWT/build/"
file = current_path + "/" + directory_path + "tfm.dat"

# Чтение параметров из первой строки
with open(file, 'r') as fl:
    first_line = fl.readline().strip()
    params = first_line.split()
    n, fn, f0, f1 = map(float, params)

# Загрузка данных, предполагая, что каждая пара чисел - это реальная и мнимая части комплексного числа
data = np.loadtxt(file, dtype=np.complex64, skiprows=1, converters={0: lambda s: complex(float(s.decode().split(',')[0]), float(s.decode().split(',')[1]))})

# Вывод параметров и первых нескольких комплексных чисел для проверки
print(f"Параметры: n={n}, fn={fn}, f0={f0}, f1={f1}")

# Переформатирование данных в форму 1000x10
data_complex = data.reshape((int(n), -1), order='F')  # Используем порядок Fortran для правильной перестановки данных

# Проверка формы данных
assert data_complex.shape == (int(n), int(fn)), "Форма массива не соответствует ожидаемой 1000x10"

data_complex_transposed = np.transpose(data_complex)


# Визуализация абсолютных значений комплексных чисел после транспонирования
plt.imshow(np.abs(data_complex_transposed), aspect='auto')
plt.colorbar()  # Добавление шкалы цвета для интерпретации значений

# Установка меток по оси Y от 0.1 до 20 с общим количеством 10 значений
ytick_labels = np.around(np.linspace(int(f1), int(f0), int(fn)), 2)
plt.yticks(fontsize=8)
ytick_positions = np.linspace(0, data_complex_transposed.shape[0]-1, int(fn))  # Распределение позиций меток по оси Y
plt.yticks(ytick_positions, ytick_labels)


plt.title("Visualization of CWT absolute values")
plt.ylabel('Frequency')  # Меняем местами метки
plt.xlabel('Time')  # Меняем местами метки
plt.show()