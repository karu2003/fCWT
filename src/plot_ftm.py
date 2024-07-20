import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap


current_path = os.getcwd()
directory_path = "lib/fCWT/build/"
file = current_path + "/" + directory_path + "tfm.cwt"

with open(file, "r") as fl:
    first_line = fl.readline().strip()
    params = first_line.split()
    n, fn, f0, f1 = map(float, params)

data = np.loadtxt(
    file,
    dtype=np.complex64,
    skiprows=1,
    converters={
        0: lambda s: complex(
            float(s.decode().split(",")[0]), float(s.decode().split(",")[1])
        )
    },
)

print(f"Signal Options: n={n}, fn={fn}, f0={f0}, f1={f1}")
data_complex = data.reshape((int(n), -1), order="F")

data_complex_transposed = np.transpose(data_complex)
print(np.argmax(np.abs(data_complex_transposed[0])))

max_indices_per_column = np.argmax(data_complex_transposed, axis=1)
# print(max_indices_per_column)
# print(len(max_indices_per_column))

# Находим индексы максимальных значений для каждой строки
max_indices = np.argmax(data_complex_transposed, axis=1)

# Создаем новый массив, заполненный нулями той же формы, что и исходный
data_max_only = np.zeros_like(data_complex_transposed)

# Заполняем новый массив максимальными значениями
for row, max_index in enumerate(max_indices):
    data_max_only[row, max_index] = data_complex_transposed[row, max_index]

print(data_max_only.shape)

colors = ["#ffffff", "#ff0000"]
cmap_name = "white_red"
cm = LinearSegmentedColormap.from_list(
    cmap_name, colors, N=100
)  # N - количество уровней

# ytick_labels = np.around(np.linspace(int(f1), int(f0), int(fn)), 2)
# plt.yticks(fontsize=8)
# ytick_positions = np.linspace(0, data_complex_transposed.shape[0] - 1, int(fn))
# plt.yticks(ytick_positions, ytick_labels)

plt.subplot(1, 2, 1)  # 1 строка, 2 колонки, первый подграфик
plt.imshow(np.abs(data_complex_transposed), aspect="auto")
plt.colorbar()
plt.title("Original CWT absolute values")
plt.ylabel("Frequency")
plt.xlabel("Time")

plt.subplot(1, 2, 2)  # 1 строка, 2 колонки, второй подграфик
plt.imshow(np.abs(data_max_only), aspect="auto", cmap=cm)
plt.colorbar()
plt.title("Max values with custom cmap")
plt.ylabel("Frequency")
plt.xlabel("Time")

plt.tight_layout()  # Автоматическая корректировка подграфиков, чтобы они не перекрывались
plt.show()
