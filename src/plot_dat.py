import numpy as np
import matplotlib.pyplot as plt
import glob
import os

current_path = os.getcwd()
print(current_path)

# Находим все файлы .dat в каталоге ./build/ и его подкаталогах
dat_files = glob.glob("**/*.dat", recursive=True)

# print(dat_files)

# for file in dat_files:
#     print(file)

# Путь к каталогу
directory_path = "lib/fCWT/build/"

# Находим все файлы .dat в указанном каталоге
dat_files = glob.glob(f"{directory_path}/*.dat")

# print(dat_files)


# Перебираем все найденные файлы
for file in dat_files:
    # Считываем данные из файла
    data = np.loadtxt(file)

    # Разделяем данные на x (индексы) и y (значения)
    x = data[:, 0]
    y = data[:, 1]

    # Создаем график для реальной части
    plt.figure(figsize=(10, 6))
    plt.plot(
        x, y, label=f"Реальная часть {file}", marker="o", linestyle="-", markersize=4
    )

    # Если в данных есть мнимая часть
    if data.shape[1] > 2:
        y_imag = data[
            :, 2
        ]  # Предполагаем, что мнимая часть находится в третьем столбце
        plt.plot(
            x,
            y_imag,
            label=f"Мнимая часть {file}",
            marker="x",
            linestyle="--",
            markersize=4,
        )

    # Добавляем заголовок и подписи к осям
    plt.title(f"График {file}")
    plt.xlabel("Индекс")
    plt.ylabel("Значение")

    # Добавляем легенду
    plt.legend()

    # Отображаем график
    plt.grid(True)

plt.show()