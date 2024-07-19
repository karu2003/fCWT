import numpy as np
import matplotlib.pyplot as plt
import os

current_path = os.getcwd()
directory_path = "lib/fCWT/build/"
file = current_path + "/" + directory_path + "tfm.dat"

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

plt.imshow(np.abs(data_complex_transposed), aspect="auto")
plt.colorbar()

ytick_labels = np.around(np.linspace(int(f1), int(f0), int(fn)), 2)
plt.yticks(fontsize=8)
ytick_positions = np.linspace(0, data_complex_transposed.shape[0] - 1, int(fn))
plt.yticks(ytick_positions, ytick_labels)


plt.title("Visualization of CWT absolute values")
plt.ylabel("Frequency")
plt.xlabel("Time")
plt.show()
