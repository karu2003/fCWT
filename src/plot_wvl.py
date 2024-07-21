import numpy as np
import matplotlib.pyplot as plt
import os

current_path = os.getcwd()
directory_path = "lib/fCWT/build/"
file = current_path + "/" + directory_path + "wavelet.wvl"

data = np.loadtxt(file)

plt.plot(data)
plt.title('Wavelet')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.show()