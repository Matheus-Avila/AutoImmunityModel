import numpy as np

file_data = np.loadtxt("result/matrix/activatedDC28.0.txt")

print("Final value adc {}".format(np.sum(file_data)))