import numpy as np

file_data = np.loadtxt("result/matrix/antibody28.0.txt")

print("Final value cd8 {}".format(np.sum(file_data)))