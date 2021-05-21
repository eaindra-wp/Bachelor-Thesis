import numpy as np
import matplotlib.pyplot as plt

matr_5120 = [2.2721, 5.57812, 4.01033, 3.01377, 1.5399]
matr_4096 = [1.49047, 4.03839, 2.30729, 1.48921, 0.84753]
matr_3072 = [3.15530, 2.86831, 1.83709, 0.96413, 0.54767]


# processors = [2,4,6,8,16]
# plt.plot(processors, matr_5120, label = "5120 x 5120")
# plt.plot(processors, matr_4096, label = "4096 x 4096")
# plt.plot(processors, matr_3072, label = "3072 x 3072")
# plt.xlabel('Number of Processes')
# plt.ylabel('Speedup')
# plt.title('Speedups to compute SVD on symmetric 2D fluid simulations')
# plt.legend()
# plt.savefig('2d-speedups')


matr_50000_2500 = [1.8661,2.0798,1.0889, 0.59871, 0.35655]
matr_65536_2048 = [1.461, 2.7511, 1.6261,1.1465,0.56362]
matr_64000_960 = [2.2691,1.5685,0.7151,0.43631, 0.23597]
processors = [2,4,6,8,16]
plt.plot(processors, matr_50000_2500, label = "50000 x 2500 ")
plt.plot(processors, matr_65536_2048, label = "65536 x 2048")
plt.plot(processors, matr_64000_960, label = "64000 x 960")
plt.xlabel('Number of Processes')
plt.ylabel('Speedup')
plt.title('Speedups to compute SVD on 3D fluid simulations')
plt.legend()
plt.savefig('3d-speedups')