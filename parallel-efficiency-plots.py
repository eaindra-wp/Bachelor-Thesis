import numpy as np
import matplotlib.pyplot as plt

# matr_5120 = [2.2721/2, 5.57812/4, 4.01033/6, 3.01377/8, 1.5399/16]
# matr_4096 = [1.49047/2, 4.03839/4, 2.30729/6, 1.48921/8, 0.84753/16]
# matr_3072 = [3.15530/2, 2.86831/4, 1.83709/6, 0.96413/8, 0.54767/16]
# processors = [2,4,6,8,16]
# plt.plot(processors, matr_5120, label = "5120 x 5120")
# plt.plot(processors, matr_4096, label = "4096 x 4096")
# plt.plot(processors, matr_3072, label = "3072 x 3072")
# plt.xlabel('Number of Processes')
# plt.ylabel('Parallel Efficiency')
# plt.title('Parallel Efficiency for SVD Computation on 2D Fluid Simulations Data')
# plt.legend()
# plt.savefig('2d-efficiency')

processors = np.asarray([2,4,6,8,16])
matr_50000_2500 = np.asarray([1.8661,2.0798,1.0889, 0.59871, 0.35655]) / processors
matr_65536_2048 = np.asarray([1.461, 2.7511, 1.6261,1.1465,0.56362]) / processors
matr_64000_960 = [2.2691,1.5685,0.7151,0.43631, 0.23597] / processors
processors = [2,4,6,8,16]
plt.plot(processors, matr_50000_2500, label = "50000 x 2500 ")
plt.plot(processors, matr_65536_2048, label = "65536 x 2048")
plt.plot(processors, matr_64000_960, label = "64000 x 960")
plt.xlabel('Number of Processes')
plt.ylabel('Parallel Efficiency')
plt.title('Parallel Efficiency for SVD Computation on 3D Fluid Simulations Data')
plt.legend()
plt.savefig('3d-efficiency')