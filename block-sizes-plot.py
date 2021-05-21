import numpy as np
import matplotlib.pyplot as plt


time_takens_4096 = np.asarray([436.903, 630.897, 563.877, 726.771, 632.482])
time_takens_5120 = np.asarray([865.691,932.098,970.149,1041.17,1089.97])
block_size_2d = np.asarray([64,128,256,512,1024])
plt.plot(block_size_2d, time_takens_4096, label = "4096 x 4096")
plt.plot(block_size_2d, time_takens_5120, label = "5120 x 5120")
plt.title('Elapsed Times For SVD with 4 parallel processes')


# time_takens_2500 = np.asarray([301.71, 335.358, 323.131, 338.055, 435.878])
# time_takens_65536 = np.asarray([250.36,274.403, 375.982, 576.191, 700.572])
# block_size_3d = np.asarray([32,64,128,512,1024])


# plt.plot(block_size_3d, time_takens_2500, label = "50000 x 2500")
# plt.plot(block_size_3d, time_takens_65536, label = "65536 x 2048")
# plt.title('Elapsed Times For SVD with 2 parallel processes')

plt.xlabel('Block-size (in square)')
plt.ylabel('Elapsed Times (seconds)')
plt.semilogy()

plt.legend()
# plt.show()
# plt.savefig('3d-block-sizes-times')
plt.savefig('2d-block-sizes-times')