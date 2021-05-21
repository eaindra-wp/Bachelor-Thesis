import numpy as np
import matplotlib.pyplot as plt

# matr_5120_u_s_time = np.asarray([0.121065,0.401301, 0.906444, 2.54588])
# matr_5120_final_times = np.asarray([1.0835, 1.90589,2.90117, 5.23698])
# matr_5120_sigma_count = np.asarray([512,1024,1534,2558])
# plt.plot(matr_5120_sigma_count, matr_5120_u_s_time, label = "U x Sigma")
# plt.plot(matr_5120_sigma_count, matr_5120_final_times, label = "(U x Sigma) x V^T")
# figname = '5120_single_process_multiply'

# matr_4096_u_s_time = np.asarray([0.0545816, 0.218445, 0.458416, 1.23913])
# matr_4096_final_times = np.asarray([0.495432, 1.02095, 1.493,  2.46616])
# matr_4096_sigma_count = np.asarray([410,818,1228,2046])
# plt.plot(matr_4096_sigma_count, matr_4096_u_s_time, label = "U x Sigma")
# plt.plot(matr_4096_sigma_count, matr_4096_final_times, label = "(U x Sigma) x V^T")
# figname = '4096_single_process_multiply'


# matr_64000_960_u_s_time = np.asarray([0.107829,0.765173,0.916023])
# matr_64000_960_final_times = np.asarray([0.774737, 1.88684, 2.52686])
# matr_64000_960_sigma_count = np.asarray([96,288,384])
# plt.plot(matr_64000_960_sigma_count, matr_64000_960_u_s_time, label = "U x Sigma")
# plt.plot(matr_64000_960_sigma_count, matr_64000_960_final_times, label = "(U x Sigma) x V^T")
# figname = '64000_960_single_process_multiply'

matr_65536_2048_u_s_time = np.asarray([2.07257, 0.119472, 0.281162])
matr_65536_2048_final_times = np.asarray([1.49394, 1.9835, 2.33402])
matr_65536_2048_sigma_count = np.asarray([58,114,172])
plt.plot(matr_65536_2048_sigma_count, matr_65536_2048_u_s_time, label = "U x Sigma")
plt.plot(matr_65536_2048_sigma_count, matr_65536_2048_final_times, label = "(U x Sigma) x V^T")
figname = '65536_2048_single_process_multiply'



plt.ylabel('Time taken (in seconds)')
plt.semilogy()
plt.xlabel('Number of Singular Values Taken')
plt.title('Time taken to multiply matrices on a single process')
plt.legend()
plt.savefig(figname)

