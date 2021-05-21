import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow
import seaborn as sns

###### 
# a user can self-determine the path of fluid simulation scene by changing the parameter values. 
# for instance, if the original fluid data values are saved in '.../fire', then, the scene name
# should be '.../fire'.
# the zip files should be placed in the folder: ".../scenename/frame_%04d" % t, where t: timestep_number. 
# sigma_property: header of the compressed sigma array (denArr_sigma, presArr_sigma, velXArr_sigma, or velYArr_sigma)
# header: fluid property type
#
# you may want to comment out each three-line block, which represents the 
# directory of each fluid simulation data to calculate the MSE. 
#####

scenename = './fire'
t = 18
grid_size = "52 x 52"


# scenename = "./guiding_3d_01_2048_65536"
# grid_size = "65536 x 2048"
# t = 2


sigma_property = "denArr_sigma"
header = "density"

orgPath =  scenename + "/"
compresssion_10 =  scenename + "_10/"
compresssion_20 =  scenename + "_20/"
compresssion_30 = scenename + "_30/"
compresssion_50 = scenename + "_50/"

# print(orgPath)

den_org_file = np.load(orgPath +  "frame_%04d" % t + "/density.npz")
den_decom_org = den_org_file['denArr']
pres_org_file = np.load(orgPath +  "frame_%04d" % t + "/pressure.npz")
pres_decom_org = pres_org_file['presArr']
velX_org_file = np.load(orgPath +  "frame_%04d" % t + "/velX.npz")
velX_decom_org = velX_org_file['velXArr']
velY_org_file = np.load(orgPath +  "frame_%04d" % t + "/velY.npz")
velY_decom_org = velY_org_file['velYArr']




print("Only 10 percent .....")
den_10_file = np.load(compresssion_10 +  "frame_%04d" % t + "_decompressed" +"/density.npz")
den_decom_10 = den_10_file['denArr']
pres_10_file = np.load(compresssion_10 +  "frame_%04d" % t + "_decompressed" + "/pressure.npz")
pres_decom_10 = pres_10_file['presArr']
velX_10_file = np.load(compresssion_10 +  "frame_%04d" % t + "_decompressed" + "/velX.npz")
velX_decom_10 = velX_10_file['velXArr']
velY_10_file = np.load(compresssion_10 +  "frame_%04d" % t + "_decompressed" + "/velY.npz")
velY_decom_10 = velY_10_file['velYArr']

den_diffs_10 = np.mean(np.power(den_decom_org - den_decom_10, 2))
pres_diffs_10 = np.mean(np.power(pres_decom_org - pres_decom_10, 2))
velX_diffs_10 = np.mean(np.power(velX_decom_org - velX_decom_10, 2))
velY_diffs_10 = np.mean(np.power(velY_decom_org - velY_decom_10, 2))
print(den_diffs_10, pres_diffs_10, velX_diffs_10, velY_diffs_10)


print("Only 20 percent .....")
den_20_file = np.load(compresssion_20 +  "frame_%04d" % t + "_decompressed" +"/density.npz")
den_decom_20 = den_20_file['denArr']
pres_20_file = np.load(compresssion_20 +  "frame_%04d" % t + "_decompressed" + "/pressure.npz")
pres_decom_20 = pres_20_file['presArr']
velX_20_file = np.load(compresssion_20 +  "frame_%04d" % t + "_decompressed" + "/velX.npz")
velX_decom_20 = velX_20_file['velXArr']
velY_20_file = np.load(compresssion_20 +  "frame_%04d" % t + "_decompressed" + "/velY.npz")
velY_decom_20 = velY_20_file['velYArr']

den_diffs_20 = np.mean(np.power(den_decom_org - den_decom_20, 2))
pres_diffs_20 = np.mean(np.power(pres_decom_org - pres_decom_20, 2))
velX_diffs_20 = np.mean(np.power(velX_decom_org - velX_decom_20, 2))
velY_diffs_20 = np.mean(np.power(velY_decom_org - velY_decom_20, 2))
print(den_diffs_20, pres_diffs_20, velX_diffs_20, velY_diffs_20)


print("Only 30 percent .....")
den_30_file = np.load(compresssion_30 +  "frame_%04d" % t + "_decompressed" +"/density.npz")
den_decom_30 = den_30_file['denArr']
pres_30_file = np.load(compresssion_30 +  "frame_%04d" % t + "_decompressed" + "/pressure.npz")
pres_decom_30 = pres_30_file['presArr']
velX_30_file = np.load(compresssion_30 +  "frame_%04d" % t + "_decompressed" + "/velX.npz")
velX_decom_30 = velX_30_file['velXArr']
velY_30_file = np.load(compresssion_30 +  "frame_%04d" % t + "_decompressed" + "/velY.npz")
velY_decom_30 = velY_30_file['velYArr']

den_diffs_30 = np.mean(np.power(den_decom_org - den_decom_30, 2))
pres_diffs_30 = np.mean(np.power(pres_decom_org - pres_decom_30, 2))
velX_diffs_30 = np.mean(np.power(velX_decom_org - velX_decom_30, 2))
velY_diffs_30 = np.mean(np.power(velY_decom_org - velY_decom_30, 2))
print(den_diffs_30, pres_diffs_30, velX_diffs_30, velY_diffs_30)


print("Only 50 percent .....")
den_50_file = np.load(compresssion_50 +  "frame_%04d" % t + "_decompressed" +"/density.npz")
den_decom_50 = den_50_file['denArr']
pres_50_file = np.load(compresssion_50 +  "frame_%04d" % t + "_decompressed" + "/pressure.npz")
pres_decom_50 = pres_50_file['presArr']
velX_50_file = np.load(compresssion_50 +  "frame_%04d" % t + "_decompressed" + "/velX.npz")
velX_decom_50 = velX_50_file['velXArr']
velY_50_file = np.load(compresssion_50 +  "frame_%04d" % t + "_decompressed" + "/velY.npz")
velY_decom_50 = velY_50_file['velYArr']

den_diffs_50 = np.mean(np.power(den_decom_org - den_decom_50, 2))
pres_diffs_50 = np.mean(np.power(pres_decom_org - pres_decom_50, 2))
velX_diffs_50 = np.mean(np.power(velX_decom_org - velX_decom_50, 2))
velY_diffs_50 = np.mean(np.power(velY_decom_org - velY_decom_50, 2))
print(den_diffs_50, pres_diffs_50, velX_diffs_50, velY_diffs_50)




# sigma_10 = np.load(compresssion_10 +  "frame_%04d/" % t + header +".npz")[sigma_property].shape[0]
# sigma_20 = np.load(compresssion_20 +  "frame_%04d/" % t + header +".npz")[sigma_property].shape[0]
# sigma_30 = np.load(compresssion_30 +  "frame_%04d/" % t + header +".npz")[sigma_property].shape[0]
# sigma_50 = np.load(compresssion_50 +  "frame_%04d/" % t + header +".npz")[sigma_property].shape[0]

# x_axis = [sigma_10, sigma_20, sigma_30, sigma_50]

# den_diffs = [den_diffs_10, den_diffs_20, den_diffs_30, den_diffs_50]
# pres_diffs = [pres_diffs_10, pres_diffs_20, pres_diffs_30, pres_diffs_50]
# velX_diffs = [velX_diffs_10, velX_diffs_20, velX_diffs_30, velX_diffs_50]
# velY_diffs = [velY_diffs_10, velY_diffs_20, velY_diffs_30, velY_diffs_50]

# if header == 'velX':
#     y_axis = np.array(velX_diffs)
# elif header == 'velY':
#     y_axis = np.array(velY_diffs)
# elif header == 'density':
#     y_axis = np.array(den_diffs)
# elif header == 'pressure':
#     y_axis == np.array(pres_diffs)


# plt.plot(x_axis, y_axis, label = header)

# plt.xlabel('Number of Principal Components')
# plt.ylabel('Mean Squared Error')
# plt.semilogy()
# plt.title('Lossy Compression Performance on a matrix of size ' + grid_size)
# plt.legend()
# # plt.savefig('/home/ewunpyae/scalapack-playground/bin_tests/figures/guiding_3d01_65536_2048_new' + sigma_property)
# plt.show()