import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow
import seaborn as sns
from scipy.interpolate import make_interp_spline
from scipy.interpolate import interp1d


###### 
# a user can self-determine the path of fluid simulation scene by changing the parameter values. 
# for instance, if the original fluid data values are saved in '.../fire', then, the scene name
# should be '.../fire'.
# the zip file name should be in the format: ".../scenenameframe_%04d" % t +  ".npz". 
# sigma_property: header of the compressed sigma array (denArr_sigma, presArr_sigma, velXArr_sigma, or velYArr_sigma)
# header: fluid property type
# 
# you may want to comment out each three-line block, which represents the 
# directory of each fluid simulation data to calculate the MSE. 
#  
#####

# grid_size = "4096"
# t = 7
# basePath = '/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/fire_'

# grid_size = "5120"
# t = 1
# basePath = '/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/fire_'


# grid_size = "2500"
# t = 1
# basePath = '/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/guiding_3d01_'


# grid_size = "64000_960"
# t = 2
# basePath = '/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/guiding_3d01_'

basePath = './fire_combined'
grid_size = ""
t = 22 



label = 'velY'
sigma_property =  label + 'Arr_sigma'

svd_10 = basePath + str(grid_size) + "_10/"
svd_20 = basePath + str(grid_size) + "_20/"
svd_30 = basePath + str(grid_size) + "_30/"
svd_50 = basePath + str(grid_size) + "_50/"


sigma_10 = np.load(svd_10 + "frame_%04d" % t + ".npz")[sigma_property].shape[0]
sigma_20 = np.load(svd_20 + "frame_%04d" % t + ".npz")[sigma_property].shape[0]
sigma_30 = np.load(svd_30 + "frame_%04d" % t + ".npz")[sigma_property].shape[0]
sigma_50 = np.load(svd_50 + "frame_%04d" % t + ".npz")[sigma_property].shape[0]

orgPath = basePath + str(grid_size) + "/"
compresssion_10 = basePath + str(grid_size) + "_10_decompressed/"
compresssion_20 = basePath + str(grid_size) + "_20_decompressed/"
compresssion_30 = basePath + str(grid_size) + "_30_decompressed/"
compresssion_50 = basePath + str(grid_size) + "_50_decompressed/"

print(orgPath)

orgFrame = np.load(orgPath + "frame_%04d" % t + ".npz")
decom_10 = np.load(compresssion_10 + "frame_%04d" % t + ".npz")
decom_20 = np.load(compresssion_20 + "frame_%04d" % t + ".npz")
decom_30 = np.load(compresssion_30 + "frame_%04d" % t + ".npz")
decom_50 = np.load(compresssion_50 + "frame_%04d" % t + ".npz")

density_npy_org = orgFrame["denArr"]
pressure_npy_org = orgFrame["presArr"]
velX_npy_org = orgFrame["velXArr"]
velY_npy_org = orgFrame["velYArr"]

print("Only 10 percent .....")
den_decom_10 = decom_10["denArr"]
pres_decom_10 = decom_10["presArr"]
velX_decom_10 = decom_10["velXArr"]
velY_decom_10 = decom_10["velYArr"]
den_diffs_10 = np.mean(np.power(density_npy_org - den_decom_10, 2))
pres_diffs_10 = np.mean(np.power(pressure_npy_org - pres_decom_10, 2))
velX_diffs_10 = np.mean(np.power(velX_npy_org - velX_decom_10, 2))
velY_diffs_10 = np.mean(np.power(velY_npy_org - velY_decom_10, 2))
print(den_diffs_10, pres_diffs_10, velX_diffs_10, velY_diffs_10)


print("Only 20 percent .....")
den_decom_20 = decom_20["denArr"]
pres_decom_20 = decom_20["presArr"]
velX_decom_20 = decom_20["velXArr"]
velY_decom_20 = decom_20["velYArr"]
den_diffs_20 = np.mean(np.power(density_npy_org - den_decom_20, 2))
pres_diffs_20 = np.mean(np.power(pressure_npy_org - pres_decom_20, 2))
velX_diffs_20 = np.mean(np.power(velX_npy_org - velX_decom_20, 2))
velY_diffs_20 = np.mean(np.power(velY_npy_org - velY_decom_20, 2))
print(den_diffs_20, pres_diffs_20, velX_diffs_20, velY_diffs_20)


print("Only 30 percent .....")
den_decom_30 = decom_30["denArr"]
pres_decom_30 = decom_30["presArr"]
velX_decom_30 = decom_30["velXArr"]
velY_decom_30 = decom_30["velYArr"]
den_diffs_30 = np.mean(np.power(density_npy_org - den_decom_30, 2))
pres_diffs_30 = np.mean(np.power(pressure_npy_org - pres_decom_30, 2))
velX_diffs_30 = np.mean(np.power(velX_npy_org - velX_decom_30, 2))
velY_diffs_30 = np.mean(np.power(velY_npy_org - velY_decom_30, 2))
print(den_diffs_30, pres_diffs_30, velX_diffs_30, velY_diffs_30)



print("Only 50 percent .....")
den_decom_50 = decom_50["denArr"]
pres_decom_50 = decom_50["presArr"]
velX_decom_50 = decom_50["velXArr"]
velY_decom_50 = decom_50["velYArr"]
den_diffs_50 = np.mean(np.power(density_npy_org - den_decom_50, 2))
pres_diffs_50 = np.mean(np.power(pressure_npy_org - pres_decom_50, 2))
velX_diffs_50 = np.mean(np.power(velX_npy_org - velX_decom_50, 2))
velY_diffs_50 = np.mean(np.power(velY_npy_org - velY_decom_50, 2))
print(den_diffs_50, pres_diffs_50, velX_diffs_50, velY_diffs_50)


den_diffs = [den_diffs_10, den_diffs_20, den_diffs_30, den_diffs_50]
pres_diffs = [pres_diffs_10, pres_diffs_20, pres_diffs_30, pres_diffs_50]
velX_diffs = [velX_diffs_10, velX_diffs_20, velX_diffs_30, velX_diffs_50]
velY_diffs = [velY_diffs_10, velY_diffs_20, velY_diffs_30, velY_diffs_50]



if label == 'velX':
    y_axis = np.array(velX_diffs)
elif label == 'velY':
    y_axis = np.array(velY_diffs)
elif label == 'density':
    y_axis = np.array(den_diffs)
elif label == 'pressure':
    y_axis == np.array(pres_diffs)


x_axis = np.array([sigma_10, sigma_20, sigma_30, sigma_50])

plt.semilogy()
plt.plot(x_axis, y_axis, label = label)



plt.xlabel('Number of Principal Components')
plt.ylabel('Mean Squared Error')
plt.title('Lossy Compression Performance on a matrix of size 64000 x 960')
plt.legend()
plt.savefig('/home/ewunpyae/scalapack-playground/bin_tests/figures/' + grid_size + "_" + label + "logscale")
