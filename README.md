# Guidline for Implementations
This guideline indicates a summary of how to run source codes for the lossy compression via SVD on ScaLAPACK, and decompress the data back. 

## Table of Contents 
* [Fluid Property Values Saved together in the Same Numpy Zip Archive](*same-npy)
* [Fluid Values for Each Fluid Property Saved As Separate Zip Files](#separate-npy)
* [Compile Programs](#compile)
* [Run Compression Programs](#run-compress)
* [Run Decompression Programs](#run-decompress)
* [Testing the Lossy Compression Performance with a sample Fluid Dataset saved in separate Numpy zips](*test-separate-npy)



--------- 
The main implementation for this thesis work are implemented in the following C++ codess. 
* `decompose-svd.cpp`
* `decompress-values.cpp`
* `decompose-svd-each-npy.cpp`
* `decompress-values-each-npy.cpp`

These C++ files can also be found on Titlis in the directory below.
```
/home/ewunpyae/scalapack-playground/bin_tests
```

# <a name="same-npy"></a>
## Fluid Values Saved together in the Same Numpy Zip Archive

The **density, pressure, velocity in x,y directions** values can be saved together as 4 Numpy arrays in a Numpy zip archive. If so, 4 Numpy arrays can be compressed iteratively with `decompose-svd.cpp`. The compressed data will be saved in a new Numpy zip archive, which contains 3 Numpy arrays (U,Sigma,Vt) for each fluid property data, resulting **12** Numpy arrays in total in the compressed archive.

The compressed fluid data for all 4 properties, can be decompressed iteratively with `decompress-values.cpp`. Afterwards, the 4 decompressed matrices will be saved in a new Numpy zip. 


# <a name="separate-npy"></a>
## Fluid Values for Each Fluid Property Saved As Separate Zip Files

The values for each fluid property can also be saved in separate Numpy zip archives. `decompose-svd-each-npy.cpp` is implemented in a way that each separate fluid data can be compressed successfully. Three arrays: U,Sigma,Vt will be saved in a new Numpy Zip Archive. 

After the lossy compression, the compressed data for just one property can be decompressed by running the program implemented in `decompress-values-each-npy.cpp`.

---------------------

# <a name="compile"></a>
## Compile Programs 
`decompose-svd.cpp` and `decompress-values.cpp` can be compiled with the command below. 
```
make combined-npz
``` 

`decompose-svd-each-npy.cpp` and `decompress-values-each-npy.cpp` can be compiled with the command below. 
```
make separate-npz
``` 
--------------

Fluid simulation datasets from two small 2D fluid simulations are added in this folder, so that you may want to test the compression and decompression performance. 

* **fire/frame_0018/** contains the data for each fluid property in separate Numpy zip archives
* **fire/** contains the data for all fluid property saved in one Numpy zip archive. 

------------------------------

# <a name="run-compress"></a>
## Run Compression Programs
### `decompose-svd`
The example command to run the executable from `decompose-svd` is as follows. 

```
mpirun -np $(np) ./decompose-svd originalNumpyZip procRow procCol blockRow blockCol compressionRatio
```
An example command is as follows. 

```
mpirun -np 4 ./decompose-svd ./fire_combined/frame_0022.npz 2 2 10 10 50
```

### `decompose-svd-each-npy`
The example command to run the executable from `decompose-svd-each-npy` is as follows. 

```
mpirun -np $(np) ./decompose-svd-each-npy originalNumpyZip procRow procCol blockRow blockCol compressionRatio
```
An example command is as follows. 

```
 mpirun -np 4 ./decompose-svd-each-npy ./fire/frame_0018/pressure.npz 2 2 10 10 50
```
If the program is executed successfully, you will see some new lines printed in the console. The following is an example. 
```
compression_ratio 50 0.5
N: 50000  p: 508  M: 2500
334.369
```
The first line indicates the compression ratio, the second line indicates the row and column size of the matrix (N,M), and the number of singular values and vectors as the compressed data. The third line indicates the maximum elapsed time to compute SVD on parallel processes. 

----------------------------

# <a name="run-decompress"></a>
## Run Decompression Programs 

### `decompress-values`
The example command to run the executable from `decompress-values` is as follows. 

```
mpirun -np $(np) ./decompress-values compressedNumpyZip procRow procCol
```
An example command is as follows. 

```
mpirun -np 1 ./decompress-values ./fire_combined_50/frame_0022.npz 1 1
```

### `decompress-values-each-npy`
The example command to run the executable from `decompress-values-each-npy` is as follows. 

```
mpirun -np $(np) ./decompress-values-each-npy compressedNumpyZip procRow procCol
```
An example command is as follows. 

```
mpirun -np 1 ./decompress-values ./fire_50/frame_0018.npz 1 1
```

If the program is executed successfully, you will see the commands similar to the example below. 
```
conpressedFileName: denArr
0.000113547 1.53463e-05
```
The first line indicates the name of the Numpy array header to be saved. In the second line,  the left number indicates the time to multiple U and Sigma, and the right number is the time to multiply U x Sigma) and V^T. 


--------------------
# <a name="test-separate-npy"></a>
## Testing the Lossy Compression Performance with a sample Fluid Dataset saved in separate Numpy zips


To compress and decompress the fluid datasets in **fire/frame_0018/**, you may run 
```
make test-data
```
to compile the test files. \
You may want to compress the original test fluid data by running 
```
./runCmd-decompose
```
and decompress that compressed test fluid data by running 
```
./runCmd-decompress
```
----------
## Running Compression and Decompression on Titlis: 
You may run these data compression and decompression files on Titlis Machine. You will need to activate 
```
conda activate ewp-project-env
```
to compile and run the data compression and decompression C++ codes. All libraries required are installed in this environment. 

`cnpy` and `zlib` libraries are locally installed and built at 
`/home/ewunpyae/scalapack-playground/cnpy/build_2` and `/home/ewunpyae/zlib-1.2.11/confg/` respectively.

----------------------------------

### 2D fluid simulation data sets used for this thesis work
```
/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/fire_4096/frame_0007.npz

/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/fire_5120/frame_0001.npz
```

### 3D fluid simulation data sets used for this thesis work
```
/home/ewunpyae/scalapack-playground/bin_tests/guiding_3d_01_2048_65536/frame_0002/

/home/ewunpyae/scalapack-playground/generate-manta-data-lossy-compression/guiding_3d01_64000_960/frame_0002.npz

```

----------

## Calculating and Plotting Lossy Compression 
These two files can also be found on Titlis in the directory.
```
/home/ewunpyae/scalapack-playground/bin_tests
```
* `lossy-diff-calculator.py`
* `lossy-diffs-calc-all-props-in-zip.py`

You may check the scripts in `lossy-diff-calculator.py` if the original and decompressed fluid data for density,pressure,velocity in X,Y directions are saved in the same Numpy zip archive. 

You may check the scrips in `lossy-diffs-calc-all-props-in-zip.py` if the original and decompressed fluid data for each property are saved in separate Numpy zip archives. 
 
---------------

## Plotting Speedups and Efficiency
These scripts were used to plot graphs for the speedups and parallel-efficiency results shown in the thesis. 
* `speedup-plots.py`
* `parallel-efficiency-plots.py`
* `decompression-time-takens.py`
* `block-sizes-plot.py`
-------------

## Elapsed times to execute SVD computation
The elapsed times measured for the SVD computation for different fluid simulation data can be checked in `time_takens` folder. 

The values in each line in the text-files represent the following format, divided by a single space.



| No.of Sigmas > 0  | No. of Sigmas and Vectors Saved | ProcessNumber |ProcGridRow | ProcGridCol| Max. Time among parallel proccesses|
| ------------- | ------------- | ------------- |------------- | ------------- | ------------- |
|   |   |    |    |   |    |    |    |


The initial elapsed time measurements may not contain the results for **No.of Sigmas > 0** and **No. of Sigmas and Vectors Saved** since the C++ source codes were fixed after some SVD computations were made. 