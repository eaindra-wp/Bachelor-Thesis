combined-npz:
	mpic++ decompose-svd.cpp -g -o decompose-svd -lscalapack -llapack -lblas -lcblas -lgfortran -I/home/ewunpyae/scalapack-playground/cnpy/build_2/include -L/home/ewunpyae/scalapack-playground/cnpy/build_2/lib -lcnpy -I/home/ewunpyae/zlib-1.2.11/confg/include -L/home/ewunpyae/zlib-1.2.11/confg/lib -lz -std=c++11 -lstdc++fs
	mpic++ decompress-values.cpp -g -o decompress-values -lscalapack -llapack -lblas -lcblas -lgfortran -I/home/ewunpyae/scalapack-playground/cnpy/build_2/include -L/home/ewunpyae/scalapack-playground/cnpy/build_2/lib -lcnpy -I/home/ewunpyae/zlib-1.2.11/confg/include -L/home/ewunpyae/zlib-1.2.11/confg/lib -lz -std=c++11 -lstdc++fs
	
separate-npz:
	mpic++ decompose-svd-each-npy.cpp -g -o decompose-svd-each-npy -lscalapack -llapack -lblas -lcblas -lgfortran -I/home/ewunpyae/scalapack-playground/cnpy/build_2/include -L/home/ewunpyae/scalapack-playground/cnpy/build_2/lib -lcnpy -I/home/ewunpyae/zlib-1.2.11/confg/include -L/home/ewunpyae/zlib-1.2.11/confg/lib -lz -std=c++11 -lstdc++fs
	mpic++ decompress-values-each-npy.cpp -g -o decompress-values-each-npy -lscalapack -llapack -lblas -lcblas -lgfortran -I/home/ewunpyae/scalapack-playground/cnpy/build_2/include -L/home/ewunpyae/scalapack-playground/cnpy/build_2/lib -lcnpy -I/home/ewunpyae/zlib-1.2.11/confg/include -L/home/ewunpyae/zlib-1.2.11/confg/lib -lz -std=c++11 -lstdc++fs
	
test-data:
	g++ -g -o runCmd-decompose runCmd-decompose.cpp -std=c++11 -lstdc++fs
	g++ -g -o runCmd-decompress runCmd-decompress.cpp -std=c++11 -lstdc++fs

