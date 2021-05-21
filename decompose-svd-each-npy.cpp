
/*
 * This program is to compress the original fluid data via SVD on parallel processes. 
 * The input Numpy zip file needs to contain only 1 Numpy array, either the values of the density, pressure, 
 * or velocity in x,y directions. 
 * 
 *
 * C++ array <---> Numpy zip archive convertion is done by using the 'cnpy' library,
 * whose repository is at https://github.com/rogersce/cnpy.
 * This 'cnpy' library is already installed and built at /home/ewunpyae/scalapack-playground/cnpy/build_2/ on Titlis.
 * 
 * The implmentation of SVD computation is inspired by: 
 * https://www.ibm.com/docs/en/pessl/5.3.0?topic=analysis-pdgesvd-pzgesvd-singular-value-decomposition-general-matrix#lgesvd__gdgesvd
 * https://github.com/saleemayman/scalapack_svd
 * 
 * Implementation of CBLACS grids and block cyclic distribution is inspired by: 
 * https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
 * 
 * The 'cnpy' library is already installed and built at /home/ewunpyae/scalapack-playground/cnpy/build_2/ on Titlis.
 * Please activate ``conda activate ewp-project-env`` on conda to run this file. 
 * 
 */

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <experimental/filesystem>

#include "/home/ewunpyae/scalapack-playground/cnpy/build_2/include/cnpy.h"

namespace fs = std::experimental::filesystem;
using namespace std;

/*
 * ScaLAPACK functions which are used for the SVD computation
 */
extern "C"
{
    void Cblacs_pinfo(int *, int *);
    void Cblacs_get(int, int, int *);
    void Cblacs_gridinit(int *, const char *, int, int);
    void Cblacs_pcoord(int, int, int *, int *);
    void Cblacs_gridexit(int);
    void Cblacs_barrier(int, const char *);
    void Cdgerv2d(int, int, int, double *, int, int, int);
    void Cdgesd2d(int, int, int, double *, int, int, int);

    int numroc_(int *, int *, int *, int *, int *);

    void pdgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
                  int *ia, int *ja, int *desc_a, double *s, double *u, int *iu,
                  int *ju, int *descu, double *vt, int *ivt, int *jvt,
                  int *desc_vt, double *work, int *lwork, int *info);

    void descinit_(int *, int *, int *, int *, int *, int *, int *,
                   int *, int *, int *);
}

int main(int argc, char **argv)
{
    int mpirank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    string fname(argv[1]);
    cnpy::npz_t my_npz = cnpy::npz_load(fname);
    /*
     * read the Numpy array inside the Numpy zip archieve
     */
    for (auto &iter : my_npz)
    {
        bool mpiroot = (mpirank == 0);
        /* Helping vars */
        int iZERO = 0;

        if (argc < 7)
        {
            if (mpiroot)
                cerr << "Usage: decompose-svd .npzfile procRow procCol Nb Mb compression_ratio" << endl;
            MPI_Finalize();
            return 1;
        }

        int N, M, Nb, Mb, procrows, proccols, p, compression_ratio;
        int info = 0;

        /*
         * initialize double heap arrays for the global matrix and local matrices,
         * and array descriptions required
         */
        double *A_glob = NULL, *A_loc = NULL;
        double *U_loc = NULL, *Vt_loc = NULL, *sigma = NULL, *work = NULL;
        double *U_glob = NULL, *Vt_glob = NULL;
        int descA[9], descU[9], descVT[9];

        /* 
         * accept the grid size, block size, and the compression ratio from the command
         */
        stringstream stream;
        stream << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6];
        stream >> procrows >> proccols >> Nb >> Mb >> compression_ratio;

        /* 
         * Initialize BLACS context to initialize grids size and coordinates of the parallel processes 
         */
        int ctxt, myid, myrow, mycol, numproc;

        Cblacs_pinfo(&myid, &numproc);

        Cblacs_get(0, 0, &ctxt);

        Cblacs_gridinit(&ctxt, "Row-major", procrows, proccols);

        Cblacs_pcoord(ctxt, myid, &myrow, &mycol);

        if (mpiroot)
        {
            /*
         * copy the values from the Numpy array to the C++ double heap array
         * Here, the elements inside the Numpy array are placed in the row-major order format, 
         * whereas Fortran access the matrix in the column-major order format. So, the elements will be copied
         * to the C++ heap array in the column-major order format for the easier data distribution and SVD
         * computation with ScaLAPACK functions. 
         */
            cnpy::NpyArray my_npy = my_npz[iter.first];

            N = static_cast<int>(my_npy.shape[0]);
            M = static_cast<int>(my_npy.shape[1]);
            p = min(N, M);

            double *my_arr = my_npy.data<double>();
            A_glob = new double[N * M];

            U_glob = new double[N * p];
            Vt_glob = new double[p * M];

            for (int r = 0; r < N; r++)
            {
                for (int c = 0; c < M; c++)
                {
                    A_glob[N * c + r] = my_arr[M * r + c];
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
        double start = MPI_Wtime();

        /*
         * broadcast the original matrix size, grid size, block size and compression_ratio 
         * to all parallel processes. 
         */
        int dimensions[8];
        if (mpiroot)
        {
            dimensions[0] = procrows;
            dimensions[1] = proccols;
            dimensions[2] = N;
            dimensions[3] = M;
            dimensions[4] = Nb;
            dimensions[5] = Mb;
            dimensions[6] = p;
            dimensions[7] = compression_ratio;
        }
        // cout << procrows << proccols << endl;
        MPI_Bcast(dimensions, 8, MPI_INT, 0, MPI_COMM_WORLD);
        procrows = dimensions[0];
        proccols = dimensions[1];
        N = dimensions[2];
        M = dimensions[3];
        Nb = dimensions[4];
        Mb = dimensions[5];
        p = dimensions[6];
        compression_ratio = dimensions[7];

        /* 
         * Determine the row and column size of the local input matrix A_loc, U_loc, Vt_loc 
         * on each parallel process
         */
        int nrows_a = numroc_(&N, &Nb, &myrow, &iZERO, &procrows);
        int ncols_a = numroc_(&M, &Mb, &mycol, &iZERO, &proccols);
        int nrows_u = numroc_(&N, &Nb, &myrow, &iZERO, &procrows);
        int ncols_u = numroc_(&p, &Mb, &mycol, &iZERO, &proccols);
        int nrows_vt = numroc_(&p, &Nb, &myrow, &iZERO, &procrows);
        int ncols_vt = numroc_(&M, &Mb, &mycol, &iZERO, &proccols);

        for (int id = 0; id < numproc; ++id)
        {
            Cblacs_barrier(ctxt, "All");
        }

        A_loc = new double[nrows_a * ncols_a];

        U_loc = new double[nrows_u * ncols_u];

        Vt_loc = new double[nrows_vt * ncols_vt];

        sigma = new double[min(N, M)];

        /*
         * Distribute the global matrix A to parallel processes via block-cyclic distribution
         */

        int sendrow = 0, sendcol = 0, recvrow = 0, recvcol = 0;
        for (int r = 0; r < N; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Determine the number of rows to be sent to the parallel process
            int nr = Nb;
            if (N - r < Nb)
                nr = N - r;

            for (int c = 0; c < M; c += Mb, sendcol = (sendcol + 1) % proccols)
            {
                // Determine the number of columns to be sent to the parallel process
                int nc = Mb;
                if (M - c < Mb)
                    nc = M - c;

                if (mpiroot)
                {
                    // Send a a sub-matrix of size (nr x nc) to the process whose coordinate is
                    // (sendrow, sendcol) on the processor grid
                    Cdgesd2d(ctxt, nr, nc, A_glob + N * c + r, N, sendrow, sendcol);
                }

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Receive the sub-matrix sent from the root process
                    Cdgerv2d(ctxt, nr, nc, A_loc + nrows_a * recvcol + recvrow, nrows_a, 0, 0);

                    recvcol = (recvcol + nc) % ncols_a;
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % nrows_a;
        }
        /*
         * Define array descriptors of the global matrices A, U, and Vt.
         */
        int llda = max(1, nrows_a);
        int lldu = max(1, nrows_u);
        int lldvt = max(1, nrows_vt);

        // DescInit
        info = 0;
        descinit_(descA, &N, &M, &Nb, &Mb, &iZERO, &iZERO, &ctxt, &llda, &info);

        if (mpiroot)
        {
            if (info < 0)
            {
                cout << "'descA' Error Info < 0: if INFO = -i, the i-th argument had an illegal value" << endl
                     << "Info = " << info << endl;
                exit(1);
            }
        }

        info = 0;
        descinit_(descU, &N, &p, &Nb, &Mb, &iZERO, &iZERO, &ctxt, &lldu, &info);

        if (mpiroot)
        {
            if (info < 0)
            {
                cout << "descU Error Info < 0: if INFO = -i, the i-th argument had an illegal value" << endl
                     << "Info = " << info << endl;
                exit(1);
            }
        }

        info = 0;

        descinit_(descVT, &p, &M, &Nb, &Mb, &iZERO, &iZERO, &ctxt, &lldvt, &info);

        if (mpiroot)
        {
            if (info < 0)
            {
                cout << "'descVT' Error Info < 0: if INFO = -i, the i-th argument had an illegal value" << endl
                     << "Info = " << info << endl;
                exit(1);
            }
        }

        work = new double[nrows_a * Nb];

        info = 0;
        int ia = 1;
        int ja = 1;
        int iu = 1;
        int ju = 1;
        int ivt = 1;
        int jvt = 1;
        int lwork = -1; //myRows * blockSize; // ?
        char jobu = 'V';
        char jobvt = 'V';

        /*
         * Firstly, lwork = -1 is defined; no SVD computation will be made yet, 
         * but the the minimum size of work required for the actual SVD computation will be returned.
         * The minimum work size will be returned in work[0].
         */
        pdgesvd_(&jobu, &jobvt, &N, &M, A_loc, &ia, &ja, descA,
                 sigma, U_loc, &iu, &ju, descU,
                 Vt_loc, &ivt, &jvt, descVT,
                 work, &lwork, &info);

        if (info != 0)
        {
            printf("first pdgesvd:: rank: %d, info: %d\n", mpirank, info);
            exit(1);
        }

        /*
         * the actual SVD is calculated after getting the work size from the first pdgesvd_() call.
         */
        lwork = work[0];
        work = new double[lwork];

        /*
         * Making sure that all processes are ready for the SVD computation 
         * The elapsed time for the SVD computation on each process is also measured.  
         */
        MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
        double start_svd = MPI_Wtime();
        pdgesvd_(&jobu, &jobvt, &N, &M, A_loc, &ia, &ja, descA,
                 sigma, U_loc, &iu, &ju, descU,
                 Vt_loc, &ivt, &jvt, descVT,
                 work, &lwork, &info);

        if (info != 0)
        {
            printf("second pdgesvd:: rank: %d, info: %d\n", mpirank, info);
            exit(1);
        }

        MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
        double end_svd = MPI_Wtime();
        double elapsed_time_svd = end_svd - start_svd;

        for (int id = 0; id < numproc; ++id)
        {
            Cblacs_barrier(ctxt, "All");
        }

        /*
         * The elements of U and Vt are on local on each process, so they will be gathered back 
         * on the root process according to the coordinate of the processe on the grid and the location of 
         * each element on the local matrix. 
         */

        /*
         * Send the elements from the local arrays of U on parallel processes to the 
         * global array U on the root process
         */
        sendrow = 0;
        for (int r = 0; r < N; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Check whether this is the last row block
            int nr = Nb;
            if (N - r < Nb)
                nr = N - r;

            for (int c = 0; c < p; c += Mb, sendcol = (sendcol + 1) % proccols)
            {
                // Check whether this is the last column block
                int nc = Mb;
                if (p - c < Mb)
                    nc = p - c;

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Send a submatrix of size (nr x nc) to the root process
                    Cdgesd2d(ctxt, nr, nc, U_loc + nrows_u * recvcol + recvrow, nrows_u, 0, 0);
                    recvcol = (recvcol + nc) % ncols_u;
                }

                if (mpiroot)
                {
                    // Receive the submatrix send from parallel processes
                    Cdgerv2d(ctxt, nr, nc, U_glob + N * c + r, N, sendrow, sendcol);
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % nrows_u;
        }

        for (int r = 0; r < p; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Check whether this is the last row block
            int nr = Nb;
            if (p - r < Nb)
                nr = p - r;

            for (int c = 0; c < M; c += Mb, sendcol = (sendcol + 1) % proccols)
            {
                // Check whether this is the last column block
                int nc = Mb;
                if (M - c < Mb)
                    nc = M - c;

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Send a submatrix of size (nr x nc) to the root process
                    Cdgesd2d(ctxt, nr, nc, Vt_loc + nrows_vt * recvcol + recvrow, nrows_vt, 0, 0);
                    recvcol = (recvcol + nc) % ncols_vt;
                }

                if (mpiroot)
                {
                    // Receive the submatrix send from parallel processes
                    Cdgerv2d(ctxt, nr, nc, Vt_glob + p * c + r, p, sendrow, sendcol);
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % nrows_vt;
        }

        /*
         * Select the maximum elapsed time for the SVD computation among paralle proccesses
         */
        double svd;
        MPI_Reduce(&elapsed_time_svd, &svd, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (mpiroot)
        {
            /*
            * count the number of total singular values greater than zero. 
            */
            int org_count = 0;
            int count;
            for (int r = 0; r < p; ++r)
            {
                if (sigma[r] > 0)
                {
                    org_count++;
                }
            }

            /*
             * Retrieve the compression ratio defined by the user in the command line, 
             * and the total number of singular values and vectors to be saved as the compressed data
             * will be determined. 
             * 'count' represents the total number of singular values and vectors to be saved.
             */
            double ratio = (1.0 * compression_ratio) / (1.0 * 100);
            cout << "compression_ratio " << compression_ratio << " " << ratio << endl;
            count = org_count * ratio;

            
            if (count % 2 == 1)
                count += 1;

            double *finalU = new double[N * count];
            double *finalSigma = new double[count];
            double *finalVt = new double[count * M];

            ulong NLong = static_cast<ulong>(N);
            ulong MLong = static_cast<ulong>(M);
            ulong countLong = static_cast<ulong>(count);

             /*
             * Copy the first 'count' singular values and vectors to new double heaps to be saved as 
             * compressed Numpy zips
             * The values inside the decomposed matrices are saved in the column major format. 
             */
            for (int r = 0; r < N; r++)
            {
                for (int c = 0; c < count; c++)
                {
                    finalU[N * c + r] = *(U_glob + N * c + r);
                }
            }
            for (int r = 0; r < count; r++)
            {
                finalSigma[r] = *(sigma + r);
            }
            for (int r = 0; r < count; r++)
            {
                for (int c = 0; c < M; c++)
                {
                    // column-wise calculation
                    finalVt[count * c + r] = *(Vt_glob + p * c + r);
                }
            }

            cout << "N: " << N << "  p: " << count << "  M: " << M << endl;
            /*
             * a new directory with the name $(original-data-folder-name)_$(compression_ratio)$
             * will be created. 
             * the name of the compressed zip file will be the same as the name of the zip for the 
             * original data. 
             */

            string compressedFileName = fname.substr(fname.find_last_of("/\\") + 1);
            string orgPath = fname.substr(0, fname.find_last_of("/\\"));

            string sceneName = orgPath.substr(0, orgPath.find_last_of("/\\"));
            string frameName = orgPath.substr(orgPath.find_last_of("/\\") + 1);

            string parentCompressedPath = sceneName + "_" + to_string(compression_ratio) + "/";
            string compressedPath = parentCompressedPath + frameName + "/";
            string zipName = compressedPath + compressedFileName;

            cout << compressedPath << endl;

            if (!fs::is_directory(compressedPath) || !fs::exists(compressedPath))
            {
                fs::create_directories(compressedPath);
            }
            
            cnpy::npz_save(zipName, iter.first + "_u", &finalU[0], {NLong, countLong}, "w");
            cnpy::npz_save(zipName, iter.first + "_sigma", &finalSigma[0], {countLong}, "a");
            cnpy::npz_save(zipName, iter.first + "_vt", &finalVt[0], {countLong, MLong}, "a");

            /* Deallocate arrays */
            delete[] finalU;
            delete[] finalSigma;
            delete[] finalVt;

            cout << svd << endl;
            string timeFileName = orgPath + "_timetakens";
            ofstream timeFile(timeFileName.c_str(), fstream::in | fstream::out | fstream::app);
            timeFile << org_count << " " << count << " " << numproc << " " << procrows << " " << proccols << " " << Nb << " " << Mb << " " << svd << endl;
            timeFile.close();
        }

        /* Deallocate arrays */
        delete[] A_glob;
        delete[] U_glob;
        delete[] Vt_glob;
        delete[] A_loc;
        delete[] U_loc;
        delete[] Vt_loc;
        delete[] sigma;
        delete[] work;

        Cblacs_gridexit(ctxt);
    }
    MPI_Finalize();

    return 0;
}