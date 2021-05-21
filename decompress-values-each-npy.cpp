/*
 * This program is to decompress the compressed fluid data on parallel processes. 
 * The input Numpy zip file needs to contain 3 Numpy arrays: Sigma, U, and Vt.
 * 
 * 
 * C++ array <---> Numpy zip archive convertion is done by using the 'cnpy' library,
 * whose repository is at https://github.com/rogersce/cnpy.
 * 
 * The matrix multiplication is inspired by: 
 * https://www.ibm.com/docs/en/pessl/5.3.0?topic=lp-pdgemm-pzgemm-matrix-matrix-product-general-matrix-its-transpose-its-conjugate-transpose
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
#include <experimental/filesystem>
#include "/home/ewunpyae/scalapack-playground/cnpy/build_2/include/cnpy.h"

using namespace std;
namespace fs = std::experimental::filesystem;

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

    void pdgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                 const double *alpha, const double *a, const int *ia, const int *ja, const int *descX,
                 const double *b, const int *ib, const int *jb, const int *descY, const double *beta,
                 double *c, const int *ic, const int *jc, const int *descZ);

    void descinit_(int *, int *, int *, int *, int *, int *, int *,
                   int *, int *, int *);
}

int main(int argc, char **argv)
{
    int mpirank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    bool mpiroot = (mpirank == 0);

    string fname(argv[1]);
    cnpy::npz_t my_npz = cnpy::npz_load(fname);

    auto start = my_npz.begin();

    int iZERO = 0;

    if (argc < 4)
    {
        if (mpiroot)
            cerr << "Usage: matrixTest matrixNPZ procRow procCol" << endl;
        MPI_Finalize();
        return 1;
    }

    int N, M, Nb, Mb, pb, procrows, proccols, p;
    int ione = 1;
    char notrans = 'N';
    /*
     * initialize double heap arrays for the global matrix and local matrices,
     * and array descriptions required
     */
    double *local_U = NULL, *local_Vt = NULL, *local_sigma = NULL, *work = NULL;
    double *global_U = NULL, *global_Vt = NULL, *global_sigma_1d = NULL, *global_sigma = NULL;

    double *global_u_sigma_mul = NULL, *local_u_sigma_mul = NULL;
    double *global_u_s_vt_mul = NULL, *local_u_s_vt_mul = NULL;

    stringstream stream;
    stream << argv[2] << " " << argv[3];
    stream >> procrows >> proccols;
    
    /*
     * The Numpy zip contains 3 Numpy arrays: Sigma,U,Vt. 
     * The values from these Numpy arrays will be retrieved here.
     */
    cnpy::NpyArray u_npy = my_npz[next(start, 1)->first];
    cnpy::NpyArray sigma_npy = my_npz[start->first];
    cnpy::NpyArray vt_npy = my_npz[next(start, 2)->first];

    N = u_npy.shape[0];
    p = sigma_npy.shape[0];
    M = vt_npy.shape[1];

    Nb = N / procrows;
    pb = p / proccols;
    Mb = M / proccols;

    /*
     * This is for conditional checking to make sure that there is at least one singular value 
     * greater than zero for the decompression. 
     *
     */
    if (p > 0 && pb > proccols)
    {
        if (mpiroot)
        {

            double *u_arr = u_npy.data<double>();
            double *sigma_arr = sigma_npy.data<double>();
            double *vt_arr = vt_npy.data<double>();

            /*
             * Copy elements inside the Numpy arrays to C++ heap arrays
             */

            global_U = new double[N * p];
            global_Vt = new double[p * M];
            global_sigma_1d = new double[p];
            global_sigma = new double[p * p];
            global_u_sigma_mul = new double[N * p];
            global_u_s_vt_mul = new double[N * M];

            for (int r = 0; r < N; r++)
            {
                for (int c = 0; c < p; c++)
                {
                    global_U[N * c + r] = u_arr[N * c + r];
                }
            }

            for (int r = 0; r < p; r++)
            {
                for (int c = 0; c < M; c++)
                {
                    global_Vt[p * c + r] = vt_arr[p * c + r];
                }
            }

            for (int r = 0; r < p; r++)
            {
                *(global_sigma_1d + r) = *(sigma_arr + r);
            }
            /*
             * The singular values are saved inside the Sigma Numpy array in a descending order, 
             * so these singular values will be placed in the diagonal of a new 2D array for the easier
             * multiplication.   
             */
            for (int r = 0; r < p; r++)
            {
                for (int c = 0; c < p; c++)
                {
                    if (r == c)
                    {
                        *(global_sigma + p * c + r) = *(global_sigma_1d + r);
                    }
                    else
                    {
                        *(global_sigma + p * c + r) = 0;
                    }
                }
            }
        }
        /* 
         * Initialize BLACS context to initialize grids size and coordinates of the parallel processes 
         */
        int ctxt, myid, myrow, mycol, numproc;
        Cblacs_pinfo(&myid, &numproc);
        Cblacs_get(0, 0, &ctxt);
        Cblacs_gridinit(&ctxt, "Row-major", procrows, proccols);
        Cblacs_pcoord(ctxt, myid, &myrow, &mycol);

        /*
         * broadcast the original matrix size, grid size, block size for U,Sigma and Vt
         * to all parallel processes. 
         */
        int dimensions[8];
        if (mpiroot)
        {
            dimensions[0] = procrows;
            dimensions[1] = proccols;
            dimensions[2] = N;
            dimensions[3] = p;
            dimensions[4] = M;
            dimensions[5] = Nb;
            dimensions[6] = Mb;
            dimensions[7] = pb;
        }
        MPI_Bcast(dimensions, 8, MPI_INT, 0, MPI_COMM_WORLD);
        procrows = dimensions[0];
        procrows = dimensions[1];
        N = dimensions[2];
        p = dimensions[3];
        M = dimensions[4];
        Nb = dimensions[5];
        Mb = dimensions[6];
        pb = dimensions[7];

        for (int id = 0; id < numproc; id++)
        {
            Cblacs_barrier(ctxt, "All");
        }

        /*
         * U and Sigma arrays will be firstly multiplied, and the multiplication result will be 
         * multiplied with Vt. 
         */

        /* 
         * Determine the row and column size of the local input matrix A_loc, U_loc, Vt_loc 
         * on each parallel process
        */
        int N_loc = numroc_(&N, &Nb, &myrow, &iZERO, &procrows);
        int p_loc = numroc_(&p, &pb, &mycol, &iZERO, &proccols);

        for (int id = 0; id < numproc; ++id)
        {
            Cblacs_barrier(ctxt, "All");
        }

        local_U = new double[N_loc * p_loc];
        local_sigma = new double[p_loc * p_loc];
        local_u_sigma_mul = new double[N_loc * p_loc];

        /*
         * Distribute the global matrix U to parallel processes via block-cyclic distribution
         */
        int sendrow = 0, sendcol = 0, recvrow = 0, recvcol = 0;
        for (int r = 0; r < N; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Determine the number of rows to be sent to the parallel process
            int nr = Nb;
            if (N - r < Nb)
                nr = N - r;

            for (int c = 0; c < p; c += pb, sendcol = (sendcol + 1) % proccols)
            {
                // Determine the number of columns to be sent to the parallel process
                int nc = pb;
                if (p - c < pb)
                    nc = p - c;

                if (mpiroot)
                {
                    // Send a a sub-matrix of size (nr x nc) to the process whose coordinate is
                    // (sendrow, sendcol) on the processor grid
                    Cdgesd2d(ctxt, nr, nc, global_U + N * c + r, N, sendrow, sendcol);
                }

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Receive the sub-matrix sent from the root process
                    Cdgerv2d(ctxt, nr, nc, local_U + N_loc * recvcol + recvrow, N_loc, 0, 0);
                    recvcol = (recvcol + nc) % p_loc;
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % N_loc;
        }
        /*
         * Distribute the global matrix Sigma to parallel processes via block-cyclic distribution
         */
        sendrow = 0, sendcol = 0, recvrow = 0, recvcol = 0;
        for (int r = 0; r < p; r += pb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            int nr = pb;
            if (p - r < pb)
                nr = p - r;

            for (int c = 0; c < p; c += pb, sendcol = (sendcol + 1) % proccols)
            {
                int nc = pb;
                if (p - c < pb)
                    nc = p - c;

                if (mpiroot)
                {
                    // Send a a sub-matrix of size (nr x nc) to the process whose coordinate is
                    // (sendrow, sendcol) on the processor grid
                    Cdgesd2d(ctxt, nr, nc, global_sigma + p * c + r, p, sendrow, sendcol);
                }

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Receive the sub-matrix sent from the root process
                    Cdgerv2d(ctxt, nr, nc, local_sigma + p_loc * recvcol + recvrow, p_loc, 0, 0);
                    recvcol = (recvcol + nc) % p_loc;
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % p_loc;
        }
        /*
         * Define array descriptors of the global matrices U, Sigma, and the multiplication result U_Sigma.
         */
        int descX[9];
        int descY[9];
        int descZ[9];
        int info;
        int lddX = N_loc > 1 ? N_loc : 1;
        int lddY = p_loc > 1 ? p_loc : 1;
        int lddZ = N_loc > 1 ? N_loc : 1;
        descinit_(descX, &N, &p, &Nb, &pb, &iZERO, &iZERO, &ctxt, &lddX, &info);
        if (info != 0)
        {
            printf("Error in descinit_descX, info = %d\n", info);
        }
        descinit_(descY, &p, &p, &pb, &pb, &iZERO, &iZERO, &ctxt, &lddY, &info);
        if (info != 0)
        {
            printf("Error in descinit_descY, info = %d\n", info);
        }
        descinit_(descZ, &N, &p, &Nb, &pb, &iZERO, &iZERO, &ctxt, &lddZ, &info);
        if (info != 0)
        {
            printf("Error in descinit_descZ, info = %d\n", info);
        }

        double alpha = 1.0;
        double beta = 0.0;
        /*
         * U and Sigma, which are already distributed on parallel processes, are multplied,
         * and the result is saved in local_u_sigma_mul locally on these parallel processes.  
         */
        MPI_Barrier(MPI_COMM_WORLD);
        double start_u_s_multiply = MPI_Wtime();
        pdgemm_(&notrans, &notrans, &N, &p, &p,
                &alpha,
                local_U, &ione, &ione, descX,
                local_sigma, &ione, &ione, descY,
                &beta,
                local_u_sigma_mul, &ione, &ione, descZ);

        if (info != 0)
            printf("rank: %d, info: %d\n", mpirank, info);

        MPI_Barrier(MPI_COMM_WORLD);
        double end_u_s_multiply = MPI_Wtime();
        double elapsed_time_u_s_multiply = end_u_s_multiply - start_u_s_multiply;

        for (int id = 0; id < numproc; ++id)
        {
            Cblacs_barrier(ctxt, "All");
        }
        /*
         * Gather the local sub-matrices of the multiplication result on parallel processes to the root process. 
         */
        sendrow = 0;
        for (int r = 0; r < N; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Check whether this is the last row block
            int nr = Nb;
            if (N - r < Nb)
                nr = N - r;

            for (int c = 0; c < p; c += pb, sendcol = (sendcol + 1) % proccols)
            {
                // Check whether this is the last column block
                int nc = pb;
                if (p - c < pb)
                    nc = p - c;

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Send a submatrix of size (nr x nc) to the root process
                    Cdgesd2d(ctxt, nr, nc, local_u_sigma_mul + N_loc * recvcol + recvrow, N_loc, 0, 0);
                    recvcol = (recvcol + nc) % p_loc;
                }

                if (mpiroot)
                {
                    // Receive the sent sub-matrix on the root process and save it in a global matrix
                    Cdgerv2d(ctxt, nr, nc, global_u_sigma_mul + N * c + r, N, sendrow, sendcol);
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % N_loc;
        }

        /**** U-Sigma multiplication process is finished. ****/
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        /*
         * Allocate memory required for the local arrays to multiply 
         * U-sigma and Vt on parallel processes
         */

        int M_loc = numroc_(&M, &Mb, &mycol, &iZERO, &proccols);
        local_Vt = new double[p_loc * M_loc];
        for (int i = 0; i < p_loc * M_loc; ++i)
            *(local_Vt + i) = 0.;

        local_u_s_vt_mul = new double[N_loc * M_loc];
        for (int i = 0; i < N_loc * M_loc; ++i)
            *(local_u_s_vt_mul + i) = 0.;

        for (int i = 0; i < N_loc * p_loc; ++i)
            *(local_u_sigma_mul + i) = 0.;

        /*
         * Distribute the global matrix U-Sigma to parallel processes via block-cyclic distribution
         */
        sendrow = 0, sendcol = 0, recvrow = 0, recvcol = 0;
        for (int r = 0; r < N; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Check whether this is the last row blocks
            int nr = Nb;
            if (N - r < Nb)
                nr = N - r;

            for (int c = 0; c < p; c += pb, sendcol = (sendcol + 1) % proccols)
            {
                // Check whether this is the last column block
                int nc = pb;
                if (p - c < pb)
                    nc = p - c;

                if (mpiroot)
                {
                    // Send a a sub-matrix of size (nr x nc) to the process whose coordinate is
                    // (sendrow, sendcol) on the processor grid
                    Cdgesd2d(ctxt, nr, nc, global_u_sigma_mul + N * c + r, N, sendrow, sendcol);
                }

                if (myrow == sendrow && mycol == sendcol)
                {
                    // paralell processes receive the sub-matrix sent from the root process respectively.
                    Cdgerv2d(ctxt, nr, nc, local_u_sigma_mul + N_loc * recvcol + recvrow, N_loc, 0, 0);
                    recvcol = (recvcol + nc) % p_loc;
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % N_loc;
        }
        /*
         * Distribute the global matrix Vt to parallel processes via block-cyclic distribution
         */

        sendrow = 0, sendcol = 0, recvrow = 0, recvcol = 0;
        for (int r = 0; r < p; r += pb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Check whether this is the last row block
            int nr = pb;
            if (p - r < pb)
                nr = p - r;

            for (int c = 0; c < M; c += Mb, sendcol = (sendcol + 1) % proccols)
            {
                // Check whether this is the last column block
                int nc = Mb;
                if (M - c < Mb)
                    nc = M - c;

                if (mpiroot)
                {
                    // Send a a sub-matrix of size (nr x nc) to the process whose coordinate is
                    // (sendrow, sendcol) on the processor grid
                    Cdgesd2d(ctxt, nr, nc, global_Vt + p * c + r, p, sendrow, sendcol);
                }

                if (myrow == sendrow && mycol == sendcol)
                {
                    // paralell processes receive the sub-matrix sent from the root process respectively.
                    Cdgerv2d(ctxt, nr, nc, local_Vt + p_loc * recvcol + recvrow, p_loc, 0, 0);
                    recvcol = (recvcol + nc) % M_loc;
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % p_loc;
        }

        /*
         * Define array descriptors of the global matrices U-Sigma, Vt, and the multiplication result U_S_Vt.
         */
        lddX = N_loc > 1 ? N_loc : 1;
        lddY = p_loc > 1 ? p_loc : 1;
        lddZ = N_loc > 1 ? N_loc : 1;
        descinit_(descX, &N, &p, &Nb, &pb, &iZERO, &iZERO, &ctxt, &lddX, &info);
        if (info != 0)
        {
            printf("Error in descinit_descX, info = %d\n", info);
        }
        descinit_(descY, &p, &M, &pb, &Mb, &iZERO, &iZERO, &ctxt, &lddY, &info);
        if (info != 0)
        {
            printf("Error in descinit_descY, info = %d\n", info);
        }
        descinit_(descZ, &N, &M, &Nb, &Mb, &iZERO, &iZERO, &ctxt, &lddZ, &info);
        if (info != 0)
        {
            printf("Error in descinit_descZ, info = %d\n", info);
        }

        alpha = 1.0;
        beta = 0.0;
        /*
         * U-Sigma and VT, which are already distributed on parallel processes, are multplied,
         * and the result is saved in local_u_s_vt_mul locally on these parallel processes.  
         */
        MPI_Barrier(MPI_COMM_WORLD);
        double start_u_s_vt_multiply = MPI_Wtime();
        pdgemm_(&notrans, &notrans, &N, &M, &p,
                &alpha,
                local_u_sigma_mul, &ione, &ione, descX,
                local_Vt, &ione, &ione, descY,
                &beta,
                local_u_s_vt_mul, &ione, &ione, descZ);

        if (info != 0)
            printf("rank: %d, info: %d\n", mpirank, info);

        MPI_Barrier(MPI_COMM_WORLD);
        double end_u_s_vt_multiply = MPI_Wtime();
        double elapsed_time_u_s_vt_multiply = end_u_s_vt_multiply - start_u_s_vt_multiply;

        for (int id = 0; id < numproc; ++id)
        {
            Cblacs_barrier(ctxt, "All");
        }

        /*
         * Gather the local sub-matrices of the multiplication result on parallel processes to the root process. 
         */
        sendrow = 0;
        for (int r = 0; r < N; r += Nb, sendrow = (sendrow + 1) % procrows)
        {
            sendcol = 0;
            // Number of rows to be sent
            // Is this the last row block?
            int nr = Nb;
            if (N - r < Nb)
                nr = N - r;

            for (int c = 0; c < M; c += Mb, sendcol = (sendcol + 1) % proccols)
            {
                // Number of cols to be sent
                // Is this the last col block?
                int nc = Mb;
                if (M - c < Mb)
                    nc = M - c;

                if (myrow == sendrow && mycol == sendcol)
                {
                    // Send the submatrix of size (nr x nc) from each parallel process
                    Cdgesd2d(ctxt, nr, nc, local_u_s_vt_mul + N_loc * recvcol + recvrow, N_loc, 0, 0);
                    recvcol = (recvcol + nc) % M_loc;
                }

                if (mpiroot)
                {
                    // Receive the sent sub-matrix on the root process and save it in a global matrix
                    Cdgerv2d(ctxt, nr, nc, global_u_s_vt_mul + N * c + r, N, sendrow, sendcol);
                }
            }

            if (myrow == sendrow)
                recvrow = (recvrow + nr) % N_loc;
        }
        for (int id = 0; id < numproc; ++id)
        {
            Cblacs_barrier(ctxt, "All");
        }
        double u_s_multiply, u_s_vt_multiply;
        MPI_Reduce(&elapsed_time_u_s_multiply, &u_s_multiply, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&elapsed_time_u_s_vt_multiply, &u_s_vt_multiply, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (mpiroot)
        {
            /*
             * Save the decompressed matrix in the Numpy zip
             * The elements inside the decompressed matrix will be saved in a row-major order format. 
             */
            double *final_result = new double[N * M];

            for (int r = 0; r < N; r++)
            {
                for (int c = 0; c < M; c++)
                {
                    final_result[M * r + c] = global_u_s_vt_mul[N * c + r];
                }
            }

            /*
             * Here, the decompressed matrix will be converted to a Numpy array, and saved in 
             * a new Numpy zip archive. 
             * The directory and the name of the decompressed zip file will be
             * $(compressed_file_path) + '_decompressed'/$(compressed_zipfile_name).
             * The header of the Numpy array will be 'denArr', 'presArr', 'velXArr', or 'velYArr'.
             */
            ulong NLong = static_cast<ulong>(N);
            ulong MLong = static_cast<ulong>(M);

            int idx1 = start->first.find("_");
            string compressedFileName = start->first.substr(0, idx1);

            cout << "conpressedFileName: " << compressedFileName << endl;

            string orgPath = fname.substr(0, fname.find_last_of("/\\"));
            // string sceneName = orgPath.substr(0, orgPath.find_last_of("_"));

            string decompressedPath = orgPath + "_decompressed/";

            string zipName = decompressedPath + fname.substr(fname.find_last_of("/\\") + 1);

            if (!fs::is_directory(decompressedPath) || !fs::exists(decompressedPath))
            {
                fs::create_directory(decompressedPath);
            }
            cnpy::npz_save(zipName, compressedFileName, &final_result[0], {NLong, MLong}, "w");
            /*
            * Save the maximum elapsed time for multplications in a text file
            */

            cout << u_s_multiply << " " << u_s_vt_multiply << endl;
            string timeFileName = orgPath + "_timetakens";
            ofstream timeFile(timeFileName.c_str(), fstream::in | fstream::out | fstream::app);
            timeFile << numproc << " " << procrows << " " << proccols << " " << N << " " << p << " " << M << " " << u_s_multiply << " " << u_s_vt_multiply << endl;
            timeFile.close();
            

            delete[] final_result;
            Cblacs_gridexit(ctxt);
        }
    }

    else
    {
        if (mpiroot)
        {
            double *final_result = new double[N * M];

            for (int r = 0; r < N * M; r++)
            {
                *(final_result + r) = 0.0;
            }

            ulong NLong = static_cast<ulong>(N);
            ulong MLong = static_cast<ulong>(M);

            int idx1 = start->first.find("_");
            string compressedFileName = start->first.substr(0, idx1);

            cout << "conpressedFileName: " << compressedFileName << endl;

            string orgPath = fname.substr(0, fname.find_last_of("/\\"));
            string sceneName = orgPath.substr(0, orgPath.find_last_of("_"));

            string decompressedPath = sceneName + "_decompressed/";

            string zipName = decompressedPath + fname.substr(fname.find_last_of("/\\") + 1);

            if (!fs::is_directory(decompressedPath) || !fs::exists(decompressedPath))
            {
                fs::create_directory(decompressedPath);
            }
            // string mode = k == 3 ? "w" : "a";

            cnpy::npz_save(zipName, compressedFileName, &final_result[0], {NLong, MLong}, "w");
           
            delete[] final_result;
        }
    }
    // Deallocate the heap arrays used for the decompression process
    delete[] local_sigma;
    delete[] local_U;
    delete[] local_Vt;
    delete[] local_u_sigma_mul;
    delete[] local_u_s_vt_mul;
    delete[] work;
    delete[] global_U;
    delete[] global_Vt;
    delete[] global_sigma;
    delete[] global_sigma_1d;
    delete[] global_u_s_vt_mul;
    delete[] global_u_sigma_mul;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}