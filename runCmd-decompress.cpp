/*
 * Decompress the compressed data of the sample fluid simulation (fire_combined)
 * 
 */

#include <fstream>
#include <iostream>
#include <experimental/filesystem>
#include <string>
#include <ctime>
#include <chrono>

using namespace std;
namespace fs = std::experimental::filesystem;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

// g++ -o runCmd runCmd.cpp -std=c++11 -lstdc++fs
int main(int argc, char **argv)
{
    int procRow = 1, procCol = 1, procNum = 1;

    string compressed_parent_path_10 = "./fire_10/frame_0018";
    string compressed_parent_path_20 = "./fire_20/frame_0018";
    string compressed_parent_path_30 = "./fire_30/frame_0018";
    string compressed_parent_path_50 = "./fire_50/frame_0018";

    for (auto &p : fs::directory_iterator(compressed_parent_path_10))
    {
        string path = p.path();
        string decompressCmd = "mpirun -np " + to_string(procNum) + " ./decompress-values-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol);
        system(decompressCmd.c_str());
    }

    for (auto &p : fs::directory_iterator(compressed_parent_path_20))
    {
        string path = p.path();
        string decompressCmd = "mpirun -np " + to_string(procNum) + " ./decompress-values-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol);
        system(decompressCmd.c_str());
    }

    for (auto &p : fs::directory_iterator(compressed_parent_path_30))
    {
        string path = p.path();
        string decompressCmd = "mpirun -np " + to_string(procNum) + " ./decompress-values-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol);
        system(decompressCmd.c_str());
    }

    for (auto &p : fs::directory_iterator(compressed_parent_path_50))
    {
        string path = p.path();
        string decompressCmd = "mpirun -np " + to_string(procNum) + " ./decompress-values-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol);
        system(decompressCmd.c_str());
    }

    return 0;
}