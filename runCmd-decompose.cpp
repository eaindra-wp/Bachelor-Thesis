/*
 * Decompose the small sample fluid simulation (fire_combined) via SVD with different 
 * compression ratios. 
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


int main(int argc, char **argv)
{
    int procRow = 2, procCol = 2, Nb = 10, Mb = 10, procNum = 4;
    int ratio_10 = 10, ratio_20 = 20, ratio_30 = 30, ratio_50 = 50;
    
    string org_parent_path = "./fire/frame_0018/";
    
    for (auto &p : fs::directory_iterator(org_parent_path))
    {
        string path = p.path();
        string compressCmd = "mpirun -np " + to_string(procNum) + " ./decompose-svd-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol) + " " + to_string(Nb) + " " + to_string(Mb) +  " " + to_string(ratio_10);
        system(compressCmd.c_str());
    }

    for (auto &p : fs::directory_iterator(org_parent_path))
    {
        string path = p.path();
        string compressCmd = "mpirun -np " + to_string(procNum) + " ./decompose-svd-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol) + " " + to_string(Nb) + " " + to_string(Mb) +  " " + to_string(ratio_20);
        system(compressCmd.c_str());
    }

    for (auto &p : fs::directory_iterator(org_parent_path))
    {
        string path = p.path();
        string compressCmd = "mpirun -np " + to_string(procNum) + " ./decompose-svd-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol) + " " + to_string(Nb) + " " + to_string(Mb) +  " " + to_string(ratio_30);
        system(compressCmd.c_str());
    }

    for (auto &p : fs::directory_iterator(org_parent_path))
    {
        string path = p.path();
        string compressCmd = "mpirun -np " + to_string(procNum) + " ./decompose-svd-each-npy " + path + " " + to_string(procRow) + " " + to_string(procCol) + " " + to_string(Nb) + " " + to_string(Mb) +  " " + to_string(ratio_50);
        system(compressCmd.c_str());
    }

    return 0;
}