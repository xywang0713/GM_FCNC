#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "clooptools.h"

int main() {
    ltini();
    ifstream
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    // 调用衰变分支比函数
    std::complex<double> Br_B_Meson_value = Br_B_Meson (i, j, GM_model&m);
    std::complex<double> Br_Kaon_value = Br_Kaon (i, j, GM_model&m);
    std::complex<double> Br_Semileptonic_pai_value = Br_Semileptonic_pai (i, j, GM_model&m);
    std::complex<double> Br_Semileptonic_k_value = Br_Semileptonic_k (i, j, GM_model&m);
    std::complex<double> Br_Semileptonic_D_value = Br_Semileptonic_D (i, j, GM_model&m);
    std::complex<double> Br_Semileptonic_Ds_value = Br_Semileptonic_Ds (i, j, GM_model&m);
    
    
    //调用H5衰变宽度函数
    std::complex<double> GAMMA_gamma_value = GAMMA_gamma(GM_model&m);
    std::complex<double> GAMMA_electron_value = GAMMA_electron(GM_model&m);
    std::complex<double> GAMMA_muon_value = GAMMA_muon(GM_model&m);
    std::complex<double> GAMMA_tau_value = GAMMA_tau(GM_model&m);
    
    std::complex<double> GAMMA_paipai_value

    if (MH5 < 2) {
        GAMMA_paipai_value = GAMMA_paipai(GM_model&m);
    } else {
        GAMMA_paipai_value = "wrong!";
    }

    std::complex<double> GAMMA_KK_value

    if (MH5 < 2) {
        GAMMA_KK_value = GAMMA_KK(GM_model&m);
    } else {
        GAMMA_KK_value = "wrong!";
    }
    

    std::ofstream outputFile("/home/wang/Desktop/FCNC_GM/GM_FCNC/outputresult.txt");
    // 将结果写入文件
    if (outputFile.is_open()) {
        outputFile << "thetaH = " << thetaH << std::endl;
        outputFile << "MH3 = " << MH3 << std::endl;
        outputFile << "MH5 = " << MH5 << std::endl;
        outputFile << "M1 = " << M1 << std::endl;
        outputFile << "M2 = " << M2 << std::endl;
        outputFile << "Br_B_Meson_value = " << Br_B_Meson_value << std::endl;
        outputFile << "Br_Kaon_value = " << Br_Kaon_value << std::endl;
        outputFile << "Br_Semileptonic_pai_value = " << Br_Semileptonic_pai_value << std::endl;
        outputFile << "Br_Semileptonic_k_value = " << Br_Semileptonic_k_value << std::endl;
        outputFile << "Br_Semileptonic_D_value = " << Br_Semileptonic_D_value << std::endl;
        outputFile << "Br_Semileptonic_Ds_value = " << Br_Semileptonic_Ds_value << std::endl;
        outputFile << "GAMMA_gamma_value = " << GAMMA_gamma_value << std::endl;
        outputFile << "GAMMA_electron_value = " << GAMMA_electron_value << std::endl;
        outputFile << "GAMMA_muon_value = " << GAMMA_muon_value << std::endl;
        outputFile << "GAMMA_tau_value = " << GAMMA_tau_value << std::endl;
        outputFile << "GAMMA_paipai_value = " << GAMMA_paipai_value << std::endl;
        outputFile << "GAMMA_KK_value = " << GAMMA_KK_value << std::endl;
        
        outputFile.close();
        std::cout << "Calculation has been written to outputresult.txt." << std::endl;
    } else {
        std::cerr << "There is somthing wrong to open file";
    }
    return 0;
    }
    ltexi();
