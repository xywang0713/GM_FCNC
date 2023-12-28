//直接注释
// GTRL+K选中,CTRL+C注释

// #include <iostream>：用于输入和输出的C++库，使用std::cout和std::cin进行输出和输入
// #include <complex>：C++中用于复数运算的库，创建和操作复数
// #include <cmath>：C++中的数学函数库，std::sqrt、std::sin和std::cos，平方std::pow(底数，次方数)
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "clooptools.h"


// 读取home目录下test.csv文件
void Calculation_test(const std::string& filename) { 
    // std::ifstream file("/home/test.csv");
    // if (!file.is_open()) {
    //     std::cerr << "Failed to open the file: " << filePath << std::endl;
    //     return 1;
    // }

    // // 跳过标题行
    // std::string line;
    // std::getline(file, line);

    // //输出到home目录下的try.csv文件里
    // std::ofstream outputFile("home/try.csv");
    // if (!outputFile.is_open()) {
    //     std::cerr << "Failed to create the output file" << std::endl;
    //     return;
    // }


    // while (std::getline(file, line)) {
    //     std::stringstream iss(line);
    //     std::string val;
    //     std::vector<double> values;
    //     while (std::getline(iss, val, ',')) {
    //         values.push_back(std::stod(val));
    //     }

    //     MH5 = values[0];
    //     sH = values[1];
    //     thetaH = std::asin(sH); 
    //     M1 = values[6];
    //     M2 = values[7];
    //     MH3 = values[9];

            Fa_ij =  xx

            xi_ij_L = 
            xi_ij_R = 
            outputFile << MH5 << "  " << MH3 << xi_ij_L << xi_ij_R<<std::endl;

    //     }
}
// 定义其他变量
std::complex<double> alpha(1/137, 0.0), SW(0.2223, 0.0), vacuum(246.0, 0.0);
double MW = 80.399;
double MH5 = 10;
double MH3 = 50;
double sH = 0.01;
double thetaH = std::asin(sH);
double M1 = 100;
double M2 = 200;

// 计算 LAMBDA_3()  函数
std::complex<double> LAMBDA_3() {
    return (std::cos(thetaH) * std::cos(thetaH) * (-3.0 * std::sin(thetaH) * MH3 * MH3 + std::sqrt(2.0) * M1 * vacuum) + std::sin(thetaH) * (MH5 * MH5 - 3.0 * std::sqrt(2.0) * M2 * std::sin(thetaH) * vacuum)) / (vacuum * std::sin(thetaH) * std::sin(thetaH) * std::sin(thetaH));
}

// 计算 LAMBDA_5() 函数
std::complex<double> LAMBDA_5() {
    return (2.0 * std::sin(thetaH) * MH3 * MH3 - std::sqrt(2.0) * M1 * vacuum) / (vacuum * std::sin(thetaH));
}

// 计算 yd(z) 的函数
std::complex<double> set_yd_z(int z) {
    std::complex<double> yd_z;
    if (z == 1) {
        yd_z = std::sqrt(2.0) * 5.04 * std::pow(10.0, -3) / (vacuum * std::cos(thetaH));
    } else if (z == 2) {
        yd_z = std::sqrt(2.0) * 0.101 / (vacuum * std::cos(thetaH));
    } else if (z == 3) {
        yd_z = std::sqrt(2.0) * 4.7 / (vacuum * std::cos(thetaH));
    } else {
        std::cerr << "z value is wrong!" << std::endl;
    }
    return yd_z;
}


// 计算 yu(y) 的函数
std::complex<double> set_yu_y(int y) {
    std::complex<double> yu_y;
    if (y == 1) {
        yu_y = std::sqrt(2.0) * 2.55 * pow(10.0, -3) / (vacuum * std::cos(thetaH));
    } else if (y == 2) {
        yu_y = std::sqrt(2.0) * 1.27 / (vacuum * std::cos(thetaH));
    } else if (y == 3) {
        yu_y = std::sqrt(2) * 172.0 / (vacuum * std::cos(thetaH));
    } else {
        std::cerr << "y value is wrong!" << std::endl;
    }
    return yu_y;
}

// 计算 Md(x) 的函数
double set_Md_x(int x) {
    double Md_x = 0.0;  // Initialize Md_x with a default value
    if (x == 1) {
        Md_x = 5.04 * pow(10.0, -3);
    } else if (x == 2) {
        Md_x = 0.101;
    } else if (x == 3) {
        Md_x = 4.7;
    } else {
        std::cerr << "x value is wrong!" << std::endl;
    }
    return Md_x;
}

// 计算 Mu(w) 的函数
double set_Mu_w(int w) {
    double Mu_w = 0.0;  // Initialize Mu_w with a default value
    if (w == 1) {
        Mu_w = 2.55 * pow(10.0, -3);
    } else if (w == 2) {
        Mu_w = 1.27;
    } else if (w == 3) {
        Mu_w = 172.0;
    } else {
        std::cerr << "w value is wrong!" << std::endl;
    }
    return Mu_w;
}


// 计算 CKM矩阵V[u,t] 的函数
 std::complex<double> V(int u, int t) {
     if (u == 1 && t == 1) {
         return 0.974688;
     } else if (u == 1 && t == 2) {
         return 0.225;
     } else if (u == 1 && t == 3) {
         return std::complex<double>(0.00157575, -0.00344881);
     } else if (u == 2 && t == 1) {
         return -0.225;
     } else if (u == 2 && t == 2) {
         return 0.974688;
     } else if (u == 2 && t == 3) {
         return 0.0418163;
     } else if (u == 3 && t == 1 ) {
         return std::complex<double>(0.00783291, -0.00344881); 
     } else if (u == 3 && t == 2 ) {
         return -0.0418163; 
     } else if (u == 3 && t == 3 ) { 
         return 1.0;
     } else {
         std::cerr << "Invalid indices u, t" << std::endl;
         return 0;
     }
 }

// 计算 CKMC矩阵VC[r,s] 的函数
 std::complex<double> VC(int r, int s) {
     if (r == 1 && s == 1) {
         return 0.974688;
     } else if (r == 1 && s == 2) {
         return -0.225;
     } else if (r == 1 && s == 3) {
         return std::complex<double>(0.00783291 , 0.00344881);
     } else if (r == 2 && s == 1) {
         return 0.225;
     } else if (r == 2 && s == 2) {
         return 0.974688;
     } else if (r == 2 && s == 3) {
         return -0.0418163;
     } else if (r == 3 && s == 1 ) {
         return std::complex<double>(0.00157575 , 0.00344881); 
     } else if (r == 3 && s == 2 ) {
         return 0.0418163; 
     } else if (r == 3 && s == 3 ) { 
         return 1.0;
     } else {
         std::cerr << "Invalid indices r, s" << std::endl;
         return 0;
     }
 }

// // 计算 CKM矩阵,CKMC矩阵
// std::complex<double> V[3][3] = {{0.974688, 0.225, std::complex<double>(0.00157575, -0.00344881)},
//                                {-0.225, 0.974688, 0.0418163}, 
//                                {std::complex<double>(0.00783291, -0.00344881), -0.0418163, 1}};
// std::complex<double> VC[3][3] = {{0.974688, -0.225, std::complex<double>(0.00783291 , 0.00344881)}, 
//                                 {0.225, 0.974688, -0.0418163}, 
//                                 {std::complex<double>(0.00157575 , 0.00344881), 0.0418163, 1}};

// (RR_Ba*Ba+RR_Bb*Bb+Ca系数*Ca+Cb系数*Cb+Cc系数*Cc+Cd系数*Cd)*P_R+(Ba系数*Ba+Bb系数*Bb+Ca系数*Ca+Cb系数*Cb+Cc系数*Cc+Cd系数*Cd)*P_L

// P_R右手
// （一）计算 RR_Ba 的函数
std::complex<double> RR_Ba (int i, int j){
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(i) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (V(1,i) * VC(1,j) + V(2,i) * VC(2,j) + V(3,i) * VC(3,j));
}

// （二）计算 RR_Bb 的函数
std::complex<double> RR_Bb (int i, int j){
    return (-alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(i) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (V(1,i) * VC(1,j) + V(2,i) * VC(2,j) + V(3,i) * VC(3,j));
}

// 计算 R_Ca 的函数
//  (三) 计算 R_Ca0_k1 的函数
std::complex<double> R_Ca0_k1 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1)) * V(1,i) * VC(1,j));
}

//  (四) 计算 R_Ca0_k2 的函数
std::complex<double> R_Ca0_k2 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2)) * V(2,i) * VC(2,j));
}

//  （五） 计算 R_Ca0_k3 的函数
std::complex<double> R_Ca0_k3 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3)) * V(3,i) * VC(3 ,j));
}

// // （）计算 R_Ca0 的函数
// std::complex<double> R_Ca0 (int i, int j){
//     return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1)) * V(1,i) * VC(1,j) + (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2)) * V(2,i) * VC(2,j) + (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3)) * V(3,i) * VC(3,j));
// }

// （六）计算 R_Ca1_k1 的函数
std::complex<double> R_Ca1_k1 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1) * set_yu_y(1)) * V(1,i) * VC(1,j));
}

// （七）计算 R_Ca1_k2 的函数
std::complex<double> R_Ca1_k2 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2) * set_yu_y(2)) * V(2,i) * VC(2,j));
}

// （八）计算 R_Ca1_k3 的函数
std::complex<double> R_Ca1_k3 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3) * set_yu_y(3)) * V(3,i) * VC(3,j));
}

// // （）计算 R_Ca1 的函数
// std::complex<double> R_Ca1 (int i, int j){
//     return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1) * set_yu_y(1)) * V(1,i) * VC(1,j) + (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2) * set_yu_y(2)) * V(2,i) * VC(2,j) + (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3) * set_yu_y(3)) * V(3,i) * VC(3,j));
// }

// （九）计算 R_Ca2_k1 的函数
std::complex<double> R_Ca2_k1 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *V(1,i) * VC(1,j));
}

// （十）计算 R_Ca2_k2 的函数
std::complex<double> R_Ca2_k2 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *V(2,i) * VC(2,j));
}

// （十一）计算 R_Ca2_k3 的函数
std::complex<double> R_Ca2_k3 (int i, int j){
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *V(3,i) * VC(3,j));
}
// // （）计算 R_Ca2 的函数
// std::complex<double> R_Ca2 (int i, int j){
//     return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *V(1,i) * VC(1,j) + set_Md_x(j) * set_yd_z(i) * set_yd_z(j) * V(2,i) * VC(2,j) + set_Md_x(j) * set_yd_z(i) * set_yd_z(j) * V(3,i) * VC(3,j));
// }

// 计算 R_Cb 的函数
// 计算 R_Cb0 的函数
//  （十二）计算 R_Cb0_k1 的函数
std::complex<double> R_Cb0_k1 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1))\
    -2 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i) + set_Mu_w(1) *set_Mu_w(1) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
           * V(1,i) * VC(1,j);
}

//  （十三）计算 R_Cb0_k2 的函数
std::complex<double> R_Cb0_k2 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2))\
    -2 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i) + set_Mu_w(2) *set_Mu_w(2) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
           * V(2,i) * VC(2,j);
}

//  （十四）计算 R_Cb0_k3 的函数
std::complex<double> R_Cb0_k3 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3))\
    -2 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i) + set_Mu_w(3) *set_Mu_w(3) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
           * V(3,i) * VC(3,j);
}

// // （）计算 R_Cb0 总和
// std::complex<double> R_Cb0 (int i, int j){
//     return R_Cb0_k1(i,j)+R_Cb0_k2(i,j)+R_Cb0_k3(i,j);
// }

// 计算 R_Cb1 的函数
//  （十五）计算 R_Cb1_k1 的函数
std::complex<double> R_Cb1_k1 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Md_x(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) - 2.0 * set_Md_x(i) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1) * set_yu_y(1)))\
           * V(1,i) * VC(1,j);
} 

//  （十六）计算 R_Cb1_k2 的函数
std::complex<double> R_Cb1_k2 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    *(-16.0 * alpha * alpha * set_Md_x(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) - 2.0 * set_Md_x(i) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2) * set_yu_y(2)))\
          * V(2,i) * VC(2,j); 
}

//  （十七）计算 R_Cb1_k3 的函数
std::complex<double> R_Cb1_k3 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Md_x(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) - 2.0 * set_Md_x(i) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3) * set_yu_y(3)))\
          * V(3,i) * VC(3,j); 
}

// // （）计算 R_Cb1 总和
// std::complex<double> R_Cb1 (int i, int j){
//     return R_Cb1_k1(i,j)+R_Cb1_k2(i,j)+R_Cb1_k3(i,j);
// }

// 计算 R_Cb2 = R_Cb2_k1+R_Cb2_k2+R_Cb2_k3 的函数
//  （十八）计算 R_Cb2_k1 的函数
std::complex<double> R_Cb2_k1 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i)))\
          * V(1,i) * VC(1,j);
}

//  （十九）计算 R_Cb2_k2 的函数
std::complex<double> R_Cb2_k2 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i)))\
          * V(2,i) * VC(2,j);
}

//  （二十）计算 R_Cb2_k3 的函数
std::complex<double> R_Cb2_k3 (int i, int j){
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i)))\
           * V(3,i) * VC(3,j); 
}

// // （）计算 R_Cb2 总和
// std::complex<double> R_Cb2 (int i, int j){
//     return R_Cb2_k1(i,j)+R_Cb2_k2(i,j)+R_Cb2_k3(i,j);
// }

// 计算 R_Cc 的函数
// 计算 R_Cc0 = R_Cc0_k1+R_Cc0_k2+R_Cc0_k3 的函数
//  （二十一）计算 R_Cc0_k1 的函数
std::complex<double> R_Cc0_k1 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW *SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - 2.0 * set_Md_x(i) * set_Md_x(i) + set_Mu_w(1) * set_Mu_w(1)) * set_yd_z(i) + 2.0 * set_Md_x(i) * set_Mu_w(1) *set_yu_y(1)))\
           * V(1,i) * VC(1,j);
}

// （二十二）计算 R_Cc0_k2 的函数
std::complex<double> R_Cc0_k2 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW *SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - 2.0 * set_Md_x(i) * set_Md_x(i) + set_Mu_w(2) * set_Mu_w(2)) * set_yd_z(i) + 2.0 * set_Md_x(i) * set_Mu_w(2) *set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

// （二十三）计算 R_Cc0_k3 的函数
std::complex<double> R_Cc0_k3 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW *SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 +(-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - 2.0 * set_Md_x(i) * set_Md_x(i) + set_Mu_w(3) * set_Mu_w(3)) * set_yd_z(i) + 2.0 * set_Md_x(i) * set_Mu_w(3) *set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 R_Cc0 总和
// std::complex<double> R_Cc0 (int i, int j){
//     return R_Cc0_k1(i,j)+R_Cc0_k2(i,j)+R_Cc0_k3(i,j);
// }

// 计算 R_Cc1 = R_Cc1_k1+R_Cc1_k2+R_Cc1_k3 的函数
//  （二十四）计算 R_Cc1_k1 的函数
std::complex<double> R_Cc1_k1 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1)*set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

//  （二十五）计算 R_Cc1_k2 的函数
std::complex<double> R_Cc1_k2 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2)*set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

//  （二十六）计算 R_Cc1_k3 的函数
std::complex<double> R_Cc1_k3 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3)*set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 R_Cc1 总和
// std::complex<double> R_Cc1 (int i, int j){
//     return R_Cc1_k1(i,j)+R_Cc1_k2(i,j)+R_Cc1_k3(i,j);
// }

// 计算 R_Cc2 = R_Cc2_k1+R_Cc2_k2+R_Cc2_k3 的函数
// （二十七）计算 R_Cc2_k1 的函数
std::complex<double> R_Cc2_k1 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j)))\
        * V(1,i) * VC(1,j);
}

//  （二十八）计算 R_Cc1_k2 的函数
std::complex<double> R_Cc2_k2 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j)))\
        * V(2,i) * VC(2,j);
}

//  （二十九）计算 R_Cc2_k3 的函数
std::complex<double> R_Cc2_k3 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32* std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() +LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 R_Cc2 总和
// std::complex<double> R_Cc2 (int i, int j){
//     return R_Cc2_k1(i,j)+R_Cc2_k2(i,j)+R_Cc2_k3(i,j);
// }

// 计算 R_Cd0 = R_Cd0_k1+R_Cd0_k2+R_Cd0_k3 的函数
//  （三十）计算 R_Cd0_k1 的函数
std::complex<double> R_Cd0_k1 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(1) * set_yd_z(i) * set_yu_y(1) + set_Md_x(i) * set_yu_y(1) * set_yu_y(1))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

//  （三十一）计算 R_Cd0_k2 的函数
std::complex<double> R_Cd0_k2 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(2) * set_yd_z(i) * set_yu_y(2) + set_Md_x(i) * set_yu_y(2) *set_yu_y(2))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

//  （三十二）计算 R_Cd0_k3 的函数
std::complex<double> R_Cd0_k3 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(3) * set_yd_z(i) * set_yu_y(3) + set_Md_x(i) * set_yu_y(3) *set_yu_y(3))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 R_Cd0 总和
// std::complex<double> R_Cd0 (int i, int j){
//     return R_Cd0_k1(i,j)+R_Cd0_k2(i,j)+R_Cd0_k3(i,j);
// }

// 计算 R_Cd1 = R_Cd1_k1+R_Cd1_k2+R_Cd1_k3 的函数
//  （三十三）计算 R_Cd1_k1 的函数
std::complex<double> R_Cd1_k1 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) *set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(j) * set_yd_z(i) * set_yd_z(j) + set_Md_x(i) * set_yu_y(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j); 
}

//  （三十四）计算 R_Cd1_k2 的函数
std::complex<double> R_Cd1_k2 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(j) * set_yd_z(i) * set_yd_z(j) + set_Md_x(i) * set_yu_y(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

//  （三十五）计算 R_Cd1_k3 的函数
std::complex<double> R_Cd1_k3 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(j) * set_yd_z(i) * set_yd_z(j) + set_Md_x(i) * set_yu_y(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 R_Cd1 总和
// std::complex<double> R_Cd1 (int i, int j){
//     return R_Cd1_k1(i,j)+R_Cd1_k2(i,j)+R_Cd1_k3(i,j);
// }

// 计算 R_Cd2 = R_Cd2_k1+R_Cd2_k2+R_Cd2_k3 的函数
//  （三十六）计算 R_Cd2_k1 的函数
std::complex<double> R_Cd2_k1 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yu_y(1) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

//  （三十七）计算 R_Cd2_k2 的函数
std::complex<double> R_Cd2_k2 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yu_y(2) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

//  （三十八）计算 R_Cd2_k3 的函数
std::complex<double> R_Cd2_k3 (int i, int j){
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yu_y(3) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

//  （）计算 R_Cd2 总和
std::complex<double> R_Cd2 (int i, int j){
    return R_Cd2_k1(i,j)+R_Cd2_k2(i,j)+R_Cd2_k3(i,j);
}


// P_L左手
// （壹）计算 LL_Ba 的函数
std::complex<double> LL_Ba (int i, int j){
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(j) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (V(1,i) * VC(1,j) + V(2,i) * VC(2,j) + V(3,i) * VC(3,j));
}

// （贰）计算 LL_Bb 的函数
std::complex<double> LL_Bb (int i, int j){
    return (-alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(j) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (V(1,i) * VC(1,j) + V(2,i) * VC(2,j) + V(3,i) * VC(3,j));
}

// 计算 L_Ca 的函数
// 计算 L_Ca0 = L_Ca0_k1+L_Ca0_k2+L_Ca0_k3 的函数
// （叁）计算 L_Ca0_k1 的函数
std::complex<double> L_Ca0_k1 (int i, int j){
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(1) * set_yd_z(j) * set_yu_y(1) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)) * V(1,i) * VC(1,j);
}

// （肆）计算 L_Ca0_k2 的函数
std::complex<double> L_Ca0_k2 (int i, int j){
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(2) * set_yd_z(j) * set_yu_y(2) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)) * V(2,i) * VC(2,j);
}

// （伍）计算 L_Ca0_k3 的函数
std::complex<double> L_Ca0_k3 (int i, int j){
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(3) * set_yd_z(j) * set_yu_y(3) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)) * V(3,i) * VC(3,j);
}

// // （）计算 L_Ca0 总和
// std::complex<double> L_Ca0 (int i, int j){
//     return L_Ca0_k1(i,j)+L_Ca0_k2(i,j)+L_Ca0_k3(i,j);
// }

// （陆）计算 L_Ca1_k1 的函数
std::complex<double> L_Ca1_k1 (int i, int j){
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)) * V(1,i) * VC(1,j);
}

// （柒）计算 L_Ca1_k2 的函数
std::complex<double> L_Ca1_k2 (int i, int j){
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)) * V(2,i) * VC(2,j);
}

// （捌）计算 L_Ca1_k3 的函数
std::complex<double> L_Ca1_k3 (int i, int j){
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)) * V(3,i) * VC(3,j);
}

// （）计算 L_Ca1 总和
std::complex<double> L_Ca1 (int i, int j){
    return L_Ca1_k1(i,j)+L_Ca1_k2(i,j)+L_Ca1_k3(i,j);
}

// （玖）计算 L_Ca2_k1 的函数
std::complex<double> L_Ca2_k1 (int i, int j){
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(j) * set_yu_y(1) * set_yu_y(1)) * V(1,i) * VC(1,j);
}

// （拾）计算 L_Ca2_k2 的函数
std::complex<double> L_Ca2_k2 (int i, int j){
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(j) * set_yu_y(2) * set_yu_y(2)) * V(2,i) * VC(2,j);
}

// （拾壹）计算 L_Ca2_k3 的函数
std::complex<double> L_Ca2_k3 (int i, int j){
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(j) * set_yu_y(3) * set_yu_y(3)) * V(3,i) * VC(3,j);
}

// // （）计算 L_Ca2 总和
// std::complex<double> L_Ca2 (int i, int j){
//     return L_Ca2_k1(i,j)+L_Ca2_k2(i,j)+L_Ca2_k3(i,j);
// }

// 计算 L_Cb 的函数
// 计算 L_Cb0 = L_Cb0_k1+L_Cb0_k2+L_Cb0_k3 的函数
// （拾贰）计算 L_Cb0_k1 的函数
std::complex<double> L_Cb0_k1 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Mu_w(1) * set_yd_z(j) * set_yu_y(1) + set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Mu_w(1) * set_Mu_w(1) * set_yd_z(j) - 2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

// （拾叁）计算 L_Cb0_k2 的函数
std::complex<double> L_Cb0_k2 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) *vacuum)) * (-set_Mu_w(2) * set_yd_z(j) * set_yu_y(2) + set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Mu_w(2) * set_Mu_w(2) * set_yd_z(j) - 2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

// （拾肆）计算 L_Cb0_k3 的函数
std::complex<double> L_Cb0_k3 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Mu_w(3) * set_yd_z(j) * set_yu_y(3) + set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Mu_w(3) * set_Mu_w(3) * set_yd_z(j) - 2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 L_Cb0 总和
// std::complex<double> L_Cb0 (int i, int j){
//     return L_Cb0_k1(i,j)+L_Cb0_k2(i,j)+L_Cb0_k3(i,j);
// }

// （拾伍）计算 L_Cb1_k1 的函数
std::complex<double> L_Cb1_k1 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yd_z(j) - set_Md_x(i) * set_Md_x(i) * set_yd_z(i) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + 2.0 * set_Md_x(j) * set_Mu_w(1) * set_yu_y(1))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

// （拾陆）计算 L_Cb1_k2 的函数
std::complex<double> L_Cb1_k2 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yd_z(j) - set_Md_x(i) * set_Md_x(i) * set_yd_z(i) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + 2.0 * set_Md_x(j) * set_Mu_w(2) * set_yu_y(2))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

// （拾柒）计算 L_Cb1_k3 的函数
std::complex<double> L_Cb1_k3 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yd_z(j) - set_Md_x(i) * set_Md_x(i) * set_yd_z(i) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + 2.0 * set_Md_x(j) * set_Mu_w(3) * set_yu_y(3))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 L_Cb1 总和
// std::complex<double> L_Cb1 (int i, int j){
//     return L_Cb1_k1(i,j)+L_Cb1_k2(i,j)+L_Cb1_k3(i,j);
// }

// （拾捌）计算 L_Cb2_k1 的函数
std::complex<double> L_Cb2_k1 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

// （拾玖）计算 L_Cb2_k2 的函数
std::complex<double> L_Cb2_k2 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

// （貳拾）计算 L_Cb2_k3 的函数
std::complex<double> L_Cb2_k3 (int i, int j){
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 L_Cb2 总和
// std::complex<double> L_Cb2 (int i, int j){
//     return L_Cb2_k1(i,j)+L_Cb2_k2(i,j)+L_Cb2_k3(i,j);
// }

// 计算 L_Cc 的函数
// 计算 L_Cc0 = L_Cc0_k1+L_Cc0_k2+L_Cc0_k3 的函数
// （貳拾壹）计算 L_Cc0_k1 的函数
std::complex<double> L_Cc0_k1 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(1) * set_yd_z(j) * set_yu_y(1) + set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
       * V(1,i) * VC(1,j); 
}

// （貳拾貳）计算 L_Cc0_k2 的函数
std::complex<double> L_Cc0_k2 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(2) * set_yd_z(j) * set_yu_y(2) + set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
       * V(2,i) * VC(2,j); 
}

// （貳拾叁）计算 L_Cc0_k3 的函数
std::complex<double> L_Cc0_k3 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(3) * set_yd_z(j) * set_yu_y(3) + set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
       * V(3,i) * VC(3,j); 
}

// // （）计算 L_Cc0 总和
// std::complex<double> L_Cc0 (int i, int j){
//     return L_Cc0_k1(i,j)+L_Cc0_k2(i,j)+L_Cc0_k3(i,j);
// }

// （貳拾肆）计算 L_Cc1_k1 的函数
std::complex<double> L_Cc1_k1 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(i) * set_yd_z(i) * set_yd_z(j) + set_Md_x(j) * set_yu_y(1) * set_yu_y(1)))\
         * V(1,i) * VC(1,j);
}

// （貳拾伍）计算 L_Cc1_k2 的函数
std::complex<double> L_Cc1_k2 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(i) * set_yd_z(i) * set_yd_z(j) + set_Md_x(j) * set_yu_y(2) * set_yu_y(2)))\
         * V(2,i) * VC(2,j);
}

// （貳拾陆）计算 L_Cc1_k3 的函数
std::complex<double> L_Cc1_k3 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(i) * set_yd_z(i) * set_yd_z(j) + set_Md_x(j) * set_yu_y(3) * set_yu_y(3)))\
         * V(3,i) * VC(3,j);
}

// // （）计算 L_Cc1 总和
// std::complex<double> L_Cc1 (int i, int j){
//     return L_Cc1_k1(i,j)+L_Cc1_k2(i,j)+L_Cc1_k3(i,j);
// }

// （貳拾柒）计算 L_Cc2_k1 的函数
std::complex<double> L_Cc2_k1 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

// （貳拾捌）计算 L_Cc2_k2 的函数
std::complex<double> L_Cc2_k2 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

// （貳拾玖）计算 L_Cc2_k3 的函数
std::complex<double> L_Cc2_k3 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 L_Cc2 总和
// std::complex<double> L_Cc2 (int i, int j){
//     return L_Cc2_k1(i,j)+L_Cc2_k2(i,j)+L_Cc2_k3(i,j);
// }

// 计算 L_Cd 的函数
// 计算 L_Cd0 = L_Cd0_k1+L_Cd0_k2+L_Cd0_k3 的函数
// （叁拾）计算 L_Cd0_k1 的函数
std::complex<double> L_Cd0_k1 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(j) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i) + set_Mu_w(1) * set_Mu_w(1)) * set_yd_z(j) + 2.0 * set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j); 
}

// （叁拾壹）计算 L_Cd0_k2 的函数
std::complex<double> L_Cd0_k2 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(j) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i) + set_Mu_w(2) * set_Mu_w(2)) * set_yd_z(j) + 2.0 * set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j); 
}

// （叁拾贰）计算 L_Cd0_k3 的函数
std::complex<double> L_Cd0_k3 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(j) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i) + set_Mu_w(3) * set_Mu_w(3)) * set_yd_z(j) + 2.0 * set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j); 
}

// // （）计算 L_Cd0 总和
// std::complex<double> L_Cd0 (int i, int j){
//     return L_Cd0_k1(i,j)+L_Cd0_k2(i,j)+L_Cd0_k3(i,j);
// }

// （叁拾叁）计算 L_Cd1_k1 的函数
std::complex<double> L_Cd1_k1 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i)) * set_yd_z(j) + set_Md_x(j) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)))\
        * V(1,i) * VC(1,j);
}

// （叁拾肆）计算 L_Cd1_k2 的函数
std::complex<double> L_Cd1_k2 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i)) * set_yd_z(j) + set_Md_x(j) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)))\
        * V(2,i) * VC(2,j);
}

// （叁拾伍）计算 L_Cd1_k3 的函数
std::complex<double> L_Cd1_k3 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i)) * set_yd_z(j) + set_Md_x(j) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 L_Cd1 总和
// std::complex<double> L_Cd1 (int i, int j){
//     return L_Cd1_k1(i,j)+L_Cd1_k2(i,j)+L_Cd1_k3(i,j);
// }

// （叁拾陆）计算 L_Cd2_k1 的函数
std::complex<double> L_Cd2_k1 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + set_Md_x(i) * set_Md_x(i)) * set_yd_z(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j)))\
        * V(1,i) * VC(1,j); 
}

// （叁拾柒）计算 L_Cd2_k2 的函数
std::complex<double> L_Cd2_k2 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + set_Md_x(i) * set_Md_x(i)) * set_yd_z(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j)))\
        * V(2,i) * VC(2,j);
}

// （叁拾捌）计算 L_Cd2_k3 的函数
std::complex<double> L_Cd2_k3 (int i, int j){
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + set_Md_x(i) * set_Md_x(i)) * set_yd_z(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j)))\
        * V(3,i) * VC(3,j);
}

// // （）计算 L_Cd2 总和
// std::complex<double> L_Cd2 (int i, int j){
//     return L_Cd2_k1(i,j)+L_Cd2_k2(i,j)+L_Cd2_k3(i,j);
// }

//BC函数，B函数直接导入，C函数先写成函数再导入
// F_Ca_0
std::complex<double> F_Ca_01 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Ca_02 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Ca_03 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Ca_1
std::complex<double> F_Ca_11 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Ca_12 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Ca_13 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Ca_2
std::complex<double> F_Ca_21 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Ca_22 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Ca_23 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_0
std::complex<double> F_Cb_01 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cb_02 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cb_03 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_1
std::complex<double> F_Cb_11 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cb_12 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cb_13 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_2
std::complex<double> F_Cb_21 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cb_22 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cb_23 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_0
std::complex<double> F_Cc_01 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cc_02 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cc_03 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_1
std::complex<double> F_Cc_11 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cc_12 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cc_13 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_2
std::complex<double> F_Cc_21 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cc_22 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cc_23 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_0
std::complex<double> F_Cd_01 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cd_02 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cd_03 (int i, int j){
    return C0i(0,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_1
std::complex<double> F_Cd_11 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cd_12 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cd_13 (int i, int j){
    return C0i(1,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_2
std::complex<double> F_Cd_21 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cd_22 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cd_23 (int i, int j){
    return C0i(2,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}

int main() {
    int i = 1; // 为i和j赋值
    int j = 2;
    //导入文件
    // Calculation_test("/home/test.csv");
    //导入B，C函数
    ltini();
    ifstream 
    while () {
      MH5 = xxxxx
      MH3 =xxxx 
    std::complex<double> F_Ba = B0(MH5*MH5,MH3*MH3,MW*MW);
    std::complex<double> F_Bb = B0(MH5*MH5,MW*MW,MW*MW);
    std::complex<double> F_Ca_01_value = F_Ca_01(i,j);
    std::complex<double> F_Ca_02_value = F_Ca_02(i,j);
    std::complex<double> F_Ca_03_value = F_Ca_03(i,j);
    std::complex<double> F_Ca_11_value = F_Ca_11(i,j);
    std::complex<double> F_Ca_12_value = F_Ca_12(i,j);
    std::complex<double> F_Ca_13_value = F_Ca_13(i,j);
    std::complex<double> F_Ca_21_value = F_Ca_21(i,j);
    std::complex<double> F_Ca_22_value = F_Ca_22(i,j);
    std::complex<double> F_Ca_23_value = F_Ca_23(i,j);
    std::complex<double> F_Cb_01_value = F_Cb_01(i,j);
    std::complex<double> F_Cb_02_value = F_Cb_02(i,j);
    std::complex<double> F_Cb_03_value = F_Cb_03(i,j);
    std::complex<double> F_Cb_11_value = F_Cb_11(i,j);
    std::complex<double> F_Cb_12_value = F_Cb_12(i,j);
    std::complex<double> F_Cb_13_value = F_Cb_13(i,j);
    std::complex<double> F_Cb_21_value = F_Cb_21(i,j);
    std::complex<double> F_Cb_22_value = F_Cb_22(i,j);
    std::complex<double> F_Cb_23_value = F_Cb_23(i,j);
    std::complex<double> F_Cc_01_value = F_Cc_01(i,j);
    std::complex<double> F_Cc_02_value = F_Cc_02(i,j);
    std::complex<double> F_Cc_03_value = F_Cc_03(i,j);
    std::complex<double> F_Cc_11_value = F_Cc_11(i,j);
    std::complex<double> F_Cc_12_value = F_Cc_12(i,j);
    std::complex<double> F_Cc_13_value = F_Cc_13(i,j);
    std::complex<double> F_Cc_21_value = F_Cc_21(i,j);
    std::complex<double> F_Cc_22_value = F_Cc_22(i,j);
    std::complex<double> F_Cc_23_value = F_Cc_23(i,j);
    std::complex<double> F_Cd_01_value = F_Cd_01(i,j);
    std::complex<double> F_Cd_02_value = F_Cd_02(i,j);
    std::complex<double> F_Cd_03_value = F_Cd_03(i,j);
    std::complex<double> F_Cd_11_value = F_Cd_11(i,j);
    std::complex<double> F_Cd_12_value = F_Cd_12(i,j);
    std::complex<double> F_Cd_13_value = F_Cd_13(i,j);
    std::complex<double> F_Cd_21_value = F_Cd_21(i,j);
    std::complex<double> F_Cd_22_value = F_Cd_22(i,j);
    std::complex<double> F_Cd_23_value = F_Cd_23(i,j);
    // 计算 RR_Ba
    std::complex<double> RR_Ba_value = RR_Ba(i, j);
    // 计算 RR_Bb
    std::complex<double> RR_Bb_value = RR_Bb(i, j);
    // 计算 R_Ca0
    std::complex<double> R_Ca0_k1_value = R_Ca0_k1(i, j);
    std::complex<double> R_Ca0_k2_value = R_Ca0_k2(i, j);
    std::complex<double> R_Ca0_k3_value = R_Ca0_k3(i, j);
    // 计算 R_Ca1
    std::complex<double> R_Ca1_k1_value = R_Ca1_k1(i, j);
    std::complex<double> R_Ca1_k2_value = R_Ca1_k2(i, j);
    std::complex<double> R_Ca1_k3_value = R_Ca1_k3(i, j);
    // 计算 R_Ca2
    std::complex<double> R_Ca2_k1_value = R_Ca2_k1(i, j);
    std::complex<double> R_Ca2_k2_value = R_Ca2_k2(i, j);
    std::complex<double> R_Ca2_k3_value = R_Ca2_k3(i, j);
    // 计算 R_Cb0
    std::complex<double> R_Cb0_k1_value = R_Cb0_k1(i, j);
    std::complex<double> R_Cb0_k2_value = R_Cb0_k2(i, j);
    std::complex<double> R_Cb0_k3_value = R_Cb0_k3(i, j);
    // 计算 R_Cb1
    std::complex<double> R_Cb1_k1_value = R_Cb1_k1(i, j);
    std::complex<double> R_Cb1_k2_value = R_Cb1_k2(i, j);
    std::complex<double> R_Cb1_k3_value = R_Cb1_k3(i, j);
    // 计算 R_Cb2
    std::complex<double> R_Cb2_k1_value = R_Cb2_k1(i, j);
    std::complex<double> R_Cb2_k2_value = R_Cb2_k2(i, j);
    std::complex<double> R_Cb2_k3_value = R_Cb2_k3(i, j);
    // 计算 R_Cc0
    std::complex<double> R_Cc0_k1_value = R_Cc0_k1(i, j);
    std::complex<double> R_Cc0_k2_value = R_Cc0_k2(i, j);
    std::complex<double> R_Cc0_k3_value = R_Cc0_k3(i, j);
    // 计算 R_Cc1
    std::complex<double> R_Cc1_k1_value = R_Cc1_k1(i, j);
    std::complex<double> R_Cc1_k2_value = R_Cc1_k2(i, j);
    std::complex<double> R_Cc1_k3_value = R_Cc1_k3(i, j);
    // 计算 R_Cc2
    std::complex<double> R_Cc2_k1_value = R_Cc2_k1(i, j);
    std::complex<double> R_Cc2_k2_value = R_Cc2_k2(i, j);
    std::complex<double> R_Cc2_k3_value = R_Cc2_k3(i, j);
    // 计算 R_Cd0
    std::complex<double> R_Cd0_k1_value = R_Cd0_k1(i, j);
    std::complex<double> R_Cd0_k2_value = R_Cd0_k2(i, j);
    std::complex<double> R_Cd0_k3_value = R_Cd0_k3(i, j);
    // 计算 R_Cd1
    std::complex<double> R_Cd1_k1_value = R_Cd1_k1(i, j);
    std::complex<double> R_Cd1_k2_value = R_Cd1_k2(i, j);
    std::complex<double> R_Cd1_k3_value = R_Cd1_k3(i, j);
    // 计算 R_Cd2
    std::complex<double> R_Cd2_k1_value = R_Cd2_k1(i, j);
    std::complex<double> R_Cd2_k2_value = R_Cd2_k2(i, j);
    std::complex<double> R_Cd2_k3_value = R_Cd2_k3(i, j);
    // 计算 LL_Ba
    std::complex<double> LL_Ba_value = LL_Ba(i, j);
    // 计算 LL_Bb
    std::complex<double> LL_Bb_value = LL_Bb(i, j);
    // 计算 L_Ca0
    std::complex<double> L_Ca0_k1_value = L_Ca0_k1(i, j);
    std::complex<double> L_Ca0_k2_value = L_Ca0_k2(i, j);
    std::complex<double> L_Ca0_k3_value = L_Ca0_k3(i, j);
    // 计算 L_Ca1
    std::complex<double> L_Ca1_k1_value = L_Ca1_k1(i, j);
    std::complex<double> L_Ca1_k2_value = L_Ca1_k2(i, j);
    std::complex<double> L_Ca1_k3_value = L_Ca1_k3(i, j);
    // 计算 L_Ca2
    std::complex<double> L_Ca2_k1_value = L_Ca2_k1(i, j);
    std::complex<double> L_Ca2_k2_value = L_Ca2_k2(i, j);
    std::complex<double> L_Ca2_k3_value = L_Ca2_k3(i, j);
    // 计算 L_Cb0
    std::complex<double> L_Cb0_k1_value = L_Cb0_k1(i, j);
    std::complex<double> L_Cb0_k2_value = L_Cb0_k2(i, j);
    std::complex<double> L_Cb0_k3_value = L_Cb0_k3(i, j);
    // 计算 L_Cb1
    std::complex<double> L_Cb1_k1_value = L_Cb1_k1(i, j);
    std::complex<double> L_Cb1_k2_value = L_Cb1_k2(i, j);
    std::complex<double> L_Cb1_k3_value = L_Cb1_k3(i, j);
    // 计算 L_Cb2
    std::complex<double> L_Cb2_k1_value = L_Cb2_k1(i, j);
    std::complex<double> L_Cb2_k2_value = L_Cb2_k2(i, j);
    std::complex<double> L_Cb2_k3_value = L_Cb2_k3(i, j);
    // 计算 L_Cc0
    std::complex<double> L_Cc0_k1_value = L_Cc0_k1(i, j);
    std::complex<double> L_Cc0_k2_value = L_Cc0_k2(i, j);
    std::complex<double> L_Cc0_k3_value = L_Cc0_k3(i, j);
    // 计算 L_Cc1
    std::complex<double> L_Cc1_k1_value = L_Cc1_k1(i, j);
    std::complex<double> L_Cc1_k2_value = L_Cc1_k2(i, j);
    std::complex<double> L_Cc1_k3_value = L_Cc1_k3(i, j);
    // 计算 L_Cc2
    std::complex<double> L_Cc2_k1_value = L_Cc2_k1(i, j);
    std::complex<double> L_Cc2_k2_value = L_Cc2_k2(i, j);
    std::complex<double> L_Cc2_k3_value = L_Cc2_k3(i, j);
    // 计算 L_Cd0
    std::complex<double> L_Cd0_k1_value = L_Cd0_k1(i, j);
    std::complex<double> L_Cd0_k2_value = L_Cd0_k2(i, j);
    std::complex<double> L_Cd0_k3_value = L_Cd0_k3(i, j);
    // 计算 L_Cd1
    std::complex<double> L_Cd1_k1_value = L_Cd1_k1(i, j);
    std::complex<double> L_Cd1_k2_value = L_Cd1_k2(i, j);
    std::complex<double> L_Cd1_k3_value = L_Cd1_k3(i, j);
    // 计算 L_Cd2
    std::complex<double> L_Cd2_k1_value = L_Cd2_k1(i, j);
    std::complex<double> L_Cd2_k2_value = L_Cd2_k2(i, j);
    std::complex<double> L_Cd2_k3_value = L_Cd2_k3(i, j);

    // 计算xi_ij_R(i,j)
    std::complex<double> xi_ij_R = RR_Ba_value * F_Ba + RR_Bb_value * F_Bb\
                                + R_Ca0_k1_value * F_Ca_01_value + R_Ca0_k2_value * F_Ca_02_value + R_Ca0_k3_value * F_Ca_03_value\
                                + R_Ca1_k1_value * F_Ca_11_value + R_Ca1_k2_value * F_Ca_12_value + R_Ca1_k3_value * F_Ca_13_value\
                                + R_Ca2_k1_value * F_Ca_21_value + R_Ca2_k2_value * F_Ca_22_value + R_Ca2_k3_value * F_Ca_23_value\
                                + R_Cb0_k1_value * F_Cb_01_value + R_Cb0_k2_value * F_Cb_02_value + R_Cb0_k3_value * F_Cb_03_value\
                                + R_Cb1_k1_value * F_Cb_11_value + R_Cb1_k2_value * F_Cb_12_value + R_Cb1_k3_value * F_Cb_13_value\
                                + R_Cb2_k1_value * F_Cb_21_value + R_Cb2_k2_value * F_Cb_22_value + R_Cb2_k3_value * F_Cb_23_value\
                                + R_Cc0_k1_value * F_Cc_01_value + R_Cc0_k2_value * F_Cc_02_value + R_Cc0_k3_value * F_Cc_03_value\
                                + R_Cc1_k1_value * F_Cc_11_value + R_Cc1_k2_value * F_Cc_12_value + R_Cc1_k3_value * F_Cc_13_value\
                                + R_Cc2_k1_value * F_Cc_21_value + R_Cc2_k2_value * F_Cc_22_value + R_Cc2_k3_value * F_Cc_23_value\
                                + R_Cd0_k1_value * F_Cd_01_value + R_Cd0_k2_value * F_Cd_02_value + R_Cd0_k3_value * F_Cd_03_value\
                                + R_Cd1_k1_value * F_Cd_11_value + R_Cd1_k2_value * F_Cd_12_value + R_Cd1_k3_value * F_Cd_13_value\
                                + R_Cd2_k1_value * F_Cd_21_value + R_Cd2_k2_value * F_Cd_22_value + R_Cd2_k3_value * F_Cd_23_value;

    //计算xi_ij_L(i,j)
    std::complex<double> xi_ij_L = LL_Ba_value * F_Ba + LL_Bb_value * F_Bb\
                                + L_Ca0_k1_value * F_Ca_01_value + L_Ca0_k2_value * F_Ca_02_value + L_Ca0_k3_value * F_Ca_03_value\
                                + L_Ca1_k1_value * F_Ca_11_value + L_Ca1_k2_value * F_Ca_12_value + L_Ca1_k3_value * F_Ca_13_value\
                                + L_Ca2_k1_value * F_Ca_21_value + L_Ca2_k2_value * F_Ca_22_value + L_Ca2_k3_value * F_Ca_23_value\
                                + L_Cb0_k1_value * F_Cb_01_value + L_Cb0_k2_value * F_Cb_02_value + L_Cb0_k3_value * F_Cb_03_value\
                                + L_Cb1_k1_value * F_Cb_11_value + L_Cb1_k2_value * F_Cb_12_value + L_Cb1_k3_value * F_Cb_13_value\
                                + L_Cb2_k1_value * F_Cb_21_value + L_Cb2_k2_value * F_Cb_22_value + L_Cb2_k3_value * F_Cb_23_value\
                                + L_Cc0_k1_value * F_Cc_01_value + L_Cc0_k2_value * F_Cc_02_value + L_Cc0_k3_value * F_Cc_03_value\
                                + L_Cc1_k1_value * F_Cc_11_value + L_Cc1_k2_value * F_Cc_12_value + L_Cc1_k3_value * F_Cc_13_value\
                                + L_Cc2_k1_value * F_Cc_21_value + L_Cc2_k2_value * F_Cc_22_value + L_Cc2_k3_value * F_Cc_23_value\
                                + L_Cd0_k1_value * F_Cd_01_value + L_Cd0_k2_value * F_Cd_02_value + L_Cd0_k3_value * F_Cd_03_value\
                                + L_Cd1_k1_value * F_Cd_11_value + L_Cd1_k2_value * F_Cd_12_value + L_Cd1_k3_value * F_Cd_13_value\
                                + L_Cd2_k1_value * F_Cd_21_value + L_Cd2_k2_value * F_Cd_22_value + L_Cd2_k3_value * F_Cd_23_value; 
    // 输出所有变量和计算结果
    std::cout<<xi_ij_L<<",  "<<xi_ij_R<<std::endl;
    }
    // outputFile << alpha << "," << SW << "," << vacuum << "," << thetaH << "," << M1 << "," << M2 << "," << MH5 << "," << MH3 << "," << xi_ij_R << "," << xi_ij_L << std::endl;
    
    // file.close();
    // outputFile.close();
    ltexi();
}