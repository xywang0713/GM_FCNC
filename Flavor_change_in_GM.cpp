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


// 定义全局变量
std::complex<double> alpha(1/137, 0.0), SW(0.2223, 0.0), vacuum(246.0, 0.0);
double MW = 80.399;
double Mb = 4.18;
double Ms = 93.4e-3;
double Md = 4.67e-3;
double Mu = 2.16e-3;
double Mc = 1.27;
double Mk = 0.493667; 
double Mpai=0.13957039;
double Meta=0.547862;
double Mmu =0.1056583755;
double MD=1.86966;
double MDs=1.96835;
double MB=5.27934;
double Gamma_k=5.3167366721e-17;
double Gamma_eta=1.31e-6; 
double GF=1.166364e-5;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ;
double EL =;


//写class
class GM_model{
public:
    double MH5;
    double MH3;
    double sH;
    double thetaH;
    double M1;
    double M2;


    void read_Data(){                           // 定义一个名为read_Data的void类型的函数
        std::ifstream file("home/test.csv")     // 用ifstream类打开home/test.csv
        if (file.is_open()){                    // 如果文件成功打开       
        std::string line;                       // 定义一个名为line的string类型变量
        int lineCount = 0;                      // 定义一个名为lineCount的int类型变量并赋值为0
        while (std::getline(file,line)){        // 当从file中读取一行内容存入line中时
            if (lineCount == 1){                // 如果行数为1
                std::istringstream iss(line);   // 使用istringstream类根据line创建一个名为iss的对象
                std:vector<double> values;      // 定义一个名为values的vector类型变量
                double value;                   // 定义一个名为value的double类型变量
                while (iss >> value) {          // 当从iss中读取值并成功时
                    values.push_back(value);    // 将value添加到values的末尾
                }
                MH5 = values[0];
                sH = values[1];
                thetaH = std::asin(sH); 
                M1 = values[6];
                M2 = values[7];
                MH3 = values[9];                // 将values中的值赋给变量
            }
            lineCount++;                        // 行数加1
        }
        file.close();
        }       
    }


    // 计算 LAMBDA_3()  函数
    std::complex<double> LAMBDA_3() {
         return (std::cos(thetaH) * std::cos(thetaH) * (-3.0 * std::sin(thetaH) * MH3 * MH3 + std::sqrt(2.0) * M1 * vacuum) + std::sin(thetaH) * (MH5 * MH5 - 3.0 * std::sqrt(2.0) * M2 * std::sin(thetaH) * vacuum)) / (vacuum * std::sin(thetaH) * std::sin(thetaH) * std::sin(thetaH));
    }

    // 计算 LAMBDA_5() 函数
    std::complex<double> LAMBDA_5() {
         return (2.0 * std::sin(thetaH) * MH3 * MH3 - std::sqrt(2.0) * M1 * vacuum) / (vacuum * std::sin(thetaH));
    }


//Potential里的参数
public:
    double LAMBDA_1;
    double LAMBDA_2;
    double LAMBDA_3();
    double LAMBDA_4;
    double LAMBDA_5();
    double MU12;
    double MU22;
    double M1;
    double M2;
}

//计算k—>pai H5的 p_H5^0
std::complex<double> P_0_PHI(GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return sqrt((Mk * Mk + Mpai * Mpai - MH5 * MH5) * (Mk * Mk + Mpai * Mpai - MH5 * MH5) / (4.0 * Mk * Mk) - Mpai * Mpai);
}

// 计算 yd(z) 的函数
std::complex<double> set_yd_z(int z, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

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
std::complex<double> set_yu_y(int y, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
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
double set_Md_x(int x, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
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
double set_Mu_w(int w, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
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
 std::complex<double> VCKM(int u, int t) {
     if (u == 1 && t == 1) {
         return 0.974352;
     } else if (u == 1 && t == 2) {
         return 0.224998;
     } else if (u == 1 && t == 3) {
         return std::complex<double>(0.0015275, -0.00335899);
     } else if (u == 2 && t == 1) {
         return std::complex<double>(-0.224865, -0.000136871);
     } else if (u == 2 && t == 2) {
         return std::complex<double>(0.973492, -0.0000316065 );
     } else if (u == 2 && t == 3) {
         return 0.0418197;
     } else if (u == 3 && t == 1 ) {
         return std::complex<double>(0.00792247, -0.00327); 
     } else if (u == 3 && t == 2 ) {
         return std::complex<double>(-0.0410911, -0.000755113);  
     } else if (u == 3 && t == 3 ) { 
         return 0.999118;
     } else {
         std::cerr << "Invalid indices u, t" << std::endl;
         return 0;
     }
 }

// 计算 CKMC矩阵VC[r,s] 的函数
 std::complex<double> VCKMC(int r, int s) {
     if (r == 1 && s == 1) {
         return 0.974352;
     } else if (r == 1 && s == 2) {
         return 0.224998;
     } else if (r == 1 && s == 3) {
         return std::complex<double>(0.0015275 , 0.00335899);
     } else if (r == 2 && s == 1) {
         return std::complex<double>(-0.224865 , 0.000136871);
     } else if (r == 2 && s == 2) {
         return std::complex<double>(0.973492 , 0.0000316065);
     } else if (r == 2 && s == 3) {
         return 0.0418197;
     } else if (r == 3 && s == 1 ) {
         return std::complex<double>(0.00792247 , 0.00327); 
     } else if (r == 3 && s == 2 ) {
         return std::complex<double>(-0.0410911 , 0.000755113); 
     } else if (r == 3 && s == 3 ) { 
         return 0.999118;
     } else {
         std::cerr << "Invalid indices r, s" << std::endl;
         return 0;
     }
 }


// P_R右手
// （一）计算 RR_Ba 的函数
std::complex<double> RR_Ba (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(i) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(1,i) * VCKMC(1,j) + VCKM(2,i) * VCKMC(2,j) + VCKM(3,i) * VCKMC(3,j));
}

// （二）计算 RR_Bb 的函数
std::complex<double> RR_Bb (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (-alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(i) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(1,i) * VCKMC(1,j) + VCKM(2,i) * VCKMC(2,j) + VCKM(3,i) * VCKMC(3,j));
}

// 计算 R_Ca 的函数
//  (三) 计算 R_Ca0_k1 的函数
std::complex<double> R_Ca0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1)) * VCKM(1,i) * VCKMC(1,j));
}

//  (四) 计算 R_Ca0_k2 的函数
std::complex<double> R_Ca0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2)) * VCKM(2,i) * VCKMC(2,j));
}

//  （五） 计算 R_Ca0_k3 的函数
std::complex<double> R_Ca0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3)) * VCKM(3,i) * VCKMC(3 ,j));
}

// （六）计算 R_Ca1_k1 的函数
std::complex<double> R_Ca1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1) * set_yu_y(1)) * VCKM(1,i) * VCKMC(1,j));
}

// （七）计算 R_Ca1_k2 的函数
std::complex<double> R_Ca1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2) * set_yu_y(2)) * VCKM(2,i) * VCKMC(2,j));
}

// （八）计算 R_Ca1_k3 的函数
std::complex<double> R_Ca1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3) * set_yu_y(3)) * VCKM(3,i) * VCKMC(3,j));
}

// （九）计算 R_Ca2_k1 的函数
std::complex<double> R_Ca2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *VCKM(1,i) * VCKMC(1,j));
}

// （十）计算 R_Ca2_k2 的函数
std::complex<double> R_Ca2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *VCKM(2,i) * VCKMC(2,j));
}

// （十一）计算 R_Ca2_k3 的函数
std::complex<double> R_Ca2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) *VCKM(3,i) * VCKMC(3,j));
}

// 计算 R_Cb 的函数
// 计算 R_Cb0 的函数
//  （十二）计算 R_Cb0_k1 的函数
std::complex<double> R_Cb0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1))\
    -2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i) + set_Mu_w(1) *set_Mu_w(1) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
           * VCKM(1,i) * VCKMC(1,j);
}

//  （十三）计算 R_Cb0_k2 的函数
std::complex<double> R_Cb0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2))\
    -2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i) + set_Mu_w(2) *set_Mu_w(2) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
           * VCKM(2,i) * VCKMC(2,j);
}

//  （十四）计算 R_Cb0_k3 的函数
std::complex<double> R_Cb0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3))\
    -2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i) + set_Mu_w(3) *set_Mu_w(3) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
           * VCKM(3,i) * VCKMC(3,j);
}

// 计算 R_Cb1 的函数
//  （十五）计算 R_Cb1_k1 的函数
std::complex<double> R_Cb1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Md_x(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) - 2.0 * set_Md_x(i) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1) * set_yu_y(1)))\
           * VCKM(1,i) * VCKMC(1,j);
} 

//  （十六）计算 R_Cb1_k2 的函数
std::complex<double> R_Cb1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    *(-16.0 * alpha * alpha * set_Md_x(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) - 2.0 * set_Md_x(i) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2) * set_yu_y(2)))\
          * VCKM(2,i) * VCKMC(2,j); 
}

//  （十七）计算 R_Cb1_k3 的函数
std::complex<double> R_Cb1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Md_x(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) - 2.0 * set_Md_x(i) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3) * set_yu_y(3)))\
          * VCKM(3,i) * VCKMC(3,j); 
}

// 计算 R_Cb2 = R_Cb2_k1+R_Cb2_k2+R_Cb2_k3 的函数
//  （十八）计算 R_Cb2_k1 的函数
std::complex<double> R_Cb2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i)))\
          * VCKM(1,i) * VCKMC(1,j);
}

//  （十九）计算 R_Cb2_k2 的函数
std::complex<double> R_Cb2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i)))\
          * VCKM(2,i) * VCKMC(2,j);
}

//  （二十）计算 R_Cb2_k3 的函数
std::complex<double> R_Cb2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Md_x(j) * set_Md_x(i) * set_yd_z(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i)) * set_yd_z(i)))\
           * VCKM(3,i) * VCKMC(3,j); 
}

// 计算 R_Cc 的函数
// 计算 R_Cc0 = R_Cc0_k1+R_Cc0_k2+R_Cc0_k3 的函数
//  （二十一）计算 R_Cc0_k1 的函数
std::complex<double> R_Cc0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW *SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(i) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - 2.0 * set_Md_x(i) * set_Md_x(i) + set_Mu_w(1) * set_Mu_w(1)) * set_yd_z(i) + 2.0 * set_Md_x(i) * set_Mu_w(1) *set_yu_y(1)))\
           * VCKM(1,i) * VCKMC(1,j);
}

// （二十二）计算 R_Cc0_k2 的函数
std::complex<double> R_Cc0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW *SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(i) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - 2.0 * set_Md_x(i) * set_Md_x(i) + set_Mu_w(2) * set_Mu_w(2)) * set_yd_z(i) + 2.0 * set_Md_x(i) * set_Mu_w(2) *set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （二十三）计算 R_Cc0_k3 的函数
std::complex<double> R_Cc0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW *SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 +(-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(i) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - 2.0 * set_Md_x(i) * set_Md_x(i) + set_Mu_w(3) * set_Mu_w(3)) * set_yd_z(i) + 2.0 * set_Md_x(i) * set_Mu_w(3) *set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 R_Cc1 = R_Cc1_k1+R_Cc1_k2+R_Cc1_k3 的函数
//  （二十四）计算 R_Cc1_k1 的函数
std::complex<double> R_Cc1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(1)*set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

//  （二十五）计算 R_Cc1_k2 的函数
std::complex<double> R_Cc1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(2)*set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

//  （二十六）计算 R_Cc1_k3 的函数
std::complex<double> R_Cc1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5) * set_yd_z(i) + set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j) - set_Md_x(i) * set_yu_y(3)*set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 R_Cc2 = R_Cc2_k1+R_Cc2_k2+R_Cc2_k3 的函数
// （二十七）计算 R_Cc2_k1 的函数
std::complex<double> R_Cc2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j)))\
        * VCKM(1,i) * VCKMC(1,j);
}

//  （二十八）计算 R_Cc1_k2 的函数
std::complex<double> R_Cc2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j)))\
        * VCKM(2,i) * VCKMC(2,j);
}

//  （二十九）计算 R_Cc2_k3 的函数
std::complex<double> R_Cc2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32* std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Md_x(i) * set_Md_x(i))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() +LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yd_z(i) * set_yd_z(j)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 R_Cd0 = R_Cd0_k1+R_Cd0_k2+R_Cd0_k3 的函数
//  （三十）计算 R_Cd0_k1 的函数
std::complex<double> R_Cd0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(1) * set_yd_z(i) * set_yu_y(1) + set_Md_x(i) * set_yu_y(1) * set_yu_y(1))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

//  （三十一）计算 R_Cd0_k2 的函数
std::complex<double> R_Cd0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(2) * set_yd_z(i) * set_yu_y(2) + set_Md_x(i) * set_yu_y(2) *set_yu_y(2))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

//  （三十二）计算 R_Cd0_k3 的函数
std::complex<double> R_Cd0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(3) * set_yd_z(i) * set_yu_y(3) + set_Md_x(i) * set_yu_y(3) *set_yu_y(3))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 R_Cd1 = R_Cd1_k1+R_Cd1_k2+R_Cd1_k3 的函数
//  （三十三）计算 R_Cd1_k1 的函数
std::complex<double> R_Cd1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) *set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(j) * set_yd_z(i) * set_yd_z(j) + set_Md_x(i) * set_yu_y(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j); 
}

//  （三十四）计算 R_Cd1_k2 的函数
std::complex<double> R_Cd1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(j) * set_yd_z(i) * set_yd_z(j) + set_Md_x(i) * set_yu_y(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

//  （三十五）计算 R_Cd1_k3 的函数
std::complex<double> R_Cd1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(j) * set_yd_z(i) * set_yd_z(j) + set_Md_x(i) * set_yu_y(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 R_Cd2 = R_Cd2_k1+R_Cd2_k2+R_Cd2_k3 的函数
//  （三十六）计算 R_Cd2_k1 的函数
std::complex<double> R_Cd2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yu_y(1) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

//  （三十七）计算 R_Cd2_k2 的函数
std::complex<double> R_Cd2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yu_y(2) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

//  （三十八）计算 R_Cd2_k3 的函数
std::complex<double> R_Cd2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yu_y(3) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Mu_w(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// P_L左手
// （壹）计算 LL_Ba 的函数
std::complex<double> LL_Ba (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(j) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(1,i) * VCKMC(1,j) + VCKM(2,i) * VCKMC(2,j) + VCKM(3,i) * VCKMC(3,j));
}

// （贰）计算 LL_Bb 的函数
std::complex<double> LL_Bb (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (-alpha * std::cos(thetaH) * std::sin(thetaH) * set_yd_z(j) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(1,i) * VCKMC(1,j) + VCKM(2,i) * VCKMC(2,j) + VCKM(3,i) * VCKMC(3,j));
}

// 计算 L_Ca 的函数
// 计算 L_Ca0 = L_Ca0_k1+L_Ca0_k2+L_Ca0_k3 的函数
// （叁）计算 L_Ca0_k1 的函数
std::complex<double> L_Ca0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(1) * set_yd_z(j) * set_yu_y(1) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)) * VCKM(1,i) * VCKMC(1,j);
}

// （肆）计算 L_Ca0_k2 的函数
std::complex<double> L_Ca0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(2) * set_yd_z(j) * set_yu_y(2) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)) * VCKM(2,i) * VCKMC(2,j);
}

// （伍）计算 L_Ca0_k3 的函数
std::complex<double> L_Ca0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(3) * set_yd_z(j) * set_yu_y(3) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)) * VCKM(3,i) * VCKMC(3,j);
}

// （陆）计算 L_Ca1_k1 的函数
std::complex<double> L_Ca1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)) * VCKM(1,i) * VCKMC(1,j);
}

// （柒）计算 L_Ca1_k2 的函数
std::complex<double> L_Ca1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)) * VCKM(2,i) * VCKMC(2,j);
}

// （捌）计算 L_Ca1_k3 的函数
std::complex<double> L_Ca1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)) * VCKM(3,i) * VCKMC(3,j);
}

// （玖）计算 L_Ca2_k1 的函数
std::complex<double> L_Ca2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(j) * set_yu_y(1) * set_yu_y(1)) * VCKM(1,i) * VCKMC(1,j);
}

// （拾）计算 L_Ca2_k2 的函数
std::complex<double> L_Ca2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(j) * set_yu_y(2) * set_yu_y(2)) * VCKM(2,i) * VCKMC(2,j);
}

// （拾壹）计算 L_Ca2_k3 的函数
std::complex<double> L_Ca2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(j) * set_yu_y(3) * set_yu_y(3)) * VCKM(3,i) * VCKMC(3,j);
}


// 计算 L_Cb 的函数
// 计算 L_Cb0 = L_Cb0_k1+L_Cb0_k2+L_Cb0_k3 的函数
// （拾贰）计算 L_Cb0_k1 的函数
std::complex<double> L_Cb0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Mu_w(1) * set_yd_z(j) * set_yu_y(1) + set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Mu_w(1) * set_Mu_w(1) * set_yd_z(j) - 2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

// （拾叁）计算 L_Cb0_k2 的函数
std::complex<double> L_Cb0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) *vacuum)) * (-set_Mu_w(2) * set_yd_z(j) * set_yu_y(2) + set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Mu_w(2) * set_Mu_w(2) * set_yd_z(j) - 2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （拾肆）计算 L_Cb0_k3 的函数
std::complex<double> L_Cb0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Mu_w(3) * set_yd_z(j) * set_yu_y(3) + set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Mu_w(3) * set_Mu_w(3) * set_yd_z(j) - 2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// （拾伍）计算 L_Cb1_k1 的函数
std::complex<double> L_Cb1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yd_z(j) - set_Md_x(i) * set_Md_x(i) * set_yd_z(i) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + 2.0 * set_Md_x(j) * set_Mu_w(1) * set_yu_y(1))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

// （拾陆）计算 L_Cb1_k2 的函数
std::complex<double> L_Cb1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yd_z(j) - set_Md_x(i) * set_Md_x(i) * set_yd_z(i) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + 2.0 * set_Md_x(j) * set_Mu_w(2) * set_yu_y(2))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （拾柒）计算 L_Cb1_k3 的函数
std::complex<double> L_Cb1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yd_z(j) - set_Md_x(i) * set_Md_x(i) * set_yd_z(i) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + 2.0 * set_Md_x(j) * set_Mu_w(3) * set_yu_y(3))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// （拾捌）计算 L_Cb2_k1 的函数
std::complex<double> L_Cb2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

// （拾玖）计算 L_Cb2_k2 的函数
std::complex<double> L_Cb2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （貳拾）计算 L_Cb2_k3 的函数
std::complex<double> L_Cb2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Md_x(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Md_x(j) * set_Md_x(j) * set_yd_z(j) - set_Md_x(i) * set_Md_x(j) * set_yd_z(i) + set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 L_Cc 的函数
// 计算 L_Cc0 = L_Cc0_k1+L_Cc0_k2+L_Cc0_k3 的函数
// （貳拾壹）计算 L_Cc0_k1 的函数
std::complex<double> L_Cc0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(1) * set_yd_z(j) * set_yu_y(1) + set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
       * VCKM(1,i) * VCKMC(1,j); 
}

// （貳拾貳）计算 L_Cc0_k2 的函数
std::complex<double> L_Cc0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(2) * set_yd_z(j) * set_yu_y(2) + set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
       * VCKM(2,i) * VCKMC(2,j); 
}

// （貳拾叁）计算 L_Cc0_k3 的函数
std::complex<double> L_Cc0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(3) * set_yd_z(j) * set_yu_y(3) + set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
       * VCKM(3,i) * VCKMC(3,j); 
}

// （貳拾肆）计算 L_Cc1_k1 的函数
std::complex<double> L_Cc1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(i) * set_yd_z(i) * set_yd_z(j) + set_Md_x(j) * set_yu_y(1) * set_yu_y(1)))\
         * VCKM(1,i) * VCKMC(1,j);
}

// （貳拾伍）计算 L_Cc1_k2 的函数
std::complex<double> L_Cc1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(i) * set_yd_z(i) * set_yd_z(j) + set_Md_x(j) * set_yu_y(2) * set_yu_y(2)))\
         * VCKM(2,i) * VCKMC(2,j);
}

// （貳拾陆）计算 L_Cc1_k3 的函数
std::complex<double> L_Cc1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(i) * set_yd_z(i) * set_yd_z(j) + set_Md_x(j) * set_yu_y(3) * set_yu_y(3)))\
         * VCKM(3,i) * VCKMC(3,j);
}

// （貳拾柒）计算 L_Cc2_k1 的函数
std::complex<double> L_Cc2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(1) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

// （貳拾捌）计算 L_Cc2_k2 的函数
std::complex<double> L_Cc2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(2) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （貳拾玖）计算 L_Cc2_k3 的函数
std::complex<double> L_Cc2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(j) * set_yu_y(3) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Md_x(i) * set_Md_x(j) * set_yd_z(i) - set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// 计算 L_Cd 的函数
// 计算 L_Cd0 = L_Cd0_k1+L_Cd0_k2+L_Cd0_k3 的函数
// （叁拾）计算 L_Cd0_k1 的函数
std::complex<double> L_Cd0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(1) * set_yd_z(j) * set_yu_y(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i) + set_Mu_w(1) * set_Mu_w(1)) * set_yd_z(j) + 2.0 * set_Md_x(j) * set_Mu_w(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j); 
}

// （叁拾壹）计算 L_Cd0_k2 的函数
std::complex<double> L_Cd0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(2) * set_yd_z(j) * set_yu_y(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i) + set_Mu_w(2) * set_Mu_w(2)) * set_yd_z(j) + 2.0 * set_Md_x(j) * set_Mu_w(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j); 
}

// （叁拾贰）计算 L_Cd0_k3 的函数
std::complex<double> L_Cd0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Mu_w(3) * set_yd_z(j) * set_yu_y(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Md_x(j) * set_Md_x(j) + 2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i) + set_Mu_w(3) * set_Mu_w(3)) * set_yd_z(j) + 2.0 * set_Md_x(j) * set_Mu_w(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j); 
}

// （叁拾叁）计算 L_Cd1_k1 的函数
std::complex<double> L_Cd1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i)) * set_yd_z(j) + set_Md_x(j) * set_Mu_w(1) * set_yu_y(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(1) * set_yu_y(1)))\
        * VCKM(1,i) * VCKMC(1,j);
}

// （叁拾肆）计算 L_Cd1_k2 的函数
std::complex<double> L_Cd1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i)) * set_yd_z(j) + set_Md_x(j) * set_Mu_w(2) * set_yu_y(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(2) * set_yu_y(2)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （叁拾伍）计算 L_Cd1_k3 的函数
std::complex<double> L_Cd1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Md_x(i) * set_Md_x(i)) * set_yd_z(j) + set_Md_x(j) * set_Mu_w(3) * set_yu_y(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j) - set_Md_x(j) * set_yu_y(3) * set_yu_y(3)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// （叁拾陆）计算 L_Cd2_k1 的函数
std::complex<double> L_Cd2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + set_Md_x(i) * set_Md_x(i)) * set_yd_z(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j)))\
        * VCKM(1,i) * VCKMC(1,j); 
}

// （叁拾柒）计算 L_Cd2_k2 的函数
std::complex<double> L_Cd2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + set_Md_x(i) * set_Md_x(i)) * set_yd_z(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j)))\
        * VCKM(2,i) * VCKMC(2,j);
}

// （叁拾捌）计算 L_Cd2_k3 的函数
std::complex<double> L_Cd2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Md_x(j) * set_Md_x(j) - 2.0 * MH5 * MH5 + set_Md_x(i) * set_Md_x(i)) * set_yd_z(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Md_x(i) * set_yd_z(i) * set_yd_z(j)))\
        * VCKM(3,i) * VCKMC(3,j);
}

// F_Ca_0
std::complex<double> F_Ca_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Ca_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Ca_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Ca_1
std::complex<double> F_Ca_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Ca_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Ca_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Ca_2
std::complex<double> F_Ca_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Ca_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Ca_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_0
std::complex<double> F_Cb_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cb_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cb_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_1
std::complex<double> F_Cb_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cb_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cb_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_2
std::complex<double> F_Cb_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cb_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cb_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_0
std::complex<double> F_Cc_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cc_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cc_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_1
std::complex<double> F_Cc_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cc_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cc_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_2
std::complex<double> F_Cc_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cc_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cc_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(i)*set_Md_x(i),set_Md_x(j)*set_Md_x(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_0
std::complex<double> F_Cd_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cd_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cd_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_1
std::complex<double> F_Cd_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cd_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cd_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_2
std::complex<double> F_Cd_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(1)*set_Mu_w(1));
}
std::complex<double> F_Cd_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> F_Cd_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Md_x(j)*set_Md_x(j),set_Md_x(i)*set_Md_x(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}

// 计算 xi_ij_R(i,j) 函数
std::complex<double> xi_ij_R (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    //调用一大堆函数
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

    return RR_Ba_value * F_Ba + RR_Bb_value * F_Bb\
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
}

//计算xi_ij_L(i,j)的函数
std::complex<double> xi_ij_L (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    //调用一大堆函数
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

    return LL_Ba_value * F_Ba + LL_Bb_value * F_Bb\
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
}

//下面涉及到的质量Mb，Mc，
//B介子衰变产生函数Br(B -> Xs phi)/Br(B -> Xc e nu) 
//注意：其中最后abs(xi_ij_R (int i, int j, GM_model&m) / VCKM[c,b])中xi是xi_bs
std::complex<double> Br_B_Meson (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (12.0 * M_PI * M_PI * vacuum * vacuum / (Mb * Mb)) * (1.0 - (MH5 * MH5) / (Mb * Mb)) * (1.0 / (((1.0 - 8.0 * ((Mc * Mc) / (Mb * Mb)) + ((Mc * Mc ) / (Mb * Mb)) * ((Mc * Mc )/ (Mb * Mb)))) * (1.0 - ((Mc * Mc) / (Mb * Mb)) * ((Mc * Mc ) / (Mb * Mb))) - 12.0 * log(((Mc * Mc ) / (Mb * Mb)) * ((Mc * Mc ) / (Mb * Mb))))) * std::norm(xi_ij_R (3, 2, GM_model&m) / VCKM[2,3]) * std::norm(xi_ij_R (3, 2, GM_model&m) / VCKM[2,3]);
}

//K介子产生Ms，Mk，Mpai，Md,GF,Gamma_k,P_0_PHI,EL
std::complex<double> Br_Kaon (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return (1.0 / Gamma_k) * (2.0 * P_0_PHI / Mk) * ( std::norm(sqrt(GF) * 4.0 * cbrt(2.0) * (EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW)) * (7.0 * 3.1 * pow (10,-7) *(Mk* Mk +Mpai * Mpai-MH5 * MH5) / 18.0) + (xi_ij_R (1, 2, GM_model&m) * Ms *0.96 /(2.0*vacuum)) *((Mk* Mk -Mpai * Mpai) / (Ms - Md))) * std::norm(sqrt(GF) * 4.0 * cbrt(2.0) * (EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW)) * (7.0 * 3.1 * pow (10,-7) *(Mk* Mk +Mpai * Mpai-MH5 * MH5) / 18.0) + (xi_ij_R (1, 2, GM_model&m) * Ms *0.96 /(2.0*vacuum)) *((Mk* Mk -Mpai * Mpai) / (Ms - Md))) / (16.0 * M_PI * Mk))
}


//Semileptonic Decay of Mesons Mx,Mmu
//BR（X->munu)
//1.BR(pai+ -> munu)
std::complex<double> Br_Semileptonic_pai(int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((sqrt(2.0) * GF * Mpai * Mpai * Mpai * Mpai * std::norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW * MW)) * std::norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW * MW))) / (96.0 * M_PI *M_PI * Mmu * Mmu * (1.0 - (Mmu * Mmu) / (Mx *Mx)) * (1.0 - (Mmu * Mmu) / (Mpai * Mpai)))) * 0.999877 * (((1.0 - 8.0 * ((MH5 * MH5 ) / (Mpai * Mpai)) + ((MH5 * MH5 )/ (Mpai * Mpai)) * ((MH5 * MH5 )/ (Mpai * Mpai)))) * (1.0 - ((MH5 * MH5 ) / (Mpai * Mpai)) * ((MH5 * MH5 ) / (Mpai * Mpai))) - 12.0 * log(((MH5 * MH5 ) / (Mpai * Mpai)) * ((MH5 * MH5 ) / (Mpai * Mpai)))) * (7.0 / 9.0) *(7.0 / 9.0);
}

//2.BR(K+ -> munu)
std::complex<double> Br_Semileptonic_k(int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((sqrt(2.0) * GF * Mk * Mk * Mk * Mk * std:norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW * MW)) * std::norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW))) / (96.0 * M_PI *M_PI * Mmu * Mmu * (1.0 - (Mmu * Mmu) / (Mk *Mk)) * (1.0 - (Mmu * Mmu) / (Mk *Mk)))) * 0.6356 * (((1.0 - 8.0 * ((MH5 * MH5 ) / (Mk * Mk)) + ((MH5 * MH5 )/ (Mk * Mk)) * ((MH5 * MH5 )/ (Mk * Mk)))) * (1.0 - ((MH5 * MH5 ) / (Mk * Mk)) * ((MH5 * MH5 ) / (Mk * Mk))) - 12.0 * log(((MH5 * MH5 ) / (Mk * Mk)) * ((MH5 * MH5 ) / (Mk * Mk)))) * (7.0 / 9.0) *(7.0 / 9.0);
}

//3.BR(D+ -> munu)
std::complex<double> Br_Semileptonic_D(int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();star
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((sqrt(2.0) * GF * MD * MD * MD * MD * std::norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW * MW)) * std::norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW))) / (96.0 * M_PI *M_PI * Mmu * Mmu * (1.0 - (Mmu * Mmu) / (MD *MD)) * (1.0 - (Mmu * Mmu) / (MD *MD)))) * std::pow(3.74,-4) * (((1.0 - 8.0 * ((MH5 * MH5 ) / (MD * MD)) + ((MH5 * MH5 )/ (MD * MD)) * ((MH5 * MH5 )/ (MD * MD)))) * (1.0 - ((MH5 * MH5 ) / (MD * MD)) * ((MH5 * MH5 ) / (MD * MD))) - 12.0 * log(((MH5 * MH5 ) / (MD * MD)) * ((MH5 * MH5 ) / (MD * MD)))) * (7.0 / 9.0) *(7.0 / 9.0);
}

//4.BR(Ds+ -> munu)
std::complex<double> Br_Semileptonic_Ds(int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((sqrt(2.0) * GF * MDs * MDs *  MDs * MDs * std:norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW * MW)) * std:norm(EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW))) / (96.0 * M_PI *M_PI * Mmu * Mmu * (1.0 - (Mmu * Mmu) / ( MDs * MDs)) * (1.0 - (Mmu * Mmu) / ( MDs * MDs)))) * std::pow(5.43,-3) * (((1.0 - 8.0 * ((MH5 * MH5 ) / ( MDs * MDs)) + ((MH5 * MH5 )/ ( MDs * MDs)) * ((MH5 * MH5 )/ ( MDs * MDs)))) * (1.0 - ((MH5 * MH5 ) / ( MDs * MDs)) * ((MH5 * MH5 ) / ( MDs * MDs))) - 12.0 * log(((MH5 * MH5 ) / ( MDs * MDs)) * ((MH5 * MH5 ) / ( MDs * MDs)))) * (7.0 / 9.0) *(7.0 / 9.0);
}

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
        outputFile.close();
        std::cout << "Calculation has been written to outputresult.txt." << std::endl;
    } else {
        std::cerr << "There is somthing wrong to open file";
    }
    return 0;
    }
    ltexi();