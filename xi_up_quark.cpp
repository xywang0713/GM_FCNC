// P_R右手
// （一）计算 up_RR_Ba 的函数
std::complex<double> up_RR_Ba (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yu_y(i) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(j,1) * VCKMC(i,1) + VCKM(j,2) * VCKMC(i,2) + VCKM(j,3) * VCKMC(i,3));
}

// （二）计算 up_RR_Bb 的函数
std::complex<double> up_RR_Bb (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (-alpha * std::cos(thetaH) * std::sin(thetaH) * set_yu_y(i) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(j,1) * VCKMC(i,1) + VCKM(j,2) * VCKMC(i,2) + VCKM(j,3) * VCKMC(i,3));
}

// 计算 up_R_Ca 的函数
//  (三) 计算 up_R_Ca0_k1 的函数，因为涉及到Md（k）的求和，这里分别对k=1,2,3进行计算，命名后缀加上了k1,k2,k3
std::complex<double> up_R_Ca0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(1) * set_yu_y(i) * set_yd_z(1)) * VCKM(j,1) * VCKMC(i,1));
}

//  (四) 计算 up_R_Ca0_k2 的函数
std::complex<double> up_R_Ca0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * ((set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(2) * set_yu_y(i) * set_yd_z(2)) * VCKM(j,2) * VCKMC(i,2));
}

//  （五） 计算 up_R_Ca0_k3 的函数
std::complex<double> up_R_Ca0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * ((set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(3) * set_yu_y(i) * set_yd_z(3)) * VCKM(j,3) * VCKMC(i,3));
}

// （六）计算 up_R_Ca1_k1 的函数
std::complex<double> up_R_Ca1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * ((set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(1) * set_yd_z(1)) * VCKM(j,1) * VCKMC(i,1));
}

// （七）计算 up_R_Ca1_k2 的函数
std::complex<double> up_R_Ca1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(2) * set_yd_z(2)) * VCKM(j,2) * VCKMC(i,2));
}

// （八）计算 up_R_Ca1_k3 的函数
std::complex<double> up_R_Ca1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * ((set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(3) * set_yd_z(3)) * VCKM(j,3) * VCKMC(i,3));
}

// （九）计算 up_R_Ca2_k1 的函数
std::complex<double> up_R_Ca2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) *VCKM(j,1) * VCKMC(i,1));
}

// （十）计算 up_R_Ca2_k2 的函数
std::complex<double> up_R_Ca2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) *VCKM(j,2) * VCKMC(i,2));
}

// （十一）计算 up_R_Ca2_k3 的函数
std::complex<double> up_R_Ca2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) / (32.0 * std::sqrt(3.0) * M_PI *M_PI)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) *VCKM(j,3) * VCKMC(i,3));
}

// 计算 R_Cb 的函数
// 计算 R_Cb0 的函数
//  （十二）计算 up_R_Cb0_k1 的函数
std::complex<double> up_R_Cb0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(1) * set_yu_y(i) * set_yd_z(1))\
    -2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Mu(i) * set_Md_x(1) * set_yd_z(1) - (set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i) + set_Md_x(1) * set_Md_x(1) * set_yu_y(i) + set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j)))\
           * VCKM(j,1) * VCKMC(i,1);
}

//  （十三）计算 up_R_Cb0_k2 的函数
std::complex<double> up_R_Cb0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(2) * set_yu_y(i) * set_yd_z(2))\
    -2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Mu(i) * set_Md_x(2) * set_yd_z(2) - (set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i) + set_Md_x(2) * set_Md_x(2) * set_yu_y(i) + set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j)))\
           * VCKM(j,2) * VCKMC(i,2);
}

//  （十四）计算 up_R_Cb0_k3 的函数
std::complex<double> up_R_Cb0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(3) * set_yu_y(i) * set_yd_z(3))\
    -2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (set_Mu(i) * set_Md_x(3) * set_yd_z(3) - (set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i) + set_Md_x(3) * set_Md_x(3) * set_yu_y(i) + set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j)))\
           * VCKM(j,3) * VCKMC(i,3);
}

// 计算 up_R_Cb1 的函数
//  （十五）计算 up_R_Cb1_k1 的函数
std::complex<double> R_Cb1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Mu_w(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j) + (set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5) * set_yu_y(i) - 2.0 * set_Mu_w(i) * set_Md_x(1) * set_yd_z(1))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(1) * set_yd_z(1)))\
           * VCKM(j,1) * VCKMC(i,1);
} 

//  （十六）计算 up_R_Cb1_k2 的函数
std::complex<double> up_R_Cb1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Mu_w(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j) + (set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5) * set_yu_y(i) - 2.0 * set_Mu_w(i) * set_Md_x(2) * set_yd_z(2))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(2) * set_yd_z(2)))\
          * VCKM(j,2) * VCKMC(i,2); 
}

//  （十七）计算 up_R_Cb1_k3 的函数
std::complex<double> up_R_Cb1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW * SW))\
    * (-16.0 * alpha * alpha * set_Mu_w(i) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::cos(thetaH) * std::sin(thetaH) * SW * SW * (set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j) + (set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5) * set_yu_y(i) - 2.0 * set_Mu_w(i) * set_Md_x(3) * set_yd_z(3))\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(3) * set_yd_z(3)))\
          * VCKM(j,3) * VCKMC(i,3); 
}

// 计算 R_Cb2 = R_Cb2_k1+R_Cb2_k2+R_Cb2_k3 的函数
//  （十八）计算 up_R_Cb2_k1 的函数
std::complex<double> up_R_Cb2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i)))\
          * VCKM(j,1) * VCKMC(i,1);
}

//  （十九）计算 up_R_Cb2_k2 的函数
std::complex<double> up_R_Cb2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i)))\
          * VCKM(j,2) * VCKMC(i,2);
}

//  （二十）计算 up_R_Cb2_k3 的函数
std::complex<double> up_R_Cb2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * std::sin(thetaH) * (-set_Mu_w(j) * set_Mu_w(i) * set_yu_y(j) + (set_Md_x(j) *set_Md_x(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i)))\
           * VCKM(j,3) * VCKMC(i,3); 
}

// 计算 R_Cc 的函数
// 计算 R_Cc0 = R_Cc0_k1+R_Cc0_k2+R_Cc0_k3 的函数
//  （二十一）计算 up_R_Cc0_k1 的函数
std::complex<double> up_R_Cc0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(1) * set_yu_y(i) * set_yd_z(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Mu_w(j) * set_Mu_w(j) + 2.0 * MH5 * MH5 - 2.0 * set_Mu_w(i) * set_Mu_w(i) + set_Md_x(1) * set_Md_x(1)) * set_yu_y(i) + 2.0 * set_Mu_w(i) * set_Md_x(1) *set_yd_z(1)))\
           * VCKM(j,1) * VCKMC(i,1);
}

// （二十二）计算 up_R_Cc0_k2 的函数
std::complex<double> up_R_Cc0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(2) * set_yu_y(i) * set_yd_z(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Mu_w(j) * set_Mu_w(j) + 2.0 * MH5 * MH5 - 2.0 * set_Mu_w(i) * set_Mu_w(i) + set_Md_x(2) * set_Md_x(2)) * set_yu_y(i) + 2.0 * set_Mu_w(i) * set_Md_x(2) *set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （二十三）计算 up_R_Cc0_k3 的函数
std::complex<double> up_R_Cc0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW *SW))\
    *(std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Md_x(3) * set_yu_y(i) * set_yd_z(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * ((-set_Mu_w(j) * set_Mu_w(j) + 2.0 * MH5 * MH5 - 2.0 * set_Mu_w(i) * set_Mu_w(i) + set_Md_x(3) * set_Md_x(3)) * set_yu_y(i) + 2.0 * set_Mu_w(i) * set_Md_x(3) *set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// 计算 R_Cc1 = R_Cc1_k1+R_Cc1_k2+R_Cc1_k3 的函数
//  （二十四）计算 up_R_Cc1_k1 的函数
std::complex<double> up_R_Cc1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5) * set_yu_y(i) + set_Mu_w(i) * set_Md_x(1) * set_yd_z(1)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

//  （二十五）计算 up_R_Cc1_k2 的函数
std::complex<double> up_R_Cc1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5) * set_yu_y(i) + set_Mu_w(i) * set_Md_x(2) * set_yd_z(2)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

//  （二十六）计算 up_R_Cc1_k3 的函数
std::complex<double> up_R_Cc1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (-((set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5) * set_yu_y(i) + set_Mu_w(i) * set_Md_x(3) * set_yd_z(3)))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(i) * set_yd_z(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// 计算 R_Cc2 = R_Cc2_k1+R_Cc2_k2+R_Cc2_k3 的函数
// （二十七）计算 up_R_Cc2_k1 的函数(这里使用的是down的表达式，down的表达式提取了1/sw^2)
std::complex<double> up_R_Cc2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i)\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j)))\
        * VCKM(j,1) * VCKMC(i,1);
}

//  （二十八）计算 up_R_Cc2_k2 的函数
std::complex<double> up_R_Cc2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i)\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j)))\
        * VCKM(j,2) * VCKMC(i,2);
}

//  （二十九）计算 up_R_Cc2_k3 的函数
std::complex<double> up_R_Cc2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI *(set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + 2.0 * set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(i)\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yu_y(i) * set_yu_y(j)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// 计算 R_Cd0 = R_Cd0_k1+R_Cd0_k2+R_Cd0_k3 的函数
//  （三十）计算 up_R_Cd0_k1 的函数
std::complex<double> up_R_Cd0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(1) * set_yu_y(i) * set_yd_z(1) + set_Mu_w(i) * set_yd_z(1) * set_yd_z(1))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

//  （三十一）计算 up_R_Cd0_k2 的函数
std::complex<double> up_R_Cd0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(2) * set_yu_y(i) * set_yd_z(2) + set_Mu_w(i) * set_yd_z(2) * set_yd_z(2))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

//  （三十二）计算 up_R_Cd0_k3 的函数
std::complex<double> up_R_Cd0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(3) * set_yu_y(i) * set_yd_z(3) + set_Mu_w(i) * set_yd_z(3) * set_yd_z(3))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// 计算 R_Cd1 = R_Cd1_k1+R_Cd1_k2+R_Cd1_k3 的函数
//  （三十三）计算 up_R_Cd1_k1 的函数
std::complex<double> up_R_Cd1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) *set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(1) * set_yd_z(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) + set_Mu_w(i) * set_yd_z(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1); 
}

//  （三十四）计算 up_R_Cd1_k2 的函数
std::complex<double> up_R_Cd1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) *set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(2) * set_yd_z(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) + set_Mu_w(i) * set_yd_z(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

//  （三十五）计算 up_R_Cd1_k3 的函数
std::complex<double> up_R_Cd1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) *set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(3) * set_yd_z(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(j) * set_yu_y(i) * set_yu_y(j) + set_Mu_w(i) * set_yd_z(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// 计算 R_Cd2 = R_Cd2_k1+R_Cd2_k2+R_Cd2_k3 的函数
//  （三十六）计算 up_R_Cd2_k1 的函数
std::complex<double> up_R_Cd2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yd_z(1) * set_yd_z(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

//  （三十七）计算 up_R_Cd2_k2 的函数
std::complex<double> up_R_Cd2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yd_z(2) * set_yd_z(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

//  （三十八）计算 up_R_Cd2_k3 的函数
std::complex<double> up_R_Cd2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::cos(thetaH) * std::sin(thetaH) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yd_z(3) * set_yd_z(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Md_x(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// P_L左手
// （壹）计算 up_LL_Ba 的函数
std::complex<double> up_LL_Ba (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yu_y(j) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(i,1) * VCKMC(j,1) + VCKM(i,2) * VCKMC(j,2) + VCKM(i,3) * VCKMC(j,3));
}

// （贰）计算 up_LL_Bb 的函数
std::complex<double> up_LL_Bb (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (-alpha * std::cos(thetaH) * std::sin(thetaH) * set_yu_y(j) / (8.0 * std::sqrt(6.0) * M_PI * std::pow(SW,2))) * (VCKM(i,1) * VCKMC(j,1) + VCKM(i,2) * VCKMC(j,2) + VCKM(i,3) * VCKMC(j,3));
}

// （叁）计算 up_L_Ca0_k1 的函数
std::complex<double> up_L_Ca0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(1) * set_yu_y(j) * set_yd_z(1) - set_Mu_w(j) * set_yd_z(1) * set_yd_z(1)) * VCKM(j,1) * VCKMC(i,1);
}

// （肆）计算 up_L_Ca0_k2 的函数
std::complex<double> up_L_Ca0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(2) * set_yu_y(j) * set_yd_z(2) - set_Mu_w(j) * set_yd_z(2) * set_yd_z(2)) * VCKM(j,2) * VCKMC(i,2);
}

// （伍）计算 up_L_Ca0_k3 的函数
std::complex<double> up_L_Ca0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Md_x(3) * set_yu_y(j) * set_yd_z(3) - set_Mu_w(j) * set_yd_z(3) * set_yd_z(3))  * VCKM(j,3) * VCKMC(i,3);
}

// （陆）计算 up_L_Ca1_k1 的函数
std::complex<double> up_L_Ca1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(1) * set_yd_z(1)) * VCKM(j,1) * VCKMC(i,1);
}

// （柒）计算 up_L_Ca1_k2 的函数
std::complex<double> up_L_Ca1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(2) * set_yd_z(2)) * VCKM(j,2) * VCKMC(i,2);
}

// （捌）计算 up_L_Ca1_k3 的函数
std::complex<double> up_L_Ca1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((-std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 *LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(3) * set_yd_z(3)) * VCKM(j,3) * VCKMC(i,3);
}

// （玖）计算 up_L_Ca2_k1 的函数
std::complex<double> up_L_Ca2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(j) * set_yd_z(1) * set_yd_z(1)) * VCKM(j,1) * VCKMC(i,1);
}

// （拾）计算 up_L_Ca2_k2 的函数
std::complex<double> up_L_Ca2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(j) * set_yd_z(2) * set_yd_z(2)) * VCKM(j,2) * VCKMC(i,2);
}

// （拾壹）计算 up_L_Ca2_k3 的函数
std::complex<double> up_L_Ca2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::sin(thetaH) * std::sin(thetaH) * (2.0 * std::cos(thetaH) * std::cos(thetaH) * (3.0 * std::sqrt(2.0) * M2 + (LAMBDA_3() - 2.0 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (std::sqrt(2.0)  * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum))) / (32.0 * std::sqrt(3.0) * M_PI * M_PI)) * (set_Mu_w(j) * set_yd_z(3) * set_yd_z(3)) * VCKM(j,3) * VCKMC(i,3);
}


// （拾贰）计算 up_L_Cb0_k1 的函数
std::complex<double> up_L_Cb0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Md_x(1) * set_yu_y(j) * set_yd_z(1) + set_Mu_w(j) * set_yd_z(1) * set_yd_z(1))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Mu_w(j) * set_Mu_w(j) * set_yu_y(j) - set_Md_x(1) * set_Md_x(1) * set_yu_y(j) - 2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + set_Mu_w(j) * set_Md_x(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

// （拾叁）计算 up_L_Cb0_k2 的函数
std::complex<double> up_L_Cb0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Md_x(2) * set_yu_y(j) * set_yd_z(2) + set_Mu_w(j) * set_yd_z(2) * set_yd_z(2))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Mu_w(j) * set_Mu_w(j) * set_yu_y(j) - set_Md_x(2) * set_Md_x(2) * set_yu_y(j) - 2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + set_Mu_w(j) * set_Md_x(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （拾肆）计算 up_L_Cb0_k3 的函数
std::complex<double> up_L_Cb0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    * (16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW *(2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0)  * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (-set_Md_x(3) * set_yu_y(j) * set_yd_z(3) + set_Mu_w(j) * set_yd_z(3) * set_yd_z(3))\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * set_Mu_w(j) * set_Mu_w(j) * set_yu_y(j) - set_Md_x(3) * set_Md_x(3) * set_yu_y(j) - 2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + set_Mu_w(j) * set_Md_x(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// （拾伍）计算 up_L_Cb1_k1 的函数
std::complex<double> up_L_Cb1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(i) * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + 2.0 * set_Mu_w(j) * set_Md_x(1) * set_yd_z(1))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

// （拾陆）计算 up_L_Cb1_k2 的函数
std::complex<double> up_L_Cb1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(i) * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + 2.0 * set_Mu_w(j) * set_Md_x(2) * set_yd_z(2))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （拾柒）计算 up_L_Cb1_k3 的函数
std::complex<double> up_L_Cb1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + 2.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (2.0 * MH5 * MH5 * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(i) * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + 2.0 * set_Mu_w(j) * set_Md_x(3) * set_yd_z(3))\
    - std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// （拾捌）计算 up_L_Cb2_k1 的函数
std::complex<double> up_L_Cb2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yd_z(1) * set_yd_z(1))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Mu_w(j) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + set_Mu_w(j) * set_Md_x(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

// （拾玖）计算 up_L_Cb2_k2 的函数
std::complex<double> up_L_Cb2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yd_z(2) * set_yd_z(2))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Mu_w(j) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + set_Mu_w(j) * set_Md_x(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （貳拾）计算 up_L_Cb2_k3 的函数
std::complex<double> up_L_Cb2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1.0 / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW * SW *SW))\
    *(16.0 * alpha * alpha * set_Mu_w(j) * M_PI * M_PI * std::sin(thetaH) * vacuum\
    + std::cos(thetaH) * std::cos(thetaH) * SW * SW * SW * SW * (2.0 * std::sin(thetaH) * std::sin(thetaH) * (3.0 * std::sqrt(2.0) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (std::sqrt(2.0) * M1 + 3.0 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yd_z(3) * set_yd_z(3))\
    + 4.0 * std::sqrt(2.0) * alpha * std::cos(thetaH) * M_PI * std::sin(thetaH) * SW * SW * (set_Mu_w(j) * set_Mu_w(j) * set_yu_y(j) - set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) + set_Mu_w(j) * set_Md_x(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}


// （貳拾壹）计算 up_L_Cc0_k1 的函数
std::complex<double> up_L_Cc0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(1) * set_yu_y(j) * set_yd_z(1) + set_Mu_w(j) * set_yd_z(1) * set_yd_z(1))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(1) * set_yd_z(1)))\
       * VCKM(j,1) * VCKMC(i,1); 
}

// （貳拾貳）计算 up_L_Cc0_k2 的函数
std::complex<double> up_L_Cc0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(2) * set_yu_y(j) * set_yd_z(2) + set_Mu_w(j) * set_yd_z(2) * set_yd_z(2))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(2) * set_yd_z(2)))\
       * VCKM(j,2) * VCKMC(i,2); 
}

// （貳拾叁）计算 up_L_Cc0_k3 的函数
std::complex<double> up_L_Cc0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH))/(32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW *(-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Md_x(3) * set_yu_y(j) * set_yd_z(3) + set_Mu_w(j) * set_yd_z(3) * set_yd_z(3))\
    + 4.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(3) * set_yd_z(3)))\
       * VCKM(j,3) * VCKMC(i,3); 
}

// （貳拾肆）计算 up_L_Cc1_k1 的函数
std::complex<double> up_L_Cc1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(1) * set_yd_z(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) + set_Mu_w(j) * set_yd_z(1) * set_yd_z(1)))\
         * VCKM(j,1) * VCKMC(i,1);
}

// （貳拾伍）计算 up_L_Cc1_k2 的函数
std::complex<double> up_L_Cc1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(2) * set_yd_z(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) + set_Mu_w(j) * set_yd_z(2) * set_yd_z(2)))\
         * VCKM(j,2) * VCKMC(i,2);
}

// （貳拾陆）计算 up_L_Cc1_k3 的函数
std::complex<double> up_L_Cc1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * (set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(3) * set_yd_z(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (-set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) + set_Mu_w(j) * set_yd_z(3) * set_yd_z(3)))\
         * VCKM(j,3) * VCKMC(i,3);
}

// （貳拾柒）计算 up_L_Cc2_k1 的函数
std::complex<double> up_L_Cc2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yd_z(1) * set_yd_z(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

// （貳拾捌）计算 up_L_Cc2_k2 的函数
std::complex<double> up_L_Cc2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yd_z(2) * set_yd_z(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （貳拾玖）计算 up_L_Cc2_k3 的函数
std::complex<double> up_L_Cc2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(j) * set_yd_z(3) * set_yd_z(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI * (2.0 * set_Mu_w(i) * set_Mu_w(j) * set_yu_y(i) - set_Mu_w(j) * set_Md_x(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}


// （叁拾）计算 up_L_Cd0_k1 的函数
std::complex<double> up_L_Cd0_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Md_x(1) * set_yu_y(j) * set_yd_z(1))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Mu_w(j) * set_Mu_w(j) + 2.0 * MH5 * MH5 - set_Mu_w(i) * set_Mu_w(i) + set_Md_x(1) * set_Md_x(1)) * set_yu_y(j) + 2.0 * set_Mu_w(j) * set_Md_x(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1); 
}

// （叁拾壹）计算 up_L_Cd0_k2 的函数
std::complex<double> up_L_Cd0_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Md_x(2) * set_yu_y(j) * set_yd_z(2))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Mu_w(j) * set_Mu_w(j) + 2.0 * MH5 * MH5 - set_Mu_w(i) * set_Mu_w(i) + set_Md_x(2) * set_Md_x(2)) * set_yu_y(j) + 2.0 * set_Mu_w(j) * set_Md_x(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2); 
}

// （叁拾贰）计算 up_L_Cd0_k3 的函数
std::complex<double> up_L_Cd0_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Md_x(3) * set_yu_y(j) * set_yd_z(3))\
    + 2.0 * std::sqrt(2.0) * alpha * M_PI *((-2.0 * set_Mu_w(j) * set_Mu_w(j) + 2.0 * MH5 * MH5 - set_Mu_w(i) * set_Mu_w(i) + set_Md_x(3) * set_Md_x(3)) * set_yu_y(j) + 2.0 * set_Mu_w(j) * set_Md_x(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3); 
}

// （叁拾叁）计算 up_L_Cd1_k1 的函数
std::complex<double> up_L_Cd1_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(j) + set_Mu_w(j) * set_Md_x(1) * set_yd_z(1))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(1) * set_yd_z(1)))\
        * VCKM(j,1) * VCKMC(i,1);
}

// （叁拾肆）计算 up_L_Cd1_k2 的函数
std::complex<double> up_L_Cd1_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(j) + set_Mu_w(j) * set_Md_x(2) * set_yd_z(2))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(2) * set_yd_z(2)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （叁拾伍）计算 up_L_Cd1_k3 的函数
std::complex<double> up_L_Cd1_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * MH5 * MH5 - set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(j) + set_Mu_w(j) * set_Md_x(3) * set_yd_z(3))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j) - set_Mu_w(j) * set_yd_z(3) * set_yd_z(3)))\
        * VCKM(j,3) * VCKMC(i,3);
}

// （叁拾陆）计算 up_L_Cd2_k1 的函数
std::complex<double> up_L_Cd2_k1 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j)))\
        * VCKM(j,1) * VCKMC(i,1); 
}

// （叁拾柒）计算 up_L_Cd2_k2 的函数
std::complex<double> up_L_Cd2_k2 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j)))\
        * VCKM(j,2) * VCKMC(i,2);
}

// （叁拾捌）计算 up_L_Cd2_k3 的函数
std::complex<double> up_L_Cd2_k3 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH)) / (32.0 * std::sqrt(3.0) * M_PI * M_PI * SW * SW))\
    * (-2.0 * std::sqrt(2.0) * alpha * M_PI * ((2.0 * set_Mu_w(j) * set_Mu_w(j) - 2.0 * MH5 * MH5 + set_Mu_w(i) * set_Mu_w(i)) * set_yu_y(j))\
    + std::cos(thetaH) * SW * SW * (-2.0 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (std::sqrt(2.0) * M1 - 6.0 * std::sqrt(2.0) * M2 + (-2.0  *LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * (set_Mu_w(i) * set_yu_y(i) * set_yu_y(j)))\
        * VCKM(j,3) * VCKMC(i,3);
}

//计算c函数
// F_Ca_0
//注意这里的F_Ca_01，第一个0代表C0函数的0,第二个1代表Md（k）的k取1
std::complex<double> up_F_Ca_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Ca_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Md_x(2)*set_Md_x(2));
}
std::complex<double> up_F_Ca_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Md_x(3)*set_Md_x(3));
}
// F_Ca_1
std::complex<double> up_F_Ca_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Ca_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Ca_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Ca_2
std::complex<double> up_F_Ca_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Ca_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Ca_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MH3*MH3,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_0
std::complex<double> up_F_Cb_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cb_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cb_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_1
std::complex<double> up_F_Cb_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cb_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cb_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cb_2
std::complex<double> up_F_Cb_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cb_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cb_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MW*MW,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_0
std::complex<double> up_F_Cc_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cc_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cc_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_1
std::complex<double> up_F_Cc_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cc_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cc_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cc_2
std::complex<double> up_F_Cc_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cc_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cc_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(i)*set_Mu_w(i),set_Mu_w(j)*set_Mu_w(j),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_0
std::complex<double> up_F_Cd_01 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cd_02 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cd_03 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_1
std::complex<double> up_F_Cd_11 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cd_12 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cd_13 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}
// F_Cd_2
std::complex<double> up_F_Cd_21 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Md_x(1)*set_Md_x(1));
}
std::complex<double> up_F_Cd_22 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Mu_w(2)*set_Mu_w(2));
}
std::complex<double> up_F_Cd_23 (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2,MH5*MH5,set_Mu_w(j)*set_Mu_w(j),set_Mu_w(i)*set_Mu_w(i),MW*MW,MH3*MH3,set_Mu_w(3)*set_Mu_w(3));
}

// 计算 up_alpha_ij_R(i,j) 函数
std::complex<double> up_alpha_ij_R (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    //调用一大堆函数
    std::complex<double> up_F_Ba = B0(MH5*MH5,MH3*MH3,MW*MW);
    std::complex<double> up_F_Bb = B0(MH5*MH5,MW*MW,MW*MW);
    std::complex<double> up_F_Ca_01_value = up_F_Ca_01(i,j);
    std::complex<double> up_F_Ca_02_value = up_F_Ca_02(i,j);
    std::complex<double> up_F_Ca_03_value = up_F_Ca_03(i,j);
    std::complex<double> up_F_Ca_11_value = up_F_Ca_11(i,j);
    std::complex<double> up_F_Ca_12_value = up_F_Ca_12(i,j);
    std::complex<double> up_F_Ca_13_value = up_F_Ca_13(i,j);
    std::complex<double> up_F_Ca_21_value = up_F_Ca_21(i,j);
    std::complex<double> up_F_Ca_22_value = up_F_Ca_22(i,j);
    std::complex<double> up_F_Ca_23_value = up_F_Ca_23(i,j);
    std::complex<double> up_F_Cb_01_value = up_F_Cb_01(i,j);
    std::complex<double> up_F_Cb_02_value = up_F_Cb_02(i,j);
    std::complex<double> up_F_Cb_03_value = up_F_Cb_03(i,j);
    std::complex<double> up_F_Cb_11_value = up_F_Cb_11(i,j);
    std::complex<double> up_F_Cb_12_value = up_F_Cb_12(i,j);
    std::complex<double> up_F_Cb_13_value = up_F_Cb_13(i,j);
    std::complex<double> up_F_Cb_21_value = up_F_Cb_21(i,j);
    std::complex<double> up_F_Cb_22_value = up_F_Cb_22(i,j);
    std::complex<double> up_F_Cb_23_value = up_F_Cb_23(i,j);
    std::complex<double> up_F_Cc_01_value = up_F_Cc_01(i,j);
    std::complex<double> up_F_Cc_02_value = up_F_Cc_02(i,j);
    std::complex<double> up_F_Cc_03_value = up_F_Cc_03(i,j);
    std::complex<double> up_F_Cc_11_value = up_F_Cc_11(i,j);
    std::complex<double> up_F_Cc_12_value = up_F_Cc_12(i,j);
    std::complex<double> up_F_Cc_13_value = up_F_Cc_13(i,j);
    std::complex<double> up_F_Cc_21_value = up_F_Cc_21(i,j);
    std::complex<double> up_F_Cc_22_value = up_F_Cc_22(i,j);
    std::complex<double> up_F_Cc_23_value = up_F_Cc_23(i,j);
    std::complex<double> up_F_Cd_01_value = up_F_Cd_01(i,j);
    std::complex<double> up_F_Cd_02_value = up_F_Cd_02(i,j);
    std::complex<double> up_F_Cd_03_value = up_F_Cd_03(i,j);
    std::complex<double> up_F_Cd_11_value = up_F_Cd_11(i,j);
    std::complex<double> up_F_Cd_12_value = up_F_Cd_12(i,j);
    std::complex<double> up_F_Cd_13_value = up_F_Cd_13(i,j);
    std::complex<double> up_F_Cd_21_value = up_F_Cd_21(i,j);
    std::complex<double> up_F_Cd_22_value = up_F_Cd_22(i,j);
    std::complex<double> up_F_Cd_23_value = up_F_Cd_23(i,j);
    // 计算 up_RR_Ba
    std::complex<double> up_RR_Ba_value = up_RR_Ba(i, j);
    // 计算 up_RR_Bb
    std::complex<double> up_RR_Bb_value = up_RR_Bb(i, j);
    // 计算 up_R_Ca0
    std::complex<double> up_R_Ca0_k1_value = up_R_Ca0_k1(i, j);
    std::complex<double> up_R_Ca0_k2_value = up_R_Ca0_k2(i, j);
    std::complex<double> up_R_Ca0_k3_value = up_R_Ca0_k3(i, j);
    // 计算 up_R_Ca1
    std::complex<double> up_R_Ca1_k1_value = up_R_Ca1_k1(i, j);
    std::complex<double> up_R_Ca1_k2_value = up_R_Ca1_k2(i, j);
    std::complex<double> up_R_Ca1_k3_value = up_R_Ca1_k3(i, j);
    // 计算 up_R_Ca2
    std::complex<double> up_R_Ca2_k1_value = up_R_Ca2_k1(i, j);
    std::complex<double> up_R_Ca2_k2_value = up_R_Ca2_k2(i, j);
    std::complex<double> up_R_Ca2_k3_value = up_R_Ca2_k3(i, j);
    // 计算 up_R_Cb0
    std::complex<double> up_R_Cb0_k1_value = up_R_Cb0_k1(i, j);
    std::complex<double> up_R_Cb0_k2_value = up_R_Cb0_k2(i, j);
    std::complex<double> up_R_Cb0_k3_value = up_R_Cb0_k3(i, j);
    // 计算 up_R_Cb1
    std::complex<double> up_R_Cb1_k1_value = up_R_Cb1_k1(i, j);
    std::complex<double> up_R_Cb1_k2_value = up_R_Cb1_k2(i, j);
    std::complex<double> up_R_Cb1_k3_value = up_R_Cb1_k3(i, j);
    // 计算 up_R_Cb2
    std::complex<double> up_R_Cb2_k1_value = up_R_Cb2_k1(i, j);
    std::complex<double> up_R_Cb2_k2_value = up_R_Cb2_k2(i, j);
    std::complex<double> up_R_Cb2_k3_value = up_R_Cb2_k3(i, j);
    // 计算 up_R_Cc0
    std::complex<double> up_R_Cc0_k1_value = up_R_Cc0_k1(i, j);
    std::complex<double> up_R_Cc0_k2_value = up_R_Cc0_k2(i, j);
    std::complex<double> up_R_Cc0_k3_value = up_R_Cc0_k3(i, j);
    // 计算 up_R_Cc1
    std::complex<double> up_R_Cc1_k1_value = up_R_Cc1_k1(i, j);
    std::complex<double> up_R_Cc1_k2_value = up_R_Cc1_k2(i, j);
    std::complex<double> up_R_Cc1_k3_value = up_R_Cc1_k3(i, j);
    // 计算 up_R_Cc2
    std::complex<double> up_R_Cc2_k1_value = up_R_Cc2_k1(i, j);
    std::complex<double> up_R_Cc2_k2_value = up_R_Cc2_k2(i, j);
    std::complex<double> up_R_Cc2_k3_value = up_R_Cc2_k3(i, j);
    // 计算 up_R_Cd0
    std::complex<double> up_R_Cd0_k1_value = up_R_Cd0_k1(i, j);
    std::complex<double> up_R_Cd0_k2_value = up_R_Cd0_k2(i, j);
    std::complex<double> up_R_Cd0_k3_value = up_R_Cd0_k3(i, j);
    // 计算up_R_Cd1
    std::complex<double> up_R_Cd1_k1_value = up_R_Cd1_k1(i, j);
    std::complex<double> up_R_Cd1_k2_value = up_R_Cd1_k2(i, j);
    std::complex<double> up_R_Cd1_k3_value = up_R_Cd1_k3(i, j);
    // 计算 up_R_Cd2
    std::complex<double> up_R_Cd2_k1_value = up_R_Cd2_k1(i, j);
    std::complex<double> up_R_Cd2_k2_value = up_R_Cd2_k2(i, j);
    std::complex<double> up_R_Cd2_k3_value = up_R_Cd2_k3(i, j);

    return up_RR_Ba_value * up_F_Ba + up_RR_Bb_value * up_F_Bb\
           + up_R_Ca0_k1_value * up_F_Ca_01_value + up_R_Ca0_k2_value * up_F_Ca_02_value + up_R_Ca0_k3_value * up_F_Ca_03_value\
           + up_R_Ca1_k1_value * up_F_Ca_11_value + up_R_Ca1_k2_value * up_F_Ca_12_value + up_R_Ca1_k3_value * up_F_Ca_13_value\
           + up_R_Ca2_k1_value * up_F_Ca_21_value + up_R_Ca2_k2_value * up_F_Ca_22_value + up_R_Ca2_k3_value * up_F_Ca_23_value\
           + up_R_Cb0_k1_value * up_F_Cb_01_value + up_R_Cb0_k2_value * up_F_Cb_02_value + up_R_Cb0_k3_value * up_F_Cb_03_value\
           + up_R_Cb1_k1_value * up_F_Cb_11_value + up_R_Cb1_k2_value * up_F_Cb_12_value + up_R_Cb1_k3_value * up_F_Cb_13_value\
           + up_R_Cb2_k1_value * up_F_Cb_21_value + up_R_Cb2_k2_value * up_F_Cb_22_value + up_R_Cb2_k3_value * up_F_Cb_23_value\
           + up_R_Cc0_k1_value * up_F_Cc_01_value + up_R_Cc0_k2_value * up_F_Cc_02_value + up_R_Cc0_k3_value * up_F_Cc_03_value\
           + up_R_Cc1_k1_value * up_F_Cc_11_value + up_R_Cc1_k2_value * up_F_Cc_12_value + up_R_Cc1_k3_value * up_F_Cc_13_value\
           + up_R_Cc2_k1_value * up_F_Cc_21_value + up_R_Cc2_k2_value * up_F_Cc_22_value + up_R_Cc2_k3_value * up_F_Cc_23_value\
           + up_R_Cd0_k1_value * up_F_Cd_01_value + up_R_Cd0_k2_value * up_F_Cd_02_value + up_R_Cd0_k3_value * up_F_Cd_03_value\
           + up_R_Cd1_k1_value * up_F_Cd_11_value + up_R_Cd1_k2_value * up_F_Cd_12_value + up_R_Cd1_k3_value * up_F_Cd_13_value\
           + up_R_Cd2_k1_value * up_F_Cd_21_value + up_R_Cd2_k2_value * up_F_Cd_22_value + up_R_Cd2_k3_value * up_F_Cd_23_value;
}

//计算up_beta_ij_L(i,j)的函数
std::complex<double> up_beta_ij_L (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    //调用一大堆函数
    std::complex<double> up_F_Ba = B0(MH5*MH5,MH3*MH3,MW*MW);
    std::complex<double> up_F_Bb = B0(MH5*MH5,MW*MW,MW*MW);
    std::complex<double> up_F_Ca_01_value = up_F_Ca_01(i,j);
    std::complex<double> up_F_Ca_02_value = up_F_Ca_02(i,j);
    std::complex<double> up_F_Ca_03_value = up_F_Ca_03(i,j);
    std::complex<double> up_F_Ca_11_value = up_F_Ca_11(i,j);
    std::complex<double> up_F_Ca_12_value = up_F_Ca_12(i,j);
    std::complex<double> up_F_Ca_13_value = up_F_Ca_13(i,j);
    std::complex<double> up_F_Ca_21_value = up_F_Ca_21(i,j);
    std::complex<double> up_F_Ca_22_value = up_F_Ca_22(i,j);
    std::complex<double> up_F_Ca_23_value = up_F_Ca_23(i,j);
    std::complex<double> up_F_Cb_01_value = up_F_Cb_01(i,j);
    std::complex<double> up_F_Cb_02_value = up_F_Cb_02(i,j);
    std::complex<double> up_F_Cb_03_value = up_F_Cb_03(i,j);
    std::complex<double> up_F_Cb_11_value = up_F_Cb_11(i,j);
    std::complex<double> up_F_Cb_12_value = up_F_Cb_12(i,j);
    std::complex<double> up_F_Cb_13_value = up_F_Cb_13(i,j);
    std::complex<double> up_F_Cb_21_value = up_F_Cb_21(i,j);
    std::complex<double> up_F_Cb_22_value = up_F_Cb_22(i,j);
    std::complex<double> up_F_Cb_23_value = up_F_Cb_23(i,j);
    std::complex<double> up_F_Cc_01_value = up_F_Cc_01(i,j);
    std::complex<double> up_F_Cc_02_value = up_F_Cc_02(i,j);
    std::complex<double> up_F_Cc_03_value = up_F_Cc_03(i,j);
    std::complex<double> up_F_Cc_11_value = up_F_Cc_11(i,j);
    std::complex<double> up_F_Cc_12_value = up_F_Cc_12(i,j);
    std::complex<double> up_F_Cc_13_value = up_F_Cc_13(i,j);
    std::complex<double> up_F_Cc_21_value = up_F_Cc_21(i,j);
    std::complex<double> up_F_Cc_22_value = up_F_Cc_22(i,j);
    std::complex<double> up_F_Cc_23_value = up_F_Cc_23(i,j);
    std::complex<double> up_F_Cd_01_value = up_F_Cd_01(i,j);
    std::complex<double> up_F_Cd_02_value = up_F_Cd_02(i,j);
    std::complex<double> up_F_Cd_03_value = up_F_Cd_03(i,j);
    std::complex<double> up_F_Cd_11_value = up_F_Cd_11(i,j);
    std::complex<double> up_F_Cd_12_value = up_F_Cd_12(i,j);
    std::complex<double> up_F_Cd_13_value = up_F_Cd_13(i,j);
    std::complex<double> up_F_Cd_21_value = up_F_Cd_21(i,j);
    std::complex<double> up_F_Cd_22_value = up_F_Cd_22(i,j);
    std::complex<double> up_F_Cd_23_value = up_F_Cd_23(i,j);
    // 计算 up_LL_Ba
    std::complex<double> up_LL_Ba_value = up_LL_Ba(i, j);
    // 计算 up_LL_Bb
    std::complex<double> up_LL_Bb_value = up_LL_Bb(i, j);
    // 计算 up_L_Ca0
    std::complex<double> up_L_Ca0_k1_value = up_L_Ca0_k1(i, j);
    std::complex<double> up_L_Ca0_k2_value = up_L_Ca0_k2(i, j);
    std::complex<double> up_L_Ca0_k3_value = up_L_Ca0_k3(i, j);
    // 计算 up_L_Ca1
    std::complex<double> up_L_Ca1_k1_value = up_L_Ca1_k1(i, j);
    std::complex<double> up_L_Ca1_k2_value = up_L_Ca1_k2(i, j);
    std::complex<double> up_L_Ca1_k3_value = up_L_Ca1_k3(i, j);
    // 计算 up_L_Ca2
    std::complex<double> up_L_Ca2_k1_value = up_L_Ca2_k1(i, j);
    std::complex<double> up_L_Ca2_k2_value = up_L_Ca2_k2(i, j);
    std::complex<double> up_L_Ca2_k3_value = up_L_Ca2_k3(i, j);
    // 计算 up_L_Cb0
    std::complex<double> up_L_Cb0_k1_value = up_L_Cb0_k1(i, j);
    std::complex<double> up_L_Cb0_k2_value = up_L_Cb0_k2(i, j);
    std::complex<double> up_L_Cb0_k3_value = up_L_Cb0_k3(i, j);
    // 计算 up_L_Cb1
    std::complex<double> up_L_Cb1_k1_value = up_L_Cb1_k1(i, j);
    std::complex<double> up_L_Cb1_k2_value = up_L_Cb1_k2(i, j);
    std::complex<double> up_L_Cb1_k3_value = up_L_Cb1_k3(i, j);
    // 计算 up_L_Cb2
    std::complex<double> up_L_Cb2_k1_value = up_L_Cb2_k1(i, j);
    std::complex<double> up_L_Cb2_k2_value = up_L_Cb2_k2(i, j);
    std::complex<double> up_L_Cb2_k3_value = up_L_Cb2_k3(i, j);
    // 计算 up_L_Cc0
    std::complex<double> up_L_Cc0_k1_value = up_L_Cc0_k1(i, j);
    std::complex<double> up_L_Cc0_k2_value = up_L_Cc0_k2(i, j);
    std::complex<double> up_L_Cc0_k3_value = up_L_Cc0_k3(i, j);
    // 计算 up_L_Cc1
    std::complex<double> up_L_Cc1_k1_value = up_L_Cc1_k1(i, j);
    std::complex<double> up_L_Cc1_k2_value = up_L_Cc1_k2(i, j);
    std::complex<double> up_L_Cc1_k3_value = up_L_Cc1_k3(i, j);
    // 计算 up_L_Cc2
    std::complex<double> up_L_Cc2_k1_value = up_L_Cc2_k1(i, j);
    std::complex<double> up_L_Cc2_k2_value = up_L_Cc2_k2(i, j);
    std::complex<double> up_L_Cc2_k3_value = up_L_Cc2_k3(i, j);
    // 计算 up_L_Cd0
    std::complex<double> up_L_Cd0_k1_value = up_L_Cd0_k1(i, j);
    std::complex<double> up_L_Cd0_k2_value = up_L_Cd0_k2(i, j);
    std::complex<double> up_L_Cd0_k3_value = up_L_Cd0_k3(i, j);
    // 计算 up_L_Cd1
    std::complex<double> up_L_Cd1_k1_value = up_L_Cd1_k1(i, j);
    std::complex<double> up_L_Cd1_k2_value = up_L_Cd1_k2(i, j);
    std::complex<double> up_L_Cd1_k3_value = up_L_Cd1_k3(i, j);
    // 计算 up_L_Cd2
    std::complex<double> up_L_Cd2_k1_value = up_L_Cd2_k1(i, j);
    std::complex<double> up_L_Cd2_k2_value = up_L_Cd2_k2(i, j);
    std::complex<double> up_L_Cd2_k3_value = up_L_Cd2_k3(i, j);

    return up_LL_Ba_value * up_F_Ba + up_LL_Bb_value * up_F_Bb\
           + up_L_Ca0_k1_value * up_F_Ca_01_value + up_L_Ca0_k2_value * up_F_Ca_02_value + up_L_Ca0_k3_value * up_F_Ca_03_value\
           + up_L_Ca1_k1_value * up_F_Ca_11_value + up_L_Ca1_k2_value * up_F_Ca_12_value + up_L_Ca1_k3_value * up_F_Ca_13_value\
           + up_L_Ca2_k1_value * up_F_Ca_21_value + up_L_Ca2_k2_value * up_F_Ca_22_value + up_L_Ca2_k3_value * up_F_Ca_23_value\
           + up_L_Cb0_k1_value * up_F_Cb_01_value + up_L_Cb0_k2_value * up_F_Cb_02_value + up_L_Cb0_k3_value * up_F_Cb_03_value\
           + up_L_Cb1_k1_value * up_F_Cb_11_value + up_L_Cb1_k2_value * up_F_Cb_12_value + up_L_Cb1_k3_value * up_F_Cb_13_value\
           + up_L_Cb2_k1_value * up_F_Cb_21_value + up_L_Cb2_k2_value * up_F_Cb_22_value + up_L_Cb2_k3_value * up_F_Cb_23_value\
           + up_L_Cc0_k1_value * up_F_Cc_01_value + up_L_Cc0_k2_value * up_F_Cc_02_value + up_L_Cc0_k3_value * up_F_Cc_03_value\
           + up_L_Cc1_k1_value * up_F_Cc_11_value + up_L_Cc1_k2_value * up_F_Cc_12_value + up_L_Cc1_k3_value * up_F_Cc_13_value\
           + up_L_Cc2_k1_value * up_F_Cc_21_value + up_L_Cc2_k2_value * up_F_Cc_22_value + up_L_Cc2_k3_value * up_F_Cc_23_value\
           + up_L_Cd0_k1_value * up_F_Cd_01_value + up_L_Cd0_k2_value * up_F_Cd_02_value + up_L_Cd0_k3_value * up_F_Cd_03_value\
           + up_L_Cd1_k1_value * up_F_Cd_11_value + up_L_Cd1_k2_value * up_F_Cd_12_value + up_L_Cd1_k3_value * up_F_Cd_13_value\
           + up_L_Cd2_k1_value * up_F_Cd_21_value + up_L_Cd2_k2_value * up_F_Cd_22_value + up_L_Cd2_k3_value * up_F_Cd_23_value; 
}
