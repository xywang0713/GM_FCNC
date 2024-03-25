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
//  (三) 计算 R_Ca0_k1 的函数，因为涉及到Mu（k）的求和，这里分别对k=1,2,3进行计算，命名后缀加上了k1,k2,k3
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

//计算c函数
// F_Ca_0
//注意这里的F_Ca_01，第一个0代表C0函数的0,第二个1代表Mu（k）的k取1
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

// 计算 alpha_ij_R(i,j) 函数
std::complex<double> alpha_ij_R (int i, int j, GM_model&m) {

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

//计算beta_ij_L(i,j)的函数
std::complex<double> beta_ij_L (int i, int j, GM_model&m) {

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
