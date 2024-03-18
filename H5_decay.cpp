//H50衰变计算
//NO.1 衰变为双光子
//计算A函数里的tau
std::complex<double> tau_W_H50(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return (MH5 * MH5) / (4.0 * MW * MW);
}

//计算A函数里的f(tau)函数
std::complex<double> f_tau(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    std::complex<double> tau_W_H50_value = tau_W_H50 (GM_model&m);

    if (tau_W_H50_value <= 1) {
        f_tau = asin(sqrt(tau_W_H50_value)) * asin(sqrt(tau_W_H50_value));
    } else if (tau_W_H50_value > 1) {
        f_tau = (-1/4) * (log((1 + std::sqrt(1 - 1 / tau_W_H50_value)) / (1 - std::sqrt(1 - 1 / tau_W_H50_value))) - std::sqrt(-1) * M_PI) * (log((1 + std::sqrt(1 - 1 / tau_W_H50_value)) / (1 - std::sqrt(1 - 1 / tau_W_H50_value))) - std::sqrt(-1) * M_PI);
    }
}

//计算A函数
std::complex<double> A_1_function(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    std::complex<double> tau_W_H50_value = tau_W_H50 (GM_model&m);
    std::complex<double> f_tau_value = f_tau (GM_model&m);

    return (-1/(2 * tau_W_H50_value * tau_W_H50_value)) * (2 * tau_W_H50_value * tau_W_H50_value + 3 * tau_W_H50_value + 3 * (2 * tau_W_H50_value - 1) * f_tau_value );

}

std::complex<double> xi_phi_gamma(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    std::complex<double> A_1_function_value = A_1_function (GM_model&m);

    return ((EL * EL * std::sin(thetaH) * vacuum * vacuum) / (4.0 * sqrt(3.0) * SW * MW * MW)) * A_1_function_value;
}

//衰变为光子的衰变width
std::complex<double> GAMMA_gamma(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    std::complex<double> xi_phi_gamma_value = xi_phi_gamma (GM_model&m);

    return ((GF * alpha_ew * alpha_ew * MH5 * MH5 * MH5)/(32 * sqrt(2.0) * M_PI * M_PI * M_PI)) * std::norm(xi_phi_gamma_value) * std::norm(xi_phi_gamma_value);
}


//衰变系数xi_phi_lepton
//用到的质量以及yl
std::complex<double> set_M_lepton_i(int i) {

    std::complex<double> M_lepton_i;
    if (i == 1) {
        M_lepton_i = 0.510998955555e-3;
    } else if (i == 2) {
        M_lepton_i = 105.6583755e-3;
    } else if (i == 3) {
        M_lepton_i = 1776.86e-3;
    } else {
        std::cerr << "M_lepton_i value is wrong!" << std::endl;
    }
    return M_lepton_i;
}

std::complex<double> set_yl_j(int j,GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    std::complex<double> yl_j;
    if (j == 1) {
        yl_j = (0.510998955555e-3*sqrt(2)) / (vacuum * std::cos(thetaH));
    } else if (j == 2) {
        yl_j = (105.6583755e-3*sqrt(2)) / (vacuum * std::cos(thetaH));
    } else if (j == 3) {
        yl_j = (1776.86e-3*sqrt(2)) / (vacuum * std::cos(thetaH));
    } else {
        std::cerr << "yl_j  value is wrong!" << std::endl;
    }
    return yl_j ;
}

//lepton right
//B0i[bb0, Mlepton(i)2, 0, MH32]
std::complex<double> lepton_R_Ba (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -(alpha * std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW *SW);
}

//B0i[bb0, Mlepton(i)2, 0, MW2]
std::complex<double> lepton_R_Bb (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW *SW);
}

//B0i[bb0, Mlepton(i)2, Mlepton(i)2, MH32]
std::complex<double> lepton_R_Bc (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * (CW * CW * CW * CW -SW *SW *SW *SW)* set_yl_j(i)) / (8 * sqrt(6) CW * CW * M_PI * SW *SW);
}

//B0i[bb0, Mlepton(i)2, Mlepton(i)2, MZ2]
std::complex<double> lepton_R_Bd (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -(alpha * std::cos(thetaH) * std::sin(thetaH) * (CW * CW * CW * CW -SW *SW *SW *SW)* set_yl_j(i)) / (8 * sqrt(6) CW * CW * M_PI * SW *SW);
}

//C0i[cc0, Mlepton(i)2, Mlepton(i)2, MH52, MZ2, Mlepton(i)2, MH32]
std::complex<double> lepton_R_Ca0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * (CW * CW *(2 * MH5 * MH5 + MZ * MZ) + (2 * MH3 * MH3 - MZ * MZ) * SW * SW)* set_yl_j(i)) / (8 * sqrt(6) CW * CW * M_PI * SW *SW);
}

//C0i[cc1, Mlepton(i)2, Mlepton(i)2, MH52, MZ2, Mlepton(i)2, MH32]
std::complex<double> lepton_R_Ca1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1/(16 * sqrt(3) * M_PI * M_PI * SW *SW)) * (std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i) * (2 * sqrt(2) * alpha * MH5 * MH5 * M_PI + std::cos(thetaH) * set_M_lepton_i(i) * SW *SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH)*( sqrt(2)*M1-6*sqrt(2)* M2+(-2* LAMBDA_3()+LAMBDA_5())std::sin(thetaH) * vacuum))* set_yl_j (i)));
}

//C0i[cc2, Mlepton(i)2, Mlepton(i)2, M H52, M Z2, Mlepton(i)2, M H32]
std::complex<double> lepton_R_Ca2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * MH5 * MH5 * std::sin(thetaH) * (3* CW * CW  -SW *SW )* set_yl_j(i)) / (8 * sqrt(6) CW * CW * M_PI * SW *SW);
}

//C0i[cc0, Mlepton(i)2, MH52, Mlepton(i)2, 0, MW2, MH32]
std::complex<double> lepton_R_Cb0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * (set_M_lepton_i(i) * set_M_lepton_i(i) - MW * MW) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW *SW);
}

//C0i[cc1, Mlepton(i)2, M H52, Mlepton(i)2, 0, MW2, MH32]
std::complex<double> lepton_R_Cb1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i))  / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (2 * sqrt(2) * (set_M_lepton_i(i) * set_M_lepton_i(i) + 2* MH5 * MH5) * M_PI + std::cos(thetaH) *set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 -6 *sqrt(2) * M2 + (-2 *LAMBDA_3() + LAMBDA_5()) *std::sin(thetaH) * vacuum )) * set_yl_j(i));
}

//C0i[cc2, Mlepton(i)2, MH52, Mlepton(i)2, 0, MW2, MH32]
std::complex<double> lepton_R_Cb2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i))  / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (2 * sqrt(2) * (set_M_lepton_i(i) * set_M_lepton_i(i) - MH5 * MH5) * M_PI + std::cos(thetaH) *set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 -6 *sqrt(2) * M2 + (-2 *LAMBDA_3() + LAMBDA_5()) *std::sin(thetaH) * vacuum )) * set_yl_j(i));
}

//C0i[cc0, Mlepton(i)2, M H52, Mlepton(i)2, 0, MW2, MW2]
std::complex<double> lepton_R_Cc0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * (-set_M_lepton_i(i) * set_M_lepton_i(i) + MW * MW) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW *SW);
}

//C0i[cc1, Mlepton(i)2, M H52, Mlepton(i)2, 0, MW2, MW2]
std::complex<double> lepton_R_Cc1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * set_yl_j(i)) / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (-4 * sqrt(2) *alpha * MH5 * MH5 * M_PI * std::sin(thetaH) + std::cos(thetaH) * set_M_lepton_i(i) * SW *SW *(2 * std::sin(thetaH) * std::sin(thetaH) *(3 *sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std:cos(thetaH) * (sqrt(2) * M1 +3 *LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i));
}


//C0i[cc2, Mlepton(i)2, MH52, Mlepton(i)2, 0, MW2, MW2]
std::complex<double> lepton_R_Cc2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1/(16 * sqrt(3) * M_PI * SW * SW *SW * SW)) * alpha * std::sin(thetaH) * (8 * alpha * set_M_lepton_i(i) * M_PI * vacuum + sqrt(2) * std::cos(thetaH) * (-2 * set_M_lepton_i(i) * set_M_lepton_i(i) + MH5 * MH5) * SW * SW * set_yl_j(i));
}

//C0i[cc0, Mlepton(i)2, MH52, Mlepton(i)2, Mlepton(i)2, MH32, MH32]
std::complex<double> lepton_R_Cd0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) *std::cos(thetaH) *(3* sqrt(2) * M2 +(LAMBDA_3()-2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc1, Mlepton(i)2, MH52, Mlepton(i)2, Mlepton(i)2, MH32, MH32]
std::complex<double> lepton_R_Cd1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) *std::cos(thetaH) *(3* sqrt(2) * M2 +(LAMBDA_3()-2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc2, Mlepton(i)2, MH52, Mlepton(i)2, Mlepton(i)2, MH32, MH32]
std::complex<double> lepton_R_Cd2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) *std::cos(thetaH) *(3* sqrt(2) * M2 +(LAMBDA_3()-2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc1, Mlepton(i)2, MH52, Mlepton(i)2, Mlepton(i)2, MH32, MH32]
std::complex<double> lepton_R_Ce1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) *std::cos(thetaH) *(3* sqrt(2) * M2 +(LAMBDA_3()-2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc0, M H52, Mlepton(i)2, Mlepton(i)2, MZ2, M Z2, Mlepton(i)2]
std::complex<double> lepton_R_Cf0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1/(16 * sqrt(3) * CW * CW * M_PI * SW * SW)) * (std::sin(thetaH) * (32 * alpha * alpha * set_M_lepton_i(i) * M_PI * vacuum - sqrt(2) * alpha * std::cos(thetaH) * (CW * CW * (2 * MH5 * MH5 - MZ * MZ) + (-2 * MH5 * MH5 + MZ * MZ) * SW * SW) * set_yl_j(i)));
}

//C0i[cc1, M H52, Mlepton(i)2, Mlepton(i)2, MZ2, M Z2, Mlepton(i)2]
std::complex<double> lepton_R_Cf1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -(1/(16 * sqrt(3) * CW * CW * CW * CW * M_PI * SW * SW * SW * SW)) * (alpha * std::sin(thetaH) * (8 * alpha * set_M_lepton_i(i) * M_PI * (CW * CW - 3 * SW *SW) * vacuum + sqrt(2) * std::cos(thetaH) * CW * CW * MH5 * MH5 * SW * SW *(3 * CW * CW -5 * SW * SW) * set_yl_j(i)));
}

//C0i[cc2, M H52, Mlepton(i)2, Mlepton(i)2, MZ2, M Z2, Mlepton(i)2]
std::complex<double> lepton_R_Cf2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1/(32 * sqrt(3) * CW * CW * CW * CW)) * (64 * alpha * alpha * set_M_lepton_i(i) * std::sin(thetaH) * vacuum + ((std::cos(thetaH) * CW * CW * set_yl_j(i))/(M_PI * M_PI * SW * SW)) * (4 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) * (- CW * CW * CW * CW + SW *SW * SW *SW) + std::cos(thetaH) * CW *CW * set_M_lepton_i(i) * SW * SW * (2 * std::sin(thetaH) * std::sin(thetaH) * (3 * sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum)+ std::cos(thetaH) * std::cos(thetaH) *(sqrt(2) * M1 +3 * LAMBDA_5() * std ::sin(thetaH) * vacuum)) * set_yl_j(i)));
}

//左手
//B0i[bb0, Mlepton(i)2, Mlepton(i)2, MH32]
std::complex<double> lepton_L_Bc (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (4 * sqrt(6) * M_PI * CW * CW);
}

//B0i[bb0, Mlepton(i)2, Mlepton(i)2, MZ2]
std::complex<double> lepton_L_Bd (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -(alpha * std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (4 * sqrt(6) * M_PI * CW * CW);
}

//C0i[cc0, Mlepton(i)2, Mlepton(i)2, M H52, M Z2, Mlepton(i)2, M H32]
std::complex<double> lepton_L_Ca0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((alpha * std::cos(thetaH) * std::sin(thetaH)) / (8 * sqrt(6) * CW * CW * M_PI * SW *SW)) * (CW * CW *(MH3 * MH3 + MH5 * MH5) + (- MH3 * MH3 + 3 * MH5 * MH5 + 2 * MZ * MZ) * SW * SW) * set_yl_j(i);
}

//C0i[cc1, Mlepton(i)2, Mlepton(i)2, M H52, M Z2, Mlepton(i)2, M H32]
std::complex<double> lepton_L_Ca1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (16 * sqrt(3) * CW * CW * M_PI * SW *SW)) * (sqrt(2) * alpha * MH5 * MH5 * M_PI * (CW * CW + 3 * SW * SW)+ std::cos(thetaH) * CW * CW * set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 - 6 * sqrt(2) * M2 + (-2 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * set_yl_j(i));
}

//C0i[cc2, Mlepton(i)2, Mlepton(i)2, M H52, M Z2, Mlepton(i)2, M H32]
std::complex<double> lepton_L_Ca2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1/(8 * sqrt(6) * CW * CW * M_PI * SW * SW)) * alpha * std::cos(thetaH) * MH5 * MH5 * std::sin(thetaH) * (CW * CW + 5 * SW * SW) * set_yl_j(i);
}


//C0i[cc0, Mlepton(i)2, M H52, Mlepton(i)2, 0, M W 2, M H32]
std::complex<double> lepton_L_Cb0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1/(8 * sqrt(6) * M_PI * SW * SW)) * alpha * std::cos(thetaH) * (set_M_lepton_i(i) * set_M_lepton_i(i) - MH3 * MH3) * std::sin(thetaH) * set_yl_j(i);
}

//C0i[cc1, Mlepton(i)2, M H52, Mlepton(i)2, 0, M W 2, M H32]
std::complex<double> lepton_L_Cb1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (2 * sqrt(2) * alpha * (set_M_lepton_i(i) * set_M_lepton_i(i) + MH5 * MH5) * M_PI + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 - 6 * sqrt(2) * M2 + (-2 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * set_yl_j(i));
}

//C0i[cc2, Mlepton(i)2, M H52, Mlepton(i)2, 0, M W 2, M H32]
std::complex<double> lepton_L_Cb2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * set_M_lepton_i(i) * std::sin(thetaH) * set_yl_j(i)) / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (2 * sqrt(2) * alpha * set_M_lepton_i(i) * M_PI + std::cos(thetaH) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 - 6 * sqrt(2) * M2 + (-2 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * set_yl_j(i));
}

//C0i[cc0, Mlepton(i)2, M H52, Mlepton(i)2, 0, M W 2, M W 2]
std::complex<double> lepton_L_Cc0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((alpha * std::cos(thetaH) * (-set_M_lepton_i(i) * set_M_lepton_i(i) + MW * MW) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW *SW));
}

//C0i[cc1, Mlepton(i)2, M H52, Mlepton(i)2, 0, MW2, MW2]
std::complex<double> lepton_L_Cc1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((set_M_lepton_i(i) * std::sin(thetaH) * (4 * alpha * alpha * M_PI * vacuum - sqrt(2) * alpha * std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * set_yl_j(i))) / (8 * sqrt(3) * M_PI * SW * SW * SW * SW));
}

//C0i[cc2, Mlepton(i)2, MH52, Mlepton(i)2, 0, MW2, MW2]
std::complex<double> lepton_L_Cc2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((std::cos(thetaH) * set_yl_j(i)) / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (-2 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW (2 * std::sin(thetaH) *std::sin(thetaH) * (3 * sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) *vaccum) + std::cos(thetaH) * std::cos(thetaH) * (sqrt(2) * M1 + 3 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i));
}

//C0i[cc0, Mlepton(i)2, M H52, Mlepton(i)2, Mlepton(i)2, M H32, M H32]
std::complex<double> lepton_L_Cd0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH))/(32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() -2 * LAMBDA_5()) * std::sin(thetaH))  + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() *std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc1, Mlepton(i)2, M H52, Mlepton(i)2, Mlepton(i)2, M H32, M H32]
std::complex<double> lepton_L_Cd1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH))/(32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() -2 * LAMBDA_5()) * std::sin(thetaH))  + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() *std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc2, Mlepton(i)2, M H52, Mlepton(i)2, Mlepton(i)2, M H32, M H32]
std::complex<double> lepton_L_Cd2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH))/(32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() -2 * LAMBDA_5()) * std::sin(thetaH))  + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() *std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc2, Mlepton(i)2, M H52, Mlepton(i)2, Mlepton(i)2, M H32, M H32]
std::complex<double> lepton_L_Ce2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH))/(32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() -2 * LAMBDA_5()) * std::sin(thetaH))  + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() *std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
}

//C0i[cc0, M H52, Mlepton(i)2, Mlepton(i)2, M2Z , M2Z , Mlepton(i)2]
std::complex<double> lepton_L_Cf0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (std::sin(thetaH) /(16 * sqrt(3) * CW * CW * CW * CW * M_PI * SW * SW *SW * SW)) * (8 * alpha * alpha * set_M_lepton_i(i) * M_PI * (CW * CW * CW * CW - SW * SW * SW *SW) * (CW * CW + 3 *SW * SW) * vacuum -sqrt(2) * alpha * std::cos(thetaH) * CW * CW *SW *SW *(CW * CW * MZ *MZ + (4 *MH5 * MH5 + MZ *MZ) * SW *SW)* set_yl_j(i));
}

//C0i[cc1, M H52, Mlepton(i)2, Mlepton(i)2, M2Z , M2Z , Mlepton(i)2]
std::complex<double> lepton_L_Cf1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return ((alpha * std::sin(thetaH)) /(16 * sqrt(3) * CW * CW * CW * CW * M_PI * SW * SW *SW * SW)) * (8 * alpha * set_M_lepton_i(i) * M_PI * (CW * CW - 3 *SW * SW) * vacuum - sqrt(2) * std::cos(thetaH) * CW * CW * MH5 * MH5 * SW *SW *(CW * CW - 7 * SW * SW)* set_yl_j(i));
}

//C0i[cc2, M H52, Mlepton(i)2, Mlepton(i)2, M2Z , M2Z , Mlepton(i)2]
std::complex<double> lepton_L_Cf2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (1 / (32 * sqrt(3) * CW * CW * CW *CW * M_PI * M_PI * SW * SW * SW * SW)) * (16 * alpha * alpha * set_M_lepton_i(i) * M_PI * M_PI * std::sin(thetaH) * (CW * CW * CW * CW - SW * SW * SW * SW) * (CW * CW * CW * CW - SW * SW * SW * SW) * vacuum + std::cos(thetaH) * CW * CW * SW * SW * SW *SW * set_yl_j(i) * (-8 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) + std::cos(thetaH) * CW * CW * set_M_lepton_i(i) * (2 * std::sin(thetaH) * std::sin(thetaH) *(3 * sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * (sqrt(2) * M1 + 3 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i)));
}

//B,C函数
std::complex<double> F_lepton_Ba (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return B0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MH3 * MH3);
}
std::complex<double> F_lepton_Bb (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return B0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW);
}

std::complex<double> F_lepton_Bc (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return B0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3);
}

std::complex<double> F_lepton_Bd (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return B0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MZ * MZ);
}

std::complex<double> F_lepton_Ca_0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5 , MZ * MZ ,set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3);
}

std::complex<double> F_lepton_Ca_1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1 , set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5 , MZ * MZ ,set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3);
}

std::complex<double> F_lepton_Ca_2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2 , set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5 , MZ * MZ ,set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3);
}

std::complex<double> F_lepton_Cb_0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW  , MH3 * MH3);
}

std::complex<double> F_lepton_Cb_1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW  , MH3 * MH3);
}

std::complex<double> F_lepton_Cb_2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW  , MH3 * MH3);
}

std::complex<double> F_lepton_Cc_0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW  , MW * MW);
}

std::complex<double> F_lepton_Cc_1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW  , MW * MW);
}

std::complex<double> F_lepton_Cc_2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MW * MW  , MW * MW);
}

std::complex<double> F_lepton_Cd_0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3  , MH3 * MH3);
}

std::complex<double> F_lepton_Cd_1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3  , MH3 * MH3);
}

std::complex<double> F_lepton_Cd_2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , set_M_lepton_i(i) * set_M_lepton_i(i) , MH3 * MH3  , MH3 * MH3);
}

std::complex<double> F_lepton_Ce_0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MH3 * MH3  , MH3 * MH3);
}

std::complex<double> F_lepton_Ce_1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MH3 * MH3  , MH3 * MH3);
}

std::complex<double> F_lepton_Ce_2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2 , set_M_lepton_i(i) * set_M_lepton_i(i) , MH5 * MH5, set_M_lepton_i(i) * set_M_lepton_i(i) , 0 , MH3 * MH3  , MH3 * MH3);
}

std::complex<double> F_lepton_Cf_0 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(0 , MH5 * MH5 , set_M_lepton_i(i) * set_M_lepton_i(i), set_M_lepton_i(i) * set_M_lepton_i(i), MZ * MZ, MZ * MZ ,set_M_lepton_i(i) * set_M_lepton_i(i));
}

std::complex<double> F_lepton_Cf_1 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(1 , MH5 * MH5 , set_M_lepton_i(i) * set_M_lepton_i(i), set_M_lepton_i(i) * set_M_lepton_i(i), MZ * MZ, MZ * MZ ,set_M_lepton_i(i) * set_M_lepton_i(i));
}

std::complex<double> F_lepton_Cf_2 (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return C0i(2 , MH5 * MH5 , set_M_lepton_i(i) * set_M_lepton_i(i), set_M_lepton_i(i) * set_M_lepton_i(i), MZ * MZ, MZ * MZ ,set_M_lepton_i(i) * set_M_lepton_i(i));
}

// 计算 alpha_lepton_i_R(i) 函数
std::complex<double> alpha_lepton_i_R (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    //调用一大堆函数
    std::complex<double> F_lepton_Ba_value = F_lepton_Ba(i);
    std::complex<double> F_lepton_Bb_value = F_lepton_Bb(i);
    std::complex<double> F_lepton_Bc_value = F_lepton_Bc(i);
    std::complex<double> F_lepton_Bd_value = F_lepton_Bd(i);
    std::complex<double> F_lepton_Ca_0_value = F_lepton_Ca_0(i);
    std::complex<double> F_lepton_Ca_1_value = F_lepton_Ca_1(i);
    std::complex<double> F_lepton_Ca_2_value = F_lepton_Ca_2(i);
    std::complex<double> F_lepton_Cb_0_value = F_lepton_Cb_0(i);
    std::complex<double> F_lepton_Cb_1_value = F_lepton_Cb_1(i);
    std::complex<double> F_lepton_Cb_2_value = F_lepton_Cb_2(i);
    std::complex<double> F_lepton_Cc_0_value = F_lepton_Cc_0(i);
    std::complex<double> F_lepton_Cc_1_value = F_lepton_Cc_1(i);
    std::complex<double> F_lepton_Cc_2_value = F_lepton_Cc_2(i);
    std::complex<double> F_lepton_Cd_0_value = F_lepton_Cd_0(i);
    std::complex<double> F_lepton_Cd_1_value = F_lepton_Cd_1(i);
    std::complex<double> F_lepton_Cd_2_value = F_lepton_Cd_2(i);
    std::complex<double> F_lepton_Ce_0_value = F_lepton_Ce_0(i);
    std::complex<double> F_lepton_Ce_1_value = F_lepton_Ce_1(i);
    std::complex<double> F_lepton_Ce_2_value = F_lepton_Ce_2(i);
    std::complex<double> F_lepton_Cf_0_value = F_lepton_Cf_0(i);
    std::complex<double> F_lepton_Cf_1_value = F_lepton_Cf_1(i);
    std::complex<double> F_lepton_Cf_2_value = F_lepton_Cf_2(i);
    //调用右手系数
    std::complex<double> lepton_R_Ba_value = lepton_R_Ba(i);
    std::complex<double> lepton_R_Bb_value = lepton_R_Bb(i);
    std::complex<double> lepton_R_Bc_value = lepton_R_Bc(i);
    std::complex<double> lepton_R_Bd_value = lepton_R_Bd(i);
    std::complex<double> lepton_R_Ca0_value = lepton_R_Ca0(i);
    std::complex<double> lepton_R_Ca1_value = lepton_R_Ca1(i);
    std::complex<double> lepton_R_Ca2_value = lepton_R_Ca2(i);
    std::complex<double> lepton_R_Cb0_value = lepton_R_Cb0(i);
    std::complex<double> lepton_R_Cb1_value = lepton_R_Cb1(i);
    std::complex<double> lepton_R_Cb2_value = lepton_R_Cb2(i);
    std::complex<double> lepton_R_Cc0_value = lepton_R_Cc0(i);
    std::complex<double> lepton_R_Cc1_value = lepton_R_Cc1(i);
    std::complex<double> lepton_R_Cc2_value = lepton_R_Cc2(i);
    std::complex<double> lepton_R_Cd0_value = lepton_R_Cd0(i);
    std::complex<double> lepton_R_Cd1_value = lepton_R_Cd1(i);
    std::complex<double> lepton_R_Cd2_value = lepton_R_Cd2(i);
    std::complex<double> lepton_R_Ce1_value = lepton_R_Ce1(i);
    std::complex<double> lepton_R_Cf0_value = lepton_R_Cf0(i);
    std::complex<double> lepton_R_Cf1_value = lepton_R_Cf1(i);
    std::complex<double> lepton_R_Cf2_value = lepton_R_Cf2(i);

    return lepton_R_Ba_value * F_lepton_Ba_value + lepton_R_Bb_value * F_lepton_Bb_value + lepton_R_Bc_value * F_lepton_Bc_value + lepton_R_Bd_value * F_lepton_Bd_value \
           + lepton_R_Ca0_value * F_lepton_Ca_0_value + lepton_R_Ca1_value * F_lepton_Ca_1_value + lepton_R_Ca2_value * F_lepton_Ca_2_value \
           + lepton_R_Cb0_value * F_lepton_Cb_0_value + lepton_R_Cb1_value * F_lepton_Cb_1_value + lepton_R_Cb2_value * F_lepton_Cb_2_value \
           + lepton_R_Cc0_value * F_lepton_Cc_0_value + lepton_R_Cc1_value * F_lepton_Cc_1_value + lepton_R_Cc2_value * F_lepton_Cc_2_value \
           + lepton_R_Cd0_value * F_lepton_Cd_0_value + lepton_R_Cd1_value * F_lepton_Cd_1_value + lepton_R_Cd2_value * F_lepton_Cd_2_value \
           + lepton_R_Ce1_value * F_lepton_Ce_1_value \
           + lepton_R_Cf0_value * F_lepton_Cf_0_value + lepton_R_Cf1_value * F_lepton_Cf_1_value + lepton_R_Cf2_value * F_lepton_Cf_2_value;
}

// 计算 beta_lepton_i_L(i) 函数
std::complex<double> beta_lepton_i_L (int i, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    //调用一大堆函数
    std::complex<double> F_lepton_Ba_value = F_lepton_Ba(i);
    std::complex<double> F_lepton_Bb_value = F_lepton_Bb(i);
    std::complex<double> F_lepton_Bc_value = F_lepton_Bc(i);
    std::complex<double> F_lepton_Bd_value = F_lepton_Bd(i);
    std::complex<double> F_lepton_Ca_0_value = F_lepton_Ca_0(i);
    std::complex<double> F_lepton_Ca_1_value = F_lepton_Ca_1(i);
    std::complex<double> F_lepton_Ca_2_value = F_lepton_Ca_2(i);
    std::complex<double> F_lepton_Cb_0_value = F_lepton_Cb_0(i);
    std::complex<double> F_lepton_Cb_1_value = F_lepton_Cb_1(i);
    std::complex<double> F_lepton_Cb_2_value = F_lepton_Cb_2(i);
    std::complex<double> F_lepton_Cc_0_value = F_lepton_Cc_0(i);
    std::complex<double> F_lepton_Cc_1_value = F_lepton_Cc_1(i);
    std::complex<double> F_lepton_Cc_2_value = F_lepton_Cc_2(i);
    std::complex<double> F_lepton_Cd_0_value = F_lepton_Cd_0(i);
    std::complex<double> F_lepton_Cd_1_value = F_lepton_Cd_1(i);
    std::complex<double> F_lepton_Cd_2_value = F_lepton_Cd_2(i);
    std::complex<double> F_lepton_Ce_0_value = F_lepton_Ce_0(i);
    std::complex<double> F_lepton_Ce_1_value = F_lepton_Ce_1(i);
    std::complex<double> F_lepton_Ce_2_value = F_lepton_Ce_2(i);
    std::complex<double> F_lepton_Cf_0_value = F_lepton_Cf_0(i);
    std::complex<double> F_lepton_Cf_1_value = F_lepton_Cf_1(i);
    std::complex<double> F_lepton_Cf_2_value = F_lepton_Cf_2(i);
    //调用左手系数
    std::complex<double> lepton_L_Ba_value = lepton_L_Ba(i);
    std::complex<double> lepton_L_Bb_value = lepton_L_Bb(i);
    std::complex<double> lepton_L_Bc_value = lepton_L_Bc(i);
    std::complex<double> lepton_L_Bd_value = lepton_L_Bd(i);
    std::complex<double> lepton_L_Ca0_value = lepton_L_Ca0(i);
    std::complex<double> lepton_L_Ca1_value = lepton_L_Ca1(i);
    std::complex<double> lepton_L_Ca2_value = lepton_L_Ca2(i);
    std::complex<double> lepton_L_Cb0_value = lepton_L_Cb0(i);
    std::complex<double> lepton_L_Cb1_value = lepton_L_Cb1(i);
    std::complex<double> lepton_L_Cb2_value = lepton_L_Cb2(i);
    std::complex<double> lepton_L_Cc0_value = lepton_L_Cc0(i);
    std::complex<double> lepton_L_Cc1_value = lepton_L_Cc1(i);
    std::complex<double> lepton_L_Cc2_value = lepton_L_Cc2(i);
    std::complex<double> lepton_L_Cd0_value = lepton_L_Cd0(i);
    std::complex<double> lepton_L_Cd1_value = lepton_L_Cd1(i);
    std::complex<double> lepton_L_Cd2_value = lepton_L_Cd2(i);
    std::complex<double> lepton_L_Ce2_value = lepton_L_Ce2(i);
    std::complex<double> lepton_L_Cf0_value = lepton_L_Cf0(i);
    std::complex<double> lepton_L_Cf1_value = lepton_L_Cf1(i);
    std::complex<double> lepton_L_Cf2_value = lepton_L_Cf2(i);

    return lepton_L_Ba_value * F_lepton_Ba_value + lepton_L_Bb_value * F_lepton_Bb_value + lepton_L_Bc_value * F_lepton_Bc_value + lepton_L_Bd_value * F_lepton_Bd_value \
           + lepton_L_Ca0_value * F_lepton_Ca_0_value + lepton_L_Ca1_value * F_lepton_Ca_1_value + lepton_L_Ca2_value * F_lepton_Ca_2_value \
           + lepton_L_Cb0_value * F_lepton_Cb_0_value + lepton_L_Cb1_value * F_lepton_Cb_1_value + lepton_L_Cb2_value * F_lepton_Cb_2_value \
           + lepton_L_Cc0_value * F_lepton_Cc_0_value + lepton_L_Cc1_value * F_lepton_Cc_1_value + lepton_L_Cc2_value * F_lepton_Cc_2_value \
           + lepton_L_Cd0_value * F_lepton_Cd_0_value + lepton_L_Cd1_value * F_lepton_Cd_1_value + lepton_L_Cd2_value * F_lepton_Cd_2_value \
           + lepton_L_Ce2_value * F_lepton_Ce_2_value \
           + lepton_L_Cf0_value * F_lepton_Cf_0_value + lepton_L_Cf1_value * F_lepton_Cf_1_value + lepton_L_Cf2_value * F_lepton_Cf_2_value;
}


//NO.2衰变为轻子
std::complex<double> GAMMA_electron(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((GF * MH5 * set_M_lepton_i(1) * set_M_lepton_i(1) * (sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1)/  (MH5 * MH5))) *  (sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1)/  (MH5 * MH5))) * (sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1)/  (MH5 * MH5)))) / (4 * sqrt(2) * M_PI)) * norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * norm(alpha_lepton_i_R(1) / set_M_lepton_i(1));
}

std::complex<double> GAMMA_muon(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((GF * MH5 * set_M_lepton_i(2) * set_M_lepton_i(2) * (sqrt(1 - 4 * set_M_lepton_i(2) * set_M_lepton_i(2)/  (MH5 * MH5))) *  (sqrt(1 - 4 * set_M_lepton_i(2) * set_M_lepton_i(2)/  (MH5 * MH5))) * (sqrt(1 - 4 * set_M_lepton_i(2) * set_M_lepton_i(2)/  (MH5 * MH5)))) / (4 * sqrt(2) * M_PI)) * norm(alpha_lepton_i_R(2) / set_M_lepton_i(2)) * norm(alpha_lepton_i_R(2) / set_M_lepton_i(2));
}

std::complex<double> GAMMA_tau(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((GF * MH5 * set_M_lepton_i(3) * set_M_lepton_i(3) * (sqrt(1 - 4 * set_M_lepton_i(3) * set_M_lepton_i(3)/  (MH5 * MH5))) *  (sqrt(1 - 4 * set_M_lepton_i(3) * set_M_lepton_i(3)/  (MH5 * MH5))) * (sqrt(1 - 4 * set_M_lepton_i(3) * set_M_lepton_i(3)/  (MH5 * MH5)))) / (4 * sqrt(3) * M_PI)) * norm(alpha_lepton_i_R(3) / set_M_lepton_i(3)) * norm(alpha_lepton_i_R(3) / set_M_lepton_i(3));
}

//NO.3 mφ ≲ 2GeV 强子衰变为介子和介子
std::complex<double> GAMMA_paipai(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return ((3 * GF) / (16 * sqrt(2) * M_PI * MH5)) * sqrt(1 - 4 * Mpai * Mpai / (MH5 * MH5) ) * norm(((Mu * alpha_ij_R (1, 1) / Mu + Md * alpha_ij_R (1, 1) / Md)/(Mu + Md)) * Mpai * Mpai ) * norm(((Mu * alpha_ij_R (1, 1) / Mu + Md * alpha_ij_R (1, 1) / Md)/(Mu + Md)) * Mpai * Mpai );
}

std::complex<double> GAMMA_KK(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return (GF / (4 * sqrt(2) * M_PI * MH5)) * sqrt(1 - 4 * Mpai * Mpai / (MH5 * MH5) ) * norm(((Mu * alpha_ij_R (1, 1) / Mu + Md * alpha_ij_R (1, 1) / Md)/(Mu + Md)) * 0.5 * Mpai * Mpai + (alpha_ij_R (2, 2) / Ms) * (Mk * Mk - 0.5 * Mpai * Mpai)) * norm(((Mu * alpha_ij_R (1, 1) / Mu + Md * alpha_ij_R (1, 1) / Md)/(Mu + Md)) * 0.5 * Mpai * Mpai + (alpha_ij_R (2, 2) / Ms) * (Mk * Mk - 0.5 * Mpai * Mpai));
}

//N0.4 mφ ≳ 2GeV 的轻标量衰变为夸克
GAMMA_ll:GAMMA_ss:GAMMA_cc:GAMMA_bb=norm(alpha_lepton_i_R(i) / set_M_lepton_i(i)) * norm(alpha_lepton_i_R(i) / set_M_lepton_i(i)) * set_M_lepton_i(i) * set_M_lepton_i(i) * sqrt(1 - 4 * set_M_lepton_i(i) * set_M_lepton_i(i) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(i) * set_M_lepton_i(i) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(i) * set_M_lepton_i(i) / (MH5 * MH5)) \
                                   :3 * norm(alpha_ij_R (2, 2) / Ms) * norm(alpha_ij_R (2, 2) / Ms) * Ms * Ms * sqrt(1 - 4 * Mk * Mk / (MH5 * MH5)) * sqrt(1 - 4 * Mk * Mk / (MH5 * MH5)) * sqrt(1 - 4 * Mk * Mk / (MH5 * MH5)) \
                                   :3 * norm(alpha_ij_R (2, 2) / Mc) * norm(alpha_ij_R (2, 2) / Mc) * Mc * Mc * sqrt(1 - 4 * MD * MD / (MH5 * MH5)) * sqrt(1 - 4 * MD * MD / (MH5 * MH5)) * sqrt(1 - 4 * MD * MD / (MH5 * MH5)) \
                                   :3 * norm(alpha_ij_R (3, 3) / Mb) * norm(alpha_ij_R (3, 3) / Mb) * Mb * Mb * sqrt(1 - 4 * MB * MB / (MH5 * MH5)) * sqrt(1 - 4 * MB * MB / (MH5 * MH5)) * sqrt(1 - 4 * MB * MB / (MH5 * MH5))

