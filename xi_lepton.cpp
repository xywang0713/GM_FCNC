//衰变系数xi_phi_lepton
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
    
    return -(alpha * std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW * SW);
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
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW * SW);
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
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * (CW * CW * CW * CW -SW * SW * SW * SW) * set_yl_j(i)) / (8 * sqrt(6) * CW * CW * M_PI * SW * SW);
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
    
    return -(alpha * std::cos(thetaH) * std::sin(thetaH) * (CW * CW * CW * CW - SW * SW * SW * SW) * set_yl_j(i)) / (8 * sqrt(6) * CW * CW * M_PI * SW * SW);
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
    
    return (alpha * std::cos(thetaH) * std::sin(thetaH) * (CW * CW * (2 * MH5 * MH5 + MZ * MZ) + (2 * MH3 * MH3 - MZ * MZ) * SW * SW) * set_yl_j(i)) / (8 * sqrt(6)* CW * CW * M_PI * SW * SW);
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
    
    return (1/(16 * sqrt(3) * M_PI * M_PI * SW * SW)) * (std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i) * (2 * sqrt(2) * alpha * MH5 * MH5 * M_PI + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 - 6 * sqrt(2) * M2 + (-2 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum))* set_yl_j(i)));
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
    
    return (alpha * std::cos(thetaH) * MH5 * MH5 * std::sin(thetaH) * (3 * CW * CW  - SW * SW ) * set_yl_j(i)) / (8 * sqrt(6) * CW * CW * M_PI * SW * SW);
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
    
    return (alpha * std::cos(thetaH) * (set_M_lepton_i(i) * set_M_lepton_i(i) - MW * MW) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW * SW);
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
    
    return ((std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i))  / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (2 * sqrt(2) * alpha * (set_M_lepton_i(i) * set_M_lepton_i(i) + 2 * MH5 * MH5) * M_PI + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 - 6 * sqrt(2) * M2 + (-2 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum )) * set_yl_j(i));
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
    
    return ((std::cos(thetaH) * std::sin(thetaH) * set_yl_j(i))  / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (2 * sqrt(2) * alpha * (set_M_lepton_i(i) * set_M_lepton_i(i) - MH5 * MH5) * M_PI + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * (-2 * std::cos(thetaH) * std::cos(thetaH) * LAMBDA_5() * vacuum + std::sin(thetaH) * (sqrt(2) * M1 - 6 * sqrt(2) * M2 + (-2 * LAMBDA_3() + LAMBDA_5()) * std::sin(thetaH) * vacuum)) * set_yl_j(i));
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
    
    return (alpha * std::cos(thetaH) * (-set_M_lepton_i(i) * set_M_lepton_i(i) + MW * MW) * std::sin(thetaH) * set_yl_j(i)) / (8 * sqrt(6) * M_PI * SW * SW);
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
    
    return ((std::cos(thetaH) * set_yl_j(i)) / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (-4 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * (2 * std::sin(thetaH) * std::sin(thetaH) * (3 *sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std:cos(thetaH) * (sqrt(2) * M1 + 3 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i));
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
    
    return (1/(16 * sqrt(3) * M_PI * SW * SW * SW * SW)) * alpha * std::sin(thetaH) * (8 * alpha * set_M_lepton_i(i) * M_PI * vacuum + sqrt(2) * std::cos(thetaH) * (-2 * set_M_lepton_i(i) * set_M_lepton_i(i) + MH5 * MH5) * SW * SW * set_yl_j(i));
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
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() - 2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
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
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() - 2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
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
    
    return -((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() - 2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
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
    
    return ((set_M_lepton_i(i) * std::sin(thetaH) * std::sin(thetaH)) / (32 * sqrt(3) * M_PI * M_PI)) * (2 * std::cos(thetaH) * std::cos(thetaH) * (3 * sqrt(2) * M2 + (LAMBDA_3() - 2 * LAMBDA_5()) * std::sin(thetaH) * vacuum) + std::sin(thetaH) * std::sin(thetaH) * (sqrt(2) * M1 - LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i) * set_yl_j(i);
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
    
    return -(1/(16 * sqrt(3) * CW * CW * CW * CW * M_PI * SW * SW * SW * SW)) * (alpha * std::sin(thetaH) * (8 * alpha * set_M_lepton_i(i) * M_PI * (CW * CW - 3 * SW * SW) * vacuum + sqrt(2) * std::cos(thetaH) * CW * CW * MH5 * MH5 * SW * SW * (3 * CW * CW - 5 * SW * SW) * set_yl_j(i)));
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
    
    return (1/(32 * sqrt(3) * CW * CW * CW * CW)) * (64 * alpha * alpha * set_M_lepton_i(i) * std::sin(thetaH) * vacuum + ((std::cos(thetaH) * CW * CW * set_yl_j(i)) / (M_PI * M_PI * SW * SW)) * (4 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) * (- CW * CW * CW * CW + SW * SW * SW *SW) + std::cos(thetaH) * CW * CW * set_M_lepton_i(i) * SW * SW * (2 * std::sin(thetaH) * std::sin(thetaH) * (3 * sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) *(sqrt(2) * M1 + 3 * LAMBDA_5() * std ::sin(thetaH) * vacuum)) * set_yl_j(i)));
}

//lepton left
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
    
    return ((std::cos(thetaH) * set_yl_j(i)) / (32 * sqrt(3) * M_PI * M_PI * SW * SW)) * (-2 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) + std::cos(thetaH) * set_M_lepton_i(i) * SW * SW * (2 * std::sin(thetaH) * std::sin(thetaH) * (3 * sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) *vaccum) + std::cos(thetaH) * std::cos(thetaH) * (sqrt(2) * M1 + 3 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i));
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
    
    return ((alpha * std::sin(thetaH)) /(16 * sqrt(3) * CW * CW * CW * CW * M_PI * SW * SW *SW * SW)) * (8 * alpha * set_M_lepton_i(i) * M_PI * (CW * CW - 3 *SW * SW) * vacuum + sqrt(2) * std::cos(thetaH) * CW * CW * MH5 * MH5 * SW *SW *(CW * CW - 7 * SW * SW)* set_yl_j(i));
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
    
    return (1 / (32 * sqrt(3) * CW * CW * CW *CW * M_PI * M_PI * SW * SW * SW * SW)) * (16 * alpha * alpha * set_M_lepton_i(i) * M_PI * M_PI * std::sin(thetaH) * (CW * CW * CW * CW - SW * SW * SW * SW) * (CW * CW * CW * CW - SW * SW * SW * SW) * vacuum + std::cos(thetaH) * CW * CW * SW * SW * SW *SW * set_yl_j(i) * (-8 * sqrt(2) * alpha * MH5 * MH5 * M_PI * std::sin(thetaH) + std::cos(thetaH) * CW * CW * set_M_lepton_i(i) * (2 * std::sin(thetaH) * std::sin(thetaH) *(3 * sqrt(2) * M2 + LAMBDA_3() * std::sin(thetaH) * vacuum) + std::cos(thetaH) * std::cos(thetaH) * (sqrt(2) * M1 + 3 * LAMBDA_5() * std::sin(thetaH) * vacuum)) * set_yl_j(i)));
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


