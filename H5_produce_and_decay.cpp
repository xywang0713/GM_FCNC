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

//下面涉及到的质量Mb，Mc，
//B介子衰变产生函数Br(B -> Xs phi)/Br(B -> Xc e nu) 
//注意：其中最后norm(xi_ij_R (int i, int j, GM_model&m) / VCKM[c,b])中xi是xi_bs
std::complex<double> Br_B_Meson (int i, int j, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    return (12.0 * M_PI * M_PI * vacuum * vacuum / (Mb * Mb)) * (1.0 - (MH5 * MH5) / (Mb * Mb)) * (1.0 / (((1.0 - 8.0 * ((Mc * Mc) / (Mb * Mb)) + ((Mc * Mc ) / (Mb * Mb)) * ((Mc * Mc )/ (Mb * Mb)))) * (1.0 - ((Mc * Mc) / (Mb * Mb)) * ((Mc * Mc ) / (Mb * Mb))) - 12.0 * log(((Mc * Mc ) / (Mb * Mb)) * ((Mc * Mc ) / (Mb * Mb))))) * std::norm(down_alpha_ij_R (3, 2) /( Ms *VCKM[2,3])) * std::norm(down_alpha_ij_R (3, 2) /( Ms *VCKM[2,3])) ;
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

    return (1.0 / Gamma_k) * (2.0 * P_0_PHI / Mk) * ( std::norm(sqrt(GF) * 4.0 * cbrt(2.0) * (EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW)) * (7.0 * 3.1 * pow (10,-7) *(Mk* Mk +Mpai * Mpai-MH5 * MH5) / 18.0) + (down_alpha_ij_R (1, 2, GM_model&m) * Ms *0.96 /(2.0*vacuum)) *((Mk* Mk -Mpai * Mpai) / (Ms - Md))) * std::norm(sqrt(GF) * 4.0 * cbrt(2.0) * (EL * EL * vacuum * vacuum * std::sin(thetaH) / (4.0 * sqrt(3.0) * SW * SW * MW *MW)) * (7.0 * 3.1 * pow (10,-7) *(Mk* Mk +Mpai * Mpai-MH5 * MH5) / 18.0) + (down_alpha_ij_R (1, 2, GM_model&m) * Ms *0.96 /(2.0*vacuum)) *((Mk* Mk -Mpai * Mpai) / (Ms - Md))) / (16.0 * M_PI * Mk))
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

    return ((EL * EL * std::sin(thetaH) * vacuum * vacuum) / (4.0 * sqrt(3.0) * SW * SW * MW * MW)) * A_1_function_value;
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

    return (GF / (4 * sqrt(2) * M_PI * MH5)) * sqrt(1 - 4 * Mpai * Mpai / (MH5 * MH5) ) * norm(((Mu * up_alpha_ij_R (1, 1) / Mu + Md * down_alpha_ij_R (1, 1) / Md)/(Mu + Md)) * 0.5 * Mpai * Mpai + (down_alpha_ij_R (2, 2) / Ms) * (Mk * Mk - 0.5 * Mpai * Mpai)) * norm(((Mu * up_alpha_ij_R (1, 1) / Mu + Md * down_alpha_ij_R (1, 1) / Md)/(Mu + Md)) * 0.5 * Mpai * Mpai + (down_alpha_ij_R (2, 2) / Ms) * (Mk * Mk - 0.5 * Mpai * Mpai));
}

//N0.4 mφ ≳ 2GeV 的轻标量衰变为夸克
std::complex<double> GAMMA_ss(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return (3 * norm(down_alpha_ij_R (2, 2) / Ms) * down_norm(alpha_ij_R (2, 2) / Ms) * Ms * Ms * sqrt(1 - 4 * Mk * Mk / (MH5 * MH5)) * sqrt(1 - 4 * Mk * Mk / (MH5 * MH5)) * sqrt(1 - 4 * Mk * Mk / (MH5 * MH5))) * GAMMA_electron() / (norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * set_M_lepton_i(1) * set_M_lepton_i(1) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)));
}

std::complex<double> GAMMA_cc(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return (3 * norm(up_alpha_ij_R (2, 2) / Mc) * up_norm(alpha_ij_R (2, 2) / Mc) * Mc * Mc * sqrt(1 - 4 * MD0 * MD0 / (MH5 * MH5)) * sqrt(1 - 4 * MD0 * MD0 / (MH5 * MH5)) * sqrt(1 - 4 * MD0 * MD0 / (MH5 * MH5))) * GAMMA_electron() / (norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * set_M_lepton_i(1) * set_M_lepton_i(1) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)));
}

//N0.4 mφ ≳ 2GeV 的轻标量衰变为夸克
std::complex<double> GAMMA_bb(GM_model&m) {
    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    return (3 * norm(down_alpha_ij_R (3, 3) / Mb) * down_norm(alpha_ij_R (3, 3) / Mb) * Mb * Mb * sqrt(1 - 4 * MB * MB / (MH5 * MH5)) * sqrt(1 - 4 * MB * MB / (MH5 * MH5)) * sqrt(1 - 4 * MB * MB / (MH5 * MH5))) * GAMMA_electron() / (norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * norm(alpha_lepton_i_R(1) / set_M_lepton_i(1)) * set_M_lepton_i(1) * set_M_lepton_i(1) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)) * sqrt(1 - 4 * set_M_lepton_i(1) * set_M_lepton_i(1) / (MH5 * MH5)));
}
