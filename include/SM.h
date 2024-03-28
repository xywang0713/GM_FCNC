#ifndef GM_FCNC_SM_H_
#define GM_FCNC_SM_H_

#include <complex>

#include "Constants.h"

enum Generation_t { FIRST = 0, SECOND = 1, THIRD = 2 };

class SM {
public:
    SM();

    double get_MZ() const { return MZ; }
    double get_MZ2() const { return MZ2; }
    double get_MW() const { return m_MW; }
    double get_MW2() const { return m_MW2; }
    double get_CW() const { return m_CW; }
    double get_CW2() const { return m_CW2; }
    double get_SW() const { return m_SW; }
    double get_SW2() const { return m_SW2; }
    double get_MT() const { return MT; }
    double get_MT2() const { return MT2; }
    double get_gweak() const { return m_g_weak; }
    double get_gphyper() const { return m_gp_hyper; }
    double get_vev() const { return m_vev; }

    double get_m_up(Generation_t i) const { return m_m_ups[i]; }
    double get_m_down(Generation_t i) const { return m_m_downs[i]; }
    std::complex<double> get_vckm(Generation_t i, Generation_t j) const { return m_VCKM[i][j]; }
    std::complex<double> get_vckmc(Generation_t i, Generation_t j) const { return std::conj(m_VCKM[i][j]); }

protected:
    double m_MW;
    double m_MW2;
    double m_thetaW;
    double m_SW;
    double m_SW2;
    double m_CW;
    double m_CW2;
    double m_alpha;
    double m_vev;       // * (sqrt(2)GF)^(-0.5)
    double m_EL;        // * EL = sqrt(4*Pi*alpha)
    double m_g_weak;    // * g = EL/SW
    double m_gp_hyper;  // * gp = EL/CW
    double m_yt;        // * sqrt(2)mt/vev

    double m_m_ups[3];
    double m_m_downs[3];
    double m_m_leps[3];
    std::complex<double> m_VCKM[3][3];
};

#endif  // GM_FCNC_SM_H_
