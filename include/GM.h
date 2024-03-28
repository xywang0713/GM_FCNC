#ifndef GM_FCNC_GM_H_
#define GM_FCNC_GM_H_

#include "SM.h"

class GM : public SM {
public:
    GM() { set_input3(288.268237, 304.221605, 339.748616, 0.194487374, -0.303281383, 100, 100); };

    // * We are following GMCALC Input Convention
    bool set_input3(double MHH, double MH3, double MH5, double sH, double sa, double M1, double M2);

    double get_mu22() const { return m_mu22; }
    double get_mu32() const { return m_mu32; }
    double get_lam1() const { return m_lam1; }
    double get_lam2() const { return m_lam2; }
    double get_lam3() const { return m_lam3; }
    double get_lam4() const { return m_lam4; }
    double get_lam5() const { return m_lam5; }
    double get_M1() const { return m_M1; }
    double get_M2() const { return m_M2; }

    double get_MHL() const { return m_MHL; }
    double get_MHH() const { return m_MHH; }
    double get_MH3() const { return m_MH3; }
    double get_MH5() const { return m_MH5; }

    double get_thetaH() const { return m_thetaH; }
    double get_alpha() const { return m_alpha; }
    double get_sH() const { return m_sH; }
    double get_cH() const { return m_cH; }
    double get_sa() const { return m_sa; }
    double get_ca() const { return m_ca; }

    double get_vphi() const { return m_vphi; }
    double get_vchi() const { return m_vchi; }

    double get_Yu(Generation_t);
    double get_Yd(Generation_t);
    double get_Yl(Generation_t);

private:
    double m_mu22, m_mu32, m_lam1, m_lam2, m_lam3, m_lam4, m_lam5, m_M1, m_M2;
    double m_MHL, m_MHH, m_MH3, m_MH5;
    double m_MHL2, m_MHH2, m_MH32, m_MH52;
    double m_thetaH, m_alpha;
    double m_sH, m_cH, m_sa, m_ca;
    double m_vphi, m_vchi;
};

#endif  // GM_FCNC_GM_H_
