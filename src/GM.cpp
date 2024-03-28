#include "GM.h"

bool GM::set_input3(double MHH, double MH3, double MH5, double sH, double sa, double M1, double M2) {
    // * theta_H in [0, pi/2]
    // * alpha in [-pi/2,pi/2]
    if (sH > 1.0 || sH < 0.0) return false;
    if (sa > 1.0 || sa < -1.0) return false;
    m_sH = sH;
    m_thetaH = asin(m_sH);
    m_cH = cos(m_thetaH);
    m_sa = sa;
    m_alpha = asin(sa);
    m_ca = cos(m_alpha);

    m_MHL = MHL;
    m_MHH = MHH;
    m_MH3 = MH3;
    m_MH5 = MH5;
    m_M1 = M1;
    m_M2 = M2;

    m_MHL2 = m_MHL * m_MHL;
    m_MHH2 = m_MHH * m_MHH;
    m_MH32 = m_MH3 * m_MH3;
    m_MH52 = m_MH5 * m_MH5;

    m_vchi = m_vev * m_sH / sqrt(8.0);
    m_vphi = m_vev * m_cH;

    double vphi2 = m_vphi * m_vphi;
    double vchi2 = m_vchi * m_vchi;
    double v2 = m_vev * m_vev;
    double sa2 = m_sa * m_sa;
    double ca2 = m_ca * m_ca;
    m_lam1 = (m_MHL2 * ca2 + m_MHH2 * sa2) / 8.0 / vphi2;
    m_lam2 = m_MH32 / v2 + m_sa * m_ca / 4.0 / sqrt(3.0) * (m_MHH2 - m_MHL2) / m_vphi / m_vchi - m_M1 / 8.0 / m_vchi;
    m_lam3 = m_MH52 / 8.0 / vchi2 - 3.0 * vphi2 / 8.0 / vchi2 * m_MH32 / v2 + vphi2 * m_M1 / 16.0 / pow(m_vchi, 3) -
             3.0 * m_M2 / 2.0 / m_vchi;
    m_lam4 = vphi2 / 8.0 / vchi2 * m_MH32 / v2 + (m_MHL2 * sa2 + m_MHH2 * ca2 - m_MH52) / 24.0 / vchi2 -
             vphi2 * m_M1 / 32.0 / pow(m_vchi, 3) + 3.0 * m_M2 / 4.0 / m_vchi;
    m_lam5 = 2.0 * m_MH32 / v2 - m_M1 / 2.0 / m_vchi;

    m_mu22 = -4.0 * m_lam1 * vphi2 + 3.0 * (m_lam5 - 2.0 * m_lam2) * vchi2 + 3.0 / 2.0 * m_M1 * m_vchi;
    m_mu32 = (m_lam5 - 2.0 * m_lam2) * vphi2 - 4.0 * (m_lam3 + 3.0 * m_lam4) * vchi2 + m_M1 * vphi2 / 4.0 / m_vchi +
             6.0 * m_M2 * m_vchi;

    return true;
}
