#include "SM.h"

#include <cmath>

SM::SM()
    : m_m_ups{2.55e-3, 1.27, MT},
      m_m_downs{5.04e-3, 0.101, 4.7},
      m_VCKM{{{0.974352, 0}, {0.224998, 0}, {0.0015275, -0.00335899}},
             {{-0.224865, -0.000136871}, {0.973492, -0.0000316065}, {0.0418197, 0}},
             {{0.00792247, -0.00327}, {-0.0410911, -0.000755113}, {0.999118, 0}}} {
    m_alpha = 1.0 / alpha_rev;
    m_EL = sqrt(4 * M_PI * m_alpha);
    m_vev = pow(M_SQRT2 * GF, -0.5);
    double A = sqrt(M_PI * m_alpha) * m_vev;
    m_thetaW = asin(2 * A / MZ) / 2.0;
    m_SW = sin(m_thetaW);
    m_SW2 = m_SW * m_SW;
    m_CW = cos(m_thetaW);
    m_CW2 = m_CW * m_CW;
    m_MW = A / m_SW;
    m_MW2 = m_MW * m_MW;
    m_g_weak = m_EL / m_SW;
    m_gp_hyper = m_EL / m_CW;
    m_yt = M_SQRT2 * MT / m_vev;
}
