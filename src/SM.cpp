#include "SM.h"

#include <cmath>

SM::SM() {
    alpha = 1.0 / alpha_rev;
    EL = sqrt(4 * M_PI * alpha);
    vev = pow(M_SQRT2 * GF, -0.5);
    double A = sqrt(M_PI * alpha) * vev;
    thetaW = asin(2 * A / MZ) / 2.0;
    SW = sin(thetaW);
    SW2 = SW * SW;
    CW = cos(thetaW);
    CW2 = CW * CW;
    MW = A / SW;
    MW2 = MW * MW;
    g_weak = EL / SW;
    gp_hyper = EL / CW;
    yt = sqrt(2) * MT / vev;
}
