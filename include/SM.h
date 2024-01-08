#ifndef GM_FCNC_SM_H_
#define GM_FCNC_SM_H_

#include "Constants.h"

class SM {
public:
    SM();

    double get_MZ() const { return MZ; }
    double get_MZ2() const { return MZ2; }
    double get_MW() const { return MW; }
    double get_MW2() const { return MW2; }
    double get_CW() const { return CW; }
    double get_CW2() const { return CW2; }
    double get_SW() const { return SW; }
    double get_SW2() const { return SW2; }
    double get_MT() const { return MT; }
    double get_MT2() const { return MT2; }
    double get_gweak() const { return g_weak; }
    double get_gphyper() const { return gp_hyper; }
    double get_vev() const { return vev; }

protected:
    double MW;
    double MW2;
    double thetaW;
    double SW;
    double SW2;
    double CW;
    double CW2;
    double alpha;
    double vev;       // * (sqrt(2)GF)^(-0.5)
    double EL;        // * EL = sqrt(4*Pi*alpha)
    double g_weak;    // * g = EL/SW
    double gp_hyper;  // * gp = EL/CW
    double yt;        // * sqrt(2)mt
};

#endif  // GM_FCNC_SM_H_
