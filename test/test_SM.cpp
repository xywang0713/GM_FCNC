#include <iostream>

#include "SM.h"

using namespace std;
int main(int argc, char const *argv[]) {
    SM mod;
    cout << "The Parameters in the SM: " << endl;
    cout << "MZ: " << mod.get_MZ() << endl;
    cout << "MW: " << mod.get_MW() << endl;
    cout << "CW: " << mod.get_CW() << endl;
    cout << "SW: " << mod.get_SW() << endl;
    cout << "MT: " << mod.get_MT() << endl;
    cout << "g: " << mod.get_gweak() << endl;
    cout << "g': " << mod.get_gphyper() << endl;
    cout << "vev: " << mod.get_vev() << endl;
    cout << "mu: " << mod.get_m_up(0) << endl;
    cout << "md: " << mod.get_m_down(0) << endl;
    cout << "mc: " << mod.get_m_up(1) << endl;
    cout << "ms: " << mod.get_m_down(1) << endl;
    cout << "mt: " << mod.get_m_up(2) << endl;
    cout << "mb: " << mod.get_m_down(2) << endl;
    cout << "CKM Matrix: " << endl;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            cout << "  " << mod.get_vckm(i, j);
        }
        cout << endl;
    }

    return 0;
}
