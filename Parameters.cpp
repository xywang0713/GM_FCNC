//定义变量
double alpha = 1/137;
double SW =0.2223;
double vacuum246.0;
double MW = 80.399;
double Mb = 4.18; //b quark
double Ms = 93.4e-3; //s quark
double Md = 4.67e-3; //d quark
double Mu = 2.16e-3; // u quark
double Mc = 1.27; //c quark
double Mk = 0.493667;  //K^{\pm}
double Mpai=0.13957039; //pi^{\pm}
double Mmu =0.1056583755; //mu
double MD=1.86966; //D^{+}
double MDs=1.96835; //Ds^{+}
double MB=5.27934; //B
double Gamma_k=5.3167366721e-17; ///K^{\pm} width
double MZ=91.1876;
double alpha_ew=1/137;
double MD0=1864.84e-3;//D^0

//写class
class GM_model{
public:
    double MH5;
    double MH3;
    double sH;
    double thetaH;
    double M1;
    double M2;


    void read_Data(){                           // 定义一个名为read_Data的void类型的函数
        std::ifstream file("home/test.csv")     // 用ifstream类打开home/test.csv
        if (file.is_open()){                    // 如果文件成功打开       
        std::string line;                       // 定义一个名为line的string类型变量
        int lineCount = 0;                      // 定义一个名为lineCount的int类型变量并赋值为0
        while (std::getline(file,line)){        // 当从file中读取一行内容存入line中时
            if (lineCount == 1){                // 如果行数为1
                std::istringstream iss(line);   // 使用istringstream类根据line创建一个名为iss的对象
                std:vector<double> values;      // 定义一个名为values的vector类型变量
                double value;                   // 定义一个名为value的double类型变量
                while (iss >> value) {          // 当从iss中读取值并成功时
                    values.push_back(value);    // 将value添加到values的末尾
                }
                MH5 = values[0];
                sH = values[1];
                thetaH = std::asin(sH); 
                M1 = values[6];
                M2 = values[7];
                MH3 = values[9];                // 将values中的值赋给变量
            }
            lineCount++;                        // 行数加1
        }
        file.close();
        }       
    }


    // 计算 LAMBDA_3()  函数
    std::complex<double> LAMBDA_3() {
         return (std::cos(thetaH) * std::cos(thetaH) * (-3.0 * std::sin(thetaH) * MH3 * MH3 + std::sqrt(2.0) * M1 * vacuum) + std::sin(thetaH) * (MH5 * MH5 - 3.0 * std::sqrt(2.0) * M2 * std::sin(thetaH) * vacuum)) / (vacuum * std::sin(thetaH) * std::sin(thetaH) * std::sin(thetaH));
    }

    // 计算 LAMBDA_5() 函数
    std::complex<double> LAMBDA_5() {
         return (2.0 * std::sin(thetaH) * MH3 * MH3 - std::sqrt(2.0) * M1 * vacuum) / (vacuum * std::sin(thetaH));
    }


//Potential里的参数
public:
    double LAMBDA_1;
    double LAMBDA_2;
    double LAMBDA_3();
    double LAMBDA_4;
    double LAMBDA_5();
    double MU12;
    double MU22;
    double M1;
    double M2;
}

// 计算 yd(z) 的函数
std::complex<double> set_yd_z(int z, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;

    std::complex<double> yd_z;
    if (z == 1) {
        yd_z = std::sqrt(2.0) * 5.04 * std::pow(10.0, -3) / (vacuum * std::cos(thetaH));
    } else if (z == 2) {
        yd_z = std::sqrt(2.0) * 0.101 / (vacuum * std::cos(thetaH));
    } else if (z == 3) {
        yd_z = std::sqrt(2.0) * 4.7 / (vacuum * std::cos(thetaH));
    } else {
        std::cerr << "z value is wrong!" << std::endl;
    }
    return yd_z;
}


// 计算 yu(y) 的函数
std::complex<double> set_yu_y(int y, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    std::complex<double> yu_y;
    if (y == 1) {
        yu_y = std::sqrt(2.0) * 2.55 * pow(10.0, -3) / (vacuum * std::cos(thetaH));
    } else if (y == 2) {
        yu_y = std::sqrt(2.0) * 1.27 / (vacuum * std::cos(thetaH));
    } else if (y == 3) {
        yu_y = std::sqrt(2) * 172.0 / (vacuum * std::cos(thetaH));
    } else {
        std::cerr << "y value is wrong!" << std::endl;
    }
    return yu_y;
}

// 计算 Md(x) 的函数
double set_Md_x(int x, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    double Md_x = 0.0;  // Initialize Md_x with a default value
    if (x == 1) {
        Md_x = 5.04 * pow(10.0, -3);
    } else if (x == 2) {
        Md_x = 0.101;
    } else if (x == 3) {
        Md_x = 4.7;
    } else {
        std::cerr << "x value is wrong!" << std::endl;
    }
    return Md_x;
}

// 计算 Mu(w) 的函数
double set_Mu_w(int w, GM_model&m) {

    thetaH = m.thetaH;
    MH3 = m.MH3;
    MH5 = m.MH5;
    LAMBDA_3() = m.LAMBDA_3();
    LAMBDA_5() = m.LAMBDA_5();
    M1 = m.M1;
    M2 = m.M2;
    
    double Mu_w = 0.0;  // Initialize Mu_w with a default value
    if (w == 1) {
        Mu_w = 2.55 * pow(10.0, -3);
    } else if (w == 2) {
        Mu_w = 1.27;
    } else if (w == 3) {
        Mu_w = 172.0;
    } else {
        std::cerr << "w value is wrong!" << std::endl;
    }
    return Mu_w;
}

//lepton用到的质量以及yl
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


// 计算 CKM矩阵V[u,t] 的函数
 std::complex<double> VCKM(int u, int t) {
     if (u == 1 && t == 1) {
         return 0.974352;
     } else if (u == 1 && t == 2) {
         return 0.224998;
     } else if (u == 1 && t == 3) {
         return std::complex<double>(0.0015275, -0.00335899);
     } else if (u == 2 && t == 1) {
         return std::complex<double>(-0.224865, -0.000136871);
     } else if (u == 2 && t == 2) {
         return std::complex<double>(0.973492, -0.0000316065 );
     } else if (u == 2 && t == 3) {
         return 0.0418197;
     } else if (u == 3 && t == 1 ) {
         return std::complex<double>(0.00792247, -0.00327); 
     } else if (u == 3 && t == 2 ) {
         return std::complex<double>(-0.0410911, -0.000755113);  
     } else if (u == 3 && t == 3 ) { 
         return 0.999118;
     } else {
         std::cerr << "Invalid indices u, t" << std::endl;
         return 0;
     }
 }

// 计算 CKMC矩阵VC[r,s] 的函数
 std::complex<double> VCKMC(int r, int s) {
     if (r == 1 && s == 1) {
         return 0.974352;
     } else if (r == 1 && s == 2) {
         return 0.224998;
     } else if (r == 1 && s == 3) {
         return std::complex<double>(0.0015275 , 0.00335899);
     } else if (r == 2 && s == 1) {
         return std::complex<double>(-0.224865 , 0.000136871);
     } else if (r == 2 && s == 2) {
         return std::complex<double>(0.973492 , 0.0000316065);
     } else if (r == 2 && s == 3) {
         return 0.0418197;
     } else if (r == 3 && s == 1 ) {
         return std::complex<double>(0.00792247 , 0.00327); 
     } else if (r == 3 && s == 2 ) {
         return std::complex<double>(-0.0410911 , 0.000755113); 
     } else if (r == 3 && s == 3 ) { 
         return 0.999118;
     } else {
         std::cerr << "Invalid indices r, s" << std::endl;
         return 0;
     }
 }

