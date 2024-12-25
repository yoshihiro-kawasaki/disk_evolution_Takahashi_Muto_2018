#include "bonner_ebert_sphere.hpp"
#include "constants.hpp"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <cmath>


BonnerEbertSphere::BonnerEbertSphere(SimulationData *pdata)
    : pdata_(pdata)
{
    double mg = cst::GAS_MOLECULAR_MASS;
    number_of_cloud_shells_          = pdata_->cloud_.nshells_;
    central_number_density_of_cloud_ = pdata_->cloud_.nc_;
    central_mass_density_of_cloud_   = mg * central_number_density_of_cloud_;
    cloud_temperature_               = pdata_->cloud_.T_;
    cloud_sound_speed_               = std::sqrt(cst::BOLTZMANN_CONSTANT*cloud_temperature_/mg);
    cloud_angular_velocity_          = pdata_->cloud_.Omega_c_;
    density_factor_                  = pdata_->cloud_.f_;
    total_cloud_mass_                = pdata_->cloud_.Mc_;
    initial_protostellar_mass_       = pdata_->star_.mass_init_;

    // pdata_->cloud_.cs_ = cloud_sound_speed_;
    // pdata_->cloud_.rhoc_= central_mass_density_of_cloud_;
}


BonnerEbertSphere::~BonnerEbertSphere()
{
    // delete arrays
    array::Delete1dArray<double>(rn_);
    array::Delete1dArray<double>(psi_);
    array::Delete1dArray<double>(dpsidrn_);
    array::Delete1dArray<double>(rho_bonner_);
    array::Delete1dArray<double>(r_bonner_);
    array::Delete1dArray<double>(mass_bonner_);   
}


void BonnerEbertSphere::SetBonnerEbertSphere()
{
    int ibmax = pdata_->cloud_.nshells_;
    // allocate arrays
    rn_          = array::Allocate1dArray<double>(ibmax);
    psi_         = array::Allocate1dArray<double>(ibmax);
    dpsidrn_     = array::Allocate1dArray<double>(ibmax);
    rho_bonner_  = array::Allocate1dArray<double>(ibmax);
    r_bonner_    = array::Allocate1dArray<double>(ibmax);
    mass_bonner_ = array::Allocate1dArray<double>(ibmax);

    rout_bonner_   = 0.0;
    rn_out_bonner_ = 0.0;
    ib_con_        = 0;

    double rn_in  = 1.0e-8; // inner normalized cloud radius
    double rn_out = 2.0e1;  // outer normalized cloud radius
    double dlnrn  = (std::log(rn_out/rn_in))/((double)ibmax);
    for (int i = 0; i < ibmax; i++) {
        rn_[i] = rn_in * std::exp(double(i) * dlnrn);
    }

    // ## calcurate bonner ebert sphere profile ## //
    double xn, h;
    int ibmax_m = ibmax - 1;
    double y[2];
    
    // boundary condition
    psi_[0]     = 0.0;
    dpsidrn_[0] = 0.0;
    y[0]        = psi_[0];
    y[1]        = dpsidrn_[0];

    // integration
    for (int i = 0; i < ibmax_m; ++i) {
        xn = rn_[i];
        h  = rn_[i+1] - rn_[i];
        SolveLaneEmdenEquation(2, xn, h, y);
        psi_[i+1]     = y[0];
        dpsidrn_[i+1] = y[1];
    }

    // calculate rho, r, mass of BE sphere
    double r_nor    = std::pow(4.0*M_PI*central_mass_density_of_cloud_*cst::GRAVITATIONAL_CONSTANT/SQR(cloud_sound_speed_), -0.5);
    double mass_nor = (CUB(cloud_sound_speed_) / cst::GRAVITATIONAL_CONSTANT) / std::sqrt(4.0*M_PI*cst::GRAVITATIONAL_CONSTANT*central_mass_density_of_cloud_);

    for (int i = 0; i < ibmax; ++i) {
        rho_bonner_[i]  = density_factor_ * central_mass_density_of_cloud_ * std::exp(-psi_[i]);
        r_bonner_[i]    = r_nor * rn_[i];
        mass_bonner_[i] = density_factor_ * mass_nor * SQR(rn_[i]) * dpsidrn_[i];
    }

    // calculate rout_bonner
    for (int i = 0; i < ibmax; ++i) {
        if (mass_bonner_[i] > total_cloud_mass_) {
            rout_bonner_   = Utils::LinearInterPolation(total_cloud_mass_, mass_bonner_[i-1], mass_bonner_[i], r_bonner_[i-1], r_bonner_[i]);
            rn_out_bonner_ = rout_bonner_ / r_nor;
            std::cout << "cloud radius : " << std::scientific << (rout_bonner_/cst::ASTRONOMICAL_UNIT) << " au" << std::endl;
            break;
        }
    }

    if (mass_bonner_[ibmax_m] < total_cloud_mass_) {
        std::cerr << "error : mass_bonner[ibmax] < MASS_BONNER" << std::endl;
        std::exit(1);
    }

    return;
}



double BonnerEbertSphere::CalculateInitialRinit()
/*
    計算開始時の最初に中心(円盤)に落下してくる分子雲コアのshellの位置を求める
    その位置は、初期の中心星質量よりコアの質量が大きくなる位置となる。
*/
{
    int i1 = 0, i2 = 0, ibmax = number_of_cloud_shells_;

    for (int i = 0; i < ibmax; ++i) {
        if (mass_bonner_[i] >= initial_protostellar_mass_) {
            i1 = i - 1;
            i2 = i;
            break;
            if (i2 == 0) {
                std::cerr << "BonnerEbertSphere::CalculateRinit error : i2 == 0" << std::endl;
                std::exit(1);
            }
        }
    }

    double r_init = Utils::LinearInterPolation(initial_protostellar_mass_, mass_bonner_[i1], mass_bonner_[i2], r_bonner_[i1], r_bonner_[i2]);

    return r_init;
}


void BonnerEbertSphere::CalculateMdotInfall(const double r_init, double &drinitdt, double &mdot_inf)
{
    if (r_init >= rout_bonner_) {
        drinitdt = 0.0;
        mdot_inf = 0.0;
        return;
    }

    int i1 = 0, i2 = 0, ibmax = number_of_cloud_shells_;
    // 分子雲コアの位置r_initにおける質量密度rho_initを求める。
    for (int i = ib_con_; i < ibmax; ++i) {
        if (r_bonner_[i] >= r_init) {
            i1 = i - 1;
            i2 = i;
            ib_con_ = i;
            break;
            if (i2 == 0) {
                std::cerr << "BonnerEbertSphere::CalculateRhoinit error : i2 == 0" << std::endl;
                std::exit(1);
            }
        }
    }

    double rho_init = Utils::LinearInterPolation(r_init, r_bonner_[i1], r_bonner_[i2], rho_bonner_[i1], rho_bonner_[i2]);

    // ref. Takahashi et al. 2013
    // 分子雲コアにおいて、初期にr_initの位置にあるshellが中心へ落下するまでの時間t_infはr_initの関数で与えられる、t_inf = tff(r_init)。
    // 逆に、r_initはtffの関数で書ける、r_init = r_init(t_inf)。
    // これの微分係数、dr_init/dt_infを求める。

    double mr_init  = Utils::LinearInterPolation(r_init, r_bonner_[i1], r_bonner_[i2], mass_bonner_[i1], mass_bonner_[i2]);
    double dmrdr    = 4.0*M_PI*SQR(r_init)*rho_init;
    double dtdrinit = C_TINF * sqrt(r_init/(2.0*cst::GRAVITATIONAL_CONSTANT*mr_init)) * (1.5 - (0.5*r_init/mr_init)*dmrdr);

    if (dtdrinit <= 0.0) {
        std::cerr << "BonnerEbertSphere::CalculateDrinitDt error : dtdrinit <= 0.0" << std::endl;
        std::cerr << std::scientific << (r_init/cst::ASTRONOMICAL_UNIT) << " au, " << (mr_init/cst::SOLAR_MASS) << " M_solar." << std::endl;
        std::exit(1);
    }

    drinitdt = 1.0 / dtdrinit;
    mdot_inf = 4.0*M_PI*rho_init*SQR(r_init)*drinitdt;
    return;
}


void BonnerEbertSphere::SolveLaneEmdenEquation(const int n, double &x, const double h, array::Double1D y)
{
    double c[5] = {0.0, 0.0, 0.5, 0.5, 1.0};
    double tx = 0.0;
    array::Double1D tmp = array::Allocate1dArray<double>(n);
    array::Double1D f   = array::Allocate1dArray<double>(n);
    array::Double2D K   = array::Allocate2dArray<double>(n, 5);

    for (int j = 1; j <= 4; ++j) {
        for (int i = 0; i < n; ++i) tmp[i] = y[i] + K[i][j-1]*c[j];
        tx   = x + c[j]*h;
        f[0] = y[1];
        f[1] = -2.0*y[1]/tx + std::exp(-y[0]);
        for (int i = 0; i < n; ++i) K[i][j] = h * f[i];
    }

    for (int i = 0; i < n; ++i) {
        y[i] = y[i] + (K[i][1] + K[i][4])/6.0 + (K[i][2] + K[i][3])/3.0;
    }

    array::Delete2dArray<double>(K);
    array::Delete1dArray<double>(f);
    array::Delete1dArray<double>(tmp);

    return;
}


void BonnerEbertSphere::Output(const std::string file_name)
{
    int ibmax = number_of_cloud_shells_;
    double a, b;

    std::ofstream file(file_name, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, file_name);

    file << std::scientific;
    file << (rout_bonner_/cst::ASTRONOMICAL_UNIT) << " " << rn_out_bonner_ << std::endl;

    for (int i = 0; i < ibmax; ++i) {

        a = SQR(rn_[i])*std::exp(-psi_[i]*0.5)*dpsidrn_[i];
        b = (1.0/std::sqrt(4.0*M_PI*CUB(cst::GRAVITATIONAL_CONSTANT)*central_mass_density_of_cloud_)) * CUB(cloud_sound_speed_) + SQR(rn_[i])*dpsidrn_[i]/cst::SOLAR_MASS;

        file << rn_[i] << " " 
             << r_bonner_[i] << " "
             << (r_bonner_[i]/cst::ASTRONOMICAL_UNIT) << " "
             << psi_[i] << " " 
             << dpsidrn_[i] << " "
             << a << " " 
             << rho_bonner_[i] << " " 
             << (rho_bonner_[i]/central_mass_density_of_cloud_) << " " 
             << mass_bonner_[i] << " "
             << (mass_bonner_[i]/cst::SOLAR_MASS) << " "
             << b << " "
             << (SQR(rn_[i])*dpsidrn_[i]) << " " 
             << std::endl;
    }

    file.close();
    return;
}