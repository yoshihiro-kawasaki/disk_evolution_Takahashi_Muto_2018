#include "bonner_ebert_sphere.hpp"
#include "physical_constants.hpp"
#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <cmath>


BonnerEbertSphere::BonnerEbertSphere(SimulationData *pdata)
    : pdata_(pdata)
{
    double mg = pc::GAS_MOLECULAR_MASS;

    number_of_cloud_shells_ = pdata_->cloud_.nshells_;
    central_number_density_of_cloud_ = pdata_->cloud_.nc_;
    central_mass_density_of_cloud_ = mg * central_number_density_of_cloud_;
    cloud_temperature_ = pdata_->cloud_.T_;
    cloud_sound_speed_ = std::sqrt(pc::BOLTZMANN_CONSTANT*cloud_temperature_/mg);
    cloud_angular_velocity_ = pdata_->cloud_.Omega_c_;
    density_factor_ = pdata_->cloud_.f_;
    total_cloud_mass_ = pdata_->cloud_.Mc_;
    initial_protostellar_mass_ = pdata_->star_.mass_init_;

    // pdata_->cloud_.cs_ = cloud_sound_speed_;
    // pdata_->cloud_.rhoc_= central_mass_density_of_cloud_;
}


BonnerEbertSphere::~BonnerEbertSphere()
{

}


void BonnerEbertSphere::SetBonnerEbertSphere()
{
    int i, ibmax;
    double rn_in = 1.0e-8, rn_out = 2.0e1;
    double dlnrn;

    ibmax = pdata_->cloud_.nshells_;

    rn_.Resize(ibmax, 0.0);
    psi_.Resize(ibmax, 0.0);
    dpsidrn_.Resize(ibmax, 0.0);
    rho_bonner_.Resize(ibmax, 0.0);
    r_bonner_.Resize(ibmax, 0.0);
    mass_bonner_.Resize(ibmax, 0.0);

    rout_bonner_ = 0.0;
    rn_out_bonner_ = 0.0;
    ib_con_ = 0;

    dlnrn = (std::log(rn_out/rn_in))/((double)ibmax);
    for (i = 0; i < ibmax; i++) {
        rn_(i) = rn_in * exp(double(i) * dlnrn);
    }

    // ## calcurate bonner ebert sphere profile ## //
    double xn, h;
    double r_nor, mass_nor;
    int ibmax_m = ibmax - 1;
    ExplicitOdeSolver ode("rk4");
    Array1D<double> y(2, 0.0);
    
    // boundary condition
    psi_(0) = 0.0;
    dpsidrn_(0) = 0.0;
    y(0) = psi_(0);
    y(1) = dpsidrn_(0);

    // integration
    for (i = 0; i < ibmax_m; ++i) {
        xn = rn_(i);
        h = rn_(i+1) - rn_(i);
        ode.RungeKutta(Func, xn, h, y, nullptr);
        psi_(i+1) = y(0);
        dpsidrn_(i+1) = y(1);
    }

    // calculate rho, r, mass of BE sphere
    r_nor = std::pow(4.0*M_PI*central_mass_density_of_cloud_*pc::GRAVITATIONAL_CONSTANT/SQR(cloud_sound_speed_), -0.5);
    mass_nor = (CUB(cloud_sound_speed_) / pc::GRAVITATIONAL_CONSTANT) / std::sqrt(4.0*M_PI*pc::GRAVITATIONAL_CONSTANT*central_mass_density_of_cloud_);

    for (i = 0; i < ibmax; ++i) {
        rho_bonner_(i) = density_factor_ * central_mass_density_of_cloud_ * std::exp(-psi_(i));
        r_bonner_(i) = r_nor*rn_(i);
        mass_bonner_(i) = density_factor_ * mass_nor * SQR(rn_(i)) * dpsidrn_(i);
    }

    // calculate rout_bonner
    for (i = 0; i < ibmax; ++i) {
        if (mass_bonner_(i) > total_cloud_mass_) {
            rout_bonner_ = Utils::LinearInterPolation(total_cloud_mass_, mass_bonner_(i-1), mass_bonner_(i), r_bonner_(i-1), r_bonner_(i));
            rn_out_bonner_ = rout_bonner_ / r_nor;
            std::cout << "cloud radius : " << std::scientific << (rout_bonner_/pc::ASTRONOMICAL_UNIT) << " au" << std::endl;
            break;
        }
    }

    if (mass_bonner_(ibmax_m) < total_cloud_mass_) {
        std::cerr << "error : mass_bonner[ibmax] < MASS_BONNER" << std::endl;
        std::exit(1);
    }

    return;
}


void BonnerEbertSphere::Func(double x, Array1D<double> &y, Array1D<double> &dydx, void *data)
{
    dydx(0) = y(1);
    dydx(1) = -2.0*y(1)/x + exp(-y(0));
    return;
}


double BonnerEbertSphere::CalculateInitialRinit()
{
    // 計算開始時の最初に中心(円盤)に落下してくる分子雲コアのshellの位置を求める
    // その位置は、初期の中心星質量よりコアの質量が大きくなる位置となる。

    int i;
    int i1 = 0, i2 = 0, ibmax = number_of_cloud_shells_;
    double r_init;

    for (i = 0; i < ibmax; ++i) {

        if (mass_bonner_(i) >= initial_protostellar_mass_) {
            i1 = i - 1;
            i2 = i;
            break;
            if (i2 == 0) {
                std::cerr << "BonnerEbertSphere::CalculateRinit error : i2 == 0" << std::endl;
                std::exit(1);
            }
        }

    }

    r_init = Utils::LinearInterPolation(initial_protostellar_mass_, mass_bonner_(i1), mass_bonner_(i2), r_bonner_(i1), r_bonner_(i2));

    return r_init;
}

void BonnerEbertSphere::CalculateMdotInfall(const double r_init, double &drinitdt, double &mdot_inf)
{
    double rho_init, dtdrinit, dmrdr, mr_init, f;
    int i;
    int i1 = 0, i2 = 0, ibmax = number_of_cloud_shells_;

    // 分子雲コアの位置r_initにおける質量密度rho_initを求める。

    for (i = ib_con_; i < ibmax; ++i) {

        if (r_bonner_(i) >= r_init) {
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

    rho_init = Utils::LinearInterPolation(r_init, r_bonner_(i1), r_bonner_(i2), rho_bonner_(i1), rho_bonner_(i2));

    // ref. Takahashi et al. 2013
    // 分子雲コアにおいて、初期にr_initの位置にあるshellが中心へ落下するまでの時間t_infはr_initの関数で与えられる、t_inf = tff(r_init)。
    // 逆に、r_initはtffの関数で書ける、r_init = r_init(t_inf)。
    // これの微分係数、dr_init/dt_infを求める。

    mr_init = Utils::LinearInterPolation(r_init, r_bonner_(i1), r_bonner_(i2), mass_bonner_(i1), mass_bonner_(i2));
    dmrdr = 4.0*M_PI*SQR(r_init)*rho_init;
    f = 2.56717;
    dtdrinit = f * sqrt(r_init/(2.0*pc::GRAVITATIONAL_CONSTANT*mr_init)) * (1.5 - (0.5*r_init/mr_init)*dmrdr);

    if (dtdrinit <= 0.0) {
        std::cerr << "BonnerEbertSphere::CalculateDrinitDt error : dtdrinit <= 0.0" << std::endl;
        std::cerr << std::scientific << (r_init/pc::ASTRONOMICAL_UNIT) << " au, " << (mr_init/pc::SOLAR_MASS) << " M_solar." << std::endl;
        std::exit(1);
    }

    drinitdt = 1.0 / dtdrinit;

    mdot_inf = 4.0*M_PI*rho_init*SQR(r_init)*drinitdt;

    return;
}


void BonnerEbertSphere::Output(std::string file_name)
{
    int i;
    int ibmax = number_of_cloud_shells_;
    double a, b;

    std::ofstream file(file_name, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, file_name);

    file << std::scientific;
    file << (rout_bonner_/pc::ASTRONOMICAL_UNIT) << " " << rn_out_bonner_ << std::endl;

    for (i = 0; i < ibmax; ++i) {

        a = SQR(rn_(i))*exp(-psi_(i)*0.5)*dpsidrn_(i);
        b = (1.0/sqrt(4.0*M_PI*CUB(pc::GRAVITATIONAL_CONSTANT)*central_mass_density_of_cloud_)) * CUB(cloud_sound_speed_)
            + SQR(rn_(i))*dpsidrn_(i)/pc::SOLAR_MASS;

        file << rn_(i) << " " 
             << r_bonner_(i) << " "
             << (r_bonner_(i)/pc::ASTRONOMICAL_UNIT) << " "
             << psi_(i) << " " 
             << dpsidrn_(i) << " "
             << a << " " 
             << rho_bonner_(i) << " " 
             << (rho_bonner_(i)/central_mass_density_of_cloud_) << " " 
             << mass_bonner_(i) << " "
             << (mass_bonner_(i)/pc::SOLAR_MASS) << " "
             << b << " "
             << (SQR(rn_(i))*dpsidrn_(i)) << " " 
             << std::endl;
    }

    file.close();
}