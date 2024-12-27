#include "../sim_data.hpp"

void Dust::InputParameters(Utils::InputConfigure &input)
{
    vfrag_  = input.GetDouble("vfrag") * 1.0e2;
    ad0_    = input.GetDouble("ad0");
    rho_di_ = input.GetDouble("rho_di");
    fdg_in_ = input.GetDouble("fdg");
    md0_    = (4.0*M_PI/3.0) * rho_di_ * CUB(ad0_);

    if (vfrag_ == 0.0) {
        std::cerr << "#input error : vfrag == 0.0" << std::endl;
        std::exit(1);
    }

    if (ad0_ == 0.0) {
        std::cerr << "#input error : ad0 == 0.0" << std::endl;
        std::exit(1);
    }

    if (rho_di_ == 0.0) {
        std::cerr << "#input error : rho_di == 0.0" << std::endl;
        std::exit(1);
    }

    if (fdg_in_ == 0.0) {
        std::cerr << "#input error : fdg == 0.0" << std::endl;
        std::exit(1);
    }

    return;
}