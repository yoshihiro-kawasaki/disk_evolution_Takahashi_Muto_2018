#include "../sim_data.hpp"

void Dust::InputParameters(Utils::InputConfigure &input)
{
    ad0_    = input.GetDouble("ad");
    rho_di_ = input.GetDouble("rho_di");
    fdg_in_ = input.GetDouble("fdg");

    if (ad0_ == 0.0) {
        std::cerr << "#input error : ad == 0.0" << std::endl;
        std::exit(1);
    }

    if (rho_di_ == 0.0) {
        std::cerr << "#input error : rho_di == 0.0" << std::endl;
        std::exit(1);
    }

    if (fdg_in_ == 0.0) {
        std::cerr << "#input error : fdg_in_ == 0.0" << std::endl;
        std::exit(1);
    }
    return;
}