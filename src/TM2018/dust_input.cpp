#include "../sim_data.hpp"

void Dust::InputParameters(Utils::InputConfigure &input)
{
    St_in_  = input.GetDouble("St"); 
    fdg_in_ = input.GetDouble("fdg");

    if (St_in_ == 0.0) {
        std::cerr << "#input error : St_in_ == 0.0" << std::endl;
        std::exit(1);
    }

    if (fdg_in_ == 0.0) {
        std::cerr << "#input error : fdg_in_ == 0.0" << std::endl;
        std::exit(1);
    }

    return;
}