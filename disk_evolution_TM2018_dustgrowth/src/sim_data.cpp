#include "sim_data.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>

SimulationData::SimulationData(Utils::InputConfigure &input)
{
    // parameter input

    std::string filename, inputfile;
    outdir_   = input.Get("outdir");
    mkdir(outdir_.c_str(), 0777);
    filename  = input.GetInputFileName();
    inputfile = outdir_ + "/input";
    Utils::FileCopy(filename, inputfile);

    grid_.nr_      = std::stoi(input.Get("nr"));
    grid_.rmin_    =  std::stod(input.Get("rmin"));
    grid_.rmax_    = std::stod(input.Get("rmax"));
    grid_.spacing_ = input.Get("spacing");

    gas_.alpha_turb_ = std::stod(input.Get("alpha_turb"));
    gas_.c_wind_     = std::stod(input.Get("c_wind"));
    dust_.vfrag_     = std::stod(input.Get("vfrag")) * 1e2;
    dust_.fdg_       = std::stod(input.Get("fdg"));
    dust_.a0_        = std::stod(input.Get("a0"));
    dust_.rho_int_   = std::stod(input.Get("rho_int"));
    dust_.m0_        = (4.0*M_PI / 3.0) * dust_.rho_int_ * CUB(dust_.a0_);

    cloud_.nshells_ = std::stoi(input.Get("nshells"));
    cloud_.nc_      = std::stod(input.Get("n_center"));
    cloud_.T_       = std::stod(input.Get("Tc"));
    cloud_.Omega_c_ = std::stod(input.Get("omega_c"));
    cloud_.f_       = std::stod(input.Get("fc"));
    cloud_.Mc_      = std::stod(input.Get("Mc")) * pc::SOLAR_MASS;

    star_.mass_init_ = std::stod(input.Get("Ms_init")) * pc::SOLAR_MASS;
    star_.mass_      = star_.mass_init_;

}


SimulationData::~SimulationData()
{

}


void SimulationData::Initialization()
{
    int ngh = 2; // number of ghost cells
    grid_.nrtotc_ = grid_.nr_ + 2*ngh;
    grid_.nrtotb_ = grid_.nr_ + 2*ngh + 1;
    // 計算(disk)領域は、is_ <= n < ie_
    grid_.is_     = ngh;
    grid_.ie_     = grid_.nr_ + ngh;

    grid_.r_cen_.Resize(grid_.nrtotc_);
    grid_.r_bnd_.Resize(grid_.nrtotb_);
    grid_.r_vol_.Resize(grid_.nrtotc_);
    grid_.dr_cen_.Resize(grid_.nrtotc_);
    grid_.dr_bnd_.Resize(grid_.nrtotc_);
    grid_.inv_r_cen_.Resize(grid_.nrtotc_);
    grid_.inv_r_vol_.Resize(grid_.nrtotc_);
    grid_.inv_dr_cen_.Resize(grid_.nrtotc_);
    grid_.cdiff_.Resize(grid_.nrtotc_, 3);
    if (grid_.spacing_ == "log") {
        SetLogSpacingGrid();
    } else if (grid_.spacing_ == "root") {
        SetRootSpacingGrid();
    }
    SetNumericalDeviative();
    OutputGrid();

    gas_.Sigma_.Resize(grid_.nrtotc_);
    gas_.Sigma_floor_.Resize(grid_.nrtotc_);
    gas_.T_.Resize(grid_.nrtotc_);
    gas_.cs_.Resize(grid_.nrtotc_);
    gas_.Omega_.Resize(grid_.nrtotc_);
    gas_.j_.Resize(grid_.nrtotc_);
    gas_.Hg_.Resize(grid_.nrtotc_);
    gas_.P_.Resize(grid_.nrtotc_);
    gas_.P_integ_.Resize(grid_.nrtotc_);
    gas_.Qt_.Resize(grid_.nrtotc_);
    gas_.Mr_.Resize(grid_.nrtotc_);
    gas_.alpha_.Resize(grid_.nrtotc_);
    gas_.vr_.Resize(grid_.nrtotc_);
    gas_.vr_vis_.Resize(grid_.nrtotc_);
    gas_.vr_eta_.Resize(grid_.nrtotc_);
    gas_.vr_src_.Resize(grid_.nrtotc_);
    gas_.cvg_.Resize(grid_.nrtotc_, 2);
    gas_.Sigma_dot_inf_.Resize(grid_.nrtotc_);
    gas_.Sigma_dot_wind_.Resize(grid_.nrtotc_);

    dust_.Sigma_.Resize(grid_.nrtotc_);
    dust_.Sigma_floor_.Resize(grid_.nrtotc_);
    dust_.m_.Resize(grid_.nrtotc_);
    dust_.a_.Resize(grid_.nrtotc_);
    dust_.mSigma_.Resize(grid_.nrtotc_);
    dust_.St_.Resize(grid_.nrtotc_);
    dust_.Hd_.Resize(grid_.nrtotc_);
    dust_.vr_.Resize(grid_.nrtotc_);
    dust_.cvd_.Resize(grid_.nrtotc_, 2);

    time_                 = 0.0;
    mdot_acc_disk_        = 0.0;
    mdot_acc_env_         = 0.0;
    total_disk_mass_      = 0.0;
    total_mass_           = star_.mass_;
    mdot_inf_             = 0.0;
    total_infall_mass_    = star_.mass_;
    mdot_wind_            = 0.0;
    mdot_wind_sweep_      = 0.0;
    wind_loss_mass_       = 0.0;
    wind_loss_mass_sweep_ = 0.0;
    total_wind_loss_mass_ = 0.0;

    return;
}


void SimulationData::SetLogSpacingGrid()
{
    int is = grid_.is_, ie = grid_.ie_;
    double log_rmin, log_rmax, dlogr;

    // log_rmin = std::log10(grid_.rmin_);
    // log_rmax = std::log10(grid_.rmax_);
    // dlogr    = (log_rmax - log_rmin) / (double(grid_.nr_));

    // for (int i = 0; i < grid_.nrtotb_; ++i) {
    //     grid_.r_bnd_(i) = std::pow(10.0, log_rmin + dlogr*(double(i) - 1.0)) * pc::ASTRONOMICAL_UNIT;
    // }

    log_rmin = std::log(grid_.rmin_);
    log_rmax = std::log(grid_.rmax_);
    dlogr    = (log_rmax - log_rmin) / (double(grid_.nr_));

    for (int i = 0; i < grid_.nrtotb_; ++i) {
        grid_.r_bnd_(i) = std::exp(log_rmin + dlogr*(double(i) - 1.0)) * pc::ASTRONOMICAL_UNIT;
    }

    for (int i = 0; i < grid_.nrtotc_; ++i) {
        grid_.r_cen_(i)     = std::sqrt(grid_.r_bnd_(i) * grid_.r_bnd_(i+1));
        grid_.r_vol_(i)     = M_PI * (SQR(grid_.r_bnd_(i+1)) - SQR(grid_.r_bnd_(i)));
        grid_.inv_r_cen_(i) = 1.0 / grid_.r_cen_(i);
        grid_.inv_r_vol_(i) = 1.0 / grid_.r_vol_(i);
        grid_.dr_bnd_(i)    = grid_.r_bnd_(i+1) -  grid_.r_bnd_(i);
    }

    for (int i = 0; i < grid_.nrtotc_-1; ++i) {
        grid_.dr_cen_(i)     = grid_.r_cen_(i+1) - grid_.r_cen_(i);
        grid_.inv_dr_cen_(i) = 1.0 / grid_.dr_cen_(i);
    }

    return;
}


void SimulationData::SetRootSpacingGrid()
{
    double p, pinv, r0p, r1p, drp;

    p = 0.5;
    pinv = 1.0 / p;
    r0p = std::pow(grid_.rmin_, p);
    r1p = std::pow(grid_.rmax_, p);
    drp = (r1p - r0p) / (double(grid_.nr_));

    for (int i = 0; i < grid_.nrtotb_; ++i) {
        grid_.r_bnd_(i) = r0p + drp * i;
        grid_.r_bnd_(i) = std::pow(grid_.r_bnd_(i), pinv) * pc::ASTRONOMICAL_UNIT;
    }

    for (int i = 0; i < grid_.nrtotc_; ++i) {
        grid_.r_cen_(i)     = std::sqrt(grid_.r_bnd_(i)*grid_.r_bnd_(i+1));
        grid_.r_vol_(i)     = M_PI * (SQR(grid_.r_bnd_(i+1)) - SQR(grid_.r_bnd_(i)));
        grid_.inv_r_cen_(i) = 1.0 / grid_.r_cen_(i);
        grid_.inv_r_vol_(i) = 1.0 / grid_.r_vol_(i);
        grid_.dr_bnd_(i) = grid_.r_bnd_(i+1) - grid_.r_bnd_(i);
    }

    for (int i = 0; i < grid_.nrtotc_-1; ++i) {
        grid_.dr_cen_(i)     = grid_.r_cen_(i+1) - grid_.r_cen_(i);
        grid_.inv_dr_cen_(i) = 1.0 / grid_.dr_cen_(i);
    }

    return;
}


void SimulationData::SetNumericalDeviative()
{
    int is = grid_.is_, ie = grid_.ie_;

    double s, t;
    for (int i = is; i < ie; ++i) {
        s = grid_.r_cen_(i) - grid_.r_cen_(i-1);
        t = grid_.r_cen_(i+1) - grid_.r_cen_(i);
        grid_.cdiff_(i, 0) = s*s / (s*t*(s+t));
        grid_.cdiff_(i, 1) = (t*t - s*s) / (s*t*(s+t));
        grid_.cdiff_(i, 2) = - t*t / (s*t*(s+t)); 
    }

    return;
}


void SimulationData::OutputGrid()
{
    int is = grid_.is_, ie = grid_.ie_;

    std::string filename = outdir_ + "/grid";
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, filename);

    file << grid_.nr_ << std::endl;

    file << std::scientific;

    for (int i = is; i < ie; ++i) {
        file << grid_.r_cen_(i) << " ";
    }
    file << std::endl;

    for (int i = is; i <= is; ++i) {
        file << grid_.r_bnd_(i) << " ";
    }
    file << std::endl;

    for (int i = is; i < ie; ++i) {
        file << grid_.r_vol_(i) << " ";
    }
    file << std::endl;

    file.close();

    return;
}


void SimulationData::OutputData(std::string filename)
{
    int is = grid_.is_, ie = grid_.ie_;

    std::ofstream file(filename, std::ios::trunc | std::ios::out | std::ios::binary);
    FAILED_TO_OPEN(file, filename);

    file.write((char*)&time_, sizeof(double));
    file.write((char*)&(star_.mass_), sizeof(double));
    file.write((char*)&(star_.Ls_), sizeof(double));
    file.write((char*)&(star_.Lacc_), sizeof(double));
    file.write((char*)&(gas_.total_mass_), sizeof(double));
    file.write((char*)&(dust_.total_mass_), sizeof(double));
    file.write((char*)&(total_mass_), sizeof(double));
    file.write((char*)&(mdot_acc_disk_), sizeof(double));
    file.write((char*)&(mdot_acc_env_), sizeof(double));
    file.write((char*)&(mdot_inf_), sizeof(double));
    file.write((char*)&(total_infall_mass_), sizeof(double));
    file.write((char*)&(mdot_wind_), sizeof(double));
    file.write((char*)&(mdot_wind_sweep_), sizeof(double));
    file.write((char*)&(wind_loss_mass_), sizeof(double));
    file.write((char*)&(wind_loss_mass_sweep_), sizeof(double));
    file.write((char*)&(total_wind_loss_mass_), sizeof(double));

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(grid_.r_cen_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Sigma_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.T_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.cs_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Omega_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Hg_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)(&gas_.Qt_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)(&gas_.Mr_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.alpha_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_vis_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_eta_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_src_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.cvg_(i, 0)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.cvg_(i, 1)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Sigma_dot_inf_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Sigma_dot_wind_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.Sigma_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.m_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.a_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.St_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.Hd_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.vr_(i)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.cvd_(i, 0)), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.cvd_(i, 1)), sizeof(double));
    }

    file.close();

    return;
}


void SimulationData::LogOutput(const int count, std::ofstream &file)
{
    file << std::scientific
         << std::setw(10) << count << " "
         << time_ << " "
         << star_.mass_ << " "
         << star_.Ls_ << " "
         << star_.Lacc_ << " "
         << gas_.total_mass_ << " "
         << dust_.total_mass_ << " "
         << total_mass_ << " "
         << mdot_acc_disk_ << " "
         << mdot_acc_env_ << " "
         << mdot_inf_ << " "
         << total_infall_mass_ << " "
         << mdot_wind_ << " "
         << mdot_wind_sweep_ << " "
         << wind_loss_mass_ << " "
         << wind_loss_mass_sweep_ << " "
         << total_wind_loss_mass_ << " "
         << std::endl;

    return;
}