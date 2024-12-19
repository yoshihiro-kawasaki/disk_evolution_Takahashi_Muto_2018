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
    outdir_ = input.Get("outdir");
    mkdir(outdir_.c_str(), 0777);
    filename = input.GetInputFileName();
    inputfile = outdir_ + "/input";
    Utils::FileCopy(filename, inputfile);

    grid_.nr_ = std::stoi(input.Get("nr"));
    grid_.rmin_ =  std::stod(input.Get("rmin"));
    grid_.rmax_ = std::stod(input.Get("rmax"));
    grid_.spacing_ = input.Get("spacing");

    gas_.alpha_turb_ = std::stod(input.Get("alpha_turb"));
    gas_.c_wind_ = std::stod(input.Get("c_wind"));
    dust_.St_ = std::stod(input.Get("St"));
    dust_.fdg_ = std::stod(input.Get("fdg"));

    cloud_.nshells_ = std::stoi(input.Get("nshells"));
    cloud_.nc_ = std::stod(input.Get("n_center"));
    cloud_.T_ = std::stod(input.Get("Tc"));
    cloud_.Omega_c_ = std::stod(input.Get("omega_c"));
    cloud_.f_ = std::stod(input.Get("fc"));
    cloud_.Mc_ = std::stod(input.Get("Mc")) * pc::SOLAR_MASS;

    star_.mass_init_ = (1.0 + dust_.fdg_) * std::stod(input.Get("Ms_init")) * pc::SOLAR_MASS;
    star_.mass_ =  star_.mass_init_;

}


SimulationData::~SimulationData()
{

}


void SimulationData::Initialization()
{
    int ngh = 1; // number of ghost cells
    grid_.nrtotc_ = grid_.nr_ + 2*ngh;
    grid_.nrtotb_ = grid_.nr_ + 2*ngh + 1;
    // 計算(disk)領域は、is_ <= n < ie_
    grid_.is_ = ngh;
    grid_.ie_ = grid_.nr_ + ngh;

    grid_.r_cen_.Resize(grid_.nrtotc_);
    grid_.r_bnd_.Resize(grid_.nrtotb_);
    grid_.r_vol_.Resize(grid_.nrtotc_);
    grid_.dr_.Resize(grid_.nrtotc_);
    grid_.inv_r_cen_.Resize(grid_.nrtotc_);
    grid_.inv_r_vol_.Resize(grid_.nrtotc_);
    if (grid_.spacing_ == "log") {
        SetLogSpacingGrid();
    } else if (grid_.spacing_ == "root") {
        SetRootSpacingGrid();
    }

    gas_.Sigma_.Resize(grid_.nrtotc_);
    gas_.Sigma_floor_.Resize(grid_.nrtotc_);
    gas_.T_.Resize(grid_.nrtotc_);
    gas_.cs_.Resize(grid_.nrtotc_);
    gas_.Omega_.Resize(grid_.nrtotc_);
    gas_.j_.Resize(grid_.nrtotc_);
    gas_.Hg_.Resize(grid_.nrtotc_);
    gas_.P_.Resize(grid_.nrtotc_);
    gas_.Mr_.Resize(grid_.nrtotc_);
    gas_.alpha_.Resize(grid_.nrtotc_);
    gas_.vr_.Resize(grid_.nrtotc_);
    gas_.vr_vis_.Resize(grid_.nrtotc_);
    gas_.vr_eta_.Resize(grid_.nrtotc_);
    gas_.vr_src_.Resize(grid_.nrtotc_);


    dust_.Sigma_.Resize(grid_.nrtotc_);
    dust_.vr_.Resize(grid_.nrtotc_);

    time_ = 0.0;

    total_mass_ = star_.mass_;
    total_disk_mass_ = 0.0;
    mdot_ = 0.0;
    total_wind_loss_mass_ = 0.0;

    return;
}


void SimulationData::SetLogSpacingGrid()
{
    int is = grid_.is_, ie = grid_.ie_;
    double log_rmin, log_rmax, dlogr;

    log_rmin = std::log10(grid_.rmin_);
    log_rmax = std::log10(grid_.rmax_);
    dlogr    = (log_rmax - log_rmin) / (double(grid_.nr_));

    for (int i = 0; i < grid_.nrtotb_; ++i) {
        grid_.r_bnd_(i) = std::pow(10.0, log_rmin + dlogr*(double(i) - 1.0)) * pc::ASTRONOMICAL_UNIT;
    }

    for (int i = 0; i < grid_.nrtotc_; ++i) {
        grid_.r_cen_(i) = std::sqrt(grid_.r_bnd_(i) * grid_.r_bnd_(i+1));
        grid_.r_vol_(i) = M_PI * (SQR(grid_.r_bnd_(i+1)) - SQR(grid_.r_bnd_(i)));
        grid_.inv_r_cen_(i) = 1.0 / grid_.r_cen_(i);
        grid_.inv_r_vol_(i) = 1.0 / grid_.r_vol_(i);
        grid_.dr_(i) = grid_.r_bnd_(i+1) -  grid_.r_bnd_(i);
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
        grid_.r_cen_(i) = std::sqrt(grid_.r_bnd_(i)*grid_.r_bnd_(i+1));
        grid_.r_vol_(i) = M_PI * (SQR(grid_.r_bnd_(i+1)) - SQR(grid_.r_bnd_(i)));
        grid_.inv_r_cen_(i) = 1.0 / grid_.r_cen_(i);
        grid_.inv_r_vol_(i) = 1.0 / grid_.r_vol_(i);
        grid_.dr_(i) = grid_.r_bnd_(i+1) - grid_.r_bnd_(i);
    }

    return;
}


void SimulationData::Output(std::string filename)
{
    int is = grid_.is_, ie = grid_.ie_;

    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, filename);

    file << std::scientific;
    file << "time(s)          = " << time_ << std::endl;
    file << "time(solar year) = " << (time_/pc::SOLAR_YEAR) << std::endl;
    file << "star mass (g)    = " << star_.mass_ << std::endl;
    file << "star mass (Msun) = " << (star_.mass_/pc::SOLAR_MASS) << std::endl;
    file << "total mass (g)   = " << total_mass_ << std::endl;
    file << "mdot (g/s)       = " << mdot_ << std::endl;
    file << "mdot (Msun/yr)   = " << (mdot_*pc::SOLAR_YEAR/pc::SOLAR_MASS) << std::endl;
    // file << "luminocity ph    = " << (luminosity_ph_/pc::SOLAR_LUMINOSITY) << std::endl;
    // file << "luminocity acc   = " << (luminosity_acc_/pc::SOLAR_LUMINOSITY) << std::endl;

    for (int i = is; i < ie; ++i) {
        file << grid_.r_cen_(i) << " " 
             << gas_.Sigma_(i) << " " 
             // << gas_.rhog_mid_(i) << " "
             << gas_.T_(i) << " " 
             << gas_.cs_(i) << " " 
             << gas_.Omega_(i) << " " 
             << gas_.Hg_(i) << " "
             << gas_.alpha_(i) << " " 
             // << gas_.nu_vis_(i) << " "
             // << gas_.Qt_(i) << " "
             << gas_.vr_(i) << " "
             << gas_.vr_vis_(i) << " "
             << gas_.vr_eta_(i) << " "
             << gas_.vr_src_(i) << " "
             << dust_.Sigma_(i) << " "
             << dust_.vr_(i) << " "
             << std::endl;
    }

    file.close();
}