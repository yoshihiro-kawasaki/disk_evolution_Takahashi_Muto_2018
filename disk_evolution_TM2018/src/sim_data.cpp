#include "sim_data.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

Gas::Gas()
    : is_allocate_arrays_(false)
{

}

Gas::~Gas()
{
    if (is_allocate_arrays_) {
        array::Delete1dArray<double>(Sigma_);
        array::Delete1dArray<double>(Sigma_floor_);
        array::Delete1dArray<double>(T_);
        array::Delete1dArray<double>(cs_);
        array::Delete1dArray<double>(Omega_);
        array::Delete1dArray<double>(j_);
        array::Delete1dArray<double>(Hg_);
        array::Delete1dArray<double>(rho_mid_);
        array::Delete1dArray<double>(P_);
        array::Delete1dArray<double>(P_integ_);
        array::Delete1dArray<double>(Qt_);
        array::Delete1dArray<double>(Mr_bnd_);
        array::Delete1dArray<double>(Mr_cen_);
        array::Delete1dArray<double>(alpha_);
        array::Delete1dArray<double>(vr_);
        array::Delete1dArray<double>(vr_vis_);
        array::Delete1dArray<double>(vr_eta_);
        array::Delete1dArray<double>(vr_src_);
        array::Delete2dArray<double>(cvg_);
        array::Delete1dArray<double>(Sigma_dot_inf_);
        array::Delete1dArray<double>(Sigma_dot_wind_);
        is_allocate_arrays_ = false;
    }
}


void Gas::SetGas(Utils::InputConfigure &input, Grid &grid)
{
    alpha_turb_ = input.GetDouble("alpha_turb");
    c_wind_     = input.GetDouble("c_wind");

    if (is_allocate_arrays_) return;

    int nrtotc = grid.nrtotc_,
        nrtotb = grid.nrtotb_;

    Sigma_          = array::Allocate1dArray<double>(nrtotc);
    Sigma_floor_    = array::Allocate1dArray<double>(nrtotc);
    T_              = array::Allocate1dArray<double>(nrtotc);
    cs_             = array::Allocate1dArray<double>(nrtotc);
    Omega_          = array::Allocate1dArray<double>(nrtotc);
    j_              = array::Allocate1dArray<double>(nrtotc);
    Hg_             = array::Allocate1dArray<double>(nrtotc);
    rho_mid_        = array::Allocate1dArray<double>(nrtotc);
    P_              = array::Allocate1dArray<double>(nrtotc);
    P_integ_        = array::Allocate1dArray<double>(nrtotc);
    Qt_             = array::Allocate1dArray<double>(nrtotc);
    Mr_bnd_         = array::Allocate1dArray<double>(nrtotb);
    Mr_cen_         = array::Allocate1dArray<double>(nrtotc);
    alpha_          = array::Allocate1dArray<double>(nrtotc);
    vr_             = array::Allocate1dArray<double>(nrtotc);
    vr_vis_         = array::Allocate1dArray<double>(nrtotc);
    vr_eta_         = array::Allocate1dArray<double>(nrtotc);
    vr_src_         = array::Allocate1dArray<double>(nrtotc);
    cvg_            = array::Allocate2dArray<double>(nrtotc, 2);
    Sigma_dot_inf_  = array::Allocate1dArray<double>(nrtotc);
    Sigma_dot_wind_ = array::Allocate1dArray<double>(nrtotc);

    total_mass_ = 0.0;

    is_allocate_arrays_ = true;
    return;
}

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

Dust::Dust()
    : is_allocate_arrays_(false)
{

}


Dust::~Dust()
{
    if (is_allocate_arrays_) {
        array::Delete1dArray<double>(Sigma_);
        array::Delete1dArray<double>(Sigma_floor_);
        array::Delete1dArray<double>(vr_);
        array::Delete2dArray<double>(cvd_);
        is_allocate_arrays_ = false;
    }
}


void Dust::SetDust(Utils::InputConfigure &input, Grid &grid)
{
    St_  = input.GetDouble("St");
    fdg_ = input.GetDouble("fdg");

    if (is_allocate_arrays_) return;

    int nrtotc = grid.nrtotc_,
        nrtotb = grid.nrtotb_;

    Sigma_       = array::Allocate1dArray<double>(nrtotc);
    Sigma_floor_ = array::Allocate1dArray<double>(nrtotc);
    vr_          = array::Allocate1dArray<double>(nrtotc);
    cvd_         = array::Allocate2dArray<double>(nrtotc, 2);

    total_mass_ = 0.0;

    is_allocate_arrays_ = true;
}

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

void Cloud::SetCloud(Utils::InputConfigure &input)
{
    nshells_ = input.GetInt("nshells");
    nc_      = input.GetDouble("n_center");
    T_       = input.GetDouble("Tc");
    Omega_c_ = input.GetDouble("omega_c");
    f_       = input.GetDouble("fc");
    Mc_      = input.GetDouble("Mc") * cst::SOLAR_MASS;
}


// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

void Star::SetStar(Utils::InputConfigure &input)
{
    mass_init_ = input.GetDouble("Ms_init") * cst::SOLAR_MASS;
    mass_      = mass_init_;
}


// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

void Time::SetTime(Utils::InputConfigure &input)
{
    tend_       = input.GetDouble("tend") * cst::SOLAR_YEAR;
    cfl_        = input.GetDouble("cfl");
    delta_tout_ = input.GetDouble("delta_tout") * cst::SOLAR_YEAR;
    t_          = 0.0;
}


// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

SimulationData::SimulationData(Utils::InputConfigure &input)
{
    // make output directory
    std::string filename, inputfile;
    outdir_   = input.GetString("outdir");
    mkdir(outdir_.c_str(), 0777);
    filename  = input.GetInputFileName();
    inputfile = outdir_ + "/input";
    Utils::FileCopy(filename, inputfile);

    // Set
    grid_.SetGrid(input);
    gas_.SetGas(input, grid_);
    dust_.SetDust(input, grid_);
    cloud_.SetCloud(input);
    star_.SetStar(input);
    time_.SetTime(input);

    Initialization();
}


SimulationData::~SimulationData()
{

}


void SimulationData::Initialization()
{
    time_.t_              = 0.0;
    total_mass_           = star_.mass_;
    total_disk_mass_      = 0.0;
    mdot_acc_disk_        = 0.0;
    mdot_acc_env_         = 0.0;
    mdot_inf_             = 0.0;
    total_infall_mass_    = star_.mass_;
    mdot_wind_            = 0.0;
    wind_mass_loss_       = 0.0;
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
        file << grid_.r_cen_[i] << " ";
    }
    file << std::endl;

    for (int i = is; i <= is; ++i) {
        file << grid_.r_bnd_[i] << " ";
    }
    file << std::endl;

    for (int i = is; i < ie; ++i) {
        file << grid_.r_vol_[i] << " ";
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

    file.write((char*)&(time_.t_), sizeof(double));
    file.write((char*)&(star_.mass_), sizeof(double));
    file.write((char*)&(gas_.total_mass_), sizeof(double));
    file.write((char*)&(dust_.total_mass_), sizeof(double));
    file.write((char*)&(total_mass_), sizeof(double));
    file.write((char*)&(mdot_acc_disk_), sizeof(double));
    file.write((char*)&(mdot_acc_env_), sizeof(double));
    file.write((char*)&(mdot_inf_), sizeof(double));
    file.write((char*)&(total_infall_mass_), sizeof(double));
    file.write((char*)&(mdot_wind_), sizeof(double));
    file.write((char*)&(wind_mass_loss_), sizeof(double));

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(grid_.r_cen_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Sigma_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.T_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.cs_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Omega_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Hg_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)(&gas_.Qt_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)(&gas_.Mr_cen_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.alpha_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_vis_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_eta_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.vr_src_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.cvg_[i][0]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.cvg_[i][1]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Sigma_dot_inf_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(gas_.Sigma_dot_wind_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.Sigma_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.vr_[i]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.cvd_[i][0]), sizeof(double));
    }

    for (int i = is; i < ie; ++i) {
        file.write((char*)&(dust_.cvd_[i][1]), sizeof(double));
    }

    file.close();

    return;
}


void SimulationData::LogOutput(const int count, std::ofstream &file)
{
    file << std::scientific
         << std::setw(10) << count << " "
         << time_.t_ << " "
         << star_.mass_ << " "
         << gas_.total_mass_ << " "
         << dust_.total_mass_ << " "
         << total_mass_ << " "
         << mdot_acc_disk_ << " "
         << mdot_acc_env_ << " "
         << mdot_inf_ << " "
         << total_infall_mass_ << " "
         << mdot_wind_ << " "
         << wind_mass_loss_
         << std::endl;

    return;
}


void SimulationData::Show(const long int count)
{
    std::cout << std::scientific; // << std::setprecision(15);

    std::cout << "total count          = " << count << std::endl;
    std::cout << "star mass            = " << (star_.mass_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk gas mass  = " << (gas_.total_mass_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk dust mass = " << (dust_.total_mass_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk mass      = " << (total_disk_mass_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total mass           = " << (total_mass_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total infall mass    = " << (total_infall_mass_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "wind mass loss       = " << (wind_mass_loss_ / cst::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "m                    = " << ((total_mass_ + wind_mass_loss_) / cst::SOLAR_MASS) << " solar mass" << std::endl;
    return;
}