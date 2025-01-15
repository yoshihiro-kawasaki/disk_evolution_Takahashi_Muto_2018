/*
円盤進化計算
円盤内の物理量の時間発展を計算するクラスの実装
ガスとダストの速度はTakahashi & Muto 2018 を用いている。
2024/12/25 作成。
*/


#include "driver_base.hpp"

DriverBase::DriverBase(SimulationData *pdata)
    : pdata_(pdata)
{
    int nrtotc = pdata_->grid_.nrtotc_;
    int nrtotb = pdata_->grid_.nrtotb_;

    Text_              = array::Allocate1dArray<double>(nrtotc);
    dOmega_dr_         = array::Allocate1dArray<double>(nrtotc);
    Nvis_              = array::Allocate1dArray<double>(nrtotc);
    dNvis_dr_          = array::Allocate1dArray<double>(nrtotc);
    dP_dr_             = array::Allocate1dArray<double>(nrtotc);
    dPinteg_dr_        = array::Allocate1dArray<double>(nrtotc);
    flux_gas_          = array::Allocate1dArray<double>(nrtotb);
    flux_dust_         = array::Allocate1dArray<double>(nrtotb);
    flux_dust_m_       = array::Allocate1dArray<double>(nrtotb);
    mdot_inf_r_        = array::Allocate1dArray<double>(nrtotb);
    mdot_wind_r_       = array::Allocate1dArray<double>(nrtotb);
    mdot_r_            = array::Allocate1dArray<double>(nrtotb);
    dSdt_gas_          = array::Allocate1dArray<double>(nrtotc);
    dSdt_gas_src_      = array::Allocate1dArray<double>(nrtotc);
    dSdt_dust_         = array::Allocate1dArray<double>(nrtotc);
    dSdt_dust_src_     = array::Allocate1dArray<double>(nrtotc);
    dmSdt_dust_        = array::Allocate1dArray<double>(nrtotc);
    dmSdt_dust_growth_ = array::Allocate1dArray<double>(nrtotc);

    is_allocate_arrays_ = true;

    mdot_acc_disk_  = 0.0;
    mdot_acc_env_   = 0.0;
    mdot_inf_       = 0.0;
    mdot_wind_      = 0.0;
}

DriverBase::~DriverBase()
{
    if (is_allocate_arrays_) {
        array::Delete1dArray<double>(Text_);
        array::Delete1dArray<double>(dOmega_dr_);
        array::Delete1dArray<double>(Nvis_);
        array::Delete1dArray<double>(dNvis_dr_);
        array::Delete1dArray<double>(dP_dr_);
        array::Delete1dArray<double>(dPinteg_dr_);
        array::Delete1dArray<double>(flux_gas_);
        array::Delete1dArray<double>(flux_dust_);
        array::Delete1dArray<double>(flux_dust_m_);
        array::Delete1dArray<double>(mdot_inf_r_);
        array::Delete1dArray<double>(mdot_wind_r_);
        array::Delete1dArray<double>(mdot_r_);
        array::Delete1dArray<double>(dSdt_gas_);
        array::Delete1dArray<double>(dSdt_gas_src_);
        array::Delete1dArray<double>(dSdt_dust_);
        array::Delete1dArray<double>(dSdt_dust_src_);
        array::Delete1dArray<double>(dmSdt_dust_);
        array::Delete1dArray<double>(dmSdt_dust_growth_);
        is_allocate_arrays_ = false;
    }
}



void DriverBase::CalculateEnclosedMassAndAngularVelocity()
{
    int is   = pdata_->grid_.is_,
        ie   = pdata_->grid_.ie_;

    pdata_->gas_.Mr_bnd_[is-1] = pdata_->star_.mass_;
    pdata_->gas_.Mr_bnd_[is]   = pdata_->star_.mass_;
    for (int i = is+1; i <= ie; ++i) pdata_->gas_.Mr_bnd_[i] = pdata_->gas_.Mr_bnd_[i-1] + pdata_->grid_.r_vol_[i-1] * pdata_->gas_.Sigma_[i-1];
    pdata_->gas_.Mr_bnd_[ie+1] = pdata_->gas_.Mr_bnd_[ie];

    for (int i = is-1; i <= ie; ++i) {
        pdata_->gas_.Mr_cen_[i] = std::sqrt(pdata_->gas_.Mr_bnd_[i+1] * pdata_->gas_.Mr_bnd_[i]);
        pdata_->gas_.Omega_[i]  = std::sqrt(cst::GRAVITATIONAL_CONSTANT*pdata_->gas_.Mr_cen_[i]*CUB(pdata_->grid_.inv_r_cen_[i]));
        pdata_->gas_.j_[i]      = SQR(pdata_->grid_.r_cen_[i]) * pdata_->gas_.Omega_[i];
        dOmega_dr_[i]           = pdata_->gas_.Omega_[i] * (M_PI*pdata_->grid_.r_cen_[i]*pdata_->gas_.Sigma_[i]/pdata_->gas_.Mr_cen_[i] - 1.5*pdata_->grid_.inv_r_cen_[i]);
    }

    return;
}


void DriverBase::CalculateDiskGasQuantities()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;

    double alpha_GI, nu, dOmegadr;

    for (int i = is-1; i <= ie; ++i) {

        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR) {
            pdata_->gas_.T_[i]       = 10.0;
            pdata_->gas_.cs_[i]      = 0.0;
            pdata_->gas_.Hg_[i]      = 0.0;
            pdata_->gas_.rho_mid_[i] = 0.0;
            pdata_->gas_.P_[i]       = 0.0;
            pdata_->gas_.P_integ_[i] = 0.0;
            pdata_->gas_.Qt_[i]      = 0.0;
            pdata_->gas_.alpha_[i]   = 0.0;
            Nvis_[i]                 = 0.0;
            continue;
        }

        pdata_->gas_.T_[i]       = Text_[i];
        pdata_->gas_.cs_[i]      = std::sqrt(cst::BOLTZMANN_CONSTANT * pdata_->gas_.T_[i] * cst::INV_MG);
        pdata_->gas_.Hg_[i]      = pdata_->gas_.cs_[i] / pdata_->gas_.Omega_[i];
        pdata_->gas_.rho_mid_[i] = pdata_->gas_.Sigma_[i] / (cst::SQRT_2PI * pdata_->gas_.Hg_[i]);
        pdata_->gas_.P_[i]       = SQR(pdata_->gas_.cs_[i]) * pdata_->gas_.rho_mid_[i];
        pdata_->gas_.P_integ_[i] = SQR(pdata_->gas_.cs_[i]) * pdata_->gas_.Sigma_[i];
        pdata_->gas_.Qt_[i]      = pdata_->gas_.Omega_[i] * pdata_->gas_.cs_[i] / (M_PI * cst::GRAVITATIONAL_CONSTANT * pdata_->gas_.Sigma_[i]);
        alpha_GI                 = std::exp(-SQR(pdata_->gas_.Qt_[i])*SQR(pdata_->gas_.Qt_[i]));
        pdata_->gas_.alpha_[i]   = alpha_GI + pdata_->gas_.alpha_turb_;
        nu                       = pdata_->gas_.alpha_[i] * pdata_->gas_.cs_[i] * pdata_->gas_.Hg_[i];
        Nvis_[i]                 = CUB(pdata_->grid_.r_cen_[i]) * nu * pdata_->gas_.Sigma_[i] * dOmega_dr_[i];
    }

    return;
}


void DriverBase::CalculateGasDustDragCoefficient()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;
    double fmass, fsigmag, fsigmad, St_prime, A, Ainv, B;

    for (int i = is; i < ie; ++i) {

        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR) {
            pdata_->gas_.cvg_[i][0]  = 0.0;
            pdata_->gas_.cvg_[i][1]  = 0.0;
            pdata_->dust_.cvd_[i][0] = 0.0;
            pdata_->dust_.cvd_[i][1] = 0.0;
            continue;
        }

        fsigmag                  = pdata_->gas_.Sigma_[i]  / (pdata_->gas_.Sigma_[i] + pdata_->dust_.Sigma_[i]);
        fsigmad                  = pdata_->dust_.Sigma_[i] / (pdata_->gas_.Sigma_[i] + pdata_->dust_.Sigma_[i]);
        fmass                    = 2.0 * M_PI * SQR(pdata_->grid_.r_cen_[i]) * pdata_->gas_.Sigma_[i] / pdata_->gas_.Mr_cen_[i];
        St_prime                 = fsigmag * pdata_->dust_.St_[i];
        A                        = 1.0 + fmass;
        Ainv                     = 1.0 / A;
        B                        = 1.0 / (1.0 + A*SQR(St_prime));
        pdata_->gas_.cvg_[i][0]  = 1.0 - fsigmad * B;
        pdata_->gas_.cvg_[i][1]  = 2.0 * fsigmad * A * St_prime * B;
        pdata_->dust_.cvd_[i][0] = (fsigmag * B + fmass*pdata_->gas_.cvg_[i][0]) * Ainv;
        pdata_->dust_.cvd_[i][1] = (fmass * pdata_->gas_.cvg_[i][1] - 2.0*fsigmag*A*St_prime*B) * Ainv;
    }
}


void DriverBase::CalculateInfallAndWindRate()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;

    // calculate infall rate
    double j_core, inv_j_core, j;
    j_core = pdata_->cloud_.Omega_c_ * SQR(pdata_->cloud_.r_ini_);
    inv_j_core = 1.0 / j_core;
    mdot_inf_ = pdata_->cloud_.mdot_inf_;

    if (mdot_inf_ == 0.0) {
        for (int i = is; i <= ie; ++i) mdot_inf_r_[i] = 0.0;
    } else {
        for (int i = is; i <= ie; ++i) {
            j = std::sqrt(cst::GRAVITATIONAL_CONSTANT*pdata_->gas_.Mr_bnd_[i]*pdata_->grid_.r_bnd_[i]);
            if (j <= j_core) {
                mdot_inf_r_[i] = mdot_inf_ * (1.0 - std::sqrt(1.0 - j*inv_j_core));
            } else {
                mdot_inf_r_[i] = mdot_inf_;
            }
        }
    }

    // calculate wind mass loss rate
    if (mdot_inf_ == 0.0 && pdata_->gas_.c_wind_ != 0.0) {
        mdot_wind_r_[is] = 0.0;
        for (int i = is+1; i <= ie; ++i) {
            mdot_wind_r_[i] = mdot_wind_r_[i-1] + pdata_->grid_.r_vol_[i-1] * pdata_->gas_.c_wind_ * pdata_->gas_.Sigma_[i-1] * pdata_->gas_.Omega_[i-1];
        }
        mdot_wind_ = mdot_wind_r_[ie];
    }

    return;
}


void DriverBase::CalculateGasAndDustVelocity()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;

    // calculate gradient
    Nvis_[is-1]           = (Nvis_[is] * pdata_->grid_.inv_r_cen_[is]) * pdata_->grid_.r_cen_[is-1];
    pdata_->gas_.P_[is-1] = pdata_->gas_.P_[is] + ((pdata_->gas_.P_[is+1] - pdata_->gas_.P_[is]) * pdata_->grid_.inv_dr_cen_[is]) * (pdata_->grid_.r_cen_[is-1] - pdata_->grid_.r_cen_[is]);
    pdata_->grid_.CalculateNumericalDifferentiation(Nvis_, dNvis_dr_);
    pdata_->grid_.CalculateNumericalDifferentiation(pdata_->gas_.P_, dP_dr_);
    // pdata_->grid_.CalculateNumericalDifferentiation(pdata_->gas_.P_integ_, dPinteg_dr_);

    // calculate velocity
    double dlnP_dlnr, eta, mdot_tot_r;
    for (int i = is; i < ie; ++i) {

        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR) {
            pdata_->gas_.vr_vis_[i] = 0.0;
            pdata_->gas_.vr_eta_[i] = 0.0;
            pdata_->gas_.vr_src_[i] = 0.0;
            pdata_->gas_.vr_[i]     = 0.0;
            pdata_->dust_.vr_[i]    = 0.0;
            continue;
        }

        pdata_->gas_.vr_vis_[i] = 2.0 * dNvis_dr_[i] / (pdata_->gas_.j_[i] * pdata_->gas_.Sigma_[i]);

        dlnP_dlnr               = (pdata_->grid_.r_cen_[i] / pdata_->gas_.P_[i]) * dP_dr_[i];
        eta                     = - 0.5 * SQR(pdata_->gas_.Hg_[i] * pdata_->grid_.inv_r_cen_[i]) * dlnP_dlnr;
        pdata_->gas_.vr_eta_[i] = eta * pdata_->grid_.r_cen_[i] * pdata_->gas_.Omega_[i];

        mdot_tot_r              = 0.5 * (mdot_inf_r_[i] + mdot_inf_r_[i+1] - (mdot_wind_r_[i] + mdot_wind_r_[i+1]));
        pdata_->gas_.vr_src_[i] = - pdata_->grid_.r_cen_[i] * mdot_tot_r / pdata_->gas_.Mr_cen_[i];

        // total radial velocity
        pdata_->gas_.vr_[i]  = pdata_->gas_.cvg_[i][0]*pdata_->gas_.vr_vis_[i]  + pdata_->gas_.cvg_[i][1]*pdata_->gas_.vr_eta_[i]  + pdata_->gas_.vr_src_[i];
        pdata_->dust_.vr_[i] = pdata_->dust_.cvd_[i][0]*pdata_->gas_.vr_vis_[i] + pdata_->dust_.cvd_[i][1]*pdata_->gas_.vr_eta_[i] + pdata_->gas_.vr_src_[i];
    }

    return;
}


double DriverBase::CalculateDtAdv()
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double dt_temp, inv_dV, dSdt_gas_adv, dSdt_gas_src, dSdt_dust_adv, dSdt_dust_src;

    CalculateFlux();

    dt_temp  = 1.0e100;

    for (int i = is; i < ie; ++i) {

        inv_dV = pdata_->grid_.inv_r_vol_[i];

        // gas
        dSdt_gas_adv = - 2.0*M_PI * (pdata_->grid_.r_bnd_[i+1]*flux_gas_[i+1] - pdata_->grid_.r_bnd_[i]*flux_gas_[i]) * inv_dV;
        dSdt_gas_src = ((1.0 - pdata_->dust_.fdg_in_) * (mdot_inf_r_[i+1] - mdot_inf_r_[i]) - (mdot_wind_r_[i+1] - mdot_wind_r_[i])) * inv_dV;
        dSdt_gas_[i] = dSdt_gas_adv + dSdt_gas_src;

        if (pdata_->gas_.vr_[i] != 0.0 && pdata_->gas_.Sigma_[i] > SIGMA_GAS_FLOOR) {
            dt_temp = std::min(dt_temp, pdata_->grid_.dr_bnd_[i] / std::abs(pdata_->gas_.vr_[i]));
        }
        dt_temp = std::min(dt_temp, 1.0/pdata_->gas_.Omega_[i]);

        if (dSdt_gas_src != 0.0 && pdata_->gas_.Sigma_[i] > SIGMA_GAS_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->gas_.Sigma_[i]/dSdt_gas_src));
        }

        // dust
        dSdt_dust_adv = - 2.0*M_PI* (pdata_->grid_.r_bnd_[i+1]*flux_dust_[i+1] - pdata_->grid_.r_bnd_[i]*flux_dust_[i]) * inv_dV;
        dSdt_dust_src = pdata_->dust_.fdg_in_*(mdot_inf_r_[i+1] - mdot_inf_r_[i]) * inv_dV;
        dSdt_dust_[i] = dSdt_dust_adv + dSdt_dust_src;

        if (pdata_->dust_.vr_[i] != 0.0 && pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_bnd_[i] / pdata_->dust_.vr_[i]));
        }

        if (dSdt_dust_src != 0.0 && pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.Sigma_[i]/dSdt_dust_src));
        }
    }

    mdot_acc_disk_  = - 2.0*M_PI*pdata_->grid_.r_bnd_[is]*(flux_gas_[is] + flux_dust_[is]);
    mdot_acc_env_   = mdot_inf_r_[is];
    mdot_acc_total_ = mdot_acc_disk_ + mdot_acc_env_;

    return dt_temp;
}


bool DriverBase::DoTimeStep(double dt)
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;

    bool success =  CalculateTimeStep(dt);

    // boundary condition (ghost cell)
    // pdata_->gas_.Sigma_[is-1]  = pdata_->gas_.Sigma_[is];
    // pdata_->dust_.Sigma_[is-1] = pdata_->dust_.Sigma_[is];

    // time step
    pdata_->time_.t_      += dt;
    pdata_->cloud_.r_ini_ += pdata_->cloud_.dr_ini_dt_*dt;

    // update others
    pdata_->star_.mass_           += mdot_acc_total_ * dt;
    pdata_->gas_.total_mass_       = sigmag_tot_;
    pdata_->dust_.total_mass_      = sigmad_tot_;
    pdata_->total_disk_mass_       = sigmag_tot_ + sigmad_tot_;
    pdata_->total_mass_           += (mdot_inf_r_[ie] - mdot_wind_) * dt;
    pdata_->mdot_acc_disk_         = mdot_acc_disk_;
    pdata_->mdot_acc_env_          = mdot_acc_env_;
    pdata_->mdot_inf_              = mdot_inf_;
    pdata_->total_infall_mass_    += mdot_inf_r_[ie] * dt;
    pdata_->mdot_wind_             = mdot_wind_;
    pdata_->wind_mass_loss_       += mdot_wind_ * dt;

    return success;
}