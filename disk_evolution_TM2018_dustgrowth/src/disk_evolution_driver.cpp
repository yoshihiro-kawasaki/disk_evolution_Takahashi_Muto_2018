/*
円盤進化計算
円盤内の物理量の時間発展を計算するクラスの実装
ガスとダストの速度はTakahashi & Muto 2018 を用いている。
2024/12/25 作成。
*/


#include "disk_evolution_driver.hpp"

DiskEvolutionDriver::DiskEvolutionDriver(SimulationData *pdata)
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
    Re_                = array::Allocate1dArray<double>(nrtotc);
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

DiskEvolutionDriver::~DiskEvolutionDriver()
{
    if (is_allocate_arrays_) {
        array::Delete1dArray<double>(Text_);
        array::Delete1dArray<double>(dOmega_dr_);
        array::Delete1dArray<double>(Nvis_);
        array::Delete1dArray<double>(dNvis_dr_);
        array::Delete1dArray<double>(dP_dr_);
        array::Delete1dArray<double>(dPinteg_dr_);
        array::Delete1dArray<double>(Re_);
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


void DiskEvolutionDriver::SetInitialCondtion()
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double temp;
    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_[i]        = 0.0;
        pdata_->gas_.Sigma_floor_[i]  = SIGMA_GAS_FLOOR;
        pdata_->dust_.Sigma_[i]       = 0.0;
        pdata_->dust_.Sigma_floor_[i] = SIGMA_DUST_FLOOR;
        pdata_->dust_.md_[i]          = 0.0;
        pdata_->dust_.ad_[i]          = 0.0;
        pdata_->dust_.mSigma_[i]      = 0.0;
        temp                          = 1.5e2 * std::pow(pdata_->grid_.r_cen_[i]/cst::ASTRONOMICAL_UNIT, -0.42857142857142855);
        Text_[i]                      = std::max(temp, 10.0);
    }

    return;
}

void DiskEvolutionDriver::CalculateDiskQuantities()
{
    CalculateEnclosedMassAndAngularVelocity();
    CalculateDiskGasQuantities();
    CalculateGasDustDragCoefficient();
    CalculateInfallAndWindRate();
    CalculateGasAndDustVelocity();
    CalculateDustGrowthRate();
    return;
}


void DiskEvolutionDriver::CalculateEnclosedMassAndAngularVelocity()
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


void DiskEvolutionDriver::CalculateDiskGasQuantities()
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
            Re_[i]                   = 0.0;
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
        Re_[i]                   = 0.5 * pdata_->gas_.alpha_turb_ * pdata_->gas_.Sigma_[i] * cst::H2_CROSS_SECTION * cst::INV_MG;
        nu                       = pdata_->gas_.alpha_[i] * pdata_->gas_.cs_[i] * pdata_->gas_.Hg_[i];
        Nvis_[i]                 = CUB(pdata_->grid_.r_cen_[i]) * nu * pdata_->gas_.Sigma_[i] * dOmega_dr_[i];
    }

    return;
}


void DiskEvolutionDriver::CalculateGasDustDragCoefficient()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;
    double inv_mfp, fmass, fsigmag, fsigmad, St_prime, A, Ainv, B;

    for (int i = is; i < ie; ++i) {

        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR) {
            pdata_->gas_.cvg_[i][0]  = 0.0;
            pdata_->gas_.cvg_[i][1]  = 0.0;
            pdata_->dust_.cvd_[i][0] = 0.0;
            pdata_->dust_.cvd_[i][1] = 0.0;
            pdata_->dust_.St_[i]     = 0.0;
            continue;
        }

        inv_mfp                  = (4.0 * pdata_->gas_.rho_mid_[i] * cst::H2_CROSS_SECTION) / (9.0 * cst::GAS_MOLECULAR_MASS);
        pdata_->dust_.St_[i]     = (0.5 * M_PI * pdata_->dust_.rho_di_ * pdata_->dust_.ad_[i] / pdata_->gas_.Sigma_[i]) * std::max(1.0, pdata_->dust_.ad_[i]*inv_mfp);
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


void DiskEvolutionDriver::CalculateDustGrowthRate()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;

    double StP, StT, HdT, Stmin, StTPast, sqr_dvbr, sqr_dvturbI, sqr_dvturbII, sqr_dvturb, dv, deltam;

    for (int i = is; i < ie; ++i) {
        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR || pdata_->dust_.Sigma_[i] < SIGMA_DUST_FLOOR) {
            pdata_->dust_.Hd_[i]  = 0.0;
            dmSdt_dust_growth_[i] = 0.0;
            continue;
        }

        StP = 0.5 * pdata_->dust_.St_[i];
        StT = pdata_->dust_.St_[i];

        pdata_->dust_.Hd_[i] = pdata_->gas_.Hg_[i] / std::sqrt(1.0 + StT*(1.0 + 2.0*StT) / (pdata_->gas_.alpha_turb_ * (1.0 + StT)));
        HdT                  = pdata_->dust_.Hd_[i];

        sqr_dvbr     = 2.0 * COE_VBR * pdata_->gas_.T_[i] / pdata_->dust_.md_[i];

        Stmin        = 1.0 / Re_[i];
        StTPast      = std::max(Stmin, std::min(1.6*StT, 1.0));
        sqr_dvturbI  = (StT-StP)/(StT + StP) * (SQR(StT)/(StTPast + StT) - SQR(StT)/(1+StT) - SQR(StP)/(StTPast + StP) + SQR(StP)/(1+StP));
        sqr_dvturbII = (StTPast - Stmin) * (((StTPast + Stmin)*StT + StTPast*Stmin)/((StT+StTPast)*(StT+Stmin)) + ((StTPast + Stmin)*StP + StTPast*Stmin)/((StP+StTPast)*(StP+Stmin)));
        sqr_dvturb   = pdata_->gas_.alpha_turb_ * SQR(pdata_->gas_.cs_[i]) * (sqr_dvturbI + sqr_dvturbII);
        if (sqr_dvturb < 0.0) sqr_dvturb = 0.0;

        dv = std::sqrt(sqr_dvbr + sqr_dvturb);

        deltam = std::min(1.0, - std::log(dv/pdata_->dust_.vfrag_)/std::log(5));

        dmSdt_dust_growth_[i] = 2.0 * cst::SQRT_PI * SQR(pdata_->dust_.ad_[i]) * dv * SQR(pdata_->dust_.Sigma_[i]) * deltam / HdT;
    }

    return;
}

void DiskEvolutionDriver::CalculateInfallAndWindRate()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;

    // calculate infall rate
    double j_core, inv_j_core, j;
    j_core = pdata_->cloud_.Omega_c_ * SQR(pdata_->cloud_.rini_);
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


void DiskEvolutionDriver::CalculateGasAndDustVelocity()
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


void DiskEvolutionDriver::CalculateFlux()
{
    int is     = pdata_->grid_.is_, 
        ie     = pdata_->grid_.ie_,
        nrtotb = pdata_->grid_.nrtotb_;
    double rb_rc, m, vgr_bnd, vdr_bnd;

    // zero gradient at boudaries
    pdata_->gas_.vr_[is-1]  = pdata_->gas_.vr_[is]; 
    pdata_->dust_.vr_[is-1] = pdata_->dust_.vr_[is];
    pdata_->gas_.vr_[ie]    = pdata_->gas_.vr_[ie-1];
    pdata_->dust_.vr_[ie]   = pdata_->dust_.vr_[ie-1];

    for (int i = is; i <= ie; ++i) {

        rb_rc = pdata_->grid_.r_bnd_[i] - pdata_->grid_.r_cen_[i-1];

        // gas
        m       = (pdata_->gas_.vr_[i] -  pdata_->gas_.vr_[i-1]) * pdata_->grid_.inv_dr_cen_[i-1];
        vgr_bnd = pdata_->gas_.vr_[i-1] + m*rb_rc;

        if (vgr_bnd >= 0.0) {
            if (i == is) {
                flux_gas_[i] = 0.0;
            } else {
                flux_gas_[i] = pdata_->gas_.Sigma_[i-1] * vgr_bnd;
            }
        } else {
            if (i == ie) {
                flux_gas_[i] = 0.0;
            } else {
                flux_gas_[i] = pdata_->gas_.Sigma_[i] * vgr_bnd;
            }
        }

        // dust
        m       = (pdata_->dust_.vr_[i] -  pdata_->dust_.vr_[i-1]) * pdata_->grid_.inv_dr_cen_[i-1];
        vdr_bnd = pdata_->dust_.vr_[i-1] + m*rb_rc;

        if (vdr_bnd >= 0.0) {
            if (i == is) {
                flux_dust_[i]   = 0.0;
                flux_dust_m_[i] = 0.0;
            } else {
                flux_dust_[i]   = pdata_->dust_.Sigma_[i-1]  * vdr_bnd;
                flux_dust_m_[i] = pdata_->dust_.mSigma_[i-1] * vdr_bnd;  
            }
        } else {
            if (i == ie) {
                flux_dust_[i]   = 0.0;
                flux_dust_m_[i] = 0.0;
            } else {
                flux_dust_[i]   = pdata_->dust_.Sigma_[i]  * vdr_bnd;
                flux_dust_m_[i] = pdata_->dust_.mSigma_[i] * vdr_bnd;
            }
        }

    }

    return;
}


double DiskEvolutionDriver::CalculateDtMin()
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double dt_temp, inv_dV, dSdt_gas_adv, dSdt_gas_src, dSdt_dust_adv, dSdt_dust_src, dmSdt_dust_adv, dmSdt_dust_src, dmSdt_dust_growth, m;

    CalculateFlux();

    dt_temp  = 1.0e100;

    for (int i = is; i < ie; ++i) {

        inv_dV = pdata_->grid_.inv_r_vol_[i];

        // gas
        dSdt_gas_adv = - 2.0*M_PI * (pdata_->grid_.r_bnd_[i+1]*flux_gas_[i+1] - pdata_->grid_.r_bnd_[i]*flux_gas_[i]) * inv_dV;
        dSdt_gas_src = ((1.0 - pdata_->dust_.fdg_) * (mdot_inf_r_[i+1] - mdot_inf_r_[i]) - (mdot_wind_r_[i+1] - mdot_wind_r_[i])) * inv_dV;
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
        dSdt_dust_src = pdata_->dust_.fdg_*(mdot_inf_r_[i+1] - mdot_inf_r_[i]) * inv_dV;
        dSdt_dust_[i] = dSdt_dust_adv + dSdt_dust_src;

        dmSdt_dust_adv    = - 2.0*M_PI* (pdata_->grid_.r_bnd_[i+1]*flux_dust_m_[i+1] - pdata_->grid_.r_bnd_[i]*flux_dust_m_[i]) * inv_dV;
        dmSdt_dust_growth = dmSdt_dust_growth_[i];
        if (pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            m = pdata_->dust_.md_[i];
        } else {
            m = pdata_->dust_.md0_;
        }
        dmSdt_dust_src = m * dSdt_dust_src;
        dmSdt_dust_[i] = dmSdt_dust_adv + dmSdt_dust_growth + dmSdt_dust_src;

        if (pdata_->dust_.vr_[i] != 0.0 && pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_bnd_[i] / pdata_->dust_.vr_[i]));
        }

        if (dSdt_dust_src != 0.0 && pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.Sigma_[i]/dSdt_dust_src));
        }

        if (dmSdt_dust_growth != 0.0 & pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.mSigma_[i] / dmSdt_dust_growth));
        }
    }

    mdot_acc_disk_  = - 2.0*M_PI*pdata_->grid_.r_bnd_[is]*(flux_gas_[is] + flux_dust_[is]);
    mdot_acc_env_   = mdot_inf_r_[is];
    mdot_acc_total_ = mdot_acc_disk_ + mdot_acc_env_;

    return dt_temp;
}


bool DiskEvolutionDriver::CalculateTimeStep(double dt)
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double sigmag_tot, sigmad_tot, dV, inv_rho_di, one_third;

    sigmag_tot = 0.0;
    sigmad_tot = 0.0;
    inv_rho_di = 1.0 / pdata_->dust_.rho_di_;
    one_third  = 1.0 / 3.0;

    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_[i]   += dSdt_gas_[i] * dt;
        pdata_->dust_.Sigma_[i]  += dSdt_dust_[i] * dt;
        pdata_->dust_.mSigma_[i] += dmSdt_dust_[i] * dt;
        if (pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            pdata_->dust_.md_[i] = pdata_->dust_.mSigma_[i] / pdata_->dust_.Sigma_[i];
            pdata_->dust_.ad_[i] = std::pow(COE_INV_SPHERE * pdata_->dust_.md_[i] * inv_rho_di, one_third);
        }

        dV          = pdata_->grid_.r_vol_[i];
        sigmag_tot += pdata_->gas_.Sigma_[i] * dV;
        sigmad_tot += pdata_->dust_.Sigma_[i] * dV;
    }

    // boundary condition (ghost cell)
    // pdata_->gas_.Sigma_[is-1]  = pdata_->gas_.Sigma_[is];
    // pdata_->dust_.Sigma_[is-1] = pdata_->dust_.Sigma_[is];

    // time step
    pdata_->time_.t_     += dt;
    pdata_->cloud_.rini_ += pdata_->cloud_.drinidt_*dt;

    // update others
    pdata_->star_.mass_           += mdot_acc_total_ * dt;
    pdata_->gas_.total_mass_       = sigmag_tot;
    pdata_->dust_.total_mass_      = sigmad_tot;
    pdata_->total_disk_mass_       = sigmag_tot + sigmad_tot;
    pdata_->total_mass_           += (mdot_inf_r_[ie] - mdot_wind_) * dt;
    pdata_->mdot_acc_disk_         = mdot_acc_disk_;
    pdata_->mdot_acc_env_          = mdot_acc_env_;
    pdata_->mdot_inf_              = mdot_inf_;
    pdata_->total_infall_mass_    += mdot_inf_r_[ie] * dt;
    pdata_->mdot_wind_             = mdot_wind_;
    pdata_->wind_mass_loss_       += mdot_wind_ * dt;

    return true;
}