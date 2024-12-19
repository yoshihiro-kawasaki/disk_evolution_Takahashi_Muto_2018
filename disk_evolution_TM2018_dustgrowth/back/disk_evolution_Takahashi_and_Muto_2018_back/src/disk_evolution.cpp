#include "disk_evolution.hpp"

DiskEvolution::DiskEvolution(SimulationData *pdata)
    : pdata_(pdata)
{
    int nrtotc = pdata_->grid_.nrtotc_;
    int nrtotb = pdata_->grid_.nrtotb_;

    Text_.Resize(nrtotc);

    torque_.Resize(nrtotc);

    flux_gas_.Resize(nrtotb);
    flux_dust_.Resize(nrtotb);

    mdot_inf_r_.Resize(nrtotb);
    mdot_wind_r_.Resize(nrtotb);

    dSdt_gas_.Resize(nrtotc);
    dSdt_dust_.Resize(nrtotc);

    mdot_acc_disk_   = 0.0;
    mdot_acc_env_    = 0.0;
    mdot_inf_        = 0.0;
    mdot_wind_       = 0.0;
    mdot_wind_sweep_ = 0.0;
}

DiskEvolution::~DiskEvolution()
{

}


void DiskEvolution::SetInitialCondtion()
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;

    double temp;

    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_(i)        = 0.0;
        pdata_->gas_.Sigma_floor_(i)  = 1.0e-15;
        pdata_->dust_.Sigma_(i)       = 0.0;
        pdata_->dust_.Sigma_floor_(i) = 1.0e-15;
        temp                          = 1.5e2 * std::pow(pdata_->grid_.r_cen_(i)*inv_AU, -0.42857142857142855);
        Text_(i)                      = std::max(temp, 10.0);
    }

    return;
}


void DiskEvolution::CalculateDiskGasParamters(double mdot_inf, double r_init, double c_wind)
{
    /*
    円盤内の物理量を計算

    Parameters
    mdot_inf : total infalling gas mass rate
    r_init   : initial cloud shell radius 
    c_wind   : efficiency of disk wind

    */

    int is = pdata_->grid_.is_;
    int ie = pdata_->grid_.ie_;

    double encl_mr, inv_encl_mr, dOmegadr, rhog_mid, Qt, nu;
    double fmass, fsigmag, fsigmad, St_prime, A, Ainv, B;

    // enclosed mass at cell boundary
    pdata_->gas_.Mr_(is-1) = pdata_->star_.mass_;
    pdata_->gas_.Mr_(is)   = pdata_->star_.mass_;
    for (int i = is+1; i <= ie; ++i) pdata_->gas_.Mr_(i) = pdata_->gas_.Mr_(i-1) + pdata_->grid_.r_vol_(i-1) * pdata_->gas_.Sigma_(i-1);
    pdata_->gas_.Mr_(ie+1) = pdata_->gas_.Mr_(ie);

    // others
    for (int i = is; i < ie; ++i) {

        encl_mr                = std::sqrt(pdata_->gas_.Mr_(i+1) * pdata_->gas_.Mr_(i));
        inv_encl_mr            = 1.0 / encl_mr;
        pdata_->gas_.Omega_(i) = std::sqrt(pc::GRAVITATIONAL_CONSTANT*encl_mr * CUB(pdata_->grid_.inv_r_cen_(i)));
        dOmegadr               = pdata_->gas_.Omega_(i) * (M_PI*pdata_->grid_.r_cen_(i)*pdata_->gas_.Sigma_(i) * inv_encl_mr - 1.5*pdata_->grid_.inv_r_cen_(i));
        pdata_->gas_.j_(i)     = SQR(pdata_->grid_.r_cen_(i)) * pdata_->gas_.Omega_(i);
        
        if (pdata_->gas_.Sigma_(i) < pdata_->gas_.Sigma_floor_(i)) {

            pdata_->gas_.T_(i)       = 10.0;
            pdata_->gas_.cs_(i)      = 0.0;
            pdata_->gas_.Hg_(i)      = 0.0;
            rhog_mid                 = 0.0;
            pdata_->gas_.P_(i)       = 0.0;
            Qt                       = 0.0;
            pdata_->gas_.alpha_(i)   = 0.0;
            torque_(i)               = 0.0;

            pdata_->gas_.cvg_(i, 0)  = 0.0;
            pdata_->gas_.cvg_(i, 1)  = 0.0;

            pdata_->dust_.cvd_(i, 0) = 0.0;
            pdata_->dust_.cvd_(i, 1) = 0.0;

        } else {

            pdata_->gas_.T_(i)       = Text_(i);
            pdata_->gas_.cs_(i)      = std::sqrt(pc::BOLTZMANN_CONSTANT * pdata_->gas_.T_(i) * inv_mg);
            pdata_->gas_.Hg_(i)      = pdata_->gas_.cs_(i) / pdata_->gas_.Omega_(i);
            rhog_mid                 = pdata_->gas_.Sigma_(i) / (sqrt_2pi * pdata_->gas_.Hg_(i));
            pdata_->gas_.P_(i)       = SQR(pdata_->gas_.cs_(i)) * rhog_mid;
            Qt                       = pdata_->gas_.Omega_(i) * pdata_->gas_.cs_(i) / (M_PI * pc::GRAVITATIONAL_CONSTANT * pdata_->gas_.Sigma_(i));
            pdata_->gas_.alpha_(i)   = std::exp(-SQR(Qt)*SQR(Qt)) + pdata_->gas_.alpha_turb_;
            nu                       = pdata_->gas_.alpha_(i) * pdata_->gas_.cs_(i) * pdata_->gas_.Hg_(i);
            torque_(i)               = CUB(pdata_->grid_.r_cen_(i)) * nu * pdata_->gas_.Sigma_(i) * dOmegadr;

            fsigmag                  = pdata_->gas_.Sigma_(i)  / (pdata_->gas_.Sigma_(i) + pdata_->dust_.Sigma_(i));
            fsigmad                  = pdata_->dust_.Sigma_(i) / (pdata_->gas_.Sigma_(i) + pdata_->dust_.Sigma_(i));
            fmass                    = 2.0 * M_PI * SQR(pdata_->grid_.r_cen_(i)) * pdata_->gas_.Sigma_(i) * inv_encl_mr;
            St_prime                 = fsigmag * pdata_->dust_.St_;
            A                        = 1 + fmass;
            Ainv                     = 1.0 / A;
            B                        = 1.0 / (1.0 + A*SQR(St_prime));
            pdata_->gas_.cvg_(i, 0)  = 1.0 - fsigmad * B;
            pdata_->gas_.cvg_(i, 1)  = 2.0 * fsigmad * A * St_prime * B;
            pdata_->dust_.cvd_(i, 0) = (fsigmag * B + fmass*pdata_->gas_.cvg_(i, 0)) * Ainv;
            pdata_->dust_.cvd_(i, 1) = (fmass * pdata_->gas_.cvg_(i, 1) - 2.0*fsigmag*A*St_prime*B) * Ainv;

        } // end if

    } // end for

    /* calculate mdot_inf and flux */
    double j, j_core, inv_j_core;
    j_core = pdata_->cloud_.Omega_c_ * SQR(r_init);
    inv_j_core = 1.0 / j_core;
    mdot_inf_ = mdot_inf;

    if (mdot_inf == 0.0) {
        mdot_inf_r_.Fill(0.0);
    } else {
        for (int i = is; i <= ie; ++i) {
            j = std::sqrt(pc::GRAVITATIONAL_CONSTANT*pdata_->gas_.Mr_(i)*pdata_->grid_.r_bnd_(i));
            if (j <= j_core) {
                mdot_inf_r_(i) = mdot_inf * (1.0 - std::sqrt(1.0 - j*inv_j_core));
            } else {
                mdot_inf_r_(i) = mdot_inf;
            }
        }
    }

    // disk wind  Takahashi and Muto 2018
    if (c_wind != 0.0) {
        mdot_wind_r_(1) = 0.0;
        for (int i = is+1; i <= ie; ++i) {
            mdot_wind_r_(i) = mdot_wind_r_(i-1) + pdata_->grid_.r_vol_(i-1) * c_wind * pdata_->gas_.Sigma_(i-1) * pdata_->gas_.Omega_(i-1);
        }
    }


    // radial velocity at cell center due to viscous torque and pressure gradient
    double gradT, gradP, inv_dr, eta;
    double mdot_tot_r, vgphi;

    // 境界条件(ゴーストセルでの値) : 改善の余地あり
    torque_(is-1)         = 0.0; // zero torque at inner boundary
    pdata_->gas_.P_(is-1) = pdata_->gas_.P_(is);

    // torque_(0) = torque_(1) 
    //     + (torque_(2) - torque_(1)) / (pdata_->grid_.r_cen_(2) - pdata_->grid_.r_cen_(1)) 
    //     * (pdata_->grid_.r_cen_(0) - pdata_->grid_.r_cen_(1));
    // pdata_->gas_.P_(0)  = pdata_->gas_.P_(1); 
    //     + (pdata_->gas_.P_(2) - pdata_->gas_.P_(1)) / (pdata_->grid_.r_cen_(2) - pdata_->grid_.r_cen_(1)) 
    //     * (pdata_->grid_.r_cen_(0) - pdata_->grid_.r_cen_(1));
    
    for (int i = is; i < ie; ++i) {

        if (pdata_->gas_.Sigma_(i) < pdata_->gas_.Sigma_floor_(i)) {
            pdata_->gas_.vr_vis_(i) = 0.0;
            pdata_->gas_.vr_eta_(i) = 0.0;
            pdata_->gas_.vr_src_(i) = 0.0;
            pdata_->gas_.vr_(i)     = 0.0;
            pdata_->dust_.vr_(i)    = 0.0;
            continue;
        }

        gradT = pdata_->grid_.cdiff_(i, 0)*torque_(i+1) + pdata_->grid_.cdiff_(i, 1)*torque_(i) + pdata_->grid_.cdiff_(i, 2)*torque_(i-1);
        gradP = (pdata_->grid_.r_cen_(i) /  pdata_->gas_.P_(i)) * pdata_->grid_.cdiff_(i, 0) * pdata_->gas_.P_(i+1) + pdata_->grid_.cdiff_(i, 1)*pdata_->gas_.P_(i)
                         + pdata_->grid_.cdiff_(i, 2)*pdata_->gas_.P_(i-1);

        pdata_->gas_.vr_vis_(i)  = 2.0 * gradT / (pdata_->gas_.j_(i) * pdata_->gas_.Sigma_(i));
        eta                      = - 0.5 * SQR(pdata_->gas_.Hg_(i) * pdata_->grid_.inv_r_cen_(i)) * gradP;
        pdata_->gas_.vr_eta_(i)  = eta * pdata_->grid_.r_cen_(i) * pdata_->gas_.Omega_(i);

        encl_mr                 = std::sqrt(pdata_->gas_.Mr_(i+1) * pdata_->gas_.Mr_(i));
        mdot_tot_r              = 0.5 * (mdot_inf_r_(i) + mdot_inf_r_(i+1) - (mdot_wind_r_(i)  + mdot_wind_r_(i+1)));
        pdata_->gas_.vr_src_(i) = - pdata_->grid_.r_cen_(i) * mdot_tot_r / encl_mr;

        pdata_->gas_.vr_(i)  = pdata_->gas_.cvg_(i, 0)*pdata_->gas_.vr_vis_(i)  
                               + pdata_->gas_.cvg_(i, 1)*pdata_->gas_.vr_eta_(i)  
                               + pdata_->gas_.vr_src_(i);

        pdata_->dust_.vr_(i) = pdata_->dust_.cvd_(i, 0)*pdata_->gas_.vr_vis_(i)  
                               + pdata_->dust_.cvd_(i, 1)*pdata_->gas_.vr_eta_(i)  
                               + pdata_->gas_.vr_src_(i);

    }

    return;

}

void DiskEvolution::InterpolationVelocityAtCellBoundary(Array1D<double> &vr_gas_bnd, Array1D<double> &vr_dust_bnd)
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    double rc, rb, inv_dr, m;

    // 
    pdata_->gas_.vr_(is-1)  = pdata_->gas_.vr_(is);
    pdata_->dust_.vr_(is-1) = pdata_->dust_.vr_(is);

    pdata_->gas_.vr_(ie)  = pdata_->gas_.vr_(ie-1);
    pdata_->dust_.vr_(ie) = pdata_->dust_.vr_(ie-1);


    for (int i = is; i < ie; ++i) {
        rc = pdata_->grid_.r_cen_(i-1);
        rb = pdata_->grid_.r_bnd_(i);
        inv_dr = pdata_->grid_.inv_dr_cen_(i-1);

        // gas 
        m = (pdata_->gas_.vr_(i) -  pdata_->gas_.vr_(i-1)) * inv_dr;
        vr_gas_bnd(i) = pdata_->gas_.vr_(i-1) + m*(rb - rc);

        // dust
        m = (pdata_->dust_.vr_(i) -  pdata_->dust_.vr_(i-1)) * inv_dr;
        vr_dust_bnd(i) = pdata_->dust_.vr_(i-1) + m*(rb - rc);
    }

    return;
}

void DiskEvolution::CalculateFlux()
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    Array1D<double> vr_gas_bnd(pdata_->grid_.nrtotb_), vr_dust_bnd(pdata_->grid_.nrtotb_);

    InterpolationVelocityAtCellBoundary(vr_gas_bnd, vr_dust_bnd);

    for (int i = is; i <= ie; ++i) {

        // gas
        if (vr_gas_bnd(i) >= 0.0) {
            if (i == is) {
                flux_gas_(i) = 0.0;
            } else {
                flux_gas_(i) = pdata_->gas_.Sigma_(i-1) * vr_gas_bnd(i);
            }
        } else {
            if (i == ie) {
                flux_gas_(i) == 0.0;
            } else {
                flux_gas_(i) = pdata_->gas_.Sigma_(i) * vr_gas_bnd(i);
            }
        }

        // dust
        if (vr_dust_bnd(i) >= 0.0) {
            if (i == is) {
                flux_dust_(i) = 0.0;
            } else {
                flux_dust_(i) = pdata_->dust_.Sigma_(i-1) * vr_dust_bnd(i);
            }
        } else {
            if (i == ie) {
                flux_dust_(i) == 0.0;
            } else {
                flux_dust_(i) = pdata_->dust_.Sigma_(i) * vr_dust_bnd(i);
            }
        }
    }

    return;
}


double DiskEvolution::CalculateDtMin()
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    double dt_temp, dV_inv, dSdt_gas_adv, dSdt_gas_src, dSdt_dust_adv, dSdt_dust_src;

    dt_temp  = 1.0e100;

    for (int i = is; i < ie; ++i) {

        dV_inv = pdata_->grid_.inv_r_vol_(i);

        // gas

        dSdt_gas_adv  = - 2.0*M_PI * (pdata_->grid_.r_bnd_(i+1)*flux_gas_(i+1) - pdata_->grid_.r_bnd_(i)*flux_gas_(i)) * dV_inv;
        dSdt_gas_src = (mdot_inf_r_(i+1) - mdot_inf_r_(i) - (mdot_wind_r_(i+1) - mdot_wind_r_(i))) * dV_inv;
        dSdt_gas_(i) = dSdt_gas_adv + dSdt_gas_src;

        if (pdata_->gas_.vr_(i) != 0.0 && pdata_->gas_.Sigma_(i) > pdata_->gas_.Sigma_floor_(i)) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_bnd_(i) / pdata_->gas_.vr_(i)));
        }

        if (dSdt_gas_src != 0.0 && pdata_->gas_.Sigma_(i) > pdata_->gas_.Sigma_floor_(i)) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->gas_.Sigma_(i)/dSdt_gas_src));
        }

        // dust
        dSdt_dust_adv = - 2.0*M_PI* (pdata_->grid_.r_bnd_(i+1)*flux_dust_(i+1) - pdata_->grid_.r_bnd_(i)*flux_dust_(i)) * dV_inv;
        dSdt_dust_src = pdata_->dust_.fdg_*(mdot_inf_r_(i+1) - mdot_inf_r_(i)) * dV_inv;
        dSdt_dust_(i) = dSdt_dust_adv + dSdt_dust_src;

        if (pdata_->dust_.vr_(i) != 0.0 && pdata_->dust_.Sigma_(i) > pdata_->dust_.Sigma_floor_(i)) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_bnd_(i) / pdata_->dust_.vr_(i)));
        }

        if (dSdt_dust_src != 0.0 && pdata_->dust_.Sigma_(i) > pdata_->dust_.Sigma_floor_(i)) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.Sigma_(i)/dSdt_dust_src));
        }

    }

    mdot_acc_disk_  = - 2.0*M_PI*pdata_->grid_.r_bnd_(is)*(flux_gas_(is) + flux_dust_(is));
    mdot_acc_env_   = (1.0 + pdata_->dust_.fdg_)*mdot_inf_r_(is);
    mdot_acc_total_ = mdot_acc_disk_ + mdot_acc_env_;

    return dt_temp;
}


void DiskEvolution::CalculateTimeStep(double dt)
{
    int i, j;
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    double sigmag_tot, sigmad_tot, dV;

    sigmag_tot = 0.0;
    sigmad_tot = 0.0;

    for (int i = is; i < ie; ++i) {

        pdata_->gas_.Sigma_(i) += dSdt_gas_(i) * dt;

        pdata_->dust_.Sigma_(i) += dSdt_dust_(i) * dt;

        dV = pdata_->grid_.r_vol_(i);
        sigmag_tot += pdata_->gas_.Sigma_(i) * dV;
        sigmad_tot += pdata_->dust_.Sigma_(i) * dV;
    }

    pdata_->star_.mass_           += mdot_acc_total_ * dt;
    pdata_->gas_.total_mass_       = sigmag_tot;
    pdata_->dust_.total_mass_      = sigmad_tot;
    pdata_->total_disk_mass_       = sigmag_tot + sigmad_tot;
    pdata_->total_mass_           += (mdot_inf_r_(ie) - mdot_wind_) * dt;
    pdata_->mdot_acc_disk_         = mdot_acc_disk_;
    pdata_->mdot_acc_env_          = mdot_acc_env_;
    pdata_->mdot_inf_              = mdot_inf_;
    pdata_->total_infall_mass_    += mdot_inf_r_(ie) * dt;
    pdata_->mdot_wind_             = mdot_wind_;
    pdata_->mdot_wind_sweep_       = mdot_wind_sweep_;
    pdata_->wind_loss_mass_       += mdot_wind_ * dt;
    pdata_->wind_loss_mass_sweep_ += mdot_wind_sweep_ * dt;
    pdata_->total_wind_loss_mass_ += (mdot_wind_ + mdot_wind_sweep_) * dt;
    
    return;
}

