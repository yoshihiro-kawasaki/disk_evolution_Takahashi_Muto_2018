#include "array.hpp"
#include "constants.hpp"
#include "sim_data.hpp"
#include "disk_evolution_driver.hpp"
#include "bonner_ebert_sphere.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    std::string input_file;
    if (argc != 2) {
        std::cerr << "# ERROR : main, no input file" << std::endl;
        std::exit(1);
    }

    input_file = argv[1];
    std::cout << "input file : " << input_file << std::endl;

    /* input */
    Utils::InputConfigure input;
    input.ReadFile(input_file);

    /* simulation data */
    SimulationData *psim = new SimulationData(input);
    psim->OutputGrid();

    /* BE sphere */
    BonnerEbertSphere *pbes = new BonnerEbertSphere(psim);
    pbes->SetBonnerEbertSphere();
    pbes->Output("be_sphere.txt");

    /* disk evolution driver */
    DiskEvolutionDriver *pded = new DiskEvolutionDriver(psim);

    /* other variables */
    double dt, dt_disk, dt_collapse, tout;
    long int count = 0;
    bool is_end_infall;
    std::string outfilename;
    std::string logfile = psim->outdir_ + "/log";
    std::ofstream logofs;

    /* set initial condition */
    psim->output_count_ = 0;
    logofs.open(logfile, std::ios::out | std::ios::trunc);
    FAILED_TO_OPEN(logofs, logfile);
    psim->cloud_.rini_  = pbes->CalculateInitialRinit();
    pbes->CalculateMdotInfall(psim->cloud_.rini_, psim->cloud_.drinidt_, psim->cloud_.mdot_inf_);
    psim->time_.t_ = std::sqrt(0.5 * CUB(psim->cloud_.rini_) / (cst::GRAVITATIONAL_CONSTANT * psim->star_.mass_init_)) * C_TINF;
    tout = psim->time_.t_ + psim->time_.delta_tout_;
    pded->SetInitialCondtion();
    pded->CalculateDiskQuantities();
    psim->Show(count);

    /*  main integration */
    while (psim->time_.t_ < psim->time_.tend_) {

        count += 1;

        // determine dt
        if (psim->cloud_.drinidt_ == 0.0) {
            dt_collapse = 1.0e2 * cst::SOLAR_YEAR;
        } else {
            dt_collapse = 1.0e-2 * psim->cloud_.rini_ / psim->cloud_.drinidt_;
        }
        dt_disk = psim->time_.cfl_ * pded->CalculateDtMin();
        dt = std::min(dt_collapse, dt_disk);
        if (dt == 0.0) {
            std::cout << (dt_disk/cst::SOLAR_YEAR) << " " << (dt_collapse/cst::SOLAR_YEAR) << std::endl;
        }
        if ((psim->time_.tend_ - psim->time_.t_) < dt) dt = psim->time_.tend_ - psim->time_.t_;

        if (count % 200000 == 0) {
            std::cout << count << " " << (dt / cst::SOLAR_YEAR) << " " << (psim->time_.t_ / cst::SOLAR_YEAR) << std::endl;
        }

        // time step (n → n + 1)
        if (pded->CalculateTimeStep(dt) == false) {
            psim->output_count_++;
            outfilename = psim->outdir_ + "/disk" + std::to_string(psim->output_count_);
            psim->OutputData(outfilename);
            psim->LogOutput(psim->output_count_, logofs);
            std::cout << "error : time step" << std::endl;
            break;
        }
        pbes->CalculateMdotInfall(psim->cloud_.rini_, psim->cloud_.drinidt_, psim->cloud_.mdot_inf_);
        pded->CalculateDiskQuantities();

        // output
        if (psim->cloud_.rini_ >= pbes->GetRoutBonner() && is_end_infall == false) {
            std::cout << "end infall : time = " << (psim->time_.t_/cst::SOLAR_YEAR) << std::endl;
            is_end_infall = true;
            psim->output_count_++;
            outfilename = psim->outdir_ + "/disk" + std::to_string(psim->output_count_);
            psim->OutputData(outfilename);
            psim->LogOutput(psim->output_count_, logofs);
            // break;
        }

        if (psim->time_.t_ >= tout || psim->time_.t_ == psim->time_.tend_) {
            psim->output_count_++;
            outfilename = psim->outdir_ + "/disk" + std::to_string(psim->output_count_);
            psim->OutputData(outfilename);
            psim->LogOutput(psim->output_count_, logofs);
            tout += psim->time_.delta_tout_;
        }
    }

    psim->Show(count);

    logofs.close();

    end = std::chrono::system_clock::now();
    const double timed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "time = " << (timed * 1.0e-3) << "sec" << std::endl;

    return 0;
}