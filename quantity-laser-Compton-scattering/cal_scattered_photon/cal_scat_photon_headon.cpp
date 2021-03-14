//////*********************************************
//
// This macro is created in order to calculate the number of scattered photons in laser Compton collisions.
// Formula is referred to Eq. 15 in the paper, https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.23.031601
// Eq. 2 in the paper, PRL 114 (2015) 195003
// created by HUANG at ELI on 10/12/2020
#include <iostream>
#include <cmath>

#include "TMath.h"
 
#include "beam.h" 

double cal_x_section_Thomson(void);
double cal_Fx(beam, beam);

int main()
{
    // output
    double num_scat_photons = 0.;
    
    const double speed_light = 3.0e8;
    
    double energy_elec = 4.2e6; // eV
    double intensity_elec = 40.e-12; //C
    double duration_elec = 3.e-12; //s
    double tran_size_elec = 40.e-6; // m
    double beta_elec = 1.0;
    double emittance_elec = 0.2e-6; 
    
    beam electron_beam;
    beta_elec = electron_beam.cal_beta_elec(energy_elec);
    double gamma = 1./sqrt(1.-beta_elec*beta_elec);

    electron_beam.set_beam_size(tran_size_elec, tran_size_elec, duration_elec, beta_elec*speed_light);
    electron_beam.emittance_electron = emittance_elec; 
    
    double num_elec = 0;
    num_elec = electron_beam.cal_num_elec_in_bunch(intensity_elec);
    
    double wavelength_laser = 1030.e-9; // m
    double intensity_laser = 1e-3; //joule
    double duration_laser = 0.43e-12; // s
    double tran_size_laser = 50.e-6; // m
    
    beam laser_beam;
    laser_beam.set_beam_size(tran_size_laser, tran_size_laser, duration_laser, speed_light);
    laser_beam.wavelength_gamma = wavelength_laser;
    
    double energy_laser = 0.;
    const double Plank_const = 6.626e-34;
    energy_laser = Plank_const * speed_light / wavelength_laser;

    std::cout << "laser energy: " << energy_laser << ", laser frequency: " << (speed_light/wavelength_laser) << std::endl;
    

    double num_photons = 0.;
    num_photons = laser_beam.cal_num_photon_in_pulse(intensity_laser, wavelength_laser);

    std::cout << "number of photons in laser: " << num_photons << ", number of electrons: " << num_elec << std::endl;
    std::cout << "electron beta: " << beta_elec << ", electron gamma: " << gamma << std::endl;
    
    double Fx = cal_Fx(electron_beam, laser_beam); 

//////
// Formula 1--https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.23.031601
//

    num_scat_photons = cal_x_section_Thomson() * num_elec * num_photons * Fx / ( sqrt(2*M_PI) *  sqrt( pow(electron_beam.beam_size[2],2) + pow(laser_beam.beam_size[2],2) ) * sqrt( pow(electron_beam.beam_size[0],2) + pow(laser_beam.beam_size[0],2) ) * sqrt( pow(electron_beam.beam_size[0]/electron_beam.emittance_electron,2) + pow(laser_beam.beam_size[0]/laser_beam.wavelength_gamma, 2) ) );
    
    std::cout << "Number of scattered photon: " << num_scat_photons << std::endl;

//////
// 
// Formula 2--PRL 114, 195003 (2015)
//
    const double elec_charge = 1.6e-19;
    const double elec_mass = 9.11e-31;
    const double hbar = 1.054e-34;
    double fine_structure =  elec_charge * elec_charge / hbar / speed_light;
    std::cout << "fine sturcture const.: " << fine_structure << std::endl;

    double norm_vec_potential = elec_charge * energy_laser / elec_mass / ( speed_light/wavelength_laser) / speed_light;

    std::cout << "normalized vector potential: " << norm_vec_potential << std::endl;
    double avg_harmonic_num = 1;


    num_scat_photons = (M_PI/3.) * fine_structure * num_photons * num_elec * pow(norm_vec_potential,2) * (1./avg_harmonic_num) * ((1.-beta_elec)/(1.+beta_elec) + (pow(norm_vec_potential,2) * (1.-beta_elec)/(4.*gamma*gamma*(1.+beta_elec)*(1.+beta_elec))) );

    std::cout <<  "Number of scattered photon second formula: " << num_scat_photons << std::endl;

return 0;    
}


double cal_x_section_Thomson(void)
{
    const double radius_clasc_elec = 2.81794e-15;  
//////***************************
// For head-on collisions
    double x_section_Thomson = 0.;
    x_section_Thomson = (8.*M_PI / 3.) * radius_clasc_elec * radius_clasc_elec;
    
//////***************************

    std::cout << "x section Thomson: " << x_section_Thomson << std::endl;
    return x_section_Thomson; 
}

double cal_Fx(beam elec_beam, beam laser_beam)
{
    double x = 0;
    double sigma_long = 0.;
    sigma_long = sqrt( pow(elec_beam.beam_size[2],2) + pow(laser_beam.beam_size[2],2) );
    double Rayleigh_length = M_PI * pow(laser_beam.beam_size[0],2) / laser_beam.wavelength_gamma; 
    x = (sqrt(2.) / sigma_long)* sqrt( (pow(elec_beam.beam_size[0],2) + pow(laser_beam.beam_size[0],2) ) / ( pow(elec_beam.beam_size[0]/elec_beam.emittance_electron,2) + pow(laser_beam.beam_size[0]/laser_beam.wavelength_gamma, 2) ) );
    double Fx = exp(pow(x,2))*(1-erf(x));
    
    return Fx;
}
