//////******************************
//
// This macros calculates the differential cross section in Laser-Compton scattering.
// The formula of the cross section is referred to [Nucl. Instr. Meth. A 495 (2002) 95] 
// *** I multiply one more r0 which is the classic radius of electron in Eq. 3. It should be there according to the cited paper. 

#include <iostream>
#include <cmath>

#include "beam.h" 

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"

double angle_degree_to_rad(double );
double cal_energy_scat_photon(double, double, double, double, double);
double cal_energy_scat_photon_headon(double, double, double, double);

double cal_differential_LCS_x_section(double alpha, double theta, beam electron_beam, beam laser_beam)
{
    // This function returns the differential cross section in unit m2. 
    // This function has inputs:
    // alpha, which is the relative angle between electron beam and scattered photon;
    // theta, which is the relative angle between electron beam and incoming laser; 

    std::cout << "theta in x section: " << theta << std::endl;
    double diff_x_section = 0.; // differential cross section
    
//    double alpha = 0.172; // this is the relative angle between electron beam and scattered photon
//    double theta = 0.; // this is the relative angle between electron beam and incoming laser
    
    double alpha_in_rad = angle_degree_to_rad(alpha);
    double theta_in_rad = angle_degree_to_rad(theta);
    
    double energy_scat_photon =0.;
    
    const double radius_clasc_elec = 2.81794e-15;
    double beta_elec = 0.;
    double energy_elec = 0.; // eV
    energy_elec = electron_beam.beam_energy;
//    beam electron_beam;
    beta_elec = electron_beam.cal_beta_elec(energy_elec);    

    double gamma_elec = 1./sqrt(1.-pow(beta_elec,2));

    std::cout << "beta: " << beta_elec << "gamma: " << gamma_elec << std::endl;
    
    double wavelength_laser = 1030.e-9;
    double energy_laser = 0.;
    
//    beam laser_beam;
    
    energy_laser = laser_beam.cal_energy_laser(wavelength_laser);
//    std::cout << "photon energy of laser: " << energy_laser << std::endl;
    
//    energy_scat_photon = cal_energy_scat_photon(theta, alpha, energy_elec, beta_elec, energy_laser);
//    std::cout << "scattered photon energy: " << energy_scat_photon << ", " << energy_scat_photon/(1.6e-19) << std::endl;

//    double temp_energy = cal_energy_scat_photon_headon(gamma_elec, energy_laser, alpha_in_rad, energy_elec);
//    energy_scat_photon =  temp_energy;
//    temp_energy =  cal_energy_scat_photon_headon(34.64, energy_laser, alpha_in_rad, energy_elec);
//    std::cout << "gamma square: " << gamma_elec*gamma_elec << std::endl;
//    std::cout << "head-on scattered photon energy: " << temp_energy << ", " << temp_energy/(1.6e-19) << std::endl;
    
    const double barn_unit_conversion = 10e-31;

//******************************
//
// Calculate the theta-diffrential cross section. Theta is the collision angle between the laser pulses and electron beams
// When alpha ~ 3 mrad. 
    alpha = 0.172;

//    for(int iTheta = 0; iTheta < 10; iTheta++)
//    {
//      theta = (double)iTheta*10.0;
      alpha_in_rad = angle_degree_to_rad(alpha);

      energy_scat_photon = cal_energy_scat_photon(theta, alpha, energy_elec, beta_elec, energy_laser);
      std::cout << "scattered photon energy: " << energy_scat_photon << std::endl;
//      energy_scat_photon = cal_energy_scat_photon_headon(gamma_elec, energy_laser, alpha_in_rad, energy_elec);

      double frt_term = M_PI * radius_clasc_elec * radius_clasc_elec* sin(alpha_in_rad) * ( (1 - beta_elec)/(1+beta_elec) ) * pow((energy_scat_photon/energy_laser) ,2);
      double scd_term = ( (1.- beta_elec*cos(alpha_in_rad)) / (1.+beta_elec) ) * ( energy_scat_photon / energy_laser);
      double thd_term = ((1.+beta_elec) / (1.-beta_elec*cos(alpha_in_rad)) )*(energy_laser / energy_scat_photon);
      double furth_term =  (1.-beta_elec*beta_elec)*(1.-cos(alpha_in_rad)*cos(alpha_in_rad)) / pow( (1.- beta_elec*cos(alpha_in_rad)), 2);
      
      diff_x_section = frt_term *( scd_term + thd_term - furth_term );
      std::cout << "differential cross section (theta): " << diff_x_section / barn_unit_conversion << std::endl;

//      diff_x_section = diff_x_section / barn_unit_conversion;
//      th1_x_section[1]->Fill(theta, diff_x_section);
//      th1_x_section[1]->SetBinError(iTheta+1, 0.0000001);

//    }

//    th1_x_section[1]->Draw();
//    c2->Print("LCS_x_section_vs_theta.pdf");


//    delete c2, th1_x_section;


  return diff_x_section;   
}

double angle_degree_to_rad(double degree)
{
    double rad = degree *M_PI / 180.;
    
    return rad;
}

double cal_energy_scat_photon(double theta, double alpha, double energy_elec, double beta, double energy_laser)
{
    double energy_scat_photon = 0.;
    double theta_in_rad = angle_degree_to_rad(theta);
    double alpha_in_rad = angle_degree_to_rad(alpha);
    
    energy_scat_photon = energy_laser*(1.+ beta*cos(theta_in_rad) ) / (1.- beta*cos(alpha_in_rad) + (energy_laser*(1.+cos(theta_in_rad-alpha_in_rad))/energy_elec ));
    
    return energy_scat_photon;
}

double cal_energy_scat_photon_headon(double gamma, double energy_laser, double alpha, double energy_elec)
{
  double energy_scat_photon = 0.;

  energy_scat_photon = 4. * gamma * gamma * energy_laser / (1. + pow(gamma*alpha, 2) );

  return energy_scat_photon;
}
