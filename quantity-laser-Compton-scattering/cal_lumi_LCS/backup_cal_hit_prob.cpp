#include <iostream>
#include <cmath>

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"

int main()
{
// output: the total cross section
  float x_section = 0.;

// input: number of scattered photons, luminosity
  float num_scat_photons = 0.;
  float lumi = 0.;
  float pow_factor = 0.;
  float geom_factor = 0.;
  const float speed_light = 3.0e8;

  float num_e_per_bunch, num_photons_in_laser = 0.;

  float e_beam_intensity = 40e-12; // in C
  const float single_e_charge = 1.6e-19;
  float beta = 0.; // the velocity of electrons relative to the speed of light
  beta = 0.985; // which is corresponding to a electron with energy of 4.2 MeV

  float lorentz_factor = 1. / sqrt(1.-pow(beta,2));
  std::cout << "lorentz factor " << lorentz_factor << std::endl;
  const float e_mass = 9.11e-31;
  float e_energy = (lorentz_factor -1) * e_mass * speed_light * speed_light;

  std::cout << "one electron energy (J): " << e_energy << std::endl;

//  num_e_per_bunch = 4.2e6 * single_e_charge / e_energy ; 
  num_e_per_bunch = e_beam_intensity / single_e_charge;
//  num_e_per_bunch = e_beam_intensity * (6.24e18);

  float e_beam_duration = 3e-12;

  float E_photon =0;
  float wavelength_photon = 1030e-9;
  const float plank_const = 6.626e-34;
  E_photon = plank_const * speed_light / wavelength_photon;
//  E_photon = 1.2 * single_e_charge; // in joule

  std::cout << "one photon energy (J): " << E_photon << std::endl;

  float laser_duration = 0.43e-12;

  float laser_beam_intensity = 10e-3; // in J
  num_photons_in_laser = laser_beam_intensity / E_photon;

  std::cout << "number of electrons in a bunch: " << num_e_per_bunch << std::endl;
  std::cout << "number of photons in laser pulse: " << num_photons_in_laser << std::endl;
  pow_factor = num_e_per_bunch * num_photons_in_laser;

  float sigma[3] = {0.};  // electron bunch size in x, y , z.
  float sigma_prime[3] = {0.}; //laser pulse bunch size in x, y, z.
  float theta = 0.; // theta the angle between electron beam and laser
   
  sigma[0] = 40e-6;
  sigma[1] = 40e-6;
  sigma[2] = e_beam_duration * beta* speed_light ; // beam duration times its speed

  sigma_prime[0] = 50e-6;
  sigma_prime[1] = 50e-6;
  sigma_prime[2] = laser_duration * speed_light; // beam duration times speed of light

  TH1F *th1_lumi = new TH1F("th1_lumi",";collisions angle #theta (deg); luminosity (mb)^{-1}",9, 0, 90);
  th1_lumi->Sumw2();
  th1_lumi->SetStats(0);
  std::cout << "pi " << M_PI << std::endl;
  std::cout << "angle, power factor, geometry factor, luminosity" << std::endl;

  for(int iBin = 0; iBin < 9; iBin++)
  {
    theta = (float)iBin* 10;

    float theta_in_pi = theta*M_PI / 180.0;
    float numerator = (1. + beta*cos(theta_in_pi));
    float denominator = (2 * M_PI * sqrt(pow(sigma[1],2)+pow(sigma_prime[1],2)) * sqrt( pow(sigma[0],2)*pow(beta+cos(theta_in_pi),2) + pow(sigma_prime[0],2)*pow(1+beta*cos(theta_in_pi),2) + (pow(sigma[2],2)+pow(sigma_prime[2],2))*pow(sin(theta_in_pi),2)) );
    geom_factor = numerator / denominator;

//    pow_factor = 1;
  
    lumi = pow_factor * geom_factor;    
    lumi = lumi / (1e31);
    std::cout << theta <<", " << pow_factor << ", " << geom_factor << ", " << lumi << std::endl;

    th1_lumi->Fill(theta, lumi);
//    th1_lumi->SetBinError(iBin+1, 0.0000001);
    th1_lumi->SetBinError(iBin+1, 0.0);
  }
//  th1_lumi->GetYaxis()->SetRangeUser(0.,90000);



  TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
  c2->SetLeftMargin(0.15);

  th1_lumi->Draw();
  c2->Print("LCS_lumi_vs_beta.pdf");

return 0;
}
