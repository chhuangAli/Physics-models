//////**************************
//
// This source code calculates luminosity in electron-laser collisions. 
// The luminosity formula is from the proceeding, luminosity increase in laser-Compton scattering by Crab crossing method, MOPVA023, proceedings of IPAC2017.

#include <iostream>
#include <cmath>

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"

float cal_pow_factor(float, float);
float cal_geom_factor(float, float *, float *, float );
float cal_error_propagation_typeone(float , float , float , float, float, float );
float cal_error_propagation_geom(float, float, float, float, float, float); 
float cal_integral_x_section_Thomson(void);
float cal_integral_x_section_KN(float);

//////***************************
//
//    this class is designed for storing parameters related to electron or laser beams.
//
class beam
{
    private:
    
        float time_duration = 0.;
        float beam_intensity = 0.;
        const float single_e_charge = 1.6e-19;
        const float mass_e = 9.1e-31;
        const float speed_light = 3.0e8;
        const float plank_const = 6.626e-34;
        
    public:
        
        float beam_energy = 0.;
        float beam_size[3] = {0}; // in x, y, z
        float number_of_particles_in_beam;
                
        void set_beam_duration(float);
        float get_beam_duration(void);
        void set_beam_size(float, float, float, float);
        float cal_num_elec_in_bunch(float); 
        float cal_num_photon_in_pulse(float, float);
        double cal_beta_elec(double); // one input which is energy in eV.
    
    
};

void beam::set_beam_duration(float time)
{
    time_duration = time;    
}

float beam::get_beam_duration(void)
{
    return time_duration;
}


void beam::set_beam_size(float x_size, float y_size, float time, float speed)
{        
    beam_size[0] = x_size; 
    beam_size[1] = y_size;
    time_duration = time;
    beam_size[2] = (time_duration * speed); // for z size    
//    return longitudinal_size;
}

float beam::cal_num_elec_in_bunch(float intensity )
{    
    float number_of_elec_in_bunch = 0.;
    
    //take care of unit of the intensity
    beam_intensity = intensity;

    number_of_elec_in_bunch = beam_intensity / single_e_charge;
    return number_of_elec_in_bunch;
}

float beam::cal_num_photon_in_pulse(float intensity, float wave_length)
{
    float num_photons_in_laser = 0.;
    
    beam_intensity = intensity;
    float E_photon = 0.;
    
    E_photon = plank_const * speed_light / wave_length;
    num_photons_in_laser = beam_intensity / E_photon;
    
    return num_photons_in_laser;
}

double beam::cal_beta_elec(double energy)
{
    double gamma = 0.;
    beam_energy = energy*single_e_charge; // change to joules
    gamma = (beam_energy)/(mass_e*speed_light*speed_light);
    double beta = 0.; // output
    beta = sqrt(1.-(1./pow(gamma,2)));
    
    return beta; 
}

int main()
{
// output: the total cross section
//  float x_section = 0.;

// input: number of scattered photons, luminosity
  float num_scat_photons = 0.;
  float lumi = 0.;
  float pow_factor = 0.;
  float geom_factor = 0.;
  const float speed_light = 3.0e8;
  

  float num_e_per_bunch, num_photons_in_laser = 0.;


  double beta = 0.; // the velocity of electrons relative to the speed of light
//  beta = 0.985; // which is corresponding to a electron with energy of 4.2 MeV


  float wavelength_photon = 1030e-9;
  float laser_duration = 0.43e-12;
  float laser_beam_intensity = 1e-3; // in J

  float theta = 0.; // theta the angle between electron beam and laser
  

  const float barn_unit_conversion = 1e34;

  TH1F *th1_lumi[3];
  for(int i = 0; i<3; i++)
  {
    th1_lumi[i] = new TH1F(Form("th1_lumi_%i ", i),";collisions angle #theta (deg); luminosity (ub)^{-1}",9, 0, 90);
    th1_lumi[i]->Sumw2();
    th1_lumi[i]->SetStats(0);
  }
  //std::cout << "pi " << M_PI << std::endl;
  //std::cout << "angle, power factor, geometry factor, luminosity" << std::endl;

  
//////*****************************
// 
// use my defined class to reproduce power factor and geometry factor
//
  
  beam elec_beam;
  float e_beam_intensity = 0; // in C
  e_beam_intensity = 10.e-12;
  float e_beam_duration = 0.;
  e_beam_duration = 100.e-15;
  elec_beam.set_beam_duration(e_beam_duration);
  
  std::cout << "number of electrons in a bunch (calculated from the class): " << elec_beam.cal_num_elec_in_bunch(e_beam_intensity) << std::endl;
  num_e_per_bunch = elec_beam.cal_num_elec_in_bunch(e_beam_intensity);
  
  double elec_beam_energy = 100.e6; // eV
//  elec_beam_energy = 4.2e6;
  beta = elec_beam.cal_beta_elec(elec_beam_energy);
  std::cout << "beta of electron: " << beta << std::endl;   
 
  float gamma_elec = 1./sqrt(1.-pow(beta,2)); 
  
  beam laser_beam;
  laser_beam_intensity = 1.;//in joule
//  laser_beam_intensity = 10e-3;
  std::cout << "number of photons in laser pulse (calculated from the class): " << laser_beam.cal_num_photon_in_pulse(laser_beam_intensity, wavelength_photon) << std::endl;
  num_photons_in_laser = laser_beam.cal_num_photon_in_pulse(laser_beam_intensity, wavelength_photon);

  pow_factor = cal_pow_factor(num_e_per_bunch, num_photons_in_laser);
  
  float laser_beam_duration = 30.e-15;
//  laser_beam_duration = 0.43e-12;
  laser_beam.set_beam_duration(laser_beam_duration);
  
  const float gaussian_std_dev = 0.683;
  
  TLegend *tlgd = new TLegend(0.7, 0.7, 1.0, 0.99);  

  bool bl_beam_duration_comparison = false;
  bool bl_beam_size_comparison = true;
  bool bl_beam_size_electron = false;
  
  for(int iSize = 0; iSize < 3; iSize++)
  {
    float beam_size = 0;
    beam_size = 1e-6;

    if(bl_beam_size_comparison)
    {	    
      if(iSize==0) beam_size = 1e-6;     
      else if(iSize == 1) beam_size = 10e-6;
      else if(iSize == 2) beam_size = 100e-6;
    
      if(bl_beam_size_electron)
      {
        elec_beam.set_beam_size(beam_size, beam_size, elec_beam.get_beam_duration(), beta* speed_light); 
        laser_beam.set_beam_size(1e-6, 1e-6, laser_beam.get_beam_duration(), speed_light);
      }
      else
      {
        elec_beam.set_beam_size(1e-6, 1e-6, elec_beam.get_beam_duration(), beta* speed_light);     
        laser_beam.set_beam_size(beam_size, beam_size, laser_beam.get_beam_duration(), speed_light);
      }
       std::cout <<" laser beam size: " << laser_beam.beam_size[0] << ", " << laser_beam.beam_size[1] << ", " << laser_beam.get_beam_duration() << "x" << speed_light <<" = " << laser_beam.beam_size[2] << std::endl;

    }
    else if(bl_beam_duration_comparison)
    {	    
      float var_duration = 0.;
      if(iSize==0) 
      {	    
        var_duration = 1e-15;
      }
      else if(iSize == 1)
      { 
        var_duration = 10e-15;
      }
      else if(iSize == 2)
      {	    
        var_duration = 100e-15;    
      }

      elec_beam.set_beam_duration(var_duration);      
    
      elec_beam.set_beam_size(beam_size, beam_size, elec_beam.get_beam_duration(), beta* speed_light);
      laser_beam.set_beam_size(beam_size, beam_size, laser_beam.get_beam_duration(), speed_light);    
    }

    // For sigma y and sigma prime y
    float unc_sigma_y = 0.;
    float unc_sigma_prime_y = 0.;
    float unc_total_sigma_y = 0.;
    
    // For sigma x and sigma prime x
    float unc_sigma_x = 0.;
    float unc_sigma_prime_x = 0.;
    float unc_total_sigma_x = 0.;
    
    for(int iBin = 0; iBin < 9; iBin++)
    {  
      theta = (float)iBin* 10;

//    theta = theta*M_PI / 180.0;
      geom_factor = cal_geom_factor(beta, elec_beam.beam_size, laser_beam.beam_size, theta);
  
      lumi = pow_factor * geom_factor;
      lumi = lumi / (barn_unit_conversion);
      std::cout << lumi << std::endl;

//////******************************
//
// Consider space jitter, uncertainty on beam size, in transeverse plane.
// Jitter is included in the equation of geometry factor.
// Two uncertainties which are on x and y are considered.
 
      unc_sigma_y = gaussian_std_dev*elec_beam.beam_size[1];
      unc_sigma_prime_y = gaussian_std_dev*laser_beam.beam_size[1];
  
      unc_total_sigma_y = cal_error_propagation_typeone(1, elec_beam.beam_size[1], unc_sigma_y, 1, laser_beam.beam_size[1], unc_sigma_prime_y );
    std::cout << "uncertainty on sigma_y: " << unc_total_sigma_y << std::endl;
      
      unc_sigma_x = gaussian_std_dev*elec_beam.beam_size[0]; 
      unc_sigma_prime_x = gaussian_std_dev*laser_beam.beam_size[0];
 
      float theta_in_rad = theta*M_PI / 180.0;
      unc_total_sigma_x = cal_error_propagation_typeone(pow((beta+cos(theta_in_rad)),2), elec_beam.beam_size[0], unc_sigma_x, pow((1.+beta*cos(theta_in_rad)),2), laser_beam.beam_size[0], unc_sigma_prime_x );
  
      std::cout << "uncertainty on sigma_x: " << unc_total_sigma_x << std::endl;

      float constant = 0.;
      constant = (1+beta*cos(theta_in_rad)) / ( 2.*M_PI);
      float var1 = sqrt(pow(elec_beam.beam_size[1],2)+pow(laser_beam.beam_size[1],2));
      float var2 = sqrt( pow(elec_beam.beam_size[0],2)*pow(beta+cos(theta_in_rad),2) + pow(laser_beam.beam_size[0],2)*pow(1+beta*cos(theta_in_rad),2) + (pow(elec_beam.beam_size[2],2)+pow(laser_beam.beam_size[2],2))*pow(sin(theta_in_rad),2)); 

      float unc_geom = 0;

      unc_geom = cal_error_propagation_geom(constant, var1, unc_total_sigma_y, constant, var2, unc_total_sigma_x );
  std::cout << "uncertainty on geom: " << unc_geom << std::endl;

      float unc_lumi = 0.; 
      unc_lumi = lumi * unc_geom / geom_factor; 
      std::cout << "lumi +/- uncertainty: " << lumi << " +/- " << unc_lumi << std::endl;

      th1_lumi[iSize]->Fill(theta, lumi);
      th1_lumi[iSize]->SetBinError(iBin+1, unc_lumi);

      double scat_photons = lumi* cal_integral_x_section_Thomson() * barn_unit_conversion;
     
      std::cout << "Thomson cross section: " << cal_integral_x_section_Thomson() << std::endl;
      std::cout << "number of scattered photons (using Thomas cross section): " << scat_photons << std::endl;

      scat_photons = lumi * cal_integral_x_section_KN(gamma_elec) * barn_unit_conversion;
      std::cout << "gamma: " << gamma_elec << std::endl;
      std::cout << "integral cross section: " << cal_integral_x_section_KN(gamma_elec) << std::endl;
      std::cout << "number of scattered photons (using K-N cross section): " << scat_photons << std::endl;
    }
    if(bl_beam_size_comparison)
    {	    
      if(bl_beam_size_electron)
      {
        if(iSize==0)
	{
          tlgd->AddEntry((TObject*)0, "electron beam transverse size:", "");
	  tlgd->AddEntry(th1_lumi[iSize], Form("1 #pm %.2f #mum", unc_sigma_y*1e6 ));
	}
        else if(iSize==1) tlgd->AddEntry(th1_lumi[iSize], Form("10 #pm %.2f #mum", unc_sigma_y*1e6 ));
        else if(iSize==2)
	{
          tlgd->AddEntry(th1_lumi[iSize], Form("100 #pm %.2f #mum", unc_sigma_y*1e6 ));
          tlgd->AddEntry((TObject*)0, Form("electron beam longitudinal size"), "");
          tlgd->AddEntry((TObject*)0, Form("%.0f #pm 0 #mum", elec_beam.beam_size[2]*1e6 ), "");
	}
      }
      else
      {
        if(iSize==0)
	{
          tlgd->AddEntry((TObject*)0, "laser transverse size:", "");
	  tlgd->AddEntry(th1_lumi[iSize], Form("1 #pm %.2f #mum", unc_sigma_prime_y*1e6 ));
	}
        else if(iSize==1) tlgd->AddEntry(th1_lumi[iSize], Form("10 #pm %.2f #mum", unc_sigma_prime_y*1e6 ));      
        else if(iSize==2)
	{
	  tlgd->AddEntry(th1_lumi[iSize], Form("100 #pm %.2f #mum", unc_sigma_prime_y*1e6 ));
          tlgd->AddEntry((TObject*)0, Form("laser longitudinal size"), "");
          tlgd->AddEntry((TObject*)0, Form("%.0f #pm 0 #mum", laser_beam.beam_size[2]*1e6 ), "");
	}
      }

    }
    else if(bl_beam_duration_comparison)
    {	 
       if(iSize==0)
       {
         tlgd->AddEntry((TObject*)0, "electron beam duration:", "");
         tlgd->AddEntry(th1_lumi[iSize], Form("%.0f fs", 1e15*elec_beam.get_beam_duration() ));
       }
       else if(iSize==1) tlgd->AddEntry(th1_lumi[iSize], Form("%.0f fs", 1e15*elec_beam.get_beam_duration() ));
       else if(iSize==2) tlgd->AddEntry(th1_lumi[iSize], Form("%.0f fs", 1e15*elec_beam.get_beam_duration() ));
    }
  }
  


   TCanvas *c2 = new TCanvas("c2","",200, 10, 800, 600);
   c2->SetLogy();
   c2->SetLeftMargin(0.15);

   th1_lumi[0]->Draw();
   if(bl_beam_size_comparison) th1_lumi[0]->SetAxisRange(1e-2, 1e5, "Y");
   else if(bl_beam_duration_comparison) th1_lumi[0]->SetAxisRange(1e2, 5e3, "Y");
   
   th1_lumi[1]->SetMarkerColor(2);
   th1_lumi[1]->SetLineColor(2);
   th1_lumi[1]->Draw("SAME");
   
   th1_lumi[2]->SetMarkerColor(8);
   th1_lumi[2]->SetLineColor(8);
   th1_lumi[2]->Draw("SAME");
   
   
   tlgd->Draw("SAME"); 
   
   
   TLatex beam_info;
   beam_info.SetTextSize(0.04);
   if(bl_beam_size_comparison)
   {
     beam_info.DrawLatex(2,3e5, Form("e beam: 100 MeV, 10 pC, 100 fs"));
     beam_info.DrawLatex(2,1.5e5, "laser: 1030 nm , 1 J, 30 fs");
     //beam_info.DrawLatex(5,5e4, Form("laser beam size: 1 #pm 0.68 #mum, %.0f #pm 0 #mum", laser_beam.beam_size[2]*1e6 ));
     beam_info.DrawLatex(2,5e4, Form("electron beam size: 1 #pm 0.68 #mum, %.0f #pm 0 #mum", elec_beam.beam_size[2]*1e6 ));
   }
   else if(bl_beam_duration_comparison)	
   {
     beam_info.DrawLatex(2,5e3, Form("e beam: 100 MeV, 10 pC, %.0f #mum ",elec_beam.beam_size[0]*1e6));
     beam_info.DrawLatex(2,4.7e3, "laser: 1030 nm , 1 J, 30 fs");
//     beam_info.DrawLatex(2,4.5e3, Form("electron beam size: 1 #pm 0.68 #mum, %.0f #pm 0 #mum", elec_beam.beam_size[2]*1e6 ));
   }

   
   
   c2->Print("LCS_lumi_vs_theta.pdf");

   delete tlgd, c2, th1_lumi;
  
return 0;
}

float cal_pow_factor(float num_electrons, float num_photons)
{
    float pow_factor = 0.;
    pow_factor = num_electrons * num_photons;
    return pow_factor;    
}

float cal_geom_factor(float rela_velo_elec, float *beam_one_size, float *beam_two_size, float theta)
{
    float geom_factor = 0.;
// unit of the theta here is degree so need to change it into rad.
    float theta_in_pi = theta*M_PI / 180.0;
    
    float sigma[3] = {0};
    sigma[0] = beam_one_size[0];
    sigma[1] = beam_one_size[1];
    sigma[2] = beam_one_size[2];
    
    float sigma_prime[3] = {0};
    sigma_prime[0] = beam_two_size[0];
    sigma_prime[1] = beam_two_size[1];
    sigma_prime[2] = beam_two_size[2];
    
    float beta = rela_velo_elec;
    
    
    float numerator = (1. + beta*cos(theta_in_pi));        
    float denominator = (2 * M_PI * sqrt(pow(sigma[1],2)+pow(sigma_prime[1],2)) * sqrt( pow(sigma[0],2)*pow(beta+cos(theta_in_pi),2) + pow(sigma_prime[0],2)*pow(1+beta*cos(theta_in_pi),2) + (pow(sigma[2],2)+pow(sigma_prime[2],2))*pow(sin(theta_in_pi),2)) );
    geom_factor = numerator / denominator;
    
    return geom_factor;
} 

float cal_error_propagation_typeone(float const_a, float var_AA, float err_var_AA, float const_b, float var_BB, float err_var_BB)
{
//////    
//
//  this calculation of error propagation is for square root of a sum of square
//  f = sqrt(aA^2 + bB^2)
//  err_f^2 ~= (A/f)^2 a^2 sigma_A^2 + (B/f)^2 b^2 sigma_B^2 +/- 2ab(AB/(f^2)) sigma_AB

 float f = sqrt( (const_a * var_AA * var_AA) + (const_b * var_BB * var_BB ) ); 
 float err_f = sqrt( ( pow( (var_AA/f),2)* pow(const_a,2) * pow(err_var_AA,2) ) + ( pow(var_BB/f, 2)* pow(const_b,2) * pow(err_var_BB,2) ) );

 return err_f;
}

float cal_error_propagation_geom(float const_a, float var1, float err_var1, float const_b, float var2, float err_var2)
{

  float geom = const_a / (var1 * var2); 
  float err_geom = sqrt( pow(const_a / var2 * (-1.) * 1. / (var1*var1) * err_var1, 2) + pow( const_b / var1 * (-1.) * 1. / (var2*var2) * err_var2 ,2));

  return err_geom;
}	

float cal_integral_x_section_Thomson(void)
{
  float classic_electron_radius = 2.81794e-15;
  return (8.*M_PI/3.) * classic_electron_radius * classic_electron_radius;
}


float cal_integral_x_section_KN(float r)
{
  float classic_electron_radius = 2.81794e-15;
  float sigma_KN = 2.* M_PI * classic_electron_radius * classic_electron_radius;
  float fst_term = ((1.+r)/r/r)*( (2.*(1.+r)/(1.+2.*r)) - (log(1.+2.*r) /r) ); 
  float snd_term = log(1.+2*r) / (2.*r);
  float thr_term = (1.+3*r)/pow(1.+2.*r,2); 
  std::cout << fst_term << ", " << snd_term << ", " << thr_term << std::endl;
  sigma_KN =  sigma_KN*(fst_term + snd_term - thr_term); 
  
  return sigma_KN;
}
