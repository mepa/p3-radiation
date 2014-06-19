#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <strstream>
#include <fstream>

using namespace std;

const double h = 6.626068e-27;
const double c = 2.99792458e10;
const double k = 1.3806503e-16;

const double Mb = 1.0e-18;
const double eV = 1.60217646e-12;

// limits of integration
const double lambdaHminus = h * c / 0.755 / eV; //1.6e4 Angstroms
const double lambdaH2 = h * c / 11.2 / eV; //1.1e3
const double lambdaH = h * c / 13.6 / eV; //9.1e2
const double lambdaHe = h * c / 24.59 / eV; //5.0e2
const double lambdaHeplus = h * c / 54.42 / eV; //2.3e2
const double lambdaMax = h * c / 50000 / eV; //2.5e-2

const double E_H2 = 11.2 * eV;
const double E_H = 13.6 * eV;
const double E_He = 24.59 * eV;
const double E_Heplus = 54.42 * eV;

inline double sq(double x) { return x * x; }

// photoionization cross sections from Verner et al. 1996
double sigma_H(double lambda)
{
  if(lambda > lambdaH) return 0.0;
  const double E = h * c / lambda;
  const double sigma0 = 5.475e4 * Mb;
  const double E0 = 0.4298 * eV;
  const double yw = 0;
  const double ya = 32.88;
  const double P = 2.963;
  const double y0 = 0;
  const double y1 = 0;
  const double x = E / E0 - y0;
  const double y = sqrt(sq(x) + sq(y1));
  const double F = (sq(x - 1) + sq(yw)) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P);
  return sigma0 * F;
}

double sigma_He(double lambda)
{
  if(lambda > lambdaHe) return 0.0;
  const double E = h * c / lambda;
  const double sigma0 = 9.492e2 * Mb;
  const double E0 = 13.61 * eV;
  const double yw = 2.039;
  const double ya = 1.469;
  const double P = 3.188;
  const double y0 = 0.4434;
  const double y1 = 2.136;
  const double x = E / E0 - y0;
  const double y = sqrt(sq(x) + sq(y1));
  const double F = (sq(x - 1) + sq(yw)) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P);
  return sigma0 * F;
}

double sigma_Heplus(double lambda)
{
  if(lambda > lambdaHeplus) return 0.0;
  const double E = h * c / lambda;
  const double sigma0 = 1.369e4 * Mb;
  const double E0 = 1.72 * eV;
  const double yw = 0;
  const double ya = 32.88;
  const double P = 2.963;
  const double y0 = 0;
  const double y1 = 0;
  const double x = E / E0 - y0;
  const double y = sqrt(sq(x) + sq(y1));
  const double F = (sq(x - 1) + sq(yw)) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P);
  return sigma0 * F;
}

double sigma_Hminus(double lambda)
{
  if(lambda > lambdaHminus) return 0.0;
  const double E = h * c / lambda;
  const double E_in_eV = E / eV;
  return 2.1e-16 * pow(E_in_eV - 0.75, 1.5) / pow(E_in_eV, 3.11);
}

double sigma_H2(double lambda)
{
  if(lambda > lambdaH2) return 0.0;
  return 3.78e-18;
}

void readHeader(istream & in, double & age, size_t & lambda_points)
{
 string line, word;
 getline(cin, line); // cout << line << endl;
 if(line != string("")) { cout << "A1: " << line << endl; abort(); }
 in >> word; // cout << word;
 if(word != string("Age")) { cout << "A2: word" << endl; abort(); }
 in >> word >> age; // cout << word << " " << age;
 in.ignore(1024, '\n'); // cout << endl;
 getline(in, line); // cout << line << endl;
 getline(in, line); // cout << line << endl;
 in >> word; // cout << word;
 if(word != string("Wavelength")) { cout << "A3: " << word << endl; abort(); }
 in >> word; // cout << word;
 if(word != string("data")) { cout << "A4: " << word << endl; abort(); }
 in >> word >> lambda_points; // cout << word << " " << lambda_points;
 in.ignore(1024, '\n'); // cout << endl;
 getline(in, line); // cout << line << endl;
 getline(in, line); // cout << line << endl;
 getline(in, line); // cout << line << endl;
}

int main(int argc, char *argv[])
{

  string line;
  getline(cin, line);

  double age = 0;
  size_t lambda_points;

  vector<double> lambda, flux;
  
  if (argc != 2) 
    {
      cout << "usage: flux age_in_Myr < inputfile > outputfile" << endl;
      return 0;
    }

  const double target_age = atof(argv[1]);

  readHeader(cin, age, lambda_points);

  double spectrum_age = -1.0;

  while(age < target_age)
    {
      lambda.resize(lambda_points);
      flux.resize(lambda_points);

      //cout << "age = " << age << endl;
      //cout << "lambda_points = " << lambda_points << endl;

      for(size_t i = 0; i < lambda_points; i ++)
	{
	  cin >> lambda[i] >> flux[i];

	  lambda[i] *= 1.0e-8;
	  flux[i] *= 1.0e8;
	}
      spectrum_age = age;
      
      readHeader(cin, age, lambda_points);
    }
  cout << endl;
  cout << "target_age = " << target_age << " Myr" << endl;
  cout << "spectrum_age = " << spectrum_age << " Myr" << endl;


  const size_t bins = 5;

  const double bin_max[bins] = { lambdaHminus, lambdaH2, lambdaH, lambdaHe, lambdaHeplus };
  const double bin_min[bins] = { lambdaH2, lambdaH, lambdaHe, lambdaHeplus, lambdaMax };

  double N_dot_Hminus[bins];
  double N_dot_H[bins];
  double N_dot_He[bins];
  double N_dot_Heplus[bins];
  double N_dot_H2[bins];

  double Gamma_Hminus[bins];
  double Gamma_H[bins];
  double Gamma_He[bins];
  double Gamma_Heplus[bins];

  double GammaE_H[bins];
  double GammaE_He[bins];
  double GammaE_Heplus[bins];
  
  for(size_t bin = 0; bin < bins; bin ++)
    {
      N_dot_Hminus[bin] = 0.0;
      N_dot_H[bin] = 0.0;
      N_dot_He[bin] = 0.0;
      N_dot_Heplus[bin] = 0.0;
      N_dot_H2[bin] = 0.0;

      Gamma_Hminus[bin] = 0.0;
      Gamma_H[bin] = 0.0;
      Gamma_He[bin] = 0.0;
      Gamma_Heplus[bin] = 0.0;

      GammaE_H[bin] = 0.0;
      GammaE_He[bin] = 0.0;
      GammaE_Heplus[bin] = 0.0;
      
      //cout << endl;
      //cout << bin << ": bin_min = " << bin_min[bin] << endl;
      //cout << bin << ": bin_max = " << bin_max[bin] << endl;
    }

  for(size_t i = 1; i < lambda_points; i ++)
    {
      const double delta_lambda = lambda[i] - lambda[i - 1];
      const double lambda_1 = lambda[i - 1];
      const double lambda_2 = lambda[i];
      const double E_1 = h * c / lambda_1;
      const double E_2 = h * c / lambda_2;
      const double flux_1 = flux[i - 1];
      const double flux_2 = flux[i];
      const double N_dot_1 = flux_1 / E_1;
      const double N_dot_2 = flux_2 / E_2;

      for(size_t bin = 0; bin < bins; bin ++)
	{
	  if(lambda[i - 1] >= bin_min[bin] && lambda[i] <= bin_max[bin])
	    {
	      N_dot_Hminus[bin] += delta_lambda * 0.5 * (N_dot_1 + N_dot_2);

	      Gamma_Hminus[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_Hminus(lambda_1) + N_dot_2 * sigma_H(lambda_2));

	      N_dot_H2[bin] += delta_lambda * 0.5 * (N_dot_1 + N_dot_2); // LW photons

	      N_dot_H[bin] += delta_lambda * 0.5 * (N_dot_1 + N_dot_2);

	      Gamma_H[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_H(lambda_1) + N_dot_2 * sigma_H(lambda_2));	 

	      GammaE_H[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_H(lambda_1) * (E_1 - E_H) + N_dot_2 * sigma_H(lambda_2) * (E_2 - E_H));	

	      N_dot_He[bin] += delta_lambda * 0.5 * (N_dot_1 + N_dot_2);

	      Gamma_He[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_He(lambda_1) + N_dot_2 * sigma_He(lambda_2));
	      
	      GammaE_He[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_He(lambda_1) * (E_1 - E_He) + N_dot_2 * sigma_He(lambda_2) * (E_2 - E_He));  

	      N_dot_Heplus[bin] += delta_lambda * 0.5 * (N_dot_1 + N_dot_2);

	      Gamma_Heplus[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_Heplus(lambda_1) + N_dot_2 * sigma_Heplus(lambda_2)); 

	      GammaE_Heplus[bin] += delta_lambda * 0.5 * (N_dot_1 * sigma_Heplus(lambda_1) * (E_1 - E_Heplus) + N_dot_2 * sigma_Heplus(lambda_2) * (E_2 - E_Heplus));
	    }
     
	}
    }

  double sigma_Hminus_average[bins], sigma_H_average[bins], sigma_He_average[bins], sigma_Heplus_average[bins];
  double sigmaE_H_average[bins], sigmaE_He_average[bins], sigmaE_Heplus_average[bins];  
  double DeltaE_H_average[bins], DeltaE_He_average[bins], DeltaE_Heplus_average[bins];

      
  for(size_t bin = 0; bin < bins; bin ++)
    {
      sigma_Hminus_average[bin] = N_dot_Hminus[bin] > 0.0 ? Gamma_Hminus[bin] / N_dot_Hminus[bin] : 0.0;
      sigma_H_average[bin] = N_dot_H[bin] > 0.0 ? Gamma_H[bin] / N_dot_H[bin] : 0.0;
      sigma_He_average[bin] = N_dot_He[bin] > 0.0 ? Gamma_He[bin] / N_dot_He[bin] : 0.0;
      sigma_Heplus_average[bin] = N_dot_Heplus[bin] > 0.0 ? Gamma_Heplus[bin] / N_dot_Heplus[bin] : 0.0;

      sigmaE_H_average[bin] = N_dot_H[bin] > 0.0 ? GammaE_H[bin] / N_dot_H[bin] : 0.0;
      sigmaE_He_average[bin] = N_dot_He[bin] > 0.0 ? GammaE_He[bin] / N_dot_He[bin] : 0.0;
      sigmaE_Heplus_average[bin] = N_dot_Heplus[bin] > 0.0 ? GammaE_Heplus[bin] / N_dot_Heplus[bin] : 0.0;

      DeltaE_H_average[bin] = sigma_H_average[bin] > 0.0 ? sigmaE_H_average[bin] / sigma_H_average[bin] : 0.0;
      DeltaE_He_average[bin] = sigma_He_average[bin] > 0.0 ? sigmaE_He_average[bin] / sigma_He_average[bin] : 0.0;
      DeltaE_Heplus_average[bin] = sigma_Heplus_average[bin] > 0.0 ? sigmaE_Heplus_average[bin] / sigma_Heplus_average[bin] : 0.0;

      DeltaE_H_average[bin] /= eV;
      DeltaE_He_average[bin] /= eV;
      DeltaE_Heplus_average[bin] /= eV;
    }

  //const double Gamma_H_reference = 1.0e-12;
  const double Gamma_H_reference = 0.3e-12;

  const double Gamma_H_total_old = Gamma_H[2] + Gamma_H[3] + Gamma_H[4];

  const double factor = Gamma_H_reference / Gamma_H_total_old;

  const double N_dot_H_total = (N_dot_H[2] + N_dot_H[3] + N_dot_H[4]) * factor;
  const double N_dot_He_total = (N_dot_He[3] + N_dot_He[4]) * factor;
  const double N_dot_Heplus_total = (N_dot_Heplus[4]) * factor;
  const double N_dot_Hminus_total = (N_dot_Hminus[0] + N_dot_Hminus[1]) * factor;
  const double N_dot_H2_total = (N_dot_H2[1]) * factor;
  
  const double Gamma_H_total = (Gamma_H[2] + Gamma_H[3] + Gamma_H[4]) * factor;
  const double Gamma_He_total = (Gamma_He[3] + Gamma_He[4]) * factor;
  const double Gamma_Heplus_total = (Gamma_Heplus[4]) * factor;

  cout << endl;

  cout << "sigma_Hminus = " << sigma_Hminus_average[0] << " " << sigma_Hminus_average[1] << " " << sigma_Hminus_average[2] << " " << sigma_Hminus_average[3] << " " << sigma_Hminus_average[4] << endl;
  cout << "sigma_H = " << sigma_H_average[0] << " " << sigma_H_average[1] << " " << sigma_H_average[2] << " " << sigma_H_average[3] << " " << sigma_H_average[4] << endl;
  cout << "sigma_He = " << sigma_He_average[0] << " " << sigma_He_average[1] << " " << sigma_He_average[2] << " " << sigma_He_average[3] << " " << sigma_He_average[4] << endl;
  cout << "sigma_Heplus = " << sigma_Heplus_average[0] << " " << sigma_Heplus_average[1] << " " << sigma_Heplus_average[2] << " " << sigma_Heplus_average[3] << " " << sigma_Heplus_average[4] << endl;
   
  cout << endl;

  cout << "sigmaE_H = " << sigmaE_H_average[0] << " " << sigmaE_H_average[1] << " " << sigmaE_H_average[2] << " " << sigmaE_H_average[3] << " " << sigmaE_H_average[4] << endl;
  cout << "sigmaE_He = " << sigmaE_He_average[0] << " " << sigmaE_He_average[1] << " " << sigmaE_He_average[2] << " " << sigmaE_He_average[3] << " " << sigmaE_He_average[4] << endl;
  cout << "sigmaE_Heplus = " << sigmaE_Heplus_average[0] << " " << sigmaE_Heplus_average[1] << " " << sigmaE_Heplus_average[2] << " " << sigmaE_Heplus_average[3] << " " << sigmaE_Heplus_average[4] << endl;
  
  cout << endl;

  cout << "DeltaE_H = " << DeltaE_H_average[0] << " " << DeltaE_H_average[1] << " " << DeltaE_H_average[2] << " " << DeltaE_H_average[3] << " " << DeltaE_H_average[4] << endl;
  cout << "DeltaE_He = " << DeltaE_He_average[0] << " " << DeltaE_He_average[1] << " " << DeltaE_He_average[2] << " " << DeltaE_He_average[3] << " " << DeltaE_He_average[4] << endl;
  cout << "DeltaE_Heplus = " << DeltaE_Heplus_average[0] << " " << DeltaE_Heplus_average[1] << " " << DeltaE_Heplus_average[2] << " " << DeltaE_Heplus_average[3] << " " << DeltaE_Heplus_average[4] << endl;
  
  cout << endl;

  cout << "N_dot_Hminus_total = " << N_dot_Hminus_total << endl;
  cout << "N_dot_H_total = " << N_dot_H_total << endl;
  cout << "N_dot_He_total = " << N_dot_He_total << endl;
  cout << "N_dot_Heplus_total = " << N_dot_Heplus_total << endl;
  cout << "N_dot_H2 = " << N_dot_H2_total << endl;

  cout << endl;

  cout << "J21_LW = " << 1.0e21 * N_dot_H2_total * 0.5 * (E_H2 + E_H) / (4.0 * M_PI * (E_H - E_H2) / h) << endl;

  cout << endl;

  cout << "Gamma_H_total = " << Gamma_H_total << endl;
  cout << "Gamma_He_total = " << Gamma_He_total << endl;
  cout << "Gamma_Heplus_total = " << Gamma_Heplus_total << endl;

  cout << endl;

  //cout << N_dot_Heplus[0] << " " << N_dot_Heplus[1] << " " << N_dot_Heplus[2] << " " << N_dot_Heplus[3] << " " << N_dot_Heplus[4] << endl;

  ostrstream filename;
  filename << "spec_" << spectrum_age << ".out" << ends; 
  ofstream spec_out(filename.str());
  for(size_t i = 0; i < lambda.size(); i ++)
    {
      spec_out << lambda[i] << " " << flux[i] * factor << endl;
    }
  
  return 0;
}

