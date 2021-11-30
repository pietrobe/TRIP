#include "Rotation_matrix.hpp"

////////////////////////////////////////////////////////////////
// rotation_matrix::
// build_matrix
////////////////////////////////////////////////////////////////
void Rotation_matrix::build_matrix(double alpha, double beta, double gamma) {

  p_alpha = alpha;
  p_beta  = beta;
  p_gamma = gamma;

  p_d_mn_1.resize(3);
  for (int i = 0; i < 3; ++i) { p_d_mn_1[i].resize(3); }
  
  p_d_mn_2.resize(5);
  for (int i = 0; i < 5; ++i) { p_d_mn_2[i].resize(5); }

  double C, S;

  C = cos(beta);  // I'm using the notation proposed in the book
  S = sin(beta);

  {
    // Build matrix d_mn_1
    p_d_mn_1[1][1] = C;

    double d11     = 1.0 / 2.0 * (1.0 + C);
    p_d_mn_1[2][2] = d11;
    p_d_mn_1[0][0] = d11;

    double ds      = 1.0 / sqrt(2.0) * S;
    p_d_mn_1[1][0] = -ds;
    p_d_mn_1[0][1] = ds;

    p_d_mn_1[2][1] = -ds;
    p_d_mn_1[1][2] = ds;

    double d02     = 1.0 / 2.0 * (1.0 - C);
    p_d_mn_1[2][0] = d02;
    p_d_mn_1[0][2] = d02;
  }

  {
    // Build matrix d_mn_2
    // see eq 2.69

    double cbh = std::cos(beta / 2.0);
    double sbh = std::sin(beta / 2.0);

    double fact[32] = {1.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 6.000000000000000e+00,
                       2.400000000000000e+01, 1.200000000000000e+02, 7.200000000000000e+02, 5.040000000000000e+03,
                       4.032000000000000e+04, 3.628800000000000e+05, 3.628800000000000e+06, 3.991680000000000e+07,
                       4.790016000000000e+08, 6.227020800000000e+09, 8.717829120000000e+10, 1.307674368000000e+12,
                       2.092278988800000e+13, 3.556874280960000e+14, 6.402373705728000e+15, 1.216451004088320e+17,
                       2.432902008176640e+18, 5.109094217170944e+19, 1.124000727777608e+21, 2.585201673888498e+22,
                       6.204484017332394e+23, 1.551121004333098e+25, 4.032914611266057e+26, 1.088886945041835e+28,
                       3.048883446117138e+29, 8.841761993739701e+30, 2.652528598121911e+32, 8.222838654177924e+33};

    int J = 2;

    for (int M = -J; M <= J; ++M) {
      for (int N = -J; N <= J; ++N) {

        double sqf = std::sqrt(fact[J + M] * fact[J - M] * fact[J + N] * fact[J - N]);


        double sum = 0.0;
        
        // int tmin = MAX(0, M - N);
        // int tmax = MIN(J + M, J - N);

        // variable = (condition) ? expressionTrue : expressionFalse;
        int tmin = (0 > M - N) ? 0 : M - N;
        int tmax = (M < - N) ? J + M : J - N;

        for (int t = tmin; t <= tmax; ++t) {

          double den = fact[J + M - t] * fact[J - N - t] * fact[t] * fact[t + N - M];
          
          sum = sum + double(std::pow(-1, t)) * (std::pow(cbh, 2 * J + M - N - 2 * t) * std::pow(sbh, 2 * t - M + N)) / den;
        }

        p_d_mn_2[M + 2][N + 2] = sqf * sum;
      }
    }
  }
}

////////////////////////////////////////////////////////////////
// rotation_matrix::
// operator()
////////////////////////////////////////////////////////////////
std::complex<double> Rotation_matrix::operator()(int J, int M, int N) const {

  complex<double> D(0.0, 0.0);  // the resulting value of the rotation matrix.

  const complex<double> I(0.0, 1.0);  // imag unit.

  D = std::exp((-1.0) * I * (p_alpha * double(M) + p_gamma * double(N)));

  if (J == 1) {
    D = D * p_d_mn_1[M + 1][N + 1];
  } else if (J == 2) {
    D = D * p_d_mn_2[M + 2][N + 2];
  } else {
    D = D * complex<double>(1.0, 0.0);
  }
  
  return D;
}
