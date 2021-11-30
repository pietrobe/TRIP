#ifndef __ROTATION_MATRIX__
#define __ROTATION_MATRIX__

#include "Utilities.hpp"

typedef std::vector<std::vector<std::complex<double > > > complex_matrix;

/**
 * @class rotation_matrix
 * @brief Define the quantum eigenstates rotation matrix described in
 * chap 2 of Polarization in Spectral Lines M. [Landi Degl'Innocenti,M. Landolfi]
 * @ingroup group02_Rii_emission_coeff
 */
class Rotation_matrix {
 public:
  /**
   * @brief Default contructor
   */
  Rotation_matrix() { Rotation_matrix(0.0, 0.0, 0.0); }

  /**
   * @brief Constructor.
   *
   * @param[in] alpha : the rotation angle alpha
   * @param[in] beta : the rotation angle beta
   * @param[in] gamma : the rotation angle gamma
   */
  Rotation_matrix(double alpha, double beta, double gamma) { build_matrix(alpha, beta, gamma); }

  /**
   * @brief Given alpha, beta and gamma build the matrix
   *
   * @param[in] : the rotation angle alpha
   * @param[in] : the rotation angle beta
   * @param[in] : the rotation angle gamma
   */
  void build_matrix(double alpha, double beta, double gamma);

  /**
   * @brief return the value of the rotation matrix related to the quantum state J, M, N
   */
  std::complex<double> operator()(int J, int M, int N) const;

 private:

  double p_alpha, p_beta, p_gamma;

  /**
   * @brief Table \f$ d_{MN}^1 \f$ described in table 2.1
   */
  complex_matrix p_d_mn_1;

  /**
   * @brief Table \f$ d_{MN}^2 \f$ described in table 2.1
   */
  complex_matrix p_d_mn_2;
};


#endif  // __ROTATION_MATRIX__