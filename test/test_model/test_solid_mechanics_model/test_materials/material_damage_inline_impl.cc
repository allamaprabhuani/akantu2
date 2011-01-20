/**
 * @file   material_damage_inline_impl.cc
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Tue Jul 27 11:57:43 2010
 *
 * @brief  Implementation of the inline functions of the material damage
 *
 * @section LICENSE
 *
 * \<insert license here\>
 *
 */

/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
inline void MaterialDamage::computeStress(Real * F, Real * sigma, Real & dam) {

  Real trace = F[0] + F[4] + F[8]; /// \F_{11} + \F_{22} + \F_{33}
  /// \sigma_{ij} = \lamda * \F_{kk} * \delta_{ij} + 2 * \mu * \F_{ij}
   sigma[0] = lambda * trace + 2*mu*F[0];
  sigma[4] = lambda * trace + 2*mu*F[4];
  sigma[8] = lambda * trace + 2*mu*F[8];

  sigma[1] = sigma[3] =  mu * (F[1] + F[3]);
  sigma[2] = sigma[6] =  mu * (F[2] + F[6]);
  sigma[5] = sigma[7] =  mu * (F[5] + F[7]);

  Real Y =
    sigma[0]*F[0] +
    sigma[1]*F[1] +
    sigma[2]*F[2] +
    sigma[3]*F[3] +
    sigma[4]*F[4] +
    sigma[5]*F[5] +
    sigma[6]*F[6] +
    sigma[7]*F[7] +
    sigma[8]*F[8];

  Y *= 0.5;

  Real Fd = Y - Yd - Sd*dam;

  if (Fd > 0) dam = (Y - Yd) / Sd;
  dam = std::min(dam,1.);

  sigma[0] *= 1-dam;
  sigma[4] *= 1-dam;
  sigma[8] *= 1-dam;
  sigma[1] *= 1-dam;
  sigma[3] *= 1-dam;
  sigma[2] *= 1-dam;
  sigma[6] *= 1-dam;
  sigma[5] *= 1-dam;
  sigma[7] *= 1-dam;
}

/* -------------------------------------------------------------------------- */
inline void MaterialDamage::computePotentialEnergy(Real * F, Real * sigma, Real * epot) {
  *epot = 0.;
  for (UInt i = 0, t = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j, ++t)
      *epot += sigma[t] * (F[t] - (i == j));
  *epot *= .5;
}

/* -------------------------------------------------------------------------- */
inline Real MaterialDamage::celerity() {
  return sqrt(E/rho);
}

/* -------------------------------------------------------------------------- */
inline Real MaterialDamage::getStableTimeStep(Real h) {
  return (h/celerity());
}
