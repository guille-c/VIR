/* Mathematical and physical constants */

#define PI     3.1415926535897932384626433832795028841971693993751
#define TWOPI  (2.0*PI)
#define HALFPI (PI/2.0)
#define RPDEG  (PI/180.0)
#define RPARCM (PI/(180.0*60.0))
#define RPARCS (PI/(180.0*3600.0))
#define RPHR   (PI/12.0)
#define RPMIN  (RPHR/60.0)
#define RPSEC  (RPHR/3600.0)
#define LN2    0.69314718055994530942     /* log_e(2) */

/* Physical Constants: primary reference: The Review of Particle
 * Physics, D. E Groom et al., The European Physical Journal, C15, 1
 * (2000)
 * http://pdg.lbl.gov/
 */

#define LIGHTSPEED 2.99792458E8   /* speed of light, m/s; defined */
#define DEGCTOK    273.15         /* add to Celsius temp to get
				   * abs temp in K */
#define TCMB       2.725          /* T(cmb) in K [Mather et al.
				   * 1999, ApJ, 512, 511 */
#define PLANCK_H   6.62606876e-34 /* Planck constant in J.s */
#define BOLTZ_K    1.3806503e-23  /* Boltzmann constant in J/K */
