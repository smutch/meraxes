#ifdef CALC_MAGS

#include "meraxes.h"
#include <math.h>

// Added for dust implementation (Kauffmann)
#define tauBstar 0.8
#define beta 0.5
#define Lbstar 3.7e8 // In solar luminosity is 5.75e10
#define mu_dust (1 / 3.)

// extinction curve of Cardelli et al. (1989)
// ( http://adsabs.harvard.edu/abs/1989ApJ...345..245C )
#define Ab_Av 1.337
#define Ar_Av 0.751
#define Ai_Av 0.479
#define Ak_Av 0.114

static double slab_model(float tau, float costheta)
{
    float correction_factor;

    if (tau == 0)
        correction_factor = 1.;
    else {
        if (costheta > 0.)
            correction_factor = (costheta / tau) * (exp(tau / costheta) - 1) / exp(tau / costheta);
        else
            correction_factor = 1 / exp(tau);

        if (correction_factor > 1)
            correction_factor = 1.;
    }

    return correction_factor;
}

void apply_dust(int n_photo_bands, galaxy_t gal, double* LumDust, int outputbin)
{
    double tauB, tauV, tauR, tauI, tauK;
    double Lum, Lum_corr, incl;

    // Dust correction using the "classical implementation" of Kauffmann et al. 1999

    // note the inclination angle
    incl = gal.Cos_Inc;

    for (int ii = 0; ii < n_photo_bands; ii++)
        LumDust[ii] = 0.0;

    // Correction for B magnitude
    if (gal.Lum[0][outputbin] > 0) {
        Lum = gal.Lum[0][outputbin];
        tauB = tauBstar * pow((Lum / Lbstar), beta);
        Lum_corr = Lum * slab_model(tauB, incl);

        LumDust[0] = (float)Lum_corr;
    }
    else
        tauB = 0.0;

    // Correction for V magnitude
    tauV = tauB / Ab_Av;
    if (gal.Lum[1][outputbin] > 0) {
        Lum = gal.Lum[1][outputbin];
        Lum_corr = Lum * slab_model(tauV, incl);

        LumDust[1] = (float)Lum_corr;
    }
    else
        tauV = 0.0;

    // Correction for R magnitude
    if (gal.Lum[2][outputbin] > 0) {
        tauR = tauV * Ar_Av;
        Lum = gal.Lum[2][outputbin];
        Lum_corr = Lum * slab_model(tauR, incl);

        LumDust[2] = (float)Lum_corr;
    }

    // Correction for I magnitude
    if (gal.Lum[3][outputbin] > 0) {
        tauI = tauV * Ai_Av;
        Lum = gal.Lum[3][outputbin];
        Lum_corr = Lum * slab_model(tauI, incl);

        LumDust[3] = (float)Lum_corr;
    }

    // Correction for K magnitude
    if (gal.Lum[4][outputbin] > 0) {
        tauK = tauV * Ak_Av;
        Lum = gal.Lum[4][outputbin];
        Lum_corr = Lum * slab_model(tauK, incl);

        LumDust[4] = (float)Lum_corr;
    }
}

#endif // CALC_MAGS
