/*
 * IDL CALL_EXTERNAL interface for slitdec_ox (REDUCE).
 *
 * Entry point: slit_func_2d, a drop-in replacement for slit_func_2d in
 * slit_func_2d_xi_zeta_bd_reduce.c and in ../slitdec_idl.c: mkslitf.pro calls it
 * unchanged (pass LIBRARY_NAME= the library built by build.sh, or link this file
 * into reduce.so instead of slitdec_idl.c).
 *
 * Arguments (25 required, as mkslitf.pro passes them; IDL long = C int):
 *   0  ncols       long   scalar          swath width in pixels
 *   1  nrows       long   scalar          swath height in pixels
 *   2  osample     long   scalar          y oversampling of the slit function
 *   3  im          double [ncols,nrows]   swath image
 *   4  ycen        double [ncols]         fractional trace position (0 <= v < 1)
 *   5  ycen_offset long   [ncols]         integer part of the trace position
 *   6  y_lower_lim long   scalar          rows below the trace row in the swath
 *   7  PSF_curve   double [3,ncols]       slit curvature: shift = c1*dy + c2*dy^2
 *   8  delta_x     long   scalar          unused (computed internally)
 *   9  lambda_sP   double scalar          unused: the spectrum penalty is mu_x
 *  10  lambda_sL   double scalar          slit-function smoothing (legacy, used
 *                                         unless mu1_y or mu2_y is set)
 *  11  mask        byte   [ncols,nrows]   in/out, 0 = bad (updated by kappa clip)
 *  12  model       double [ncols,nrows]   OUTPUT model swath
 *  13  sP          double [ncols]         OUTPUT detector-pixel spectrum
 *  14  sP_old      double [ncols]         unused
 *  15  sL          double [ny]            in/out slit function, ny = (nrows+1)*osample+1
 *  16  unc         double [ncols]         OUTPUT spectrum uncertainty
 *  17-23                                  unused scratch arrays
 *  24  sL_preset   long   scalar          non-zero: sL is given and held fixed
 * Optional (each needs the preceding ones):
 *  25  maxiter     long   scalar          default 100
 *  26  osample_x   long   scalar          fine spectrum bins per column, default 8
 *  27  mu1_x       double scalar          spectrum penalties on the 1st, 2nd and
 *  28  mu2_x       double scalar            3rd derivative (physical units, see
 *  29  mu3_x       double scalar            NOTES.md); defaults below
 *  30  mu1_y       double scalar          slit-function penalties on the 1st and
 *  31  mu2_y       double scalar            2nd derivative; both 0: lambda_sL
 *  32  sP_fine     double [ncols*osample_x] OUTPUT fine spectrum
 *  33  info        double [11]            OUTPUT: success, sP change, status,
 *                                         iterations, delta_x, dropped zeta,
 *                                         sP band, damping x, resolution x [px],
 *                                         damping y, resolution y [px] (damping
 *                                         < 1 can ring; legacy lambda_sL gives
 *                                         -mu1_y equivalent as damping y)
 *
 * When an optional argument is not passed, the environment variable
 * SLITDEC_OX_<NAME> (SLITDEC_OX_MAXITER, _OSAMPLE_X, _MU1_X, _MU2_X, _MU3_X, _MU1_Y,
 * _MU2_Y, and for testing _T_OFFSET) overrides the default, so mkslitf.pro can be
 * steered without editing it:
 *   IDL> setenv, 'SLITDEC_OX_OSAMPLE_X=16'
 * osample_x = 1 with all mu_x = 0 is the historic slitdec model (apart from the
 * curvature being evaluated at the y-subpixel centre, see NOTES.md).
 *
 * Memory layout: IDL [ncols,nrows] is column-major, i.e. C [nrows][ncols], so the
 * flat index y*ncols + x is the same on both sides and no copy is needed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "slitdec_ox.h"

extern double slitdec_ox_t_offset;

/* Default spectrum penalty (see NOTES.md) */
#define SLITDEC_OX_MU1_X 0.0
#define SLITDEC_OX_MU2_X 0.0
#define SLITDEC_OX_MU3_X 1.64e-4

#ifdef BUILDING_REDUCE_WIN_DLL
#define REDUCE_DLL __declspec(dllexport)
#else
#define REDUCE_DLL
#endif

static double env_or(const char *name, double dflt)
{
    const char *s = getenv(name);
    char *end;
    double v;
    if (s == NULL || *s == '\0')
        return dflt;
    v = strtod(s, &end);
    return end == s ? dflt : v;
}

int REDUCE_DLL slit_func_2d(int argc, void *argv[])
{
    int ncols, nrows, osample, y_lower_lim, sL_preset, maxiter, osample_x;
    double lambda_sL, mu_x[3], mu_y[2];
    double *im, *PSF_curve, *sP, *sL, *unc, *model, *ycen_frac, *sP_fine = NULL, *info_out = NULL;
    int *ycen_offset;
    unsigned char *mask;
    int x, y, ny, ret;
    double info[11] = {0};

    if (argc < 25)
    {
        printf("slit_func_2d (slitdec_ox): requires at least 25 arguments, got %d\n", argc);
        printf("  ncols, nrows, osample, im, ycen_frac, ycen_offset, y_lower_lim, PSF_curve,\n");
        printf("  delta_x, lambda_sP, lambda_sL, mask, model, sP, sP_old, sL, unc,\n");
        printf("  7 x scratch, sL_preset [, maxiter, osample_x, mu1_x, mu2_x, mu3_x,\n");
        printf("  mu1_y, mu2_y, sP_fine, info]\n");
        return -1;
    }

    ncols       = *(int *)argv[0];
    nrows       = *(int *)argv[1];
    osample     = *(int *)argv[2];
    im          = (double *)argv[3];
    ycen_frac   = (double *)argv[4];
    ycen_offset = (int *)argv[5];
    y_lower_lim = *(int *)argv[6];
    PSF_curve   = (double *)argv[7];
    lambda_sL   = *(double *)argv[10];
    mask        = (unsigned char *)argv[11];
    model       = (double *)argv[12];
    sP          = (double *)argv[13];
    sL          = (double *)argv[15];
    unc         = (double *)argv[16];
    sL_preset   = *(int *)argv[24];

    maxiter   = argc > 25 ? *(int *)argv[25]    : (int)env_or("SLITDEC_OX_MAXITER", 100);
    osample_x = argc > 26 ? *(int *)argv[26]    : (int)env_or("SLITDEC_OX_OSAMPLE_X", 8);
    mu_x[0]   = argc > 27 ? *(double *)argv[27] : env_or("SLITDEC_OX_MU1_X", SLITDEC_OX_MU1_X);
    mu_x[1]   = argc > 28 ? *(double *)argv[28] : env_or("SLITDEC_OX_MU2_X", SLITDEC_OX_MU2_X);
    mu_x[2]   = argc > 29 ? *(double *)argv[29] : env_or("SLITDEC_OX_MU3_X", SLITDEC_OX_MU3_X);
    mu_y[0]   = argc > 30 ? *(double *)argv[30] : env_or("SLITDEC_OX_MU1_Y", 0.0);
    mu_y[1]   = argc > 31 ? *(double *)argv[31] : env_or("SLITDEC_OX_MU2_Y", 0.0);
    /* Testing: where along a y-subpixel the tilt is evaluated (0.5 = centre,
       0 = slitdec.c, see NOTES.md); only via the environment.             */
    slitdec_ox_t_offset = env_or("SLITDEC_OX_T_OFFSET", 0.5);
    if (argc > 32)
        sP_fine = (double *)argv[32];
    if (argc > 33)
        info_out = (double *)argv[33];
    info[4] = -11; /* ask for the extended info */

    if (lambda_sL <= 0.0)
        lambda_sL = 1.0;
    if (osample_x < 1)
        osample_x = 1;

    ny = osample * (nrows + 1) + 1;

    /* slitdec_ox takes the full trace position and splits it itself;
       it overwrites this copy with the fractional part.                    */
    double *ycen = malloc(ncols * sizeof(double));
    double *slitcurve = calloc((size_t)ncols * 6, sizeof(double));
    double *slitdeltas = calloc(ny, sizeof(double));
    double *pix_unc = calloc((size_t)nrows * ncols, sizeof(double)); /* unused: weight_mode 0 */
    if (!ycen || !slitcurve || !slitdeltas || !pix_unc)
    {
        free(ycen); free(slitcurve); free(slitdeltas); free(pix_unc);
        return -1;
    }
    for (x = 0; x < ncols; x++)
    {
        ycen[x] = ycen_frac[x] + ycen_offset[x];
        /* IDL PSF_curve[3,ncols]: element [k,x] at k + 3*x. slitdec uses
           slitcurve[x*6 + k]; the constant term k = 0 is not used.         */
        slitcurve[x * 6 + 1] = PSF_curve[1 + 3 * x];
        slitcurve[x * 6 + 2] = PSF_curve[2 + 3 * x];
    }

    /* Initial spectrum from the image, as the original slit_func_2d does */
    for (x = 0; x < ncols; x++)
    {
        double norm = 0.0;
        sP[x] = 0.0;
        for (y = 0; y < nrows; y++)
        {
            norm += mask[y * ncols + x];
            sP[x] += fabs(im[y * ncols + x] * mask[y * ncols + x]);
        }
        norm /= nrows;
        sP[x] = norm > 0.0 ? sP[x] / norm : -1000.0;
    }

    ret = slitdec_ox2(ncols, nrows, im, pix_unc, mask, ycen, slitcurve, slitdeltas,
                          osample, osample_x, mu_x, 3, mu_y, 2, lambda_sL, maxiter,
                          10.0,      /* kappa: 10-sigma outlier rejection */
                          sL_preset, /* use_preset */
                          0,         /* weight_mode: none */
                          -1.0,      /* noise_floor: unused */
                          sP, sP_fine, sL, model, unc, info, y_lower_lim);

    if (info_out != NULL)
        for (x = 0; x < 11; x++)
            info_out[x] = info[x];

    free(ycen);
    free(slitcurve);
    free(slitdeltas);
    free(pix_unc);
    return ret;
}
