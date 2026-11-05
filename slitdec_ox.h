/* slitdec with oversampling of the spectrum along x (see slitdec_ox.c). */

typedef struct
{
    int x;    /* Contributing fine spectrum bin (x * osample_x + p) */
    int iy;   /* Contributing slit-function subpixel                */
    double w; /* Contribution weight                                */
} zeta_ref_ox;

typedef struct
{
    int min_iy, max_iy;
    int min_x, max_x;   /* in fine-bin units */
} zeta_rng_ox;

int slitdec_ox(     int ncols,
                    int nrows,
                    double *im,
                    double *pix_unc,
                    unsigned char *mask,
                    double *ycen,
                    double *slitcurve,
                    double *slitdeltas,
                    int osample,
                    int osample_x,
                    double mu_x,
                    int sp_order,
                    double lambda_sL,
                    int maxiter,
                    double kappa,
                    int use_preset,
                    int weight_mode,
                    double noise_floor,
                    double *sP,
                    double *sP_fine,
                    double *sL,
                    double *model,
                    double *unc,
                    double *info);

#define SLITDEC_OX_KMAX 15 /* highest derivative order of the penalties */

/* General entry point. Penalties sum_k mu[k-1] * F * int (d^k s)^2 in physical
   units for the spectrum (mu_x, nk_x orders) and the slit function (mu_y, nk_y
   orders; none: slitdec's lambda_sL). y_lower_lim: rows below the trace row in
   the swath (IDL REDUCE), < 0 means nrows / 2. See slitdec_ox.c. */
int slitdec_ox2(int ncols, int nrows, double *im, double *pix_unc,
                unsigned char *mask, double *ycen, double *slitcurve,
                double *slitdeltas, int osample, int osample_x,
                const double *mu_x, int nk_x, const double *mu_y, int nk_y,
                double lambda_sL, int maxiter, double kappa,
                int use_preset, int weight_mode, double noise_floor,
                double *sP, double *sP_fine, double *sL, double *model,
                double *unc, double *info, int y_lower_lim);
