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

/* As slitdec_ox, with the row of the trace in the swath given explicitly
   (IDL REDUCE y_lower_lim; < 0 means nrows / 2). */
int slitdec_ox_ylim(int ncols, int nrows, double *im, double *pix_unc,
                    unsigned char *mask, double *ycen, double *slitcurve,
                    double *slitdeltas, int osample, int osample_x, double mu_x,
                    int sp_order, double lambda_sL, int maxiter, double kappa,
                    int use_preset, int weight_mode, double noise_floor,
                    double *sP, double *sP_fine, double *sL, double *model,
                    double *unc, double *info, int y_lower_lim);
