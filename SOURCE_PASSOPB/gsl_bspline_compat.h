#pragma once
// GSL 2.x moved derivative workspace storage into gsl_bspline_workspace.
// Keep the existing numerical calls unchanged on both GSL 1.x and 2.x.
#include <gsl/gsl_version.h>
#if GSL_MAJOR_VERSION >= 2
typedef gsl_bspline_workspace channel_deriv_workspace;
inline channel_deriv_workspace *channel_deriv_alloc(size_t) { return 0; }
inline void channel_deriv_free(channel_deriv_workspace *) {}
inline int channel_deriv_eval(double x, size_t n, gsl_matrix *d,
                         gsl_bspline_workspace *w, channel_deriv_workspace *) {
    return gsl_bspline_deriv_eval(x, n, d, w);
}
inline int channel_deriv_eval_nonzero(double x, size_t n, gsl_matrix *d,
                                 size_t *first, size_t *last,
                                 gsl_bspline_workspace *w, channel_deriv_workspace *) {
    return gsl_bspline_deriv_eval_nonzero(x, n, d, first, last, w);
}
#else
typedef gsl_bspline_deriv_workspace channel_deriv_workspace;
#define channel_deriv_alloc gsl_bspline_deriv_alloc
#define channel_deriv_free gsl_bspline_deriv_free
#define channel_deriv_eval gsl_bspline_deriv_eval
#define channel_deriv_eval_nonzero gsl_bspline_deriv_eval_nonzero
#endif
