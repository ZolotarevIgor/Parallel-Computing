#include <math.h>
static inline double k(const double u) {
    return 1.0 + u*u;
}
static inline double ki(const double u1, const double u2) {
    const double k1 = k(u1);
    const double k2 = k(u2);
    return 2*(k1*k2)/(k1+k2);
}

static inline double q(const double u) {
    return 0.5*M_PI*M_PI*(1-u*u);
}
static inline double u(const double x1, const double x2) {
    return cos(0.5*M_PI*(x1+x2));
}

static inline double f(const double u) {
    return M_PI*M_PI*u*u*u;
}