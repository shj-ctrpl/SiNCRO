#define EPS 1.0E-6
#define LEFT 0
#define RIGHT 1
#define CENTER 2

#define ZDIR 0
#define YDIR 1
#define XDIR 2

#define FLUXZERO 1.0E+10
#define REFLECTIVE 1.0E-10
#define VACUUM 0.5

/// @brief Pre-calculated 'constants' from integration of P(xi)*P(xi). Scalar value.
#define mu00 2.000000000000

/// @brief Pre-calculated 'constants' from integration of P(xi)*P(xi). Scalar value.
#define mu11 0.666666666667

#define rmu11 1.5

/// @brief Pre-calculated 'constants' from integration of P(xi)*P(xi). Scalar value.
#define mu22 0.400000000000

/// @brief Pre-calculated 'constants' from integration of P(xi)*P(xi). Scalar value.
#define mu33 0.285714285714

/// @brief Pre-calculated 'constants' from integration of P(xi)*P(xi). Scalar value.
#define mu44 0.222222222222

/// @brief Pre-calculated 'constants' from integration of P''(xi)*P(xi). Scalar value.
#define k20 6.00

#define rk20 0.166666666667

/// @brief Pre-calculated 'constants' from integration of P''(xi)*P(xi). Scalar value.
#define k31 10.0

/// @brief Pre-calculated 'constants' from integration of P''(xi)*P(xi). Scalar value.
#define k40 20.0

/// @brief Pre-calculated 'constants' from integration of P''(xi)*P(xi). Scalar value.
#define k42 14.0

/// @brief Loop over 3D
#define LOOP3D(nz, ny, nx)             \
    for (int z = 0; z < (nz); z++)     \
        for (int y = 0; y < (ny); y++) \
            for (int x = 0; x < (nx); x++)

#define LOOP3DInside(nz, ny, nx)           \
    for (int z = 0; z < (nz - 1); z++)     \
        for (int y = 0; y < (ny - 1); y++) \
            for (int x = 0; x < (nx - 1); x++)

