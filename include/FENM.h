#pragma once

#define EPS 0.0000000001
#define LEFT 0
#define RIGHT 1
#define CENTER 2

#define ZDIR 0
#define YDIR 1
#define XDIR 2

#include "../external/jiarray/JIArray.h"
#include "../external/jiarray/JIVector.h"
#include "Node.hpp"
#include <cmath>
#include <algorithm>

using namespace dnegri::jiarray;

class FENM
{
private:
    /// @brief The number of energy groups.
    const int _ng = 2;

    /// @brief The number of directions.
    const int _ndir = 3;

    /// @brief The number of nodes in each directions.
    int _nz, _ny, _nx;

    /// @brief The size of nodes in xy direction.
    double _hz, _hy, _hx;

    double revk;

    /// @brief Transverse leakage coefficients. (z, y, x, dir, g)
    zdouble5 tlcoeff0, tlcoeff1, tlcoeff2;

    double mu00, mu11, mu22, mu33, mu44;

    zdouble6 mu_m;

    /// @brief Pre-calculated 'constants' from integration of P''(xi)*P(xi). Scalar value.
    double k20, k31, k40, k42;

    /// @brief Pre-calculated 'variables' from integration of P''(xi)*P(xi). (z, y, x, dir, g)
    zdouble5 k51, k53, k60, k62, k64;

    zdouble5 T5, T6;

    /// @brief Diagonal matrix of diffusion coefficients, D (z, y, x, dir, g)
    zdouble5 diagD;

    /// @brief Diagonal matrix of the inverse of diffusion coefficients, D^(-1). (z, y, x, dir, g)
    zdouble5 invdiagD;

    /// @brief Full matrix of the material cross-section data. (z, y, x, from-g, to-g)
    zdouble5 matM, invmatM, matS, matF;

    /// @brief Net current coefficients. (z, y, x, dir, g)
    zdouble5 c0, c1, c2, c3, c4, c5, c6;

    zdouble4 avgFlux;



    JIArray<JIArray<double, 2>, 4> mat;

    zdouble6 jnet;
    zdouble6 surfFlux;
    
    Node _node;

    void Initialize();

    /// @brief Update constants.
    void CalcConst(const int z, const int y, const int x);

    void CalcTLCoeffs(const int z, const int y, const int x);

    void CalcMatrix(const int z, const int y, const int x);

    void CalcEvenCoeffs(const int z, const int y, const int x);

    void CalcOddCoeffsOneNode(const int z, const int y, const int x, const int d, const int face);

    void CalcOddCoeffsTwoNode(const int z, const int y, const int x);

    void CalcCurrentAndFlux(const int z, const int y, const int x);
public:
    /// @brief Basic constructor
    FENM(const int _nz, const int _ny, const int _nx, const double _hz, const double _hy, const double _hx, Node &node);

    ~FENM();
};
