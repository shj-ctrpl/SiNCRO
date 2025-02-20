#pragma once

#include "../external/jiarray/JIArray.h"
#include "../external/jiarray/JIVector.h"

#include "Constants.h"
#include "Node.hpp"
#include <cmath>
#include <algorithm>

using namespace dnegri::jiarray;

class FENM
{
private:
    Node &_node;
    zdouble4 &flux;
    zdouble5 &jnet;
    zdouble5 &sFlux;
    double &keff;

    /// @brief The number of energy groups.
    int ng;

    /// @brief The number of directions.
    int nd;

    /// @brief The number of nodes in each directions.
    int nz, ny, nx;

    double revk;

    /// @brief Transverse leakage coefficients. (z, y, x, dir, g)
    zdouble5 tlcoeff0, tlcoeff1, tlcoeff2;

    zdouble6 mu_m;

    /// @brief Pre-calculated 'variables' from integration of P''(xi)*P(xi). (z, y, x, dir, g)
    zdouble5 k51, k53, k60, k62, k64;

    zdouble5 T5, T6;

    /// @brief Diagonal matrix of diffusion coefficients, D (z, y, x, dir, g)
    zdouble5 diagD;

    /// @brief Diagonal matrix of the inverse of diffusion coefficients, D^(-1). (z, y, x, dir, g)
    zdouble5 invD;

    /// @brief Full matrix of the material cross-section data. (z, y, x, from-g, to-g)
    zdouble5 matM, invM, matS, matF;

    /// @brief Net current coefficients. (z, y, x, dir, g)
    zdouble5 c0, c1, c2, c3, c4, c5, c6;

    /// @brief Update constants.
    void CalcConst(const int z, const int y, const int x);
    
    void CalcTLCoeffs0(const int z, const int y, const int x);

    void CalcTLCoeffs12(const int z, const int y, const int x);

    void CalcMatrix(const int z, const int y, const int x);

    void CalcEvenCoeffs(const int z, const int y, const int x);

    void CalcOddCoeffsOneNode(const int z, const int y, const int x, const int d, const int face);

    void CalcOddCoeffsTwoNode(const int z, const int y, const int x);

    void CalcCurrentAndFlux(const int z, const int y, const int x);

public:
    /// @brief Basic constructor
    FENM(Node &node);

    ~FENM();

    void Drive();
};
