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
    zdouble5 matM, invM;

    /// @brief Net current coefficients. (z, y, x, dir, g)
    zdouble5 c1, c2, c3, c4, c5, c6;

    /// @brief Update constants.
    void CalcConst();
    void CalcTLCoeffs0();
    void CalcTLCoeffs12();
    void CalcMatrix();
    void CalcEvenCoeffs();
    void CalcOddCoeffsOneNode(int z, int y, int x, int d, int face);
    void CalcOddCoeffsTwoNode(int z, int y, int x, int d);
    void CalcCurrentAndFlux();

public:
    /// @brief Basic constructor
    FENM(Node &node);

    ~FENM();

    void Drive();

    void printTLCoeffs(int n, int d, int g)
    {
        for (int z = 0; z < nz; z++)
        {
            std::cout << "Z level: " << z << std::endl;
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    if (n == 0)
                    {
                        std::cout << tlcoeff0(z, y, x, d, g) << " ";
                    }
                    else if (n == 1)
                    {
                        std::cout << tlcoeff1(z, y, x, d, g) << " ";
                    }
                    else
                    {
                        std::cout << tlcoeff2(z, y, x, d, g) << " ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    void printJnetcoeffs(int n, int d, int g)
    {
        for (int z = 0; z < nz; z++)
        {
            std::cout << "Z level: " << z << std::endl;
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    if (n == 1)
                    {
                        std::cout << c1(z, y, x, d, g) << " ";
                    }
                    else if (n == 2)
                    {
                        std::cout << c2(z, y, x, d, g) << " ";
                    }
                    else if (n == 3)
                    {
                        std::cout << c3(z, y, x, d, g) << " ";
                    }
                    else if (n == 4)
                    {
                        std::cout << c4(z, y, x, d, g) << " ";
                    }
                    else if (n == 5)
                    {
                        std::cout << c5(z, y, x, d, g) << " ";
                    }
                    else if (n == 6)
                    {
                        std::cout << c6(z, y, x, d, g) << " ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
};
