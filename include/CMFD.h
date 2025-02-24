#pragma once

#include "../external/jiarray/JIArray.h"

#include "Constants.h"
#include "Node.hpp"
#include <iostream>
#include <math.h>

using namespace dnegri::jiarray;

class CMFD
{
private:
    Node &node;

    bool FDMmode = false;
    bool Samesizemode = false;

    /// @brief The nubmer of the corase mesh of each direction
    int nx, ny, nz;

    /// @brief The number of fine mesh of each direction in a single coarse mesh
    int mx, my, mz;

    /// @brief Mesh size of each direction
    double hx, hy, hz;

    /// @brief Reverse of hx, hy and hz
    double rhx, rhy, rhz;

    double vol;

    int ng;

    int nd;

    double sumFluxOld, sumFluxNew;

    ///========================================
    /// Fine mesh
    ///========================================

    /// @brief Effective multiplication factor
    double &keff;

    /// @brief Scalar flux of the fine mesh [z, y, x, g]
    zdouble4 &fine_Flux;

    /// @brief Net current of the fine mesh [z, y, x, d, g]
    zdouble5 &fine_Jnet;

    ///========================================
    /// Coarse mesh
    ///========================================

    zdouble4 coarse_Old;

    /// @brief Scalar flux of the coarse mesh [z, y, x, g]
    zdouble4 coarse_Flux;

    /// @brief Net current of the corase mesh [z, y, x, d, g]
    zdouble5 coarse_Jnet;

    /// @brief Removal cross section of the corase mesh [z, y, x, g]
    zdouble4 coarse_R;

    zdouble4 coarse_Src;

    /// @brief Fission cross section matrix of the corase mesh [z, y, x, from_g, to_g]
    zdouble5 coarse_F;

    /// @brief Scattering cross section matrix of the corase mesh [z, y, x, from_g, to_g]
    zdouble5 coarse_S;

    /// @brief Diffusion coefficient of the corase mesh
    zdouble4 coarse_D;

    /// @brief Mesh-size normalized diffusion coefficient of the corase mesh [z, y, x, d, g]
    zdouble5 coarse_Dt;

    /// @brief Correction factor of the corase mesh [z, y, x, d, g]
    zdouble5 coarse_Dh;

public:
    CMFD(Node &node, int nx, int ny, int nz, bool FDMmode);

    /// @brief Calculate D-tilde, homogenized cross-sections and homogenized diffusion coefficients
    void UpdateSource();
    void UpdateCrossSection();
    void UpdateFlux(int numIter);
    void InitFDM();
    void UpdateKeff(bool print);
    void UpdateJnet();
    void UpdateDhat();
    void UpdateDhatNodal();
    void Redistribute();
    void printDhat(int g, int dir)
    {
        for (int z = 0; z < nz; z++)
        {
            std::cout << "z = " << z << std::endl;
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    std::cout << coarse_Dh(z, y, x, dir, g) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    inline double harmonicMean(const double a, const double b)
    {
        return 2 * (a * b) / (a + b);
    }
};