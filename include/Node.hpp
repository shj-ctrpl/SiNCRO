#pragma once

#include "../external/jiarray/JIArray.h"
#include "../external/jiarray/JIVector.h"
#include "Constants.h"

#include <iostream>
#include <map>

using namespace dnegri::jiarray;

struct Material
{
    std::vector<double> xsAot;
    std::vector<double> xsDfs;
    std::vector<double> xsFis;
    std::vector<double> xsSct;
};

class Node
{
private:
public:
    /// @brief  Average flux of the mesh [z, y, x, g]
    zdouble4 flux;

    /// @brief Net current of the mesh [z, y, x, dir, g, face]
    zdouble5 jnet;

    /// @brief Surface flux of the mesh [z, y, x, dir, g, face]
    zdouble5 sFlux;

    /// @brief Eigenvalue of the system
    double keff;

    /// @brief Mesh size of the core [dir]
    zdouble1 _hMesh;

    /// @brief Albedo of boundary [dir, face]
    zdouble2 albedo;

    /// @brief Absorption cross section of the mesh [z, y, x, g]
    zdouble4 xsA;

    /// @brief Removal cross section of the mesh [z, y, x, g]
    zdouble4 xsR;

    /// @brief Diffusion coefficient of the mesh [z, y, x, g]
    zdouble4 xsD;

    /// @brief Diffusion coefficient divided by mesh size
    zdouble5 beta;

    /// @brief Mesh size normalized diffusion coefficient (D_tilde) of the mesh [z, y, x, g, dir]
    zdouble5 xsDt;

    /// @brief Nu-Fission cross section matrix with fission spectrum (χ) [z, y, x, from-g, to-g]
    zdouble5 xsF;

    /// @brief Scattering cross section matrix [z, y, x, from-g, to-g]
    zdouble5 xsS;

    std::map<int, Material> materials;
    std::map<int, int *> batches;

    int nz, ny, nx, ng, nd;
    double hx, hy, hz, vol;
    Node()
    {
    }

    Node(const int nz, const int ny, const int nx, const double hz, const double hy, const double hx, const int ng = 2, const int nd = 3)
    {
        this->nz = nz;
        this->ny = ny;
        this->nx = nx;
        this->ng = ng;
        this->hx = hx;
        this->hy = hy;
        this->hz = hz;
        this->nd = nd;

        xsA = zdouble4(nz, ny, nx, ng);
        xsR = zdouble4(nz, ny, nx, ng);
        xsD = zdouble4(nz, ny, nx, ng);
        beta = zdouble5(nz, ny, nx, ng, nd);
        xsDt = zdouble5(nz, ny, nx, ng, nd);
        xsF = zdouble5(nz, ny, nx, ng, ng);
        xsS = zdouble5(nz, ny, nx, ng, ng);

        flux = zdouble4(nz, ny, nx, ng);
        flux = 1.0;
        jnet = zdouble5(nz + 1, ny + 1, nx + 1, nd, ng);
        sFlux = zdouble5(nz + 1, ny + 1, nx + 1, nd, ng);

        keff = 1.0;

        _hMesh = zdouble1(3);
        _hMesh = {hz, hy, hx};
        vol = hx * hy * hz;
    }

    ~Node()
    {
        for (auto batch : batches)
        {
            delete[] batch.second;
        }
    }

    void AddMaterial(int matIdx, const std::vector<double> &xsAot, const std::vector<double> &xsDfs,
                     const std::vector<double> &xsFis, const std::vector<double> &xsSct)
    {
        materials.insert({matIdx, {xsAot, xsDfs, xsFis, xsSct}});
    }

    /// @brief Add a new batch, which represents z-directional fuel assembly
    /// @param batchIdx Index of the batch
    void AddBatch(int batchIdx, const std::initializer_list<std::pair<int, int>> matZmesh)
    {
        batches.insert(std::pair(batchIdx, new int[nz]));
        int z = 0;
        for (auto mat : matZmesh)
        {
            for (int i = 0; i < mat.second; i++)
            {
                batches.find(batchIdx)->second[z] = mat.first;
                z++;
            }
        }
    }

    void SetCore(zint2 core)
    {
        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    int matIdx = batches.find(core(x, y))->second[z];
                    Material &mat = materials.find(matIdx)->second;

                    for (int from_g = 0; from_g < ng; from_g++)
                    {
                        xsA(z, y, x, from_g) = mat.xsAot[from_g];
                        xsR(z, y, x, from_g) = xsA(z, y, x, from_g);
                        xsD(z, y, x, from_g) = mat.xsDfs[from_g];

                        for (int dir = 0; dir < nd; dir++)
                        {
                            beta(z, y, x, from_g, dir) = xsD(z, y, x, from_g) / _hMesh(dir);
                            xsDt(z, y, x, from_g, dir) = 0.5 * xsD(z, y, x, from_g) * _hMesh(dir);
                        }

                        for (int to_g = 0; to_g < ng; to_g++)
                        {  
                            xsF(z, y, x, from_g, to_g) = mat.xsFis[ng * to_g + from_g];
                            xsS(z, y, x, from_g, to_g) = mat.xsSct[ng * to_g + from_g];
                            xsR(z, y, x, from_g) += xsS(z, y, x, from_g, to_g);
                        }
                    }
                }
            }
        }
    }

    void SetAlbedo(double zLeft, double zRight, double yLeft, double yRight, double xLeft, double xRight)
    {
        albedo = zdouble2(3, 2);
        albedo = {zLeft, yLeft, xLeft, zRight, yRight, xRight};
    }

    /// @brief Macroscopic absorption cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param g group index
    inline const double XsA(const int z, const int y, const int x, const int g) const
    {
        return xsA(z, y, x, g);
    };

    /// @brief Macroscopic removal cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param g group index
    inline const double XsR(const int z, const int y, const int x, const int g) const
    {
        return xsR(z, y, x, g);
    };

    /// @brief Macroscopic diffusion coefficient for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param g group index
    inline const double XsD(const int z, const int y, const int x, const int g) const
    {
        return xsD(z, y, x, g);
    };

    /// @brief Macroscopic nu-fission cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param from_g from-group index
    /// @param to_g to-group index
    inline const double XsF(const int z, const int y, const int x, const int from_g, const int to_g) const
    {
        return xsF(z, y, x, from_g, to_g);
    };

    /// @brief Macroscopic scattering cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param from_g from-group index
    /// @param to_g to-group index
    inline const double XsS(const int z, const int y, const int x, const int from_g, const int to_g) const
    {
        return xsS(z, y, x, from_g, to_g);
    };

    /// @brief Mesh size of the core
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param dir
    inline const double HMesh(const int dir) const
    {
        return _hMesh(dir);
    }

    inline const double Albedo(const int dir, const int face)
    {
        return albedo(dir, face);
    }

    inline const double Beta(const int z, const int y, const int x, const int g, const int dir) const
    {
        return beta(z, y, x, g, dir);
    }

    void PrintMaterial()
    {
        std::cout << "Total cross section of group 0: " << std::endl;
        for (int z = 0; z < nz; z++)
        {
            std::cout << "zdepth = " << z << std::endl;
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    std::cout << std::scientific << xsA(z, y, x, 0) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        std::cout << "Total cross section of group 1: " << std::endl;
        for (int z = 0; z < nz; z++)
        {
            std::cout << "zdepth = " << z << std::endl;
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    std::cout << std::scientific << xsA(z, y, x, 1) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    zdouble4 &Flux()
    {
        return flux;
    }

    zdouble5 &SFlux()
    {
        return sFlux;
    }

    zdouble5 &Jnet()
    {
        return jnet;
    }

    double &Keff()
    {
        return keff;
    }
};