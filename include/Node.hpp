#pragma once

#define LEFT 1
#define RIGHT 2

#include "../external/jiarray/JIArray.h"
#include "../external/jiarray/JIVector.h"
#include <iostream>
#include <map>

using namespace dnegri::jiarray;

struct Material
{
    std::vector<double> xsTot;
    std::vector<double> xsDfs;
    std::vector<double> xsFis;
    std::vector<double> xsSct;
};

class Node
{
private:
    int _nz, _ny, _nx, _ng, _ndir;
    double _hx, _hy, _hz;

    /// @brief Mesh size of the core [z, y, x, direction]
    zdouble4 _hMesh;

    /// @brief Albedo of boundary [dir, face]
    zdouble2 albedo;

    std::map<int, Material> materials;
    std::map<int, int *> batches;

    /// @brief Total cross section of the mesh [z, y, x, g]
    zdouble4 _xsTot;

    /// @brief Diffusion coefficient of the mesh [z, y, x, g]
    zdouble4 _xsDfs;

    /// @brief Mesh size normalized diffusion coefficient (D_tilde) of the mesh [z, y, x, g, dir]
    zdouble5 _xsDtd;

    /// @brief Nu-Fission cross section matrix with fission spectrum (χ) [z, y, x, from-g, to-g]
    zdouble5 _xsFis;

    /// @brief Scattering cross section matrix [z, y, x, from-g, to-g]
    zdouble5 _xsSct;

public:
    Node()
    {

    }
    Node(const int nz, const int ny, const int nx, const double hz, const double hy, const double hx, const int ng = 2, const int ndir = 3)
    {
        _nz = nz;
        _ny = ny;
        _nx = nx;
        _ng = ng;
        _hx = hx;
        _hy = hy;
        _hz = hz;
        _ndir = ndir;

        _xsTot = zdouble4(_nz, _ny, _nx, _ng);
        _xsDfs = zdouble4(_nz, _ny, _nx, _ng);
        _xsDtd = zdouble5(_nz, _ny, _nx, _ng, _ndir);
        _xsFis = zdouble5(_nz, _ny, _nx, _ng, _ng);
        _xsSct = zdouble5(_nz, _ny, _nx, _ng, _ng);
        _hMesh = zdouble4(nz, ny, nx, ndir);
    }
    ~Node(){
        for(auto batch: batches){
            delete[] batch.second;
        }
    }
    void AddMaterial(int matIdx, const std::vector<double>& xsTot, const std::vector<double>& xsDfs,
                     const std::vector<double>& xsFis, const std::vector<double>& xsSct)
    {
        materials.insert({matIdx, {xsTot, xsDfs, xsFis, xsSct}});
    }

    /// @brief Add a new batch, which represents z-directional fuel assembly
    /// @param batchIdx Index of the batch
    void AddBatch(int batchIdx, const std::initializer_list<std::pair<int, int>> matZmesh)
    {
        batches.insert(std::pair(batchIdx, new int[_nz]));
        int z = 0;
        for (auto mat : matZmesh)
        {
            for (int i = 0; i < mat.second; i++){
                batches.find(batchIdx)->second[z] = mat.first;
                z++;
            }
        }
    }

    void SetCore(zint2 core)
    {
        for (int z = 0; z < _nz; z++)
        {
            for (int y = 0; y < _ny; y++)
            {
                for (int x = 0; x < _nx; x++)
                {
                    int matIdx = batches.find(core(x, y))->second[z - 1];
                    Material& mat = materials.find(matIdx)->second;

                    for (int from_g = 0; from_g < _ng; from_g++)
                    {
                        _xsTot(z, y, x, from_g) = mat.xsTot[from_g];
                        _xsDfs(z, y, x, from_g) = mat.xsDfs[from_g];

                        for (int dir = 0; dir < _ndir; dir++)
                        {
                            _xsDtd(z, y, x, from_g, dir) = 0.5 * _xsDfs(z, y, x, from_g) * _hMesh(z, y, x, dir);
                        }

                        for (int to_g = 0; to_g < _ng; to_g++)
                        {
                            _xsFis(z, y, x, from_g, to_g) = mat.xsFis[_ng * to_g + from_g];
                            _xsSct(z, y, x, from_g, to_g) = mat.xsSct[_ng * to_g + from_g];
                        }
                    }
                }
            }
        }
    }

    /// @brief Macroscopic total cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param g group index
    inline const double XsTot(const int z, const int y, const int x, const int g) const
    {
        return _xsTot(z, y, x, g);
    };

    /// @brief Macroscopic diffusion coefficient for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param g group index
    inline const double XsDfs(const int z, const int y, const int x, const int g) const
    {
        return _xsDfs(z, y, x, g);
    };
    /// @brief Macroscopic nu-fission cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param from_g from-group index
    /// @param to_g to-group index
    inline const double XsFis(const int z, const int y, const int x, const int from_g, const int to_g) const
    {
        return _xsFis(z, y, x, from_g, to_g);
    };

    /// @brief Macroscopic scattering cross section for the node
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param from_g from-group index
    /// @param to_g to-group index
    inline const double XsSct(const int z, const int y, const int x, const int from_g, const int to_g) const
    {
        return _xsSct(z, y, x, from_g, to_g);
    };

    /// @brief Mesh size of the core
    /// @param z z directional index
    /// @param y y directional index
    /// @param x x directional index
    /// @param dir
    inline const double HMesh(const int z, const int y, const int x, const int dir) const
    {
        return _hMesh(z, y, x, dir);
    }

    inline const double Albedo(const int dir, const int face)
    {
        return albedo(dir, face);
    }

    void PrintMaterial()
    {
        std::cout << "Total cross section of group 0: " << std::endl;
        for (int z = 0; z < _nz; z++)
        {
            std::cout << "zdepth = " << z << std::endl;
            for (int y = 0; y < _ny; y++)
            {
                for (int x = 0; x < _nx; x++)
                {
                    std::cout << std::scientific << _xsTot(z, y, x, 0) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        std::cout << "Total cross section of group 1: " << std::endl;
        for (int z = 0; z < _nz; z++)
        {
            std::cout << "zdepth = " << z << std::endl;
            for (int y = 0; y < _ny; y++)
            {
                for (int x = 0; x < _nx; x++)
                {
                    std::cout << std::scientific << _xsTot(z, y, x, 1) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
};