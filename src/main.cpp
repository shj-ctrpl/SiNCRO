#include "Node.hpp"
#include "../include/FENM.h"
#include "../include/CMFD.h"

using namespace dnegri::jiarray;

int main(void)
{
    Node node;

    int nx = 6, ny = 6, nz = 14, ng = 2, ndir = 3;
    double hx = 5.0, hy = 5.0, hz = 5.0;

    node = Node(nz, ny, nx, hz, hy, hx, ng, ndir);

    node.AddMaterial(1,
                     {0.010, 0.100},
                     {1.400, 0.200},
                     {0.005, 0.140,
                      0.000, 0.000},
                     {0.000, 0.000,
                      0.060, 0.000});

    node.AddMaterial(2,
                     {0.050, 0.500},
                     {1.100, 0.120},
                     {0.000, 0.000,
                      0.000, 0.000},
                     {0.000, 0.000,
                      0.010, 0.000});

    node.AddMaterial(3,
                     {0.010, 0.100},
                     {1.400, 0.200},
                     {0.009, 0.240,
                      0.000, 0.000},
                     {0.000, 0.000,
                      0.060, 0.000});

    node.AddMaterial(4,
                     {0.0006, 0.019},
                     {1.230, 0.160},
                     {0.000, 0.000,
                      0.000, 0.000},
                     {0.000, 0.000,
                      0.050, 0.000});

    node.AddBatch(1, {std::pair<int, int>(1, 6), std::pair<int, int>(2, 1), std::pair<int, int>(3, 3), std::pair<int, int>(4, 4)});

    zint2 core;
    core = zint2(ny, nx);
    core = {1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1};

    node.SetCore(core);

    node.SetAlbedo(REFLECTIVE, VACUUM, REFLECTIVE, REFLECTIVE, REFLECTIVE, REFLECTIVE);
    // node.PrintMaterial();

    CMFD coarse(node, 6, 6, 14, false);
    // coarse.ShowXS(0, true);
    CMFD fine(node, 6, 6, 14, true);

    FENM fenm(node);

    // FENM fenm(node);

    bool yesfine = false;

    for (int i = 0; i < 20; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            coarse.UpdateCrossSection();
            coarse.UpdateSource();
            coarse.UpdateFlux(5);
            coarse.UpdateKeff(false);
            coarse.Redistribute();
        }
        coarse.showKeff();
        coarse.UpdateJnet();
        fenm.Drive();
        coarse.UpdateDhat();
    }

    return 0;
}