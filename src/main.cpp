// #include "FENM.h"
#include "Node.hpp"

using namespace dnegri::jiarray;

int main(void)
{
    Node node;

    int nx = 4, ny = 4, nz = 4, ng = 2;
    double hx = 1.0, hy = 1.0, hz = 2.0;

    node = Node(nz, ny, nx, hz, hy, hx, 2, 3);

    node.AddMaterial(1,
                     {0.1, 0.2},
                     {0.3, 0.4},
                     {0.1, 0.2,
                      0.3, 0.4},
                     {0.1, 0.2,
                      0.3, 0.4});

    node.AddMaterial(2,
                     {1.1, 1.2},
                     {1.3, 1.4},
                     {1.1, 1.2,
                      1.3, 1.4},
                     {1.1, 1.2,
                      1.3, 1.4});

    node.AddBatch(1, {std::pair<int, int>(1, 2), std::pair<int, int> (2, 2)});
    node.AddBatch(2, {std::pair<int, int>(1, 1), std::pair<int, int> (2, 3)});

    zint2 core;
    core = zint2(ny, nx);
    core = {1, 1, 1, 1,
            1, 2, 2, 1,
            1, 2, 2, 1,
            1, 1, 1, 1};
            
    node.SetCore(core);
    node.PrintMaterial();

    return 0;
}