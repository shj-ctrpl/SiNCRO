#include "../include/CMFD.h"
#include <float.h>

CMFD::CMFD(Node &node, int nx, int ny, int nz, bool FDMmode) : node(node), keff(node.Keff()), fine_Flux(node.Flux()), fine_Jnet(node.Jnet())
{
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    ng = node.ng;
    nd = node.nd;
    this->FDMmode = FDMmode;

    // Check the size between corase-mesh and fine-mesh
    mx = node.nx / nx;
    my = node.ny / ny;
    mz = node.nz / nz;
    hx = node.hx * mx;
    hy = node.hy * my;
    hz = node.hz * mz;
    rhx = 1 / hx;
    rhy = 1 / hy;
    rhz = 1 / hz;
    vol = hx * hy * hz;

    // Memory allocations and variable initializations
    coarse_Old = zdouble4(nz, ny, nx, ng);
    coarse_Flux = zdouble4(nz, ny, nx, ng);
    coarse_Flux = 1.0;
    coarse_Jnet = zdouble5(nz + 1, ny + 1, nx + 1, nd, ng);
    coarse_R = zdouble4(nz, ny, nx, ng);
    coarse_Src = zdouble4(nz, ny, nx, ng);
    coarse_F = zdouble5(nz, ny, nx, ng, ng);
    coarse_S = zdouble5(nz, ny, nx, ng, ng);
    coarse_D = zdouble4(nz, ny, nx, ng);
    coarse_Dt = zdouble5(nz + 1, ny + 1, nx + 1, nd, ng);
    coarse_Dh = zdouble5(nz + 1, ny + 1, nx + 1, nd, ng);
    coarse_Dh = 0.0;

    sumFluxNew = 0.0;
    sumFluxOld = 0.0;

    if (mx == 1 && my == 1 && mz == 1)
    {
        Samesizemode = true;
    }

    UpdateCrossSection();

    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                for (int g = 0; g < ng; g++)
                {
                    for (int to_g = 0; to_g < ng; to_g++)
                    {
                        sumFluxOld += coarse_Flux(z, y, x, g) * coarse_F(z, y, x, g, to_g);
                    }
                }
            }
        }
    }
}

void CMFD::UpdateCrossSection()
{
    // Generate Homogenized cross-section
    if (Samesizemode)
    {
        for (int z = 0; z < node.nz; z++)
        {
            for (int y = 0; y < node.ny; y++)
            {
                for (int x = 0; x < node.nx; x++)
                {
                    for (int from_g = 0; from_g < ng; from_g++)
                    {
                        coarse_Flux(z, y, x, from_g) = node.flux(z, y, x, from_g);
                        coarse_R(z, y, x, from_g) = node.xsA(z, y, x, from_g);
                        coarse_D(z, y, x, from_g) = node.xsD(z, y, x, from_g);

                        for (int to_g = 0; to_g < ng; to_g++)
                        {
                            coarse_S(z, y, x, from_g, to_g) = node.xsS(z, y, x, from_g, to_g);
                            coarse_F(z, y, x, from_g, to_g) = node.xsF(z, y, x, from_g, to_g);
                            coarse_R(z, y, x, from_g) += coarse_S(z, y, x, from_g, to_g);
                        }
                    }
                }
            }
        }
    }
    else
    {
        double div = (mx * my * mz);
        div = 1 / div;
        coarse_Flux = 0.0;
        coarse_R = 0.0;
        coarse_D = 0.0;
        coarse_S = 0.0;
        coarse_F = 0.0;

        for (int z = 0; z < node.nz; z++)
        {
            for (int y = 0; y < node.ny; y++)
            {
                for (int x = 0; x < node.nx; x++)
                {
                    int cz = z / mz, cy = y / my, cx = x / mx;
                    for (int from_g = 0; from_g < ng; from_g++)
                    {
                        coarse_Flux(cz, cy, cx, from_g) += node.flux(z, y, x, from_g);
                        coarse_R(cz, cy, cx, from_g) += node.xsA(z, y, x, from_g) * node.flux(z, y, x, from_g);
                        coarse_D(cz, cy, cx, from_g) += node.flux(z, y, x, from_g) / node.xsD(z, y, x, from_g);

                        for (int to_g = 0; to_g < ng; to_g++)
                        {
                            coarse_S(cz, cy, cx, from_g, to_g) += node.xsS(z, y, x, from_g, to_g) * node.flux(z, y, x, from_g);
                            coarse_F(cz, cy, cx, from_g, to_g) += node.xsF(z, y, x, from_g, to_g) * node.flux(z, y, x, from_g);
                        }
                    }
                }
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    for (int from_g = 0; from_g < ng; from_g++)
                    {
                        coarse_R(z, y, x, from_g) = coarse_R(z, y, x, from_g) / coarse_Flux(z, y, x, from_g);
                        coarse_D(z, y, x, from_g) = coarse_Flux(z, y, x, from_g) / coarse_D(z, y, x, from_g);
                        for (int to_g = 0; to_g < ng; to_g++)
                        {
                            coarse_S(z, y, x, from_g, to_g) = coarse_S(z, y, x, from_g, to_g) / coarse_Flux(z, y, x, from_g);
                            coarse_F(z, y, x, from_g, to_g) = coarse_F(z, y, x, from_g, to_g) / coarse_Flux(z, y, x, from_g);
                            coarse_R(z, y, x, from_g) += coarse_S(z, y, x, from_g, to_g);
                        }
                        coarse_Flux(z, y, x, from_g) = coarse_Flux(z, y, x, from_g) * div;
                    }
                }
            }
        }
    }

    // Initialize D-tilde values
    for (int g = 0; g < ng; g++)
    {
        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                coarse_Dt(z, y, 0, XDIR, g) = harmonicMean(node.Albedo(XDIR, LEFT) * 0.5, coarse_D(z, y, 0, g) * rhx);
                for (int x = 0; x < nx - 1; x++)
                {
                    coarse_Dt(z, y, x + 1, XDIR, g) = harmonicMean(coarse_D(z, y, x, g) * rhx, coarse_D(z, y, x + 1, g) * rhx);
                }
                coarse_Dt(z, y, nx, XDIR, g) = harmonicMean(node.Albedo(XDIR, RIGHT) * 0.5, coarse_D(z, y, nx - 1, g) * rhx);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int x = 0; x < nx; x++)
            {
                coarse_Dt(z, 0, x, YDIR, g) = harmonicMean(node.Albedo(YDIR, LEFT) * 0.5, coarse_D(z, 0, x, g) * rhy);
                for (int y = 0; y < ny - 1; y++)
                {
                    coarse_Dt(z, y + 1, x, YDIR, g) = harmonicMean(coarse_D(z, y, x, g) * rhy, coarse_D(z, y + 1, x, g) * rhy);
                }
                coarse_Dt(z, ny, x, YDIR, g) = harmonicMean(node.Albedo(YDIR, RIGHT) * 0.5, coarse_D(z, ny - 1, x, g) * rhy);
            }
        }

        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                coarse_Dt(0, y, x, ZDIR, g) = harmonicMean(node.Albedo(ZDIR, LEFT) * 0.5, coarse_D(0, y, x, g) * rhz);
                for (int z = 0; z < nz - 1; z++)
                {
                    coarse_Dt(z + 1, y, x, ZDIR, g) = harmonicMean(coarse_D(z + 1, y, x, g) * rhz, coarse_D(z, y, x, g) * rhz);
                }
                coarse_Dt(nz, y, x, ZDIR, g) = harmonicMean(node.Albedo(ZDIR, RIGHT) * 0.5, coarse_D(nz - 1, y, x, g) * rhz);
            }
        }
    }
}

void CMFD::InitFDM()
{
    if (Samesizemode)
    {
        for (int z = 0; z < node.nz; z++)
        {
            for (int y = 0; y < node.ny; y++)
            {
                for (int x = 0; x < node.nx; x++)
                {
                    for (int from_g = 0; from_g < ng; from_g++)
                    {
                        coarse_Flux(z, y, x, from_g) = node.flux(z, y, x, from_g);
                    }
                }
            }
        }
    }
}

void CMFD::UpdateSource()
{
    coarse_Src = 0.0;

    // Initialize source term
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                for (int from_g = 0; from_g < ng; from_g++)
                {
                    for (int to_g = 0; to_g < ng; to_g++)
                    {
                        coarse_Src(z, y, x, to_g) += coarse_F(z, y, x, from_g, to_g) * coarse_Flux(z, y, x, from_g) / keff;
                        coarse_Src(z, y, x, to_g) += coarse_S(z, y, x, from_g, to_g) * coarse_Flux(z, y, x, from_g);
                    }
                }
            }
        }
    }
}

void CMFD::UpdateFlux(int numIter)
{
    for (int i = 0; i < numIter; i++)
    {
        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    for (int g = 0; g < ng; g++)
                    {
                        double src = 0;
                        double center = 0;
                        coarse_Old(z, y, x, g) = coarse_Flux(z, y, x, g);
                        center += (coarse_Dt(z, y, x, XDIR, g) + coarse_Dt(z, y, x + 1, XDIR, g) - coarse_Dh(z, y, x, XDIR, g) + coarse_Dh(z, y, x + 1, XDIR, g)) * rhx;
                        center += (coarse_Dt(z, y, x, YDIR, g) + coarse_Dt(z, y + 1, x, YDIR, g) - coarse_Dh(z, y, x, YDIR, g) + coarse_Dh(z, y + 1, x, YDIR, g)) * rhy;
                        center += (coarse_Dt(z, y, x, ZDIR, g) + coarse_Dt(z + 1, y, x, ZDIR, g) - coarse_Dh(z, y, x, ZDIR, g) + coarse_Dh(z + 1, y, x, ZDIR, g)) * rhz;
                        center += coarse_R(z, y, x, g);

                        src = coarse_Src(z, y, x, g);

                        // XDIR Left boudndary
                        if (x != 0)
                        {
                            src += (coarse_Dt(z, y, x, XDIR, g) + coarse_Dh(z, y, x, XDIR, g)) * coarse_Flux(z, y, x - 1, g) * rhx;
                        }

                        // XDIR Right boudndary
                        if (x != nx - 1)
                        {
                            src += (coarse_Dt(z, y, x + 1, XDIR, g) - coarse_Dh(z, y, x + 1, XDIR, g)) * coarse_Flux(z, y, x + 1, g) * rhx;
                        }

                        // YDIR Left boudndary
                        if (y != 0)
                        {
                            src += (coarse_Dt(z, y, x, YDIR, g) + coarse_Dh(z, y, x, YDIR, g)) * coarse_Flux(z, y - 1, x, g) * rhy;
                        }

                        // YDIR Right boudndary
                        if (y != ny - 1)
                        {
                            src += (coarse_Dt(z, y + 1, x, YDIR, g) - coarse_Dh(z, y + 1, x, YDIR, g)) * coarse_Flux(z, y + 1, x, g) * rhy;
                        }

                        // ZDIR Left boudndary
                        if (z != 0)
                        {
                            src += (coarse_Dt(z, y, x, ZDIR, g) + coarse_Dh(z, y, x, ZDIR, g)) * coarse_Flux(z - 1, y, x, g) * rhz;
                        }

                        // ZDIR Right boudndary
                        if (z != nz - 1)
                        {
                            src += (coarse_Dt(z + 1, y, x, ZDIR, g) - coarse_Dh(z + 1, y, x, ZDIR, g)) * coarse_Flux(z + 1, y, x, g) * rhz;
                        }

                        coarse_Flux(z, y, x, g) = src / center;
                    }
                }
            }
        }
    }
}

void CMFD::UpdateKeff(bool print)
{
    sumFluxNew = 0;

    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                for (int g = 0; g < ng; g++)
                {
                    for (int to_g = 0; to_g < ng; to_g++)
                    {
                        sumFluxNew += coarse_Flux(z, y, x, g) * coarse_F(z, y, x, g, to_g);
                    }
                }
            }
        }
    }

    keff = keff * (sumFluxNew / sumFluxOld);

    if (print)
    {
        if (FDMmode)
        {
            std::cout << "FDMmode, oldSum: " << sumFluxOld << " newSum: " << sumFluxNew << "    ";
            std::cout << "FDM keff: " << keff << std::endl;
        }
        else
        {
            std::cout << "CMFDmode, oldSum: " << sumFluxOld << " newSum: " << sumFluxNew << "    ";
            std::cout << "CMFD keff: " << keff << std::endl;
        }
    }

    sumFluxOld = sumFluxNew;
}

void CMFD::UpdateJnet()
{
    for (int g = 0; g < ng; g++)
    {
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                coarse_Jnet(0, y, x, ZDIR, g) = -(coarse_Dt(0, y, x, ZDIR, g) + coarse_Dh(0, y, x, ZDIR, g)) * coarse_Flux(0, y, x, g);
                for (int z = 1; z < nz; z++)
                {
                    coarse_Jnet(z, y, x, ZDIR, g) = -coarse_Dt(z, y, x, ZDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z - 1, y, x, g)) - coarse_Dh(z, y, x, ZDIR, g) * (coarse_Flux(z, y, x, g) + coarse_Flux(z - 1, y, x, g));
                }
                coarse_Jnet(nz, y, x, ZDIR, g) = (coarse_Dt(nz, y, x, ZDIR, g) - coarse_Dh(nz, y, x, ZDIR, g)) * coarse_Flux(nz - 1, y, x, g);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int x = 0; x < nx; x++)
            {
                coarse_Jnet(z, 0, x, YDIR, g) = -(coarse_Dt(z, 0, x, YDIR, g) + coarse_Dh(z, 0, x, YDIR, g)) * coarse_Flux(z, 0, x, g);
                for (int y = 1; y < ny; y++)
                {
                    coarse_Jnet(z, y, x, YDIR, g) = -coarse_Dt(z, y, x, YDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z, y - 1, x, g)) - coarse_Dh(z, y, x, YDIR, g) * (coarse_Flux(z, y, x, g) + coarse_Flux(z, y - 1, x, g));
                }
                coarse_Jnet(z, ny, x, YDIR, g) = (coarse_Dt(z, ny, x, YDIR, g) - coarse_Dh(z, ny, x, YDIR, g)) * coarse_Flux(z, ny - 1, x, g);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                coarse_Jnet(z, y, 0, XDIR, g) = -(coarse_Dt(z, y, 0, XDIR, g) + coarse_Dh(z, y, 0, XDIR, g)) * coarse_Flux(z, y, 0, g);
                for (int x = 1; x < nx; x++)
                {
                    coarse_Jnet(z, y, x, XDIR, g) = -coarse_Dt(z, y, x, XDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z, y, x - 1, g)) - coarse_Dh(z, y, x, XDIR, g) * (coarse_Flux(z, y, x, g) + coarse_Flux(z, y, x - 1, g));
                }
                coarse_Jnet(z, y, nx, XDIR, g) = (coarse_Dt(z, y, nx, XDIR, g) - coarse_Dh(z, y, nx, XDIR, g)) * coarse_Flux(z, y, nx - 1, g);
            }
        }
    }

    if (Samesizemode)
    {
        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    for (int g = 0; g < ng; g++)
                    {
                        fine_Flux(z, y, x, g) = coarse_Flux(z, y, x, g);
                    }
                }
            }
        }

        for (int z = 0; z <= nz; z++)
        {
            for (int y = 0; y <= ny; y++)
            {
                for (int x = 0; x <= nx; x++)
                {
                    for (int d = 0; d < nd; d++)
                    {
                        for (int g = 0; g < ng; g++)
                        {
                            fine_Jnet(z, y, x, d, g) = coarse_Jnet(z, y, x, d, g);
                        }
                    }
                }
            }
        }
    }
}

void CMFD::UpdateDhat()
{
    if (FDMmode)
    {
        return;
    }

    double jnet_fdm;

    for (int g = 0; g < ng; g++)
    {
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                jnet_fdm = -coarse_Dt(0, y, x, ZDIR, g) * coarse_Flux(0, y, x, g);
                coarse_Dh(0, y, x, ZDIR, g) = (jnet_fdm - coarse_Jnet(0, y, x, ZDIR, g)) / coarse_Flux(0, y, x, g);
                for (int z = 1; z < nz; z++)
                {
                    jnet_fdm = -coarse_Dt(z, y, x, ZDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z - 1, y, x, g));
                    coarse_Dh(z, y, x, ZDIR, g) = (jnet_fdm - coarse_Jnet(z, y, x, ZDIR, g)) / (coarse_Flux(z, y, x, g) + coarse_Flux(z - 1, y, x, g));
                }
                jnet_fdm = coarse_Dt(nz, y, x, ZDIR, g) * coarse_Flux(nz - 1, y, x, g);
                coarse_Dh(nz, y, x, ZDIR, g) = (jnet_fdm - coarse_Jnet(nz, y, x, ZDIR, g)) / coarse_Flux(nz - 1, y, x, g);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int x = 0; x < nx; x++)
            {
                jnet_fdm = -coarse_Dt(z, 0, x, YDIR, g) * coarse_Flux(z, 0, x, g);
                coarse_Dh(z, 0, x, YDIR, g) = (jnet_fdm - coarse_Jnet(z, 0, x, YDIR, g)) / coarse_Flux(z, 0, x, g);
                for (int y = 1; y < ny; y++)
                {

                    jnet_fdm = -coarse_Dt(z, y, x, YDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z, y - 1, x, g));
                    coarse_Dh(z, y, x, YDIR, g) = (jnet_fdm - coarse_Jnet(z, y, x, YDIR, g)) / (coarse_Flux(z, y, x, g) + coarse_Flux(z, y - 1, x, g));
                }
                jnet_fdm = coarse_Dt(z, ny, x, YDIR, g) * coarse_Flux(z, ny - 1, x, g);
                coarse_Dh(z, ny, x, YDIR, g) = (jnet_fdm - coarse_Jnet(z, ny, x, YDIR, g)) / coarse_Flux(z, ny - 1, x, g);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                jnet_fdm = -coarse_Dt(z, y, 0, XDIR, g) * coarse_Flux(z, y, 0, g);
                coarse_Dh(z, y, 0, XDIR, g) = (jnet_fdm - coarse_Jnet(z, y, 0, XDIR, g)) / coarse_Flux(z, y, 0, g);
                for (int x = 1; x < nx; x++)
                {
                    jnet_fdm = -coarse_Dt(z, y, x, XDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z, y, x - 1, g));
                    coarse_Dh(z, y, x, XDIR, g) = (jnet_fdm - coarse_Jnet(z, y, x, XDIR, g)) / (coarse_Flux(z, y, x, g) + coarse_Flux(z, y, x - 1, g));
                }
                jnet_fdm = coarse_Dt(z, y, nx, XDIR, g) * coarse_Flux(z, y, nx - 1, g);
                coarse_Dh(z, y, nx, XDIR, g) = (jnet_fdm - coarse_Jnet(z, y, nx, XDIR, g)) / coarse_Flux(z, y, nx - 1, g);
            }
        }
    }
}

void CMFD::UpdateDhatNodal()
{
    if (FDMmode)
    {
        return;
    }

    double jnet_fdm;

    for (int g = 0; g < ng; g++)
    {
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                jnet_fdm = -coarse_Dt(0, y, x, ZDIR, g) * coarse_Flux(0, y, x, g);
                coarse_Dh(0, y, x, ZDIR, g) = (jnet_fdm - fine_Jnet(0, y, x, ZDIR, g)) / coarse_Flux(0, y, x, g);
                for (int z = 1; z < nz; z++)
                {
                    jnet_fdm = -coarse_Dt(z, y, x, ZDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z - 1, y, x, g));
                    coarse_Dh(z, y, x, ZDIR, g) = (jnet_fdm - fine_Jnet(z, y, x, ZDIR, g)) / (coarse_Flux(z, y, x, g) + coarse_Flux(z - 1, y, x, g));
                }
                jnet_fdm = coarse_Dt(nz, y, x, ZDIR, g) * coarse_Flux(nz - 1, y, x, g);
                coarse_Dh(nz, y, x, ZDIR, g) = (jnet_fdm - fine_Jnet(nz, y, x, ZDIR, g)) / coarse_Flux(nz - 1, y, x, g);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int x = 0; x < nx; x++)
            {
                jnet_fdm = -coarse_Dt(z, 0, x, YDIR, g) * coarse_Flux(z, 0, x, g);
                coarse_Dh(z, 0, x, YDIR, g) = (jnet_fdm - fine_Jnet(z, 0, x, YDIR, g)) / coarse_Flux(z, 0, x, g);
                for (int y = 1; y < ny; y++)
                {

                    jnet_fdm = -coarse_Dt(z, y, x, YDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z, y - 1, x, g));
                    coarse_Dh(z, y, x, YDIR, g) = (jnet_fdm - fine_Jnet(z, y, x, YDIR, g)) / (coarse_Flux(z, y, x, g) + coarse_Flux(z, y - 1, x, g));
                }
                jnet_fdm = coarse_Dt(z, ny, x, YDIR, g) * coarse_Flux(z, ny - 1, x, g);
                coarse_Dh(z, ny, x, YDIR, g) = (jnet_fdm - fine_Jnet(z, ny, x, YDIR, g)) / coarse_Flux(z, ny - 1, x, g);
            }
        }

        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                jnet_fdm = -coarse_Dt(z, y, 0, XDIR, g) * coarse_Flux(z, y, 0, g);
                coarse_Dh(z, y, 0, XDIR, g) = (jnet_fdm - fine_Jnet(z, y, 0, XDIR, g)) / coarse_Flux(z, y, 0, g);
                for (int x = 1; x < nx; x++)
                {
                    jnet_fdm = -coarse_Dt(z, y, x, XDIR, g) * (coarse_Flux(z, y, x, g) - coarse_Flux(z, y, x - 1, g));
                    coarse_Dh(z, y, x, XDIR, g) = (jnet_fdm - fine_Jnet(z, y, x, XDIR, g)) / (coarse_Flux(z, y, x, g) + coarse_Flux(z, y, x - 1, g));
                }
                jnet_fdm = coarse_Dt(z, y, nx, XDIR, g) * coarse_Flux(z, y, nx - 1, g);
                coarse_Dh(z, y, nx, XDIR, g) = (jnet_fdm - fine_Jnet(z, y, nx, XDIR, g)) / coarse_Flux(z, y, nx - 1, g);
            }
        }
    }
}

void CMFD::Redistribute()
{

    if (Samesizemode)
    {
        for (int z = 0; z < node.nz; z++)
        {
            for (int y = 0; y < node.ny; y++)
            {
                for (int x = 0; x < node.nx; x++)
                {
                    for (int g = 0; g < ng; g++)
                    {
                        node.flux(z, y, x, g) = coarse_Flux(z, y, x, g);
                    }
                }
            }
        }
    }
    else
    {
        for (int z = 0; z < node.nz; z++)
        {
            for (int y = 0; y < node.ny; y++)
            {
                for (int x = 0; x < node.nx; x++)
                {
                    int cz = z / mz, cy = y / my, cx = x / mx;
                    for (int g = 0; g < ng; g++)
                    {
                        node.flux(z, y, x, g) = node.flux(z, y, x, g) * (coarse_Flux(cz, cy, cx, g) / coarse_Old(cz, cy, cx, g));
                    }
                }
            }
        }
    }
}