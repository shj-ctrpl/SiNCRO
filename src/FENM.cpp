#include "../include/FENM.h"


#define pos z,y,x
#define imesh z,y,x,d,g

FENM::FENM(const int nz, const int ny, const int nx, const double hz, const double hy, const double hx, Node &node)
{
    _nz = nz + 2;
    _ny = ny + 2;
    _nx = nx + 2;
    _hz = hz + 2;
    _hy = hy + 2;
    _hx = hx + 2;
    _node = node;

    jnet = zdouble6(_nz, _ny, _nx, _ndir, _ng, 2);
    tlcoeff0 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    tlcoeff1 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    tlcoeff2 = zdouble5(_nz, _ny, _nx, _ndir, _ng);

    mu00 = 2.0;
    mu11 = 2.0 / 3.0;
    mu22 = 2.0 / 5.0;
    mu33 = 2.0 / 7.0;
    mu44 = 2.0 / 9.0;

    mu_m = zdouble6(_nz, _ny, _nx, _ndir, _ng, _ng);

    k20 = 6.00;
    k31 = 10.0;
    k40 = 20.0;
    k42 = 14.0;

    k51 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    k53 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    k60 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    k62 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    k64 = zdouble5(_nz, _ny, _nx, _ndir, _ng);

    T5 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    T6 = zdouble5(_nz, _ny, _nx, _ndir, _ng);

    diagD = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    invdiagD = zdouble5(_nz, _ny, _nx, _ndir, _ng);

    matM = zdouble5(_nz, _ny, _nx, _ng, _ng);
    matS = zdouble5(_nz, _ny, _nx, _ng, _ng);
    matF = zdouble5(_nz, _ny, _nx, _ng, _ng);
    invmatM = zdouble5(_nz, _ny, _nx, _ng, _ng);

    c0 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    c1 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    c2 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    c3 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    c4 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    c5 = zdouble5(_nz, _ny, _nx, _ndir, _ng);
    c6 = zdouble5(_nz, _ny, _nx, _ndir, _ng);

    avgFlux = zdouble4(_nz, _ny, _nx, _ng);
}

FENM::~FENM()
{
}

void FENM::Initialize()
{
}

void FENM::CalcConst(const int z, const int y, const int x)
{
    for (int d = 0; d < _ndir; d++)
    {
        for (int g = 0; g < _ng; g++)
        {
            // Pre-calculated kappa values
            double kp2 = _node.XsTot(pos, g) * _node.HMesh(pos, d) * _node.HMesh(pos, d) / (4 * _node.XsDfs(pos, g));
            double kp = std::sqrt(kp2);
            double kp3 = kp2 * kp;
            double kp4 = kp2 * kp2;
            double kp5 = kp3 * kp2;

            // Pre-calculated 1 over kappa
            double rkp = 1 / kp;
            double rkp2 = rkp * rkp;
            double rkp3 = rkp2 * rkp;
            double rkp4 = rkp2 * rkp2;
            double rkp5 = rkp3 * rkp2;

            // Pre-calculated hyperbolic functions
            double sinhkp = std::sinh(kp);
            double coshkp = std::cosh(kp);

            // Free parameter m_gi
            double mg0 = -sinhkp * rkp;
            double mg1 = 3 * rkp2 * (sinhkp - kp * coshkp);
            double mg2 = -5 * rkp3 * (3 * sinhkp - 3 * kp * coshkp + kp2 * sinhkp);
            double mg3 = 7 * rkp4 * (15 * sinhkp - 15 * kp * coshkp + 6 * kp2 * sinhkp + kp3 * coshkp);
            double mg4 = -9 * rkp5 * (105 * sinhkp - 105 * kp * coshkp + 45 * kp2 * sinhkp - 10 * kp3 * coshkp + kp4 * sinhkp);

            // Common nominator and denominators for constants calculation
            double oddDenom = 1 / (sinhkp + mg1 + mg3);
            double evenDenom = 1 / (coshkp + mg0 + mg2 + mg4);

            T5(imesh) = (kp * coshkp + mg1 + 6 * mg3) * oddDenom;
            T6(imesh) = (kp * sinhkp + 3 * mg2 + 10 * mg4) * evenDenom;

            // Calculate the inner proudcts of the basis functions and 2nd order derivatives
            k51(imesh) = 2 * (kp * coshkp - sinhkp + 5 * mg3) * oddDenom;
            k53(imesh) = 2 * (kp * (15 + kp2) * coshkp - 3 * (5 + 2 * kp2) * sinhkp) * rkp2 * oddDenom;
            k60(imesh) = 2 * T6(imesh);
            k62(imesh) = 2 * (-3 * coshkp + (3 + kp2) * (sinhkp * rkp) + 7 * mg4) * evenDenom;
            k64(imesh) = 2 * (-5 * kp * (21 + 2 * kp2) * coshkp + (105 + 45 * kp2 + kp4) * sinhkp) * rkp3 * evenDenom;

            diagD(imesh) = 4 * _node.XsDfs(pos, g) / (_node.HMesh(pos, d) * _node.HMesh(pos, d));
            invdiagD(imesh) = 1.0 / diagD(imesh);
        }
    }
}

void FENM::CalcTLCoeffs(const int z, const int y, const int x)
{
    // Calculate 0th transverse leakage coefficient, L_x0 = (J_yr-J_yl)/h_y + (J_zr-J_zl)/h_z
    double jAvgZ, jAvgY, jAvgX;

    for (int g = 0; g < _ng; g++)
    {
        jAvgZ = jnet(pos, ZDIR, g, RIGHT) - jnet(pos, ZDIR, g, LEFT);
        jAvgY = jnet(pos, YDIR, g, RIGHT) - jnet(pos, YDIR, g, LEFT);
        jAvgX = jnet(pos, XDIR, g, RIGHT) - jnet(pos, XDIR, g, LEFT);
        tlcoeff0(pos, ZDIR, g) = jAvgY + jAvgX;
        tlcoeff0(pos, YDIR, g) = jAvgZ + jAvgX;
        tlcoeff0(pos, XDIR, g) = jAvgZ + jAvgY;
    }

    // Calculate the 1st and 2nd transverse leakage

    for (int d = 0; d < _ndir; d++)
    {
        for (int g = 0; g < _ng; g++)
        {
            double h[3];
            double TL[3];
            double tempDenom;

            // Check if LUT is better...

            h[LEFT] = _node.HMesh(z - (d == ZDIR), y - (d == YDIR), x - (d == XDIR), d);
            h[RIGHT] = _node.HMesh(z + (d == ZDIR), y + (d == YDIR), x + (d == XDIR), d);
            h[CENTER] = _node.HMesh(pos, d);

            TL[CENTER] = tlcoeff0(imesh);

            if ((x == 0 || y == 0 || z == 0) && _node.Albedo(d, LEFT) < EPS)
            {
                TL[LEFT] = TL[CENTER];
            }
            else
            {
                TL[LEFT] = tlcoeff0(z - (d == ZDIR), y - (d == YDIR), x - (d == XDIR), d, g);
            }

            TL[RIGHT] = tlcoeff0(z + (d == ZDIR), y + (d == YDIR), x + (d == XDIR), d, g);
            // Consider albedo conditions
            if (_node.Albedo(d, LEFT) < EPS)
            {
                TL[LEFT] = TL[CENTER];
            }
            else
            {
            }

            // Left boundary
            if (h[LEFT] < EPS)
            {
                tlcoeff1 = 0.125 * (5.0 * TL[CENTER] + TL[RIGHT]);
                tlcoeff2 = 0.125 * (-3.0 * TL[CENTER] + TL[RIGHT]);
            }

            // Right boundary
            else if (h[RIGHT] < EPS)
            {
                tlcoeff1 = -0.125 * (5.0 * TL[CENTER] + TL[LEFT]);
                tlcoeff2 = 0.125 * (-3.0 * TL[CENTER] + TL[LEFT]);
            }

            // Inner nodes
            else
            {
                tempDenom = 1 / ((h[CENTER] + h[LEFT]) * (h[CENTER] + h[RIGHT]) * (h[CENTER] + h[LEFT] + h[RIGHT]));
                tlcoeff1(imesh) = 0.5 * tempDenom * h[CENTER] *
                                          ((TL[CENTER] - TL[LEFT]) * (h[CENTER] + 2 * h[RIGHT]) * (h[CENTER] + h[RIGHT]) + (TL[RIGHT] - TL[CENTER]) * (2 * h[LEFT] + h[CENTER]) * (h[LEFT] + h[CENTER]));

                tlcoeff2(imesh) = 0.5 * tempDenom * (h[CENTER] * h[CENTER]) *
                                          ((TL[LEFT] - TL[CENTER]) * (h[CENTER] + h[RIGHT]) + (TL[RIGHT] - TL[CENTER]) * (h[LEFT] + h[CENTER]));
            }
        }
    }
}

void FENM::CalcMatrix(const int z, const int y, const int x)
{
    for (int to_g = 0; to_g < _ng; to_g++)
    {
        matS(pos, to_g, to_g) = _node.XsTot(pos, to_g);
        for (int from_g = 0; from_g < _ng; from_g++)
        {
            matS(pos, from_g, to_g) -= _node.XsSct(pos, from_g, to_g);
            matF(pos, from_g, to_g) = _node.XsFis(pos, from_g, to_g);
        }
    }
    matM = matS - revk * matF;

    double det = matM(pos, 0, 0) * matM(pos, 1, 1) - matM(pos, 1, 0) * matM(pos, 0, 1);
    if (std::abs(det) > EPS)
    {
        double rdet = 1 / det;
        invmatM(pos, 0, 0) = rdet * matM(pos, 1, 1);
        invmatM(pos, 0, 1) = -rdet * matM(pos, 0, 1);
        invmatM(pos, 1, 0) = -rdet * matM(pos, 1, 0);
        invmatM(pos, 1, 1) = rdet * matM(pos, 0, 0);
    }
    else
    {
        invmatM(pos, 0, 0) = 0;
        invmatM(pos, 0, 1) = 0;
        invmatM(pos, 1, 0) = 0;
        invmatM(pos, 1, 1) = 0;
    }

    for (int d = 0; d < _ndir; d++)
    {
        // H = mu_33 * invD * invk_53 * matM
        zdouble2 H(2, 2);
        H(0, 0) = mu33 * invdiagD(pos, d, 0) * matM(pos, 0, 0) / k53(pos, d, 0);
        H(0, 1) = mu33 * invdiagD(pos, d, 0) * matM(pos, 0, 1) / k53(pos, d, 0);
        H(1, 0) = mu33 * invdiagD(pos, d, 1) * matM(pos, 1, 0) / k53(pos, d, 1);
        H(1, 1) = mu33 * invdiagD(pos, d, 1) * matM(pos, 1, 1) / k53(pos, d, 1);

        // mu_m = (1/mu11) * invmatM * D (k_31 + k_51H)
        double rmu11 = 1 / mu11;
        mu_m(pos, d, 0, 0) = rmu11 * (invmatM(pos, 0, 0) * diagD(pos, d, 0) * (k31 + k51(pos, d, 0) * H(0, 0)) +
                                          (invmatM(pos, 0, 1) * diagD(pos, d, 1) * k51(pos, d, 1) * H(1, 0)));

        mu_m(pos, d, 0, 1) = rmu11 * (invmatM(pos, 0, 1) * diagD(pos, d, 0) * (k31 + k51(pos, d, 1) * H(1, 1)) +
                                          (invmatM(pos, 0, 0) * diagD(pos, d, 0) * k51(pos, d, 0) * H(0, 1)));

        mu_m(pos, d, 1, 0) = rmu11 * (invmatM(pos, 1, 0) * diagD(pos, d, 1) * (k31 + k51(pos, d, 0) * H(0, 0)) +
                                          (invmatM(pos, 1, 1) * diagD(pos, d, 1) * k51(pos, d, 1) * H(1, 0)));

        mu_m(pos, d, 1, 1) = rmu11 * (invmatM(pos, 1, 1) * diagD(pos, d, 1) * (k31 + k51(pos, d, 1) * H(1, 1)) +
                                          (invmatM(pos, 1, 0) * diagD(pos, d, 0) * k51(pos, d, 0) * H(0, 1)));
    }
}

void FENM::CalcEvenCoeffs(const int z, const int y, const int x)
{
    zdouble2 A(2, 2);
    zdouble2 T(2, 2); // mu44 * invk64 * invD * A

    zdouble1 D(2);
    zdouble1 DI(2); // inversed diag D
    zdouble1 K60(2);
    zdouble1 K62(2);

    zdouble2 B(2, 2);
    zdouble1 a(2); // l_0 + A*Phi
    zdouble1 b(2);

    for (int d = 0; d <= _ndir; d++)
    {
        A = matM.slice(pos);
        D(0) = diagD(pos, d, 0);
        D(1) = diagD(pos, d, 1);
        DI(0) = invdiagD(pos, d, 0);
        DI(1) = invdiagD(pos, d, 1);
        K60(0) = k60(pos, d, 0);
        K60(1) = k60(pos, d, 1);
        K62(0) = k62(pos, d, 0);
        K62(1) = k62(pos, d, 1);

        T(0, 0) = 4.5 * A(0, 0) * DI(0) / k64(pos, d, 0);
        T(0, 1) = 4.5 * A(0, 1) * DI(0) / k64(pos, d, 0);
        T(1, 0) = 4.5 * A(1, 0) * DI(1) / k64(pos, d, 1);
        T(1, 1) = 4.5 * A(1, 1) * DI(1) / k64(pos, d, 1);

        B(0, 0) = 15 * (A(0, 0) * (K60(0) * T(0, 0) + 20) + A(0, 1) * K60(1) * T(0, 0)) + D(0) * (14 + K62(0) * T(0, 0));
        B(1, 0) = T(0, 1) * (15 * A(0, 0) * K60(0) + D(0) * K62(0)) + 15 * A(0, 1) * (K60(1) * T(1, 1) + 20);
        B(0, 1) = T(1, 0) * (15 * A(1, 1) * K60(1) + D(1) * K62(1)) + 15 * A(1, 0) * (K60(0) * T(0, 0) + 20);
        B(1, 1) = 15 * (A(1, 0) * (K60(0) * T(0, 1) + 20) + A(1, 1) * K60(1) * T(1, 1)) + D(1) * (14 + K62(1) * T(1, 1));

        a(0) = A(0, 0) * avgFlux(pos, 0) + A(0, 1) * avgFlux(pos, 1) + tlcoeff0(pos, d, 0);
        a(1) = A(1, 0) * avgFlux(pos, 0) + A(1, 1) * avgFlux(pos, 1) + tlcoeff0(pos, d, 1);

        b(0) = (-0.833333333) * (a(0) * A(0, 0) * DI(0) + a(1) * A(0, 1) * DI(1) - 3 * tlcoeff2(pos, d, 0));
        b(1) = (-0.833333333) * (a(0) * A(1, 0) * DI(0) + a(1) * A(1, 1) * DI(1) - 3 * tlcoeff2(pos, d, 1));

        double rdet = B(0, 0) * B(1, 1) - B(1, 0) * B(0, 1);
        if (abs(rdet) > EPS)
        {
            rdet = 1 / rdet;
            c4(pos, d, 0) = rdet * B(1, 1) * b(0) - B(0, 1) * b(1);
            c4(pos, d, 1) = rdet * B(0, 0) * b(0) - B(1, 0) * b(1);
        }
        else
        {
            c4(pos, d, 0) = 0;
            c4(pos, d, 1) = 0;
        }
        c6(pos, d, 0) = T(0, 0) * c4(pos, d, 0) + T(0, 1) * c4(pos, d, 1);
        c6(pos, d, 1) = T(1, 0) * c4(pos, d, 0) + T(1, 1) * c4(pos, d, 1);

        c2(pos, d, 0) = 0.166666667 * (2 * a(0) - D(0) * (c4(pos, d, 0) * (T(0, 0) + 20) + c4(pos, d, 1) * T(0, 1)));
        c2(pos, d, 1) = 0.166666667 * (2 * a(1) - D(1) * (c4(pos, d, 0) * T(1, 0) + c4(pos, d, 1) * (T(1, 1) + 20)));

        // Old code
        /*
        // Temporary 2-by-2 material matrix
        zdouble2 M(2, 2);
        M = matM.slice(pos);

        // tmp1 = μ_22 * A * (invK_20)
        zdouble2 tmp1(2, 2);
        tmp1 = 2 * mu22 * M / k20;

        // tmp2 = μ_44 * k_60 * (invk_64) * invD * A
        zdouble2 tmp2(2, 2);
        tmp2(0, 0) = mu44 * k60(pos, d, 0) * invdiagD(pos, d, 0) * M(0, 0) / k64(pos, d, 0);
        tmp2(0, 1) = mu44 * k60(pos, d, 0) * invdiagD(pos, d, 0) * M(0, 1) / k64(pos, d, 0);
        tmp2(1, 0) = mu44 * k60(pos, d, 1) * invdiagD(pos, d, 1) * M(1, 0) / k64(pos, d, 1);
        tmp2(1, 1) = mu44 * k60(pos, d, 1) * invdiagD(pos, d, 1) * M(1, 1) / k64(pos, d, 1);

        // B = (DK_42 + mu44 * D * K_62 * invK_64 * invD * A + mu22 * A * invK_20 (K_40 + mu44 * K_60 * invK_64 * invD * A))
        zdouble2 B(2, 2);
        B(0, 0) = (diagD(pos, d, 0) * k42);
        B(0, 0) += mu44 * diagD(pos, d, 0) * k62(pos, d, 0) * invdiagD(pos, d, 0) * M(0, 0) / k64(pos, d, 0);
        B(0, 0) += k40 * tmp1(0, 0) + tmp1(0, 0) * tmp2(0, 0) + tmp1(1, 0) * tmp2(0, 1);

        B(0, 1) = (diagD(pos, d, 0) * k42);
        B(0, 1) += mu44 * diagD(pos, d, 0) * k62(pos, d, 0) * invdiagD(pos, d, 0) * M(0, 1) / k64(pos, d, 0);
        B(0, 1) += k40 * tmp1(0, 1) + tmp1(0, 1) * tmp2(0, 0) + tmp1(1, 1) * tmp2(0, 1);

        B(1, 0) = (diagD(pos, d, 1) * k42);
        B(1, 0) += mu44 * diagD(pos, d, 1) * k62(pos, d, 1) * invdiagD(pos, d, 1) * M(1, 0) / k64(pos, d, 1);
        B(1, 0) += k40 * tmp1(1, 0) + tmp1(0, 0) * tmp2(1, 0) + tmp1(0, 1) * tmp2(1, 0);

        B(1, 1) = (diagD(pos, d, 1) * k42);
        B(1, 1) += mu44 * diagD(pos, d, 1) * k62(pos, d, 1) * invdiagD(pos, d, 1) * M(1, 1) / k64(pos, d, 1);
        B(1, 1) += k40 * tmp1(1, 1) + tmp1(0, 1) * tmp2(1, 0) + tmp1(1, 1) * tmp2(1, 1);

        // tmp3 = l_0 + A * Phibar
        zdouble1 tmp3(2);
        tmp3(0) = tlcoeff0(pos, d, 0) + M(0, 0) * avgFlux(pos, 0) + M(1, 0) * avgFlux(pos, 1);
        tmp3(1) = tlcoeff0(pos, d, 0) + M(0, 1) * avgFlux(pos, 0) + M(1, 1) * avgFlux(pos, 1);

        // b = (-mu22 * L2 - 2 * mu22 * A * invK_20 * invD * (L0 + A * Phi))
        zdouble1 b(2);
        b(0) = -mu22 * tlcoeff2(pos, d, 0) - 2 * tmp1(0, 0) * invdiagD(pos, d, 0) * tmp3(0) + tmp1(1, 0) * invdiagD(pos, d, 1) * tmp3(1);
        b(1) = -mu22 * tlcoeff2(pos, d, 1) - 2 * tmp1(0, 1) * invdiagD(pos, d, 0) * tmp3(0) + tmp1(1, 1) * invdiagD(pos, d, 1) * tmp3(1);

        double rdet = 1 / (B(0, 0) * B(1, 1) - B(0, 1) * B(1, 0));
        if (abs(rdet) > EPS)
        {
            c4(pos, d, 0) = rdet * (B(1, 1) * b(0) - B(1, 0) * b(1));
            c4(pos, d, 1) = rdet * (B(0, 0) * b(0) - B(0, 1) * b(1));
        }
        else
        {
            c4(pos, d, 0) = 0;
            c4(pos, d, 1) = 0;
        }

        c6(pos, d, 0) = tmp2(0, 0) * c4(pos, d, 0) / k60(pos, d, 0) + tmp2(1, 0) * c4(pos, d, 1) / k60(pos, d, 1);
        c6(pos, d, 1) = tmp2(0, 1) * c4(pos, d, 0) / k60(pos, d, 0) + tmp2(1, 1) * c4(pos, d, 1) / k60(pos, d, 1);
        c2(pos, d, 0) = 2 * invdiagD() - ((k40 + tmp2(0, 0)) * c4(pos, d, 0) + (k40 + tmp2(1, 0)) * c4(pos, d, 1));
        c2(pos, d, 1) = 0;
        */
    }
}

void FENM::CalcOddCoeffsOneNode(const int z, const int y, const int x, const int d, const int face)
{
    zdouble2 D(2, 2);
    zdouble2 Dt(2, 2);
    zdouble2 a11(2, 2);
    zdouble1 a12(2);
    zdouble1 a13(2);
    zdouble2 a22(2, 2);
    zdouble1 a23(2);
    zdouble1 a31(2);
    zdouble1 a32(2);
    zdouble1 a33(2);
    zdouble1 b1(2);
    zdouble1 b2(2);

    // a11 = mu_11 * A
    a11 = mu11 * matM.slice(pos);

    // a22 = mu_33 * A
    a22 = mu33 * matM.slice(pos);

    for (int g = 0; g < _ng; g++)
    {
        // a12 = -DK_31
        a12(g) = -k31 * diagD(imesh);

        // a13 = -DK_51
        a13(g) = -k51(imesh) * diagD(imesh);

        // a23 = -DK_53
        a23(g) = -k53(imesh) * diagD(imesh);

        // a31 = Dtilde_m + albedo
        a31(g) = 0.5 * _node.HMesh(pos, d) * diagD(imesh) + _node.Albedo(d, face);

        // a32 = 6 * Dtilde_m + albedo
        a32(g) = 3 * _node.HMesh(pos, d) * diagD(imesh) + _node.Albedo(d, face);

        // a33 = Dtilde_m * T5m + albedo
        a33(g) = T5(imesh) * _node.HMesh(pos, d) * diagD(imesh) + _node.Albedo(d, face);

        // b1 = -m011 * l_0
        b1(g) = -mu11 * tlcoeff0(imesh);

        // b2 = D * (3 * c_2 + 10 * c_4 + T6 * c_6) + albedo * (flux + c_2 + c_4 + c_6)
        b2(g) = diagD(imesh) * (3 * c2(imesh) + 10 * c4(imesh) + T6(imesh) * c6(imesh));
        b2(g) += _node.Albedo(d, face) * (avgFlux(pos, g) + c2(imesh) + c4(imesh) + c6(imesh));
    }

    zdouble2 A(2, 2);

    a22(0, 0) = a22(0, 0) / a23(0);
    a22(0, 1) = a22(0, 1) / a23(0);
    a22(1, 0) = a22(1, 0) / a23(1);
    a22(1, 1) = a22(1, 1) / a23(1);

    a11(0, 0) = a11(0, 0) / a31(0);
    a11(0, 1) = a11(0, 1) / a31(0);
    a11(1, 0) = a11(1, 0) / a31(1);
    a11(1, 1) = a11(1, 1) / a31(1);

    A(0, 0) = a22(0, 0) * a13(0) - a11(0, 0) - a32(0) + a12(0);
    A(0, 1) = a22(0, 1) * a13(0) - a11(0, 1) - a32(0);
    A(1, 0) = a22(1, 0) * a13(1) - a11(1, 0) - a32(1);
    A(1, 1) = a22(1, 1) * a13(1) - a11(1, 1) - a32(1) + a12(1);

    A(0, 0) -= (a11(0, 0) * a22(0, 0) + a11(0, 1) * a22(1, 0)) * a33(0);
    A(0, 1) -= (a11(0, 0) * a22(0, 1) + a11(0, 1) * a22(1, 1)) * a33(0);
    A(1, 0) -= (a11(1, 0) * a22(0, 0) + a11(1, 1) * a22(1, 0)) * a33(1);
    A(1, 1) -= (a11(1, 0) * a22(0, 1) + a11(1, 1) * a22(1, 1)) * a33(1);

    b1(0) = b1(0) - a11(0, 0) * b2(0) - a11(0, 1) * b2(1);
    b1(1) = b1(1) - a11(1, 0) * b2(0) - a11(1, 1) * b2(1);

    double rdet = 1 / (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
    zdouble2 invA(2, 2);

    invA(0, 0) = rdet * A(1, 1);
    invA(0, 1) = -rdet * A(1, 0);
    invA(1, 0) = -rdet * A(0, 1);
    invA(1, 1) = rdet * A(0, 0);

    c3(pos, d, 0) = invA(0, 0) * b1(0) + invA(0, 1) * b1(1);
    c3(pos, d, 1) = invA(1, 0) * b1(0) + invA(1, 1) * b1(1);

    c5(pos, d, 0) = a22(0, 0) * c3(pos, d, 0) + a22(0, 1) * c3(pos, d, 1);
    c5(pos, d, 1) = a22(1, 0) * c3(pos, d, 0) + a22(1, 1) * c3(pos, d, 1);

    c1(pos, d, 0) = b2(0) - a32(0) * c3(pos, d, 0) - a33(0) * c5(pos, d, 0) / a31(0);
    c1(pos, d, 1) = b2(1) - a32(1) * c3(pos, d, 1) - a33(1) * c5(pos, d, 1) / a31(1);

    int sign = 0;

    for (int g = 0; g < _ng; g++)
    {
        jnet(pos, d, g, face) = _node.HMesh(pos, d) * 0.5 * D(g, g) *
                                    (c1(imesh) + 6 * c3(imesh) + T5(imesh) * c5(imesh) + sign * (3 * c2(imesh) + 10 * c4(imesh) + T6(imesh) * c6(imesh)));
        surfFlux(pos, d, g, face) = avgFlux(pos, g) + sign * (c1(imesh) + c3(imesh) + c5(imesh)) + c2(imesh) + c4(imesh) + c6(imesh);
    }
}

void FENM::CalcOddCoeffsTwoNode(const int z, const int y, const int x)
{
    double rdet;

    zdouble2 Al(2, 2);
    zdouble2 Ali(2, 2);
    zdouble2 Hl(2, 2);
    zdouble2 mul(2, 2);
    zdouble1 Dl(2);
    zdouble1 Dtl(2);

    zdouble2 Ar(2, 2);
    zdouble2 Ari(2, 2);
    zdouble2 Hr(2, 2);
    zdouble2 mur(2, 2);
    zdouble1 Dr(2);
    zdouble1 Dtr(2);

    zdouble2 tmp(2, 2);
    zdouble1 tmp2(2);
    zdouble2 Z1(2, 2);
    zdouble1 Z2(2);

    zdouble1 bp(2);
    zdouble1 bj(2);

    zdouble2 M(2, 2);
    zdouble1 vg(2);

    for (int d = 0; d <= _ndir; d++)
    {
        int zr = z + (d == ZDIR);
        int yr = y + (d == YDIR);
        int xr = x + (d == XDIR);

        Al = matM.slice(pos);
        Ali = invmatM.slice(pos);
        Ar = matM.slice(zr, yr, xr);
        Ari = invmatM.slice(zr, yr, xr);

        Dl(0) = diagD(pos, d, 0);
        Dl(1) = diagD(pos, d, 1);
        Dtl = 0.5 * _node.HMesh(pos, d) * Dl;
        Dr(0) = diagD(zr, yr, xr, d, 0);
        Dr(1) = diagD(zr, yr, xr, d, 1);
        Dtr = 0.5 * _node.HMesh(zr, yr, xr, d) * Dr;

        bp(0) = avgFlux(pos, 0) + c2(pos, d, 0) + c4(pos, d, 0) + c6(pos, d, 0);
        bp(1) = avgFlux(pos, 1) + c2(pos, d, 1) + c4(pos, d, 1) + c6(pos, d, 1);
        bp(0) += avgFlux(zr, yr, xr, 0) + c2(zr, yr, xr, d, 0) + c4(zr, yr, xr, d, 0) + c6(zr, yr, xr, d, 0);
        bp(1) += avgFlux(zr, yr, xr, 1) + c2(zr, yr, xr, d, 1) + c4(zr, yr, xr, d, 1) + c6(zr, yr, xr, d, 1);

        bj(0) = Dtl(0) * (3 * c2(pos, d, 0) + 10 * c2(pos, d, 0) + T6(pos, d, 0) * c6(pos, d, 0));
        bj(1) = Dtl(1) * (3 * c2(pos, d, 1) + 10 * c2(pos, d, 1) + T6(pos, d, 1) * c6(pos, d, 1));
        bj(0) += Dtr(0) * (3 * c2(zr, yr, xr, d, 0) + 10 * c2(zr, yr, xr, d, 0) + T6(zr, yr, xr, d, 0) * c6(zr, yr, xr, d, 0));
        bj(1) += Dtr(1) * (3 * c2(zr, yr, xr, d, 1) + 10 * c2(zr, yr, xr, d, 1) + T6(zr, yr, xr, d, 1) * c6(zr, yr, xr, d, 1));

        Hl(0, 0) = mu33 * Al(0, 0) / (Dl(0) * k53(pos, d, 0));
        Hl(0, 1) = mu33 * Al(0, 1) / (Dl(0) * k53(pos, d, 0));
        Hl(1, 0) = mu33 * Al(1, 0) / (Dl(1) * k53(pos, d, 1));
        Hl(1, 1) = mu33 * Al(1, 1) / (Dl(1) * k53(pos, d, 1));
        Hr(0, 0) = mu33 * Ar(0, 0) / (Dr(0) * k53(zr, yr, xr, d, 0));
        Hr(0, 1) = mu33 * Ar(0, 1) / (Dr(0) * k53(zr, yr, xr, d, 0));
        Hr(1, 0) = mu33 * Ar(1, 0) / (Dr(1) * k53(zr, yr, xr, d, 1));
        Hr(1, 1) = mu33 * Ar(1, 1) / (Dr(1) * k53(zr, yr, xr, d, 1));

        mul(0, 0) = -Dl(0) * (10 + Hl(0, 0) + k51(pos, d, 0));
        mul(0, 1) = -Dl(0) * Hl(0, 1) * k51(pos, d, 0);
        mul(1, 0) = -Dl(1) * Hl(1, 0) * k51(pos, d, 1);
        mul(1, 1) = -Dl(1) * (10 + Hl(1, 1) + k51(pos, d, 1));
        mur(0, 0) = -Dr(0) * (10 + Hr(0, 0) + k51(zr, yr, xr, d, 0));
        mur(0, 1) = -Dr(0) * Hr(0, 1) * k51(zr, yr, xr, d, 0);
        mur(1, 0) = -Dr(1) * Hr(1, 0) * k51(zr, yr, xr, d, 1);
        mur(1, 1) = -Dr(1) * (10 + Hl(1, 1) + k51(zr, yr, xr, d, 1));

        rdet = 1 / ((1 + Hr(0, 0) + mur(0, 0)) * (1 + Hr(1, 1) + mur(1, 1)) - (Hr(0, 1) + mur(0, 1)) * (Hr(1, 0) + mur(1, 0)));

        // tmp = inv(I + H_r + mu_r)
        tmp(0, 0) = rdet * (1 + Hr(1, 1) + mur(1, 1));
        tmp(0, 1) = -rdet * (Hr(0, 1) + mur(0, 1));
        tmp(1, 0) = -rdet * (Hr(1, 0) + mur(1, 0));
        tmp(1, 1) = rdet * (1 + Hr(0, 0) + mur(0, 0));

        Z1(0, 0) = tmp(0, 0) * (1 + Hl(0, 0) + mul(0, 0)) + tmp(0, 1) * (Hl(1, 0) + mul(1, 0));
        Z1(0, 1) = tmp(0, 0) * (Hl(0, 1) + mul(0, 1)) + tmp(0, 1) * (1 + Hl(1, 1) + mul(1, 1));
        Z1(1, 0) = tmp(1, 1) * (Hl(1, 0) + mul(1, 0)) + tmp(1, 0) * (1 + Hl(0, 0) + mul(0, 0));
        Z1(1, 1) = tmp(1, 1) * (1 + Hl(1, 1) + mul(1, 1)) + tmp(1, 0) * (Hl(0, 1) + mul(0, 1));

        // tmp2 = (bφ + inv(A_l)*L1l + inv(A_r)*L_1r)
        tmp2(0) = (bp(0) + Ali(0, 0) * tlcoeff1(pos, d, 0) + Ali(0, 1) * tlcoeff1(pos, d, 1) + Ari(0, 0) * tlcoeff1(zr, yr, xr, d, 0) + Ari(0, 1) * tlcoeff1(zr, yr, xr, d, 1));
        tmp2(1) = (bp(1) + Ali(1, 0) * tlcoeff1(pos, d, 0) + Ali(1, 1) * tlcoeff1(pos, d, 1) + Ari(1, 0) * tlcoeff1(zr, yr, xr, d, 0) + Ari(1, 1) * tlcoeff1(zr, yr, xr, d, 1));

        Z2(0) = tmp(0, 0) * tmp2(0) + tmp(0, 1) * tmp2(1);
        Z2(1) = tmp(1, 0) * tmp2(0) + tmp(1, 1) * tmp2(1);

        // tmp = D_r * (6I + mu_r + T5_r * H_r)
        tmp(0, 0) = Dr(0) * (6 + mur(0, 0) + T5(zr, yr, xr, d, 0) * Hr(0, 0));
        tmp(0, 1) = Dr(0) * (mur(0, 1) + T5(zr, yr, xr, d, 0) * Hr(0, 1));
        tmp(1, 0) = Dr(1) * (mur(1, 0) + T5(zr, yr, xr, d, 1) * Hr(1, 0));
        tmp(1, 1) = Dr(1) * (6 + mur(1, 1) + T5(zr, yr, xr, d, 1) * Hr(1, 1));

        M(0, 0) = -Dl(0) * (6 + mul(0, 0) + T5(pos, d, 0) * Hl(0, 0)) - tmp(0, 0) * Z1(0, 0) + tmp(0, 1) * Z1(1, 0);
        M(0, 1) = -Dl(0) * (mul(0, 1) + T5(pos, d, 0) * Hl(0, 1)) - tmp(0, 0) * Z1(0, 1) + tmp(0, 1) * Z1(1, 1);
        M(1, 0) = -Dl(1) * (mul(1, 0) + T5(pos, d, 1) * Hl(1, 0)) - tmp(1, 0) * Z1(0, 0) + tmp(1, 1) * Z1(1, 0);
        M(1, 1) = -Dl(1) * (6 + mul(1, 1) + T5(pos, d, 1) * Hl(1, 1)) - tmp(1, 0) * Z1(0, 1) + tmp(1, 1) * Z1(1, 1);

        vg(0) = bj(0);
        vg(0) -= Dtl(0) * Ali(0, 0) * tlcoeff1(pos, d, 0) - Dtl(1) * Ali(0, 1) * tlcoeff1(pos, d, 1);
        vg(0) += Dtr(0) * Ari(0, 0) * tlcoeff1(zr, yr, xr, d, 0) - Dtr(1) * Ari(0, 1) * tlcoeff1(zr, yr, xr, d, 1);
        vg(0) -= (tmp(0, 0) * Z2(0) + tmp(0, 1) * Z2(1));
        vg(1) = bj(1);
        vg(1) -= Dtl(0) * Ali(1, 0) * tlcoeff1(pos, d, 0) - Dtl(1) * Ali(1, 1) * tlcoeff1(pos, d, 1);
        vg(1) += Dtr(0) * Ari(1, 0) * tlcoeff1(zr, yr, xr, d, 0) - Dtr(1) * Ari(1, 1) * tlcoeff1(zr, yr, xr, d, 1);
        vg(1) -= (tmp(1, 0) * Z2(0) + tmp(1, 1) * Z2(1));

        rdet = 1 / (M(1, 1) * M(0, 0) - M(0, 1) * M(1, 0));

        c3(pos, d, 0) = rdet * (M(1, 1) * vg(0) - M(0, 1) * vg(1));
        c3(pos, d, 1) = rdet * (M(0, 0) * vg(0) - M(1, 0) * vg(1));

        c5(pos, d, 0) = Hl(0, 0) * c3(pos, d, 0) + Hl(0, 1) * c3(pos, d, 1);
        c5(pos, d, 1) = Hl(1, 0) * c3(pos, d, 0) + Hl(1, 1) * c3(pos, d, 1);

        c1(pos, d, 0) = -(Ali(0, 0) * tlcoeff1(pos, d, 0) + Ali(0, 1) * tlcoeff1(pos, d, 1));
        c1(pos, d, 0) += mur(0, 0) * c3(pos, d, 0) + mur(0, 1) * c3(pos, d, 1);
        c1(pos, d, 1) = -(Ali(1, 0) * tlcoeff1(pos, d, 0) + Ali(1, 1) * tlcoeff1(pos, d, 1));
        c1(pos, d, 1) += mur(1, 0) * c3(pos, d, 0) + mur(1, 1) * c3(pos, d, 1);

        int sign = 1;
        for (int g = 0; g < _ng; g++)
        {
            jnet(pos, d, g, RIGHT) = _node.HMesh(pos, d) * 0.5 * Dr(0) * (-1 * c1(pos, d, 0) + 3 * c2(pos, d, 0) - 6 * c3(pos, d, 0) + 10 * c4(pos, d, 0) - T5(pos, d, 0) * c5(pos, d, 0) + T6(pos, d, 0) * c6(pos, d, 0));
            surfFlux(pos, d, g, RIGHT) = avgFlux(pos, g) - (c1(imesh) + c3(imesh) + c3(imesh)) + c2(imesh) + c4(imesh) + c6(imesh);
        }
    }
}
