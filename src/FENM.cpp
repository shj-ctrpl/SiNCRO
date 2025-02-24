#include "../include/FENM.h"

#define pos z, y, x
#define posr zr, yr, xr
#define im z, y, x, d, g

FENM::FENM(Node &node) : _node(node), flux(_node.Flux()), sFlux(_node.SFlux()), jnet(_node.Jnet()), keff(_node.Keff())
{
    this->nz = node.nz;
    this->ny = node.ny;
    this->nx = node.nx;
    this->ng = node.ng;
    this->nd = node.nd;

    tlcoeff0 = zdouble5(nz, ny, nx, nd, ng);
    tlcoeff1 = zdouble5(nz, ny, nx, nd, ng);
    tlcoeff2 = zdouble5(nz, ny, nx, nd, ng);

    mu_m = zdouble6(nz, ny, nx, nd, ng, ng);

    k51 = zdouble5(nz, ny, nx, nd, ng);
    k53 = zdouble5(nz, ny, nx, nd, ng);
    k60 = zdouble5(nz, ny, nx, nd, ng);
    k62 = zdouble5(nz, ny, nx, nd, ng);
    k64 = zdouble5(nz, ny, nx, nd, ng);

    T5 = zdouble5(nz, ny, nx, nd, ng);
    T6 = zdouble5(nz, ny, nx, nd, ng);

    diagD = zdouble5(nz, ny, nx, nd, ng);

    matM = zdouble5(ng, ng, nz, ny, nx);
    invM = zdouble5(ng, ng, nz, ny, nx);

    c1 = zdouble5(nz, ny, nx, nd, ng);
    c2 = zdouble5(nz, ny, nx, nd, ng);
    c3 = zdouble5(nz, ny, nx, nd, ng);
    c4 = zdouble5(nz, ny, nx, nd, ng);
    c5 = zdouble5(nz, ny, nx, nd, ng);
    c6 = zdouble5(nz, ny, nx, nd, ng);
}

FENM::~FENM()
{
}

void FENM::CalcConst()
{
    LOOP3D(nz, ny, nx)
    {
        for (int d = 0; d < nd; d++)
        {
            for (int g = 0; g < ng; g++)
            {
                // Pre-calculated kappa values
                double kp2 = _node.XsR(pos, g) * _node.HMesh(d) * _node.HMesh(d) / (4 * _node.XsD(pos, g));
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
                double mg3 = 7 * rkp4 * (15 * sinhkp - 15 * kp * coshkp + 6 * kp2 * sinhkp - kp3 * coshkp);
                double mg4 = -9 * rkp5 * (105 * sinhkp - 105 * kp * coshkp + 45 * kp2 * sinhkp - 10 * kp3 * coshkp + kp4 * sinhkp);

                // Common nominator and denominators for constants calculation
                double oddDenom = 1 / (sinhkp + mg1 + mg3);
                double evenDenom = 1 / (coshkp + mg0 + mg2 + mg4);

                T5(im) = (kp * coshkp + mg1 + 6 * mg3) * oddDenom;
                T6(im) = (kp * sinhkp + 3 * mg2 + 10 * mg4) * evenDenom;

                // Calculate the inner proudcts of the basis functions and 2nd order derivatives
                k51(im) = 2 * (kp * coshkp - sinhkp + 5 * mg3) * oddDenom;
                k53(im) = 2 * (kp * (15 + kp2) * coshkp - 3 * (5 + 2 * kp2) * sinhkp) * rkp2 * oddDenom;
                k60(im) = 2 * T6(im);
                k62(im) = 2 * (-3 * coshkp + (3 + kp2) * (sinhkp * rkp) + 7 * mg4) * evenDenom;
                k64(im) = 2 * (-5 * kp * (21 + 2 * kp2) * coshkp + (105 + 45 * kp2 + kp4) * sinhkp) * rkp3 * evenDenom;

                diagD(im) = 4 * _node.XsD(pos, g) / (_node.HMesh(d) * _node.HMesh(d));
            }
        }
    }
}

void FENM::CalcTLCoeffs0()
{
    LOOP3D(nz, ny, nx)
    {
        // Calculate 0th transverse leakage coefficient, L_x0 = (J_yr-J_yl)/h_y + (J_zr-J_zl)/h_z
        double jAvgZ, jAvgY, jAvgX;

        for (int g = 0; g < ng; g++)
        {
            jAvgZ = (jnet(z + 1, y, x, ZDIR, g) - jnet(z, y, x, ZDIR, g)) / _node.hz;
            jAvgY = (jnet(z, y + 1, x, YDIR, g) - jnet(z, y, x, YDIR, g)) / _node.hy;
            jAvgX = (jnet(z, y, x + 1, XDIR, g) - jnet(z, y, x, XDIR, g)) / _node.hx;
            tlcoeff0(pos, ZDIR, g) = jAvgY + jAvgX;
            tlcoeff0(pos, YDIR, g) = jAvgZ + jAvgX;
            tlcoeff0(pos, XDIR, g) = jAvgZ + jAvgY;
        }
    }
}

void FENM::CalcTLCoeffs12()
{
    LOOP3D(nz, ny, nx)
    {
        // Calculate the 1st and 2nd transverse leakage
        for (int d = 0; d < nd; d++)
        {
            int lIdx[3] = {z, y, x};
            lIdx[d] -= 1;

            int rIdx[3] = {z, y, x};
            rIdx[d] += 1;

            for (int g = 0; g < ng; g++)
            {
                double hC = _node.HMesh(d);
                double hL = hC, hR = hC;
                double tlC = tlcoeff0(im);
                double tlL = tlC, tlR = tlC;

                double rh = 0.5 / ((hC + hL) * (hC + hR) * (hC + hL + hR));

                // Boundary check
                if (lIdx[d] < 0)
                {
                    tlR = tlcoeff0(rIdx[ZDIR], rIdx[YDIR], rIdx[XDIR], d, g);
                    if (_node.Albedo(d, LEFT) < EPS)
                    {
                        tlcoeff1(im) = rh * hC * ((tlC - tlL) * (hC + 2 * hR) * (hC + hR) + (tlR - tlC) * (hC + 2 * hL) * (hC + hL));
                        tlcoeff2(im) = rh * (hC * hC) * ((tlL - tlC) * (hC + hR) + (tlR - tlC) * (hC + hL));
                    }
                    else
                    {
                        tlcoeff1(im) = 0.125 * (5.0 * tlC + tlR);
                        tlcoeff2(im) = 0.125 * (-3.0 * tlC + tlR);
                    }
                }
                else if (rIdx[d] >= _node._nMesh(d))
                {
                    tlL = tlcoeff0(lIdx[ZDIR], lIdx[YDIR], lIdx[XDIR], d, g);
                    if (_node.Albedo(d, RIGHT) < EPS)
                    {
                        tlcoeff1(im) = rh * hC * ((tlC - tlL) * (hC + 2 * hR) * (hC + hR) + (tlR - tlC) * (hC + 2 * hL) * (hC + hL));
                        tlcoeff2(im) = rh * (hC * hC) * ((tlL - tlC) * (hC + hR) + (tlR - tlC) * (hC + hL));
                    }
                    else
                    {
                        tlcoeff1(im) = -0.125 * (5.0 * tlC + tlL);
                        tlcoeff2(im) = 0.125 * (-3.0 * tlC + tlL);
                    }
                }
                else
                {
                    tlL = tlcoeff0(lIdx[ZDIR], lIdx[YDIR], lIdx[XDIR], d, g);
                    tlR = tlcoeff0(rIdx[ZDIR], rIdx[YDIR], rIdx[XDIR], d, g);
                    tlcoeff1(im) = rh * hC * ((tlC - tlL) * (hC + 2 * hR) * (hC + hR) + (tlR - tlC) * (hC + 2 * hL) * (hC + hL));
                    tlcoeff2(im) = rh * (hC * hC) * ((tlL - tlC) * (hC + hR) + (tlR - tlC) * (hC + hL));
                }
            }
        }
    }
}

void FENM::CalcMatrix()
{
    matM = 0.0;
    LOOP3D(nz, ny, nx)
    {
        // M = Σr - Σs - λF
        for (int from_g = 0; from_g < ng; from_g++)
        {
            matM(from_g, from_g, pos) = _node.XsR(pos, from_g);
            for (int to_g = 0; to_g < ng; to_g++)
            {
                matM(from_g, to_g, pos) -= (_node.XsS(pos, from_g, to_g) + _node.XsF(pos, from_g, to_g) / keff);
            }
        }

        double det = matM(0, 0, pos) * matM(1, 1, pos) - matM(1, 0, pos) * matM(0, 1, pos);

        if (std::abs(det) > EPS)
        {
            det = 1 / det;
            invM(0, 0, pos) = det * matM(1, 1, pos);
            invM(0, 1, pos) = -det * matM(0, 1, pos);
            invM(1, 0, pos) = -det * matM(1, 0, pos);
            invM(1, 1, pos) = det * matM(0, 0, pos);
        }
        else
        {
            invM(0, 0, pos) = 0;
            invM(0, 1, pos) = 0;
            invM(1, 0, pos) = 0;
            invM(1, 1, pos) = 0;
        }
        /*
                for (int d = 0; d < nd; d++)
                {
                    // H = mu_33 * invD * invk_53 * matM
                    zdouble2 H(2, 2);
                    H(0, 0) = mu33 * invD(pos, d, 0) * matM(0, 0, pos) / k53(pos, d, 0);
                    H(0, 1) = mu33 * invD(pos, d, 0) * matM(0, 1, pos) / k53(pos, d, 0);
                    H(1, 0) = mu33 * invD(pos, d, 1) * matM(1, 0, pos) / k53(pos, d, 1);
                    H(1, 1) = mu33 * invD(pos, d, 1) * matM(1, 1, pos) / k53(pos, d, 1);


                    // mu_m = (1/mu11) * invM * D (k_31 + k_51H)
                    mu_m(pos, d, 0, 0) = rmu11 * (invM(0, 0, pos) * diagD(pos, d, 0) * (k31 + k51(pos, d, 0) * H(0, 0)) +
                                                  (invM(0, 1, pos) * diagD(pos, d, 1) * k51(pos, d, 1) * H(1, 0)));

                    mu_m(pos, d, 0, 1) = rmu11 * (invM(0, 1, pos) * diagD(pos, d, 0) * (k31 + k51(pos, d, 1) * H(1, 1)) +
                                                  (invM(0, 0, pos) * diagD(pos, d, 0) * k51(pos, d, 0) * H(0, 1)));

                    mu_m(pos, d, 1, 0) = rmu11 * (invM(1, 0, pos) * diagD(pos, d, 1) * (k31 + k51(pos, d, 0) * H(0, 0)) +
                                                  (invM(1, 1, pos) * diagD(pos, d, 1) * k51(pos, d, 1) * H(1, 0)));

                    mu_m(pos, d, 1, 1) = rmu11 * (invM(1, 1, pos) * diagD(pos, d, 1) * (k31 + k51(pos, d, 1) * H(1, 1)) +
                                                  (invM(1, 0, pos) * diagD(pos, d, 0) * k51(pos, d, 0) * H(0, 1)));

                }
            */
    }
}

void FENM::CalcEvenCoeffs()
{
    LOOP3D(nz, ny, nx)
    {
        zdouble2 A(2, 2), T(2, 2), U(2, 2), B(2, 2);
        zdouble1 D(2), DI(2), K60(2), K62(2), a(2), b(2);

        for (int d = 0; d < nd; d++)
        {
            A = matM.slice(pos);
            D(0) = diagD(pos, d, 0);
            D(1) = diagD(pos, d, 1);
            DI(0) = 1 / diagD(pos, d, 0);
            DI(1) = 1 / diagD(pos, d, 1);
            K60(0) = k60(pos, d, 0);
            K60(1) = k60(pos, d, 1);
            K62(0) = k62(pos, d, 0);
            K62(1) = k62(pos, d, 1);

            // T = μ₄₄K₆₄⁻¹D⁻¹A
            T(0, 0) = mu44 * DI(0) * A(0, 0) / k64(pos, d, 0);
            T(0, 1) = mu44 * DI(0) * A(0, 1) / k64(pos, d, 0);
            T(1, 0) = mu44 * DI(1) * A(1, 0) / k64(pos, d, 1);
            T(1, 1) = mu44 * DI(1) * A(1, 1) / k64(pos, d, 1);

            // U = μ₂₂AK₂₀⁻¹ = (2/5) * A / 6 = 0.0666666667 * A
            U = mu22 * rk20 * A;

            // B = DK₄₂ + DK₆₂μ₄₄K₆₄⁻¹D⁻¹A + μ₂₂AK₂₀⁻¹(K₄₀ + K₆₀μ₄₄K₆₄⁻¹D⁻¹A) → DK₄₂ + DK₆₂T + U(K₄₀ + K₆₀T)
            B(0, 0) = D(0) * k42 + D(0) * K62(0) * T(0, 0) + U(0, 0) * k40 + (U(0, 0) * T(0, 0) * K60(0) + U(0, 1) * T(1, 0) * K60(1));
            B(0, 1) = D(0) * k42 + D(0) * K62(0) * T(0, 1) + U(0, 1) * k40 + (U(0, 0) * T(0, 1) * K60(0) + U(0, 1) * T(1, 1) * K60(1));
            B(1, 0) = D(1) * k42 + D(1) * K62(1) * T(1, 0) + U(1, 0) * k40 + (U(1, 0) * T(0, 0) * K60(0) + U(1, 1) * T(1, 0) * K60(1));
            B(1, 1) = D(1) * k42 + D(1) * K62(1) * T(1, 1) + U(1, 1) * k40 + (U(1, 0) * T(0, 1) * K60(0) + U(1, 1) * T(1, 1) * K60(1));

            // b = -μ₂₂l₂ - 2μ₂₂AK₂₀⁻¹D⁻¹(l₀+AΦ) → -μ₂₂l₂ - 2UD⁻¹(l₀+AΦ)
            // First, calculate a = (l₀+AΦ)
            a(0) = tlcoeff0(pos, d, 0) + A(0, 0) * flux(pos, 0) + A(0, 1) * flux(pos, 1);
            a(1) = tlcoeff0(pos, d, 1) + A(1, 0) * flux(pos, 0) + A(1, 1) * flux(pos, 1);

            // Next, add -μ₂₂l₂ - 2UD⁻¹a
            b(0) = -mu22 * tlcoeff2(pos, d, 0) - 2 * (U(0, 0) * DI(0) * a(0) + U(0, 1) * DI(1) * a(1));
            b(1) = -mu22 * tlcoeff2(pos, d, 1) - 2 * (U(1, 0) * DI(0) * a(0) + U(1, 1) * DI(1) * a(1));

            double rdet = B(0, 0) * B(1, 1) - B(1, 0) * B(0, 1);
            if (std::abs(rdet) > EPS)
            {
                rdet = 1 / rdet;
                c4(pos, d, 0) = rdet * (B(1, 1) * b(0) - B(0, 1) * b(1));
                c4(pos, d, 1) = rdet * (B(0, 0) * b(1) - B(1, 0) * b(0));
            }
            else
            {
                c4(pos, d, 0) = 0;
                c4(pos, d, 1) = 0;
            }

            // c₆ = μ₄₄K₆₄⁻¹D⁻¹Ac₄ → Tc₄
            c6(pos, d, 0) = T(0, 0) * c4(pos, d, 0) + T(0, 1) * c4(pos, d, 1);
            c6(pos, d, 1) = T(1, 0) * c4(pos, d, 0) + T(1, 1) * c4(pos, d, 1);

            // c₂ = K₂₀⁻¹(2D⁻¹(l₀+AΦ) - K₄₀c₄ - K₆₀μ₄₄K₆₄⁻¹D⁻¹Ac₄) → K₂₀⁻¹(2 * D⁻¹a - K₄₀c₄ - K₆₀c₆)
            c2(pos, d, 0) = rk20 * (2 * DI(0) * a(0) - k40 * c4(pos, d, 0) - K60(0) * c6(pos, d, 0));
            c2(pos, d, 1) = rk20 * (2 * DI(1) * a(1) - k40 * c4(pos, d, 1) - K60(1) * c6(pos, d, 1));
        }
    }
}

void FENM::CalcOddCoeffsOneNode(int z, int y, int x, int d, int face)
{
    int jIdx[3] = {z, y, x};
    jIdx[d] += face;

    zdouble2 invA = invM.slice(pos);
    zdouble2 B(2, 2), invB(2, 2);
    zdouble1 a12(2), a13(2), a23(2), a31(2), a32(2), a33(2), bj(2), bp(2), b(2), Dt(2);
    zdouble2 a11 = mu11 * matM.slice(pos); // a₁₁ = μ₁₁A
    zdouble2 a22 = mu33 * matM.slice(pos); // a₂₂ = μ₃₃A

    for (int g = 0; g < ng; g++)
    {
        Dt(g) = 0.5 * _node.HMesh(d) * diagD(im);                          // D tilde
        a12(g) = k31 * diagD(im);                                          // a₁₂ = DK₃₁
        a13(g) = k51(im) * diagD(im);                                      // a₁₃ = DK₅₁
        a23(g) = k53(im) * diagD(im);                                      // a₂₃ = DK₅₃
        a31(g) = Dt(g) + _node.Albedo(d, face);                            // a₃₁ = (Dt + α)
        a32(g) = 6.0 * Dt(g) + _node.Albedo(d, face);                      // a₃₂ = (6Dt + α)
        a33(g) = T5(im) * Dt(g) + _node.Albedo(d, face);                   // a₃₃ = (T₅Dt + α)
        bj(g) = -Dt(g) * (3.0 * c2(im) + 10.0 * c4(im) + T6(im) * c6(im)); // bⱼ = Dt(3c₂ + 10C₄ + T₆c₆)
        bp(g) = (flux(pos, g) + c2(im) + c4(im) + c6(im));                 // bᵩ = Φ + c₂ + C₄ + c₆
    }

    // Base equation: {DK₃₁ + DK₅₁H + μ₁₁A(Dt + α)⁻¹(6Dt + α) + μ₁₁A(Dt + α)⁻¹(T₅Dt + α)H}c₃ = μ₁₁A(Dt + α)⁻¹(bⱼ - αbᵩ) + μ₁₁l₁
    // 1. H = K₅₃⁻¹D⁻¹μ₃₃A = a₂₂ / a₂₃ → a₂₂
    a22(0, 0) = a22(0, 0) / a23(0);
    a22(0, 1) = a22(0, 1) / a23(0);
    a22(1, 0) = a22(1, 0) / a23(1);
    a22(1, 1) = a22(1, 1) / a23(1);

    // 2. μ₁₁A(Dt + α)⁻¹ = a₁₁a₃₁⁻¹ → a₁₁
    a11(0, 0) = a11(0, 0) / a31(0);
    a11(0, 1) = a11(0, 1) / a31(1);
    a11(1, 0) = a11(1, 0) / a31(0);
    a11(1, 1) = a11(1, 1) / a31(1);

    // Modified equation: {DK₃₁ + DK₅₁H + μ₁₁A(Dt + α)⁻¹(6Dt + α) + μ₁₁A(Dt + α)⁻¹(T₅Dt + α)H}c₃ → {a₁₂ + a₁₃a₂₂ + a₁₁a₃₂ + a₁₁a₃₃H}c₃
    B(0, 0) = a12(0) + a13(0) * a22(0, 0) + a11(0, 0) * a32(0) + a11(0, 0) * a33(0) * a22(0, 0) + a11(0, 1) * a33(1) * a22(1, 0);
    B(0, 1) = a13(0) * a22(0, 1) + a11(0, 1) * a32(1) + a11(0, 0) * a33(0) * a22(0, 1) + a11(0, 1) * a33(1) * a22(1, 1);
    B(1, 0) = a13(1) * a22(1, 0) + a11(1, 0) * a32(0) + a11(1, 0) * a33(0) * a22(0, 0) + a11(1, 1) * a33(1) * a22(1, 0);
    B(1, 1) = a12(1) + a13(1) * a22(1, 1) + a11(1, 1) * a32(1) + a11(1, 0) * a33(0) * a22(0, 1) + a11(1, 1) * a33(1) * a22(1, 1);

    // Calculate μ₁₁A(Dt + α)⁻¹(bⱼ - αbᵩ) + μ₁₁l₁ → a₁₁(bⱼ - αbᵩ) + μ₁₁l₁
    if (face == LEFT)
    {
        bj = -1.0 * bj;
    }
    b(0) = a11(0, 0) * (bj(0) - _node.albedo(d, face) * bp(0)) + a11(0, 1) * (bj(1) - _node.albedo(d, face) * bp(1)) + mu11 * tlcoeff1(pos, d, 0);
    b(1) = a11(1, 0) * (bj(0) - _node.albedo(d, face) * bp(0)) + a11(1, 1) * (bj(1) - _node.albedo(d, face) * bp(1)) + mu11 * tlcoeff1(pos, d, 1);

    double rdet = 1 / (B(0, 0) * B(1, 1) - B(0, 1) * B(1, 0));
    invB(0, 0) = rdet * B(1, 1);
    invB(0, 1) = -rdet * B(0, 1);
    invB(1, 0) = -rdet * B(1, 0);
    invB(1, 1) = rdet * B(0, 0);

    // Solve for c₃
    c3(pos, d, 0) = invB(0, 0) * b(0) + invB(0, 1) * b(1);
    c3(pos, d, 1) = invB(1, 0) * b(0) + invB(1, 1) * b(1);

    // μ₃₃Ac₃ = DK₅₃c₅ → c₅ = D⁻¹K₅₃⁻¹μ₃₃Ac₃ = a₂₂c₃
    c5(pos, d, 0) = a22(0, 0) * c3(pos, d, 0) + a22(0, 1) * c3(pos, d, 1);
    c5(pos, d, 1) = a22(1, 0) * c3(pos, d, 0) + a22(1, 1) * c3(pos, d, 1);

    // μ₁₁Ac₁ = (DK₃₁c₃ + DK₅₁c₅) - μ₁₁l₁ → c₁ = μ₁₁⁻¹A⁻¹(a₁₂c₃ + a₁₃c₅) - A⁻¹l₁
    c1(pos, d, 0) = rmu11 * (invA(0, 0) * (a12(0) * c3(pos, d, 0) + a13(0) * c5(pos, d, 0) - mu11 * tlcoeff1(pos, d, 0)) + invA(0, 1) * (a12(1) * c3(pos, d, 1) + a13(1) * c5(pos, d, 1) - mu11 * tlcoeff1(pos, d, 1)));
    c1(pos, d, 1) = rmu11 * (invA(1, 0) * (a12(0) * c3(pos, d, 0) + a13(0) * c5(pos, d, 0) - mu11 * tlcoeff1(pos, d, 0)) + invA(1, 1) * (a12(1) * c3(pos, d, 1) + a13(1) * c5(pos, d, 1) - mu11 * tlcoeff1(pos, d, 1)));

    for (int g = 0; g < ng; g++)
    {
        jnet(jIdx[0], jIdx[1], jIdx[2], d, g) = -Dt(g) * (c1(im) + 6.0 * c3(im) + T5(im) * c5(im)) + bj(g);
        sFlux(jIdx[0], jIdx[1], jIdx[2], d, g) = (c1(im) + c3(im) + c5(im)) + bp(g);
    }
}

void FENM::CalcOddCoeffsTwoNode(int z, int y, int x, int d)
{
    double rdet;

    zdouble2 Al(2, 2), Ali(2, 2), Hl(2, 2), mul(2, 2), Ar(2, 2), Ari(2, 2), Hr(2, 2), mur(2, 2);
    zdouble1 Dl(2), Dtl(2), Dr(2), Dtr(2), bp(2), bj(2);

    zdouble2 tmp(2, 2), Z1(2, 2), M(2, 2);
    zdouble1 tmp2(2), Z2(2), v1(2);

    int zr = z + (d == ZDIR);
    int yr = y + (d == YDIR);
    int xr = x + (d == XDIR);

    Al = matM.slice(pos);
    Ali = invM.slice(pos);
    Ar = matM.slice(posr);
    Ari = invM.slice(posr);

    for (int g = 0; g < ng; g++)
    {
        Dl(g) = diagD(pos, d, g);
        Dtl(g) = 0.5 * _node.HMesh(d) * Dl(g);
        Dr(g) = diagD(posr, d, g);
        Dtr(g) = 0.5 * _node.HMesh(d) * Dr(g);
        bp(g) = (flux(posr, g) + c2(posr, d, g) + c4(posr, d, g) + c6(posr, d, g)) - (flux(pos, g) + c2(pos, d, g) + c4(pos, d, g) + c6(pos, d, g));
        bj(g) = -(Dtl(g) * (3 * c2(pos, d, g) + 10 * c4(pos, d, g) + T6(pos, d, g) * c6(pos, d, g)) + Dtr(g) * (3 * c2(posr, d, g) + 10 * c4(posr, d, g) + T6(posr, d, g) * c6(posr, d, g)));
    }

    // H = K₅₃⁻¹D⁻¹μ₃₃A
    Hl(0, 0) = mu33 * Al(0, 0) / (Dl(0) * k53(pos, d, 0));
    Hl(0, 1) = mu33 * Al(0, 1) / (Dl(0) * k53(pos, d, 0));
    Hl(1, 0) = mu33 * Al(1, 0) / (Dl(1) * k53(pos, d, 1));
    Hl(1, 1) = mu33 * Al(1, 1) / (Dl(1) * k53(pos, d, 1));
    Hr(0, 0) = mu33 * Ar(0, 0) / (Dr(0) * k53(posr, d, 0));
    Hr(0, 1) = mu33 * Ar(0, 1) / (Dr(0) * k53(posr, d, 0));
    Hr(1, 0) = mu33 * Ar(1, 0) / (Dr(1) * k53(posr, d, 1));
    Hr(1, 1) = mu33 * Ar(1, 1) / (Dr(1) * k53(posr, d, 1));

    // μ = -μ₁₁⁻¹A⁻¹D(K₃₁+K₅₁H)
    mul(0, 0) = -rmu11 * (Ali(0, 0) * Dl(0) * (k31 + Hl(0, 0) * k51(pos, d, 0)) + Ali(0, 1) * Dl(1) * Hl(1, 0) * k51(pos, d, 1));
    mul(0, 1) = -rmu11 * (Ali(0, 1) * Dl(1) * (k31 + Hl(1, 1) * k51(pos, d, 1)) + Ali(0, 0) * Dl(0) * Hl(0, 1) * k51(pos, d, 0));
    mul(1, 0) = -rmu11 * (Ali(1, 0) * Dl(0) * (k31 + Hl(0, 0) * k51(pos, d, 0)) + Ali(1, 1) * Dl(1) * Hl(1, 0) * k51(pos, d, 1));
    mul(1, 1) = -rmu11 * (Ali(1, 1) * Dl(1) * (k31 + Hl(1, 1) * k51(pos, d, 1)) + Ali(1, 0) * Dl(0) * Hl(0, 1) * k51(pos, d, 0));

    mur(0, 0) = -rmu11 * (Ari(0, 0) * Dr(0) * (k31 + Hr(0, 0) * k51(posr, d, 0)) + Ari(0, 1) * Dr(1) * Hr(1, 0) * k51(posr, d, 1));
    mur(0, 1) = -rmu11 * (Ari(0, 1) * Dr(1) * (k31 + Hr(1, 1) * k51(posr, d, 1)) + Ari(0, 0) * Dr(0) * Hr(0, 1) * k51(posr, d, 0));
    mur(1, 0) = -rmu11 * (Ari(1, 0) * Dr(0) * (k31 + Hr(0, 0) * k51(posr, d, 0)) + Ari(1, 1) * Dr(1) * Hr(1, 0) * k51(posr, d, 1));
    mur(1, 1) = -rmu11 * (Ari(1, 1) * Dr(1) * (k31 + Hr(1, 1) * k51(posr, d, 1)) + Ari(1, 0) * Dr(0) * Hr(0, 1) * k51(posr, d, 0));

    rdet = 1 / ((1 + Hr(0, 0) + mur(0, 0)) * (1 + Hr(1, 1) + mur(1, 1)) - (Hr(0, 1) + mur(0, 1)) * (Hr(1, 0) + mur(1, 0)));

    // tmp = (I + Hᵣ + μᵣ)⁻¹
    tmp(0, 0) = rdet * (1 + Hr(1, 1) + mur(1, 1));
    tmp(0, 1) = -rdet * (Hr(0, 1) + mur(0, 1));
    tmp(1, 0) = -rdet * (Hr(1, 0) + mur(1, 0));
    tmp(1, 1) = rdet * (1 + Hr(0, 0) + mur(0, 0));

    // ζ₁ = (I + Hᵣ + μᵣ)⁻¹(I + Hₗ + μₗ) = tmp * (I + Hₗ + μₗ)
    Z1(0, 0) = tmp(0, 0) * (1 + Hl(0, 0) + mul(0, 0)) + tmp(0, 1) * (Hl(1, 0) + mul(1, 0));
    Z1(0, 1) = tmp(0, 0) * (Hl(0, 1) + mul(0, 1)) + tmp(0, 1) * (1 + Hl(1, 1) + mul(1, 1));
    Z1(1, 0) = tmp(1, 0) * (1 + Hl(0, 0) + mul(0, 0)) + tmp(1, 1) * (Hl(1, 0) + mul(1, 0));
    Z1(1, 1) = tmp(1, 0) * (Hl(0, 1) + mul(0, 1)) + tmp(1, 1) * (1 + Hl(1, 1) + mul(1, 1));

    // tmp2 = (bᵩ + Aₗ⁻¹l₁ₗ + Aᵣ⁻¹l₁ᵣ)
    tmp2(0) = (bp(0) + Ali(0, 0) * tlcoeff1(pos, d, 0) + Ali(0, 1) * tlcoeff1(pos, d, 1) + Ari(0, 0) * tlcoeff1(posr, d, 0) + Ari(0, 1) * tlcoeff1(posr, d, 1));
    tmp2(1) = (bp(1) + Ali(1, 0) * tlcoeff1(pos, d, 0) + Ali(1, 1) * tlcoeff1(pos, d, 1) + Ari(1, 0) * tlcoeff1(posr, d, 0) + Ari(1, 1) * tlcoeff1(posr, d, 1));

    // ζ₂ = (I + Hᵣ + μᵣ)⁻¹(bᵩ + Aₗ⁻¹l₁ₗ + Aᵣ⁻¹l₁ᵣ) = tmp * tmp2
    Z2(0) = tmp(0, 0) * tmp2(0) + tmp(0, 1) * tmp2(1);
    Z2(1) = tmp(1, 0) * tmp2(0) + tmp(1, 1) * tmp2(1);

    // tmp = Dtᵣ(6I + μᵣ + T₅ᵣHᵣ)
    tmp(0, 0) = Dtr(0) * (6 + mur(0, 0) + T5(posr, d, 0) * Hr(0, 0));
    tmp(0, 1) = Dtr(0) * (mur(0, 1) + T5(posr, d, 0) * Hr(0, 1));
    tmp(1, 0) = Dtr(1) * (mur(1, 0) + T5(posr, d, 1) * Hr(1, 0));
    tmp(1, 1) = Dtr(1) * (6 + mur(1, 1) + T5(posr, d, 1) * Hr(1, 1));

    // M = -{Dtₗ(6I + μₗ + T₅ₗHₗ) + Dtᵣ(6I + μᵣ + T₅ᵣHᵣ)ζ₁}
    M(0, 0) = -Dtl(0) * (6 + mul(0, 0) + T5(pos, d, 0) * Hl(0, 0)) - tmp(0, 0) * Z1(0, 0) - tmp(0, 1) * Z1(1, 0);
    M(0, 1) = -Dtl(0) * (mul(0, 1) + T5(pos, d, 0) * Hl(0, 1)) - tmp(0, 0) * Z1(0, 1) - tmp(0, 1) * Z1(1, 1);
    M(1, 0) = -Dtl(1) * (mul(1, 0) + T5(pos, d, 1) * Hl(1, 0)) - tmp(1, 0) * Z1(0, 0) - tmp(1, 1) * Z1(1, 0);
    M(1, 1) = -Dtl(1) * (6 + mul(1, 1) + T5(pos, d, 1) * Hl(1, 1)) - tmp(1, 0) * Z1(0, 1) - tmp(1, 1) * Z1(1, 1);

    // v₁ = bⱼ - DtₗAₗ⁻¹l₁ₗ + DtᵣAᵣ⁻¹l₁ᵣ - Dtᵣ(6I + μᵣ + T₅ᵣHᵣ)ζ₂ = bⱼ - DtₗAₗ⁻¹l₁ₗ + DtᵣAᵣ⁻¹l₁ᵣ - tmp * ζ₂
    v1(0) = bj(0);
    v1(0) -= (Dtl(0) * Ali(0, 0) * tlcoeff1(pos, d, 0) + Dtl(0) * Ali(0, 1) * tlcoeff1(pos, d, 1));
    v1(0) += (Dtr(0) * Ari(0, 0) * tlcoeff1(posr, d, 0) + Dtr(0) * Ari(0, 1) * tlcoeff1(posr, d, 1));
    v1(0) -= (tmp(0, 0) * Z2(0) + tmp(0, 1) * Z2(1));
    v1(1) = bj(1);
    v1(1) -= (Dtl(1) * Ali(1, 0) * tlcoeff1(pos, d, 0) + Dtl(1) * Ali(1, 1) * tlcoeff1(pos, d, 1));
    v1(1) += (Dtr(1) * Ari(1, 0) * tlcoeff1(posr, d, 0) + Dtr(1) * Ari(1, 1) * tlcoeff1(posr, d, 1));
    v1(1) -= (tmp(1, 0) * Z2(0) + tmp(1, 1) * Z2(1));

    // c₃ᵣ = ζ₂ - ζ₁M⁻¹v₁
    rdet = 1 / (M(1, 1) * M(0, 0) - M(0, 1) * M(1, 0));
    c3(posr, d, 0) = Z2(0) + rdet * (v1(0) * (M(1, 1) * Z1(0, 0) - M(1, 0) * Z1(0, 1)) + v1(1) * (M(0, 0) * Z1(0, 1) - M(0, 1) * Z1(0, 0)));
    c3(posr, d, 1) = Z2(1) + rdet * (v1(0) * (M(1, 1) * Z1(1, 0) - M(1, 0) * Z1(1, 1)) + v1(1) * (M(0, 0) * Z1(1, 1) - M(0, 1) * Z1(1, 0)));

    // c₅ᵣ = Hᵣc₃ᵣ
    c5(posr, d, 0) = Hr(0, 0) * c3(posr, d, 0) + Hr(0, 1) * c3(posr, d, 1);
    c5(posr, d, 1) = Hr(1, 0) * c3(posr, d, 0) + Hr(1, 1) * c3(posr, d, 1);

    // c₁ᵣ = -Aᵣ⁻¹l₁ᵣ + μᵣc₃ᵣ
    c1(posr, d, 0) = -Ari(0, 0) * tlcoeff1(posr, d, 0) - Ari(0, 1) * tlcoeff1(posr, d, 1);
    c1(posr, d, 0) += mur(0, 0) * c3(posr, d, 0) + mur(0, 1) * c3(posr, d, 1);
    c1(posr, d, 1) = -Ari(1, 0) * tlcoeff1(posr, d, 0) - Ari(1, 1) * tlcoeff1(posr, d, 1);
    c1(posr, d, 1) += mur(1, 0) * c3(posr, d, 0) + mur(1, 1) * c3(posr, d, 1);

    for (int g = 0; g < ng; g++)
    {
        jnet(posr, d, g) = -Dtr(g) * (c1(posr, d, g) - 3.0 * c2(posr, d, g) + 6.0 * c3(posr, d, g) - 10.0 * c4(posr, d, g) + T5(posr, d, g) * c5(posr, d, g) - T6(posr, d, g) * c6(posr, d, g));
        sFlux(posr, d, g) = flux(posr, g) - (c1(posr, d, g) + c3(posr, d, g) + c5(posr, d, g)) + c2(posr, d, g) + c4(posr, d, g) + c6(posr, d, g);
    }
}

void FENM::Drive()
{
    CalcConst();
    CalcTLCoeffs0();
    CalcTLCoeffs12();
    CalcMatrix();
    CalcEvenCoeffs();

    // Boundary calculation
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            CalcOddCoeffsOneNode(z, y, 0, XDIR, LEFT);
            CalcOddCoeffsOneNode(z, y, nx - 1, XDIR, RIGHT);
        }
    }
    for (int z = 0; z < nz; z++)
    {
        for (int x = 0; x < nx; x++)
        {
            CalcOddCoeffsOneNode(z, 0, x, YDIR, LEFT);
            CalcOddCoeffsOneNode(z, ny - 1, x, YDIR, RIGHT);
        }
    }
    for (int y = 0; y < ny; y++)
    {
        for (int x = 0; x < nx; x++)
        {
            CalcOddCoeffsOneNode(0, y, x, ZDIR, LEFT);
            CalcOddCoeffsOneNode(nz - 1, y, x, ZDIR, RIGHT);
        }
    }

    // Boundary calculation
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx - 1; x++)
            {
                CalcOddCoeffsTwoNode(z, y, x, XDIR);
            }
        }
    }
    for (int z = 0; z < nz; z++)
    {
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny - 1; y++)
            {
                CalcOddCoeffsTwoNode(z, y, x, YDIR);
            }
        }
    }
    for (int y = 0; y < ny; y++)
    {
        for (int x = 0; x < nx; x++)
        {
            for (int z = 0; z < nz - 1; z++)
            {
                CalcOddCoeffsTwoNode(z, y, x, ZDIR);
            }
        }
    }
}