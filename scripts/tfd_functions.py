"""
Functions for computing bremsstrahlung cross-sections using the TFD method from Ref:

Based on B. Martinez et al, Phys. Plasmas 26, 103109 (2019), notes provided by B. Martinez, and Mathieu Lobert's thesis

Written by Michael Lavell, May 4, 2024

"""


import numpy as np
from scipy import constants
from scipy.special import lambertw

good_colors=["#db6d00","#006ddb","#920000","#52a736","#9B30FF"]

# electron rest mass
me_c2 = constants.electron_mass * constants.c**2
# Compton radius
r_c = constants.hbar / (constants.electron_mass * constants.c)
# normalized fine structure constant
alpha_f = constants.e**2 / (4.0 * constants.pi * constants.epsilon_0 * constants.hbar * constants.c)
# normalized electron radius
r_e = constants.e**2 / (4.0 * constants.pi * constants.epsilon_0 * constants.electron_mass * constants.c**2)
# unit conversion
m2_to_barn = 1.0e28

# (gamma - 1)mc^2 = E
# (1 - (v/c)^2)^-1/2 = 1 + E/(mc^2)
def calc_beam_velocity_EeV(EeV, m):
    return constants.c * (1.0 - (1.0 + constants.eV * EeV / (m * constants.c**2))**-2)**0.5


# Functions for characteristic path lengths

def lThomasFermi(Z):
    ltf = 4.0 * constants.pi * constants.epsilon_0 * constants.hbar**2 \
          / (constants.m_e * constants.e**2) * Z**(-1.0/3.0)
    # BM code has extra factor of 0.885
    return ltf / r_c * 0.885

def lDebye(TeV, ni, Zstar):
    ld = (constants.epsilon_0 * constants.e * TeV
          / (constants.e**2 * ni * Zstar
             * (Zstar + 1.0)))**0.5
    return ld / r_c

def lDebye2(Te, ne, Ti, ni, Zstar):
    ld = (constants.epsilon_0 * constants.e * Te / (constants.e**2 * ne) +
          constants.epsilon_0 * constants.e * Ti / (constants.e**2 * Zstar**2 * ni))**0.5
    return ld / r_c

def lInteratomic(ni):
    lia = (4.0 * np.pi * ni / 3.0)**(-1.0/3.0)
    return lia / r_c

def lReducedPotential(etaf, Lf, etad, Ld):
    # function from Bertrand.

    # ITFD from Eq.(19) (typo in T2 '-Lf**2')
    T1 = (etaf**2 / 2.) * ( (1. + Lf**2) * np.log(1. + Lf**2) - Lf**2) / (1. + Lf**2)
    T2 = (etad**2 / 2.) * ( (1. + Ld**2) * np.log(1. + Ld**2) - Ld**2) / (1. + Ld**2)
    T3 = etaf * etad * ( Ld**2 * np.log(1. + Lf**2) - Lf**2 * np.log(1. + Ld**2)) / (Ld**2 - Lf**2)
    ITFD = T1 + T2 + T3

    # notation
    a = ITFD

    # Eq.(22)
    Lr = np.sqrt(np.exp(lambertw(-np.exp(-(1. + 2. * a))) + 1. + 2.*a) - 1.)
    Lr = np.real(Lr)

    return Lr

def Elwert_factor(Z, k, g1):
    ec = 1.0

    if (k > 0.0) and (k < g1 - 1.0):
        g2 = g1 - k
        b1 = np.sqrt(1.0 - 1.0 / g1**2)
        b2 = np.sqrt(1.0 - 1.0 / g2**2)

        if (Z * alpha_f * (1.0 / b2 - 1.0 / b1)) < 100.0:
            ec = b1/b2 * (1.0 - np.exp(-2.0 * np.pi * Z * alpha_f / b1)) / \
                 (1.0 - np.exp(-2.0 * np.pi * Z * alpha_f / b2))
    return ec

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Differential for non-relativistic case
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Thomas-Fermi-Debye (screened) differential cross section
# Non-relativistic case

def g_func(eta, dp, dm, l):
    dpl_p1 = (dp * l)**2 + 1
    dml_p1 = (dm * l)**2 + 1
    val = 0.5 * eta**2 * (np.log(dpl_p1 / dml_p1) +
                          1.0 / dpl_p1 - 1.0 / dml_p1)
    return val

def dif_cs_tfd_nr(k, Z, Zstar, Ti, Te, ni, ne, g1, enable_ec=True):
    cs = 0.0

    if (k > 0.0) and (k < g1 - 1.0):

        g2 = g1 - k
        p1 = np.sqrt(g1**2 - 1.0)
        p2 = np.sqrt(g2**2 - 1.0)
        # b1 = np.sqrt(1.0 - 1.0 / g1**2)
        # b2 = np.sqrt(1.0 - 1.0 / g2**2)

        dp = p1 + p2 # maximum momentum (defined differently in paper..)
        dm = p1 - p2 # minimum momentum

        eta_tf = 1.0 - Zstar / Z
        eta_d = Zstar / Zax 

        ltf = lThomasFermi(Z)

        # result shown in paper has 1/b1**2 instead of 1/p1**2 in eq.8
        T_coef = 16.0 * (Z * r_e)**2 * alpha_f / (3.0 * k * p1**2)

        T_tf = g_func(eta_tf, dp, dm, ltf)

        T_d = 0.0
        T_coupling = 0.0

        if Zstar > 0:
            ri = min(lInteratomic(ne), lInteratomic(ni))
            # ld = max(lDebye(Ti, ni, Zstar), ri)
            ld = max(lDebye2(Te, ne, Ti, ni, Zstar), ri)
            T_d = g_func(eta_d, dp, dm, ld)

            T_coupling = eta_tf * eta_d / (ld**2 - ltf**2) * \
                         (ltf**2 * np.log(((dm * ld)**2 + 1.0) / ((dp * ld)**2 + 1.0)) +
                          ld**2 * np.log(((dp * ltf)**2 + 1.0) / ((dm * ltf)**2 + 1.0)))

        T_ec = 1.0
        if enable_ec:
            T_ec = Elwert_factor(Z, k, g1)

        cs = T_coef * (T_tf + T_d + T_coupling) * T_ec

    return cs


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Differential for moderately relativistic case
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# equation 16
def I1_screening(d, l, eta):
    T1 = l * d * (np.arctan(l * d) - np.arctan(l))
    T2 = - 0.5 * l**2 * (1.0 - d)**2 / (1.0 + l**2)
    T3 = 0.5 * np.log((1.0 + l**2) / (1.0 + (l * d)**2))
    return eta * (T1 + T2 + T3) # notes from Martinez had eta^2 here

# equation 17
def I2_screening(d, l, eta):
    T1 = 4.0 * (l * d)**3 * (np.arctan(d * l) - np.arctan(l))
    T2 = (1.0 + 3.0 * (l * d)**2) * np.log((1.0 + l**2) / (1.0 + (l * d)**2))
    T3 = (6.0 * l**4 * d**2) * np.log(d) / (1.0 + l**2)
    T4 = l**2 * (d - 1.0) * (d + 1.0 - 4.0 * (l * d)**2) / (1.0 + l**2)
    return 0.5 * eta * (T1 + T2 + T3 + T4)

# equation 19
def I_screening(ltf, ld, eta_tf, eta_d):
    ltf2_p1 = ltf**2 + 1.0
    ld2_p1 = ld**2 + 1.0

    Ttf = 0.5 * eta_tf**2 * (ltf2_p1 * np.log(ltf2_p1) - ltf**2) / ltf2_p1
    Td = 0.0
    Tc = 0.0
    if eta_d > 0.0:
        Td = 0.5 * eta_d**2 * (ld2_p1 * np.log(ld2_p1) - ld**2) / ld2_p1
        Tc = eta_tf * eta_d * (ld**2 * np.log(ltf2_p1) - ltf**2 * np.log(ld2_p1)) \
             / (ld**2 - ltf**2)

    return Ttf + Td + Tc

# equation 12
def dif_cs_tfd_mr(k, Z, Zstar, Ti, Te, ni, ne, g1, enable_ec=False):
    cs = 0.0

    if (k > 0.0) and (k < g1 - 1.0):

        eta_tf = 1.0 - Zstar / Z
        eta_d = Zstar / Z

        d = k / (2.0 * g1 * (g1 - k)) # momentum transfer
        ltf = lThomasFermi(Z)

        I1_tf = I1_screening(d, ltf, eta_tf)
        I2_tf = I2_screening(d, ltf, eta_tf)
        I1_d = 0.0
        I2_d = 0.0

        if Zstar > 0:
            ri = min(lInteratomic(ne), lInteratomic(ni))
            # ld = max(lDebye(Ti, ni, Zstar), ri)
            ld = max(lDebye2(Te, ne, Ti, ni, Zstar), ri)
            I1_d = I1_screening(d, ld, eta_d)
            I2_d = I2_screening(d, ld, eta_d)

        I1 = I1_tf + I1_d
        I2 = I2_tf + I2_d

        T_coef = 4.0 * (Z * r_e)**2 * alpha_f / k
        T1 = (1.0 + ((g1 - k) / g1)**2) * (I1 + 1.0)
        T2 = 2.0 / 3.0 * ((g1 - k)/ g1) * (I2 + 5.0 / 6.0)

        T_ec = 1.0
        if enable_ec:
            T_ec = Elwert_factor(Z, k, g1)

        cs = T_coef * (T1 - T2) * T_ec

    return cs

def dif_cs_r_mr(k, Z, Zstar, Ti, Te, ni, ne, g1, enable_ec=False):
    cs = 0.0

    if (k > 0.0) and (k < g1 - 1.0):

        eta_tf = 1.0 - Zstar / Z
        eta_d = Zstar / Z

        d = k / (2.0 * g1 * (g1 - k)) # momentum transfer
        ltf = lThomasFermi(Z)

        ld = 0.0
        if Zstar > 0:
            ri = min(lInteratomic(ne), lInteratomic(ni))
            # ld = max(lDebye(Ti, ni, Zstar), ri)
            ld = max(lDebye2(Te, ne, Ti, ni, Zstar), ri)

        lr = lReducedPotential(eta_tf, ltf, eta_d, ld)

        eta_r = 1.0
        I1 = I1_screening(d, lr, eta_r)
        I2 = I2_screening(d, lr, eta_r)

        T_coef = 4.0 * (Z * r_e) ** 2 * alpha_f / k
        T1 = (1.0 + ((g1 - k) / g1) ** 2) * (I1 + 1.0)
        T2 = 2.0 / 3.0 * ((g1 - k)/ g1) * (I2 + 5.0 / 6.0)

        cs = T_coef * (T1 - T2)

    return cs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Differential for ultra-relativistic case
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def zeta_Riemann(s, N):
    z = 0.0
    for i in range(N):
        z += 1.0 / (i+1)**s
    return z


# equation 24
def Coulomb_correction(Z):
    aZ2 = (alpha_f * Z)**2
    T1 = aZ2 / (1.0 + aZ2)
    T2 = 0.0
    for n in np.arange(1,4):
        z = zeta_Riemann(2.0 * n + 1.0, 100)
        T2 += (-aZ2)**n * (z - 1.0)
    return T1 * T2


# equation 23
def dif_cs_tfd_ur(k, Z, Zstar, Ti, Te, ni, ne, g1, enable_ec=False):
    cs = 0.0

    if (k > 0.0) and (k < g1 - 1.0):

        eta_tf = 1.0 - Zstar / Z
        eta_d = Zstar / Z

        d = k / (2.0 * g1 * (g1 - k)) # momentum transfer
        ltf = lThomasFermi(Z)

        I1_tf = I1_screening(d, ltf, eta_tf)
        I2_tf = I2_screening(d, ltf, eta_tf)
        I1_d = 0.0
        I2_d = 0.0

        if Zstar > 0:
            ri = min(lInteratomic(ne), lInteratomic(ni))
            # ld = max(lDebye(Ti, ni, Zstar), ri)
            ld = max(lDebye2(Te, ne, Ti, ni, Zstar), ri)
            I1_d = I1_screening(d, ld, eta_d)
            I2_d = I2_screening(d, ld, eta_d)

        I1 = I1_tf + I1_d
        I2 = I2_tf + I2_d
        fC = Coulomb_correction(Z)

        T_coef = 4.0 * (Z * r_e)**2 * alpha_f / k
        T1 = (1.0 + ((g1 - k) / g1)**2) * (I1 + 1.0 - fC)
        T2 = 2.0 / 3.0 * ((g1 - k)/ g1) * (I2 + 5.0 / 6.0 - fC)

        T_ec = 1.0
        if enable_ec:
            T_ec = Elwert_factor(Z, k, g1)

        cs = T_coef * (T1 - T2) * T_ec

    return cs


def dif_cs_r_ur(k, Z, Zstar, Ti, Te, ni, ne, g1, enable_ec=False):
    cs = 0.0

    if (k > 0.0) and (k < g1 - 1.0):

        eta_tf = 1.0 - Zstar / Z
        eta_d = Zstar / Z

        d = k / (2.0 * g1 * (g1 - k)) # momentum transfer
        ltf = lThomasFermi(Z)

        ld = 0.0
        if Zstar > 0:
            ri = min(lInteratomic(ne), lInteratomic(ni))
            # ld = max(lDebye(Ti, ni, Zstar), ri)
            ld = max(lDebye2(Te, ne, Ti, ni, Zstar), ri)

        lr = lReducedPotential(eta_tf, ltf, eta_d, ld)

        eta_r = 1.0
        I1 = I1_screening(d, lr, eta_r)
        I2 = I2_screening(d, lr, eta_r)
        fC = Coulomb_correction(Z)

        T_coef = 4.0 * (Z * r_e) ** 2 * alpha_f / k
        T1 = (1.0 + ((g1 - k) / g1) ** 2) * (I1 - fC + 1.0)
        T2 = 2.0 / 3.0 * ((g1 - k)/ g1) * (I2 - fC + 5.0 / 6.0)

        cs = T_coef * (T1 - T2)

    return cs


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Total cross-section
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# def get_dsdk(ks, Z, Zstar, T, ni, gamma):
#
#     npts = len(ks)
#     dsdk = np.zeros(npts)
#
#     if gamma > 1 and gamma <= 2.0:
#         for i in range(npts):
#             dsdk[i] = dif_cs_tfd_nr(ks[i], Z, Zstar, T, ni, gamma)
#     elif gamma > 2 and gamma <= 100.0:
#         for i in range(npts):
#             # dsdk[i] = dif_cs_tfd_mr(ks[i], Z, Zstar, T, ni, gamma)
#             dsdk[i] = dif_cs_r_mr(ks[i], Z, Zstar, T, ni, gamma)
#     elif gamma > 100:
#         for i in range(npts):
#             # dsdk[i] = dif_cs_tfd_ur(ks[i], Z, Zstar, T, ni, gamma)
#             dsdk[i] = dif_cs_r_ur(ks[i], Z, Zstar, T, ni, gamma)
#
#     return dsdk
#
#
# def total_cross_section(Z, Zstar, T, ni, gamma):
#
#     gm1 = gamma - 1.0
#
#     npts = 100
#     k_over_gm1_min = 1.0-10
#     k_over_gm1 = np.linspace(k_over_gm1_min, 1.0, npts)
#     ks = k_over_gm1 * gm1
#
#     dsdk = get_dsdk(ks, Z, Zstar, T, ni, gamma)
#
#     dk = np.zeros(npts)
#     dk[0] = ks[0]
#     for i in np.arange(1, npts):
#         dk[i] = ks[i] - ks[i - 1]
#
#     sigma_ttl = sum(dsdk * dk)
#
#     return sigma_ttl


def f_diff_cs(k, Z, Zstar, Ti, Te, ni, ne, gamma):

    dsdk = 0.0
    if 1.0 < gamma <= 2.0:
        dsdk = dif_cs_tfd_nr(k, Z, Zstar, Ti, Te, ni, ne, gamma)
    elif 2.0 < gamma <= 100.0:
        dsdk = dif_cs_r_mr(k, Z, Zstar, Ti, Te, ni, ne, gamma)
    elif gamma > 100.0:
        dsdk = dif_cs_r_ur(k, Z, Zstar, Ti, Te, ni, ne, gamma)

    return dsdk


def ttl_cs_gauss(Z, Zstar, Ti, Te, ni, ne, gamma, k_over_gm1_min=1.0e-7):

    f = lambda x: x * f_diff_cs(x, Z, Zstar, Ti, Te, ni, ne, gamma)

    gm1 = gamma - 1.0
    deg = 10

    a = np.log(k_over_gm1_min * gm1)
    b = np.log(gm1)

    # https://stackoverflow.com/questions/33457880/different-intervals-for-gauss-legendre-quadrature-in-numpy
    x, w = np.polynomial.legendre.leggauss(deg)
    ks = np.exp(0.5 * (x + 1.0) * (b - a) + a)

    integral = 0.0
    for j in range(deg):
        # kdsdk = ks[j] * f_diff_cs(ks[j], Z, Zstar, T, ni, gamma)
        integral += np.sum(w[j] * f(ks[j])) * 0.5 * (b - a)

    return integral


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Thomas-Fermi model for computing ionization state
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Thomas-Fermi ionization state Z_eff for any temperature, density, atomic number, atomic weight
# It comes from Table IV of the paper cited by B Martinez.

def ThomasFermi_Zeff(T, density, atn, atw):
    a1 = 0.003323
    a2 = 0.9718
    a3 = 9.26148e-5
    a4 = 3.10165
    b0 = -1.7630
    b1 = 1.43175
    b2 = 0.31546
    c1 = -0.366667
    c2 = 0.983333

    Temp = T # temperature in eV
    rho = density # density in g/cc
    Z = atn # atomic number
    A = atw # atomic weight, total nucleons
    alpha = 14.3139
    beta = 0.6624

    T0 = Temp / Z**(4.0/3.0)
    R = rho / (Z * A)
    Tf = T0 / (1.0 + T0)
    A = (a1 * T0**a2) + (a3 * T0**a4);
    B = -1.0*np.exp(b0 + b1*Tf + b2*Tf**7.0);
    C = c1*Tf + c2;
    Q1 = A * R**B;
    Q  = ( R**C + Q1**C )**(1.0 / C);
    x = alpha * (Q**beta);
    Zeff = Z * x / (1.0 + x + np.sqrt(1.0 + 2.0*x));
    return Zeff


def ThomasFermi_T0(rho, Z, A):
    # Calculation for T = 0
    alpha = 14.3139
    beta = 0.6624
    x = alpha*(rho/(Z*A))**beta
    f = x/(1. + x + np.sqrt(1.0+2*x))
    Zeff = f*Z
    # dZeffdrho = (Zeff/rho)*(beta/sqrt(1.+2.*x));

    return Zeff

# # Solid copper is density=8.9, atn=29, atw=63.546
# cu_rho = 8.9 # density in g/cc
# cu_Z = 29 # atomic number
# cu_A = 63.546 # atomic weight, total nucleons

# Temps = [0, 1, 100, 10e3]

# myTs = np.logspace(-5, 5, 100) #10^span(-5.,5.,601);
# myZs = ThomasFermi_Zeff(myTs, cu_rho, cu_Z, cu_A);

# fig,ax = plt.subplots(1,1,figsize=(5,4))
# ax.plot(myTs, myZs)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel("T (eV)")
# ax.set_ylabel(r"$Z_{\mathrm{eff}}$")
# ax.grid()









