"""
fit_omega.py
------------
Підбір параметрів Lambda-CDM з публічних даних BAO.
Вільні параметри: H0, Omega_m, Omega_Lambda, Omega_k = 1 - Omega_m - Omega_Lambda
Omega_r фіксовано (мале, нехтуємо).
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.optimize import curve_fit

BAO_DATA = [
    # z       D_V/r_s  sigma  тип
    (0.106,   3.047,   0.137, 'DV'),   # 6dFGS, Beutler+2011
    (0.150,   4.480,   0.168, 'DV'),   # SDSS MGS, Ross+2015
    (0.320,   8.467,   0.167, 'DV'),   # BOSS LOWZ, Alam+2017
    (0.570,  13.773,   0.134, 'DV'),   # BOSS CMASS, Alam+2017
    (0.440,  11.550,   0.400, 'DV'),   # WiggleZ, Blake+2011
    (0.600,  14.940,   0.420, 'DV'),   # WiggleZ, Blake+2011
    (0.730,  16.850,   0.600, 'DV'),   # WiggleZ, Blake+2011
    (1.520,  26.130,   0.580, 'DV'),   # eBOSS QSO, Ata+2018
]

R_S = 147.09    # Мпк (Planck 2018, фіксовано)
C   = 299792.458  # км/с

def E(z, Om, OL):
    Or = 9.4e-5
    Ok = 1.0 - Om - OL - Or
    return np.sqrt(abs(Or*(1+z)**4 + Om*(1+z)**3 + Ok*(1+z)**2 + OL))

def d_C(z, H0, Om, OL):
    integrand = lambda zp: 1.0 / E(zp, Om, OL)
    # epsabs та epsrel покращують точність чисельного інтегрування
    val, _ = quad(integrand, 0, z, epsabs=1e-8, epsrel=1e-8)
    return (C / H0) * val

def d_V(z, H0, Om, OL):
    dc  = d_C(z, H0, Om, OL)
    Hz  = H0 * E(z, Om, OL)
    return (z * dc**2 * C / Hz) ** (1.0/3.0)


def main():
    print("=" * 55)
    print("  Підбір параметрів ΛCDM через curve_fit (МНК)")
    print("=" * 55)

    # 1. Готуємо дані
    z_vals = np.array([x[0] for x in BAO_DATA])
    obs_vals = np.array([x[1] for x in BAO_DATA])
    sigmas = np.array([x[2] for x in BAO_DATA])

    # 2. Функція-обгортка для curve_fit
    def model_func(z, H0, Om, OL):
        # Повертаємо масив теоретичних значень
        return np.array([d_V(zi, H0, Om, OL) / R_S for zi in z])

    # 3. "Класика": автоматичний підбір МНК
    # curve_fit використовує алгоритм Левенберга-Марквардта
    p0 = [1, 1, 1]

    print("  Оптимізація параметрів за допомогою curve_fit...")
    popt, pcov = curve_fit(model_func, z_vals, obs_vals, p0=p0, sigma=sigmas, absolute_sigma=True)

    H0_fit, Om_fit, OL_fit = popt
    perr = np.sqrt(np.diag(pcov)) # Корінь з діагоналі = стандартні відхилення
    sH0, sOm, sOL = perr

    # Розрахунок Ω_k
    Or = 9.4e-5
    Ok_fit = 1.0 - Om_fit - OL_fit - Or
    sOk = np.sqrt(sOm**2 + sOL**2)

    # 4. Вивід результатів
    print(f"\n  Результати підбору:")
    print(f"  H₀  = {H0_fit:.2f} ± {sH0:.2f} км/с/Мпк")
    print(f"  Ω_m = {Om_fit:.4f} ± {sOm:.4f}")
    print(f"  Ω_Λ = {OL_fit:.4f} ± {sOL:.4f}")
    print(f"  Ω_k = {Ok_fit:.4f} ± {sOk:.4f}")

    # Оцінка якості (chi2)
    theory_vals = model_func(z_vals, *popt)
    chi2_val = np.sum(((obs_vals - theory_vals) / sigmas)**2)
    ndof = len(BAO_DATA) - 3
    print(f"\n  χ² = {chi2_val:.2f} | χ²/dof = {chi2_val/ndof:.2f}")

    if abs(Ok_fit) < 2*sOk:
        print(f"  ✓ Всесвіт плаский (Ω_k сумісне з 0)")

    print("=" * 55)

if __name__ == '__main__':
    main()
