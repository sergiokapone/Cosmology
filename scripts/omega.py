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

def chi2(params):
    H0, Om, OL = params
    if H0 < 50 or H0 > 90: return 1e10
    if Om < 0.05 or Om > 0.8: return 1e10
    if OL < 0 or OL > 1.2: return 1e10
    if abs(1 - Om - OL) > 0.5: return 1e10
    total = 0.0
    for (z, obs, sigma, kind) in BAO_DATA:
        if kind == 'DV':
            theory = d_V(z, H0, Om, OL) / R_S
        else:
            theory = d_C(z, H0, Om, OL) / (1+z) / R_S
        total += ((obs - theory) / sigma) ** 2
    return total


def callback_func(params):
    H0, Om, OL = params
    print(f"  Поточні спроби: H0={H0:.2f}, Om={Om:.4f}, OL={OL:.4f}")


def main():
    print("=" * 55)
    print("  Підбір параметрів ΛCDM з даних BAO")
    print("=" * 55)
    print(f"  Дані: {len(BAO_DATA)} точок BAO")
    print(f"  r_s = {R_S} Мпк (Planck 2018, фіксовано)")
    print(f"  Вільні: H0, Ω_m, Ω_Λ  |  Ω_k = 1 - Ω_m - Ω_Λ - Ω_r")
    print("-" * 55)

   # Оптимізація з правильними параметрами точності
# Оптимізація з відображенням ітерацій
    x0 = [67.4, 0.315, 0.685]
    print("  Початок оптимізації...")
    result = minimize(chi2, x0, method='BFGS',
                      callback=callback_func,
                      options={'gtol': 1e-6, 'maxiter': 10000})

    H0_fit, Om_fit, OL_fit = result.x
    Or = 9.4e-5
    Ok_fit = 1.0 - Om_fit - OL_fit - Or

    try:
        cov = np.array(result.hess_inv)
        sH0, sOm, sOL = np.sqrt(np.diag(cov))
        sOk = np.sqrt(sOm**2 + sOL**2)
    except Exception:
        sH0 = sOm = sOL = sOk = float('nan')

    chi2_min = result.fun
    ndof = len(BAO_DATA) - 3
    chi2_red = chi2_min / ndof

    print(f"\n  Результати підбору:")
    print(f"  H₀  = {H0_fit:.2f} ± {sH0:.2f} км/с/Мпк")
    print(f"  Ω_m = {Om_fit:.4f} ± {sOm:.4f}")
    print(f"  Ω_Λ = {OL_fit:.4f} ± {sOL:.4f}")
    print(f"  Ω_k = {Ok_fit:.4f} ± {sOk:.4f}")
    print(f"  Ω_r = {Or:.2e}  (фіксовано)")
    print()
    print(f"  χ² = {chi2_min:.2f}  |  dof = {ndof}  |  χ²/dof = {chi2_red:.2f}")
    print()

    if abs(Ok_fit) < 2*sOk:
        print(f"  ✓ Ω_k сумісне з нулем (|Ω_k| < 2σ) → Всесвіт плаский")
    else:
        print(f"  ✗ Ω_k відрізняється від нуля на {abs(Ok_fit)/sOk:.1f}σ")

    print()
    print("  Порівняння з Planck 2018:")
    print("  H₀  = 67.36 ± 0.54")
    print("  Ω_m = 0.3153 ± 0.0073")
    print("  Ω_Λ = 0.6847 ± 0.0073")
    print("  Ω_k = 0.0007 ± 0.0019")
    print()
    print("  Залишки (obs - theory) / sigma:")
    print(f"  {'z':>6}  {'obs':>7}  {'theory':>7}  {'резидуал':>8}")
    print("  " + "-"*38)
    for (z, obs, sigma, kind) in BAO_DATA:
        if kind == 'DV':
            th = d_V(z, H0_fit, Om_fit, OL_fit) / R_S
        else:
            th = d_C(z, H0_fit, Om_fit, OL_fit) / (1+z) / R_S
        res = (obs - th) / sigma
        print(f"  {z:>6.3f}  {obs:>7.3f}  {th:>7.3f}  {res:>+8.2f}σ")
    print("=" * 55)

if __name__ == '__main__':
    main()
