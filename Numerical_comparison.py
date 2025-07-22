import numpy as np
from scipy.integrate import quad
from mpmath import mp, sqrt, ln, mpc, nstr

mp.dps = 16

def integrand_total(x, a):
    return 1 / ((a + x**2)**1.5 * (a + (1 - x)**2))

def I1(x, a):
    return (x - 2*a) / (a * sqrt(a + x**2))

def I2(x, a):
    sqrt_a = sqrt(a)
    factor = 1 / (2 * sqrt_a) - 1j
    root = sqrt(-2j * sqrt_a - 1)
    num = 1j * (sqrt(a + x**2) - x) - (sqrt_a - 1j) - root
    den = 1j * (sqrt(a + x**2) - x) - (sqrt_a - 1j) + root
    return factor / root * ln(num / den)

def I3(x, a):
    sqrt_a = sqrt(a)
    factor = 1 / (2 * sqrt_a) + 1j
    root = sqrt(2j * sqrt_a - 1)
    num = -1j * (sqrt(a + x**2) - x) - (1j + sqrt_a) - root
    den = -1j * (sqrt(a + x**2) - x) - (1j + sqrt_a) + root
    return factor / root * ln(num / den)

a_values = np.linspace(0.1, 1.0, 10)
x0 = 0
x1 = 1.0

print("Validación de la expresión total:\n")
for a in a_values:
    num_val, _ = quad(lambda x: integrand_total(x, a), x0, x1)

    I1_val = I1(x1, a) - I1(x0, a)
    I2_val = I2(x1, a) - I2(x0, a)
    I3_val = I3(x1, a) - I3(x0, a)

    total_analytic = (I1_val + I2_val + I3_val) / (4*a + 1)

    print(f"a = {a:.2f}")
    print(f"Numérica : {nstr(num_val, 12)}")
    print(f"Analítica: {nstr(total_analytic, 12)}")
    print(f"Error    : {nstr(abs(total_analytic - num_val), 4)}\n")
