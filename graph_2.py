import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from mpmath import mp, sqrt, ln, mpc, fabs

plt.rcParams["font.family"] = "Latin Modern Math"
mp.dps = 50  

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

a_values = np.linspace(0.0001, 1.0, 10000)
sqrt_a_vals = []
percent_error_vals = []

x0, x1 = 0, 1

for a in a_values:
    I1_val = I1(x1, a) - I1(x0, a)
    I2_val = I2(x1, a) - I2(x0, a)
    I3_val = I3(x1, a) - I3(x0, a)

    combo = (I1_val + I2_val + I3_val) / (4*a + 1)
    exact = 1 / a
    error_abs = fabs(combo - exact)
    percent_error = 100 * error_abs / fabs(exact)

    sqrt_a_vals.append(float(sqrt(a)))
    percent_error_vals.append(float(percent_error))

sqrt_a_vals = np.array(sqrt_a_vals)
percent_error_vals = np.array(percent_error_vals)

# Gráfica
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(sqrt_a_vals, percent_error_vals, color="deeppink", lw=2, label="Relative error (%)")

below_10 = percent_error_vals < 10
ax.fill_between(sqrt_a_vals, percent_error_vals, where=below_10, color="pink", alpha=0.4, label="Error < 10%")

ax.axhline(10, color="gray", linestyle="--", linewidth=1)

ax.set_xlabel(r'$\sqrt{a}$', fontsize=14)
ax.set_ylabel('Relative error (%)', fontsize=14)
ax.set_title(r'Error porcentual de $\frac{I_1 + I_2 + I_3}{4a+1}$ frente a $\frac{1}{a}$', fontsize=15)

ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.01))
ax.tick_params(axis='x', which='minor', length=3)
ax.legend(fontsize=40)
ax.tick_params(axis='both', which='major', labelsize=13)
ax.tick_params(axis='both', which='minor', labelsize=10)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend()
plt.tight_layout()
plt.show()
