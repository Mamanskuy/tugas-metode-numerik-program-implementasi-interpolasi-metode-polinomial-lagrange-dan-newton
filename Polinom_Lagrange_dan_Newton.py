print("||    PROGRAM IMPLEMENTASI INTERPOLASI   ||")
print("||          LAGRANGE DAN NEWTON          ||")
print("||   ALMAN KAMAL MAHDI - 21120122120024  ||")
print("||         METODE NUMERIK KELAS B        ||")

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


class LagrangeInterp:
    def __init__(self, stress, failure_time):
        self.stress = stress
        self.failure_time = failure_time
    
    def interpolate(self, x_val):
        total_sum = 0
        n = len(self.stress)
        for i in range(n):
            term = self.failure_time[i]
            for j in range(n):
                if j != i:
                    term = term * (x_val - self.stress[j]) / (self.stress[i] - self.stress[j])
            total_sum += term
        return total_sum

    def get_poly(self):
        x_sym = sp.symbols('x_sym')
        total_sum = 0
        n = len(self.stress)
        for i in range(n):
            term = self.failure_time[i]
            for j in range(n):
                if j != i:
                    term *= (x_sym - self.stress[j]) / (self.stress[i] - self.stress[j])
            total_sum += term
        return sp.simplify(total_sum)
    
    def plot(self, x_vals, ax, plot_color):
        y_vals = [self.interpolate(xi) for xi in x_vals]
        ax.plot(self.stress, self.failure_time, 'ro', label='Data points')
        ax.plot(x_vals, y_vals, color=plot_color, label='Lagrange Interpolation')
        ax.set_xlabel('Stress (kg/mm^2)')
        ax.set_ylabel('Failure Time (hours)')
        ax.legend()
        ax.set_title('Lagrange Interpolation')


class NewtonInterp:
    def __init__(self, stress, failure_time):
        self.stress = stress
        self.failure_time = failure_time
        self.div_diff_table = self.calc_div_diff()
    
    def calc_div_diff(self):
        n = len(self.stress)
        div_diff = np.zeros([n, n])
        div_diff[:, 0] = self.failure_time

        for j in range(1, n):
            for i in range(n - j):
                div_diff[i][j] = (div_diff[i + 1][j - 1] - div_diff[i][j - 1]) / (self.stress[i + j] - self.stress[i])
        
        return div_diff[0, :]
    
    def interpolate(self, x_val):
        n = len(self.stress)
        result = self.div_diff_table[0]
        for i in range(1, n):
            term = self.div_diff_table[i]
            for j in range(i):
                term *= (x_val - self.stress[j])
            result += term
        return result

    def get_poly(self):
        x_sym = sp.symbols('x_sym')
        n = len(self.stress)
        result = 0
        for i in range(n):
            term = self.div_diff_table[i]
            for j in range(i):
                term *= (x_sym - self.stress[j])
            result += term
        return sp.simplify(result)
    
    def plot(self, x_vals, ax, plot_color):
        y_vals = [self.interpolate(xi) for xi in x_vals]
        ax.plot(self.stress, self.failure_time, 'ro', label='Data points')
        ax.plot(x_vals, y_vals, color=plot_color, label='Newton Interpolation')
        ax.set_xlabel('Stress (kg/mm^2)')
        ax.set_ylabel('Failure Time (hours)')
        ax.legend()
        ax.set_title('Newton Interpolation')

# Data
stress_vals = np.array([5, 10, 15, 20, 25, 30, 35, 40])
failure_times = np.array([40, 30, 25, 40, 18, 20, 22, 15])

# Test Interpolasi Lagrange
x_test_lagrange = np.linspace(5, 40)
lagrange_interp = LagrangeInterp(stress_vals, failure_times)

# Test Interpolasi Newton
x_test_newton = np.linspace(5, 40)
newton_interp = NewtonInterp(stress_vals, failure_times)

# Membuat subplots
fig, axs = plt.subplots(1, 2, figsize=(14, 7))

# Plot Lagrange Interpolation dengan warna biru
lagrange_interp.plot(x_test_lagrange, axs[0], 'b')

# Plot Newton Interpolation dengan warna hijau
newton_interp.plot(x_test_newton, axs[1], 'g')

# Tampilkan plot
plt.tight_layout()
plt.show()

# Menampilkan hasil fungsi polinomial lagrange
lagrange_poly = lagrange_interp.get_poly()
print(f"Fungsi Polinomial hasil Interpolasi Lagrange:\n {lagrange_poly.expand()}")

# Menampilkan hasil fungsi polinomial newton
newton_poly = newton_interp.get_poly()
print(f"Fungsi Polinomial hasil Interpolasi Newton:\n {newton_poly.expand()}")
