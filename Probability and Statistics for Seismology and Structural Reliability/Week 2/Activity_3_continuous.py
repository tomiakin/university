import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

# Define the PDF function
def pdf(x):
    return np.where((0 <= x) & (x <= 1), 2 * x, 0)

# Calculate P(1/4 <= X <= 1/2)
probability, _ = spi.quad(pdf, 1/4, 1/2)

print(f'P(1/4 <= X <= 1/2) is approximately: {probability:.4f}')

# Calculate the integral of the PDF over the entire range [0, 1]
total_integral, _ = spi.quad(pdf, 0, 1)

print(f'The integral of the PDF over the entire range [0, 1] is approximately: {total_integral:.4f}')

# Plot the PDF
x_values = np.linspace(0, 1, 100)
pdf_values = pdf(x_values)

plt.figure(figsize=(8, 4))
plt.plot(x_values, pdf_values, label='PDF: 2x')
plt.xlabel('x')
plt.ylabel('PDF')
plt.title('Probability Density Function (PDF)')
plt.legend()
plt.grid(True)
plt.show()
