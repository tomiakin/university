import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
# NOTE: you can also use import scipy.integrate as spi, but notation changes see Activity 3

def pdf(x):
    return x / 210

def cdf(x):
    return np.cumsum(pdf(x))

x_values = np.arange(1, 21)
pdf_values = pdf(x_values)
cdf_values = cdf(x_values)

plt.figure(figsize=(12, 4))

# Plot PDF
plt.subplot(1, 2, 1)
plt.bar(x_values, pdf_values, align='center', alpha=0.7)
plt.xlabel('X')
plt.ylabel('PDF')
plt.title('Probability Density Function (PDF)')

# Plot CDF
plt.subplot(1, 2, 2)
plt.plot(x_values, cdf_values, marker='o', linestyle='-')
plt.xlabel('X')
plt.ylabel('CDF')
plt.title('Cumulative Distribution Function (CDF)')

plt.tight_layout()
plt.show()

# Finding p(X>17)
cdf_17th = cdf_values[16]
prob = 1 -cdf_17th

print('p(X>17) is', prob)




