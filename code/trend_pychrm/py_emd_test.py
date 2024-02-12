import numpy as np
import matplotlib.pyplot as plt
from pyemd import EMD

# Set seed for reproducibility
np.random.seed(42)

# Create synthetic time series with secular trend and two periodic components
time = np.arange(0, 40 * 12, 1)  # 40 years of monthly data
secular_trend = 0.2 * time / 120  # 0.2 degrees per decade
periodic_1 = 0.1 * np.sin(2 * np.pi * (1 / (2 + np.random.uniform(-0.5, 0.5))) * time)
periodic_2 = 0.15 * np.sin(2 * np.pi * (1 / (16 + np.random.uniform(-1, 1))) * time)
noise = np.random.normal(0, 0.1, len(time))

temperature_anomalies = secular_trend + periodic_1 + periodic_2 + noise

# Plot the synthetic time series
plt.figure(figsize=(12, 6))
plt.plot(time, temperature_anomalies, label='Synthetic Temperature Anomalies')
plt.xlabel('Time (months)')
plt.ylabel('Temperature Anomalies')
plt.title('Synthetic Temperature Anomalies Data')
plt.legend()
plt.show()

# Apply Empirical Mode Decomposition (EMD)
emd = EMD()
IMFs = emd.emd(temperature_anomalies)

# Compute the periods of each IMF
fs = 1.0 / (time[1] - time[0])
periods = [1 / freq for freq in emd..frequencies]

# Plot the original time series and its IMFs
plt.figure(figsize=(12, 6))
plt.plot(time, temperature_anomalies, label='Original Temperature Anomalies')
for i, imf in enumerate(IMFs):
    plt.plot(time, imf, label=f'IMF {i+1}')

plt.xlabel('Time (months)')
plt.ylabel('Temperature Anomalies')
plt.title('EMD Decomposition of Temperature Anomalies')
plt.legend()
plt.show()

# Reconstruction without periodic components
trend_component = IMFs[5]

# Plot the original and reconstructed time series
plt.figure(figsize=(12, 6))
plt.plot(time, temperature_anomalies, label='Original Temperature Anomalies')
plt.plot(time, trend_component, label='Trend Component')
plt.xlabel('Time (months)')
plt.ylabel('Temperature Anomalies')
plt.title('Reconstructed Temperature Anomalies without Periodic Components')
plt.legend()
plt.show()