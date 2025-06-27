"""
使用 MARCS 模型光谱数据，计算有效宽度，并画出归一化光谱，以调整scale
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import h, c, k_B
from scipy.signal import savgol_filter

# 读取 MARCS 模型光谱（单位：erg / (s cm² Å)）
data = np.loadtxt('marcs_spectrum_teff8000_feh0_logg4.dat')
wavelength_all = data[:, 0]  # 单位：Å
flux_all = data[:, 1]        # 单位：erg / (s cm² Å)

# --- 滤波预处理 ---
flux_smoothed = savgol_filter(flux_all, window_length=51, polyorder=3)


# 参数
T_eff = 8000  # 与模型一致

# 黑体函数，输出单位 erg / (s cm² Å steradian)
def blackbody_lambda(wavelength, T):
    lam_cm = wavelength * 1e-8  # Å → cm
    h_val = h.cgs.value
    c_val = c.cgs.value
    k_val = k_B.cgs.value
    exponent = h_val * c_val / (lam_cm * k_val * T)
    B = (2.0 * h_val * c_val**2 / lam_cm**5) / (np.exp(exponent) - 1.0)
    return B*1e-8 # 单位：erg / (s cm² Å steradian)

# 黑体光谱 × π，得到通量密度（与模型匹配的单位）
B_lam = blackbody_lambda(wavelength_all, T_eff) * np.pi  # 单位：erg / (s cm² Å)

# ---------- 归一化谱 ----------
flux_norm = flux_smoothed / B_lam

# 取flux_norm <= 1.1的值
mean_value = np.mean(flux_norm)

print(f"Mean of flux_norm: {mean_value:.8e}")

plt.plot(wavelength_all, flux_norm, '--', label='normalized spectrum')
plt.ylabel('normalized spectrum' )
plt.xlim(4821.7, 4900.9)
plt.ylim(0, 1.1 * max(flux_norm))
plt.legend()
plt.savefig('normalized_spectrum.jpg')
plt.show()


wave_min = 3876.3
wave_max = 3901.7
# 选择目标波段
mask = (wavelength_all >= wave_min) & (wavelength_all <= wave_max)
wavelength = wavelength_all[mask]
flux = flux_smoothed[mask]
B_lam_1=blackbody_lambda(wavelength, T_eff) * np.pi 
scale=1.2
continuum = B_lam_1 * scale

# 归一化后的谱（模型 / 连续谱）
flux_norm_1 = flux / continuum

ew = np.trapz(1 - flux_norm_1, wavelength)

print(f'使用黑体连续谱计算的 H 等效宽度 EW = {ew:.4f} Å')

plt.figure(figsize=(8, 5))
plt.plot(wavelength, 1-flux_norm_1, label='Normalized Spectrum')
plt.xlabel('Wavelength (Å)')
plt.ylabel('Normalized Flux')
plt.title(f"Hα Equivalent Width (T={T_eff}K) = {ew:.4f} Å")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
