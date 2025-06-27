"""
1.提取同一吸收线不同温度的等效宽度数据并绘制，标上scale。
"""
import pandas as pd
import matplotlib.pyplot as plt

# 读取数据
df = pd.read_csv("balmer_line_stellar.txt", delim_whitespace=True, comment="#", header=None,
                 names=["Line", "Wavelength", "Teff", "Range1", "EW1",
                        "Range2", "EW2", "Scale", "EW3", "EW4"])

# 筛选 Hβ 数据
df_Hb = df[df["Line"] == "Hβ"].sort_values(by="Teff")

# 提取横轴和纵轴
temps = df_Hb["Teff"]
ew4 = df_Hb["EW4"]
scales = df_Hb["Scale"]

# 绘图
plt.figure(figsize=(8, 5))
plt.plot(temps, ew4, marker='o', color='blue', label="Equivalent Width")

# 标注每个点的 scale 值
for x, y, s in zip(temps, ew4, scales):
    plt.text(x, y + 0.2, f"scale={s}", ha='center', fontsize=9, color='blue')

# 美化图形
plt.title("Hβ Equivalent Width at different Temperature")
plt.xlabel("Effective Temperature (K)")
plt.ylabel("Equivalent Width(Å)")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("Hβ_EW.png", dpi=300)
plt.show()


"""
2.绘制不同方法计算7000K的Hβ吸收率谱。红线是玻尔兹曼分布方法的谱，StellarSpecModel 恒星光谱模型工具的谱。第二部分绘图是绘制两种方法的吸收率的差值绝对值，并计算了差值绝对值最大值和平均值。
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import h, c, k_B
from scipy.signal import savgol_filter

# 读取 MARCS 模型光谱（单位：erg / (s cm² Å)）
data = np.loadtxt('marcs_spectrum_teff7000_feh0_logg4.dat')
wavelength_all = data[:, 0]  # 单位：Å
flux_all = data[:, 1]        # 单位：erg / (s cm² Å)

# --- 滤波预处理 ---
flux_smoothed = savgol_filter(flux_all, window_length=51, polyorder=3)

# 恒星参数
T_eff = 7000  # 与模型一致

# 黑体函数（输出单位：erg / (s cm² Å)）
def blackbody_lambda(wavelength, T):
    lam_cm = wavelength * 1e-8  # Å → cm
    h_val = h.cgs.value
    c_val = c.cgs.value
    k_val = k_B.cgs.value
    exponent = h_val * c_val / (lam_cm * k_val * T)
    B = (2.0 * h_val * c_val**2 / lam_cm**5) / (np.exp(exponent) - 1.0)
    return B * 1e-8 * np.pi  # erg / (s cm² Å)

# 波长范围（Hζ）
wave_min = 4846.5
wave_max = 4876.1

# 截取目标波段
mask = (wavelength_all >= wave_min) & (wavelength_all <= wave_max)
wavelength = wavelength_all[mask]
flux = flux_smoothed[mask]

# 黑体连续谱计算并缩放
B_lam = blackbody_lambda(wavelength, T_eff)
scale = 1.1
continuum = B_lam * scale

# 归一化光谱
flux_norm = flux / continuum

# 吸收率：1 - 归一化通量
absorption = 1 - flux_norm

# 计算等效宽度（积分）
ew = np.trapz(absorption, wavelength)

# 读取 MARCS 模型光谱（单位：erg / (s cm² Å)）
data_1 = np.loadtxt('Hβ_7000K_profile.txt')
wavelength_profile = data_1[:, 0]         # Å
I_cont = data_1[:, 1]             # erg/s/cm²/Å/sr
I_abs = data_1[:, 2]  


mask_1 = (wavelength_profile >= wave_min) & (wavelength_profile <= wave_max)

wavelength_profile= wavelength_profile[mask_1]
I_cont = I_cont[mask_1]
I_abs = I_abs[mask_1]

# 计算吸收率谱（归一化吸收线）
absorption_1 = 1 - (I_abs / I_cont)
# 计算等效宽度（单位：Å）
ew_1 = np.trapz(absorption_1, wavelength_profile)

# --- 绘图 ---
plt.figure(figsize=(10, 5))
plt.plot(wavelength, absorption, color='blue', lw=2, label='Absorption = 1 - F / Fc')
plt.plot(wavelength_profile, absorption_1, color='red', lw=2, label='profile_Absorption = 1 - F / Fc')
plt.fill_between(wavelength, absorption, color='blue', alpha=0.2)
plt.title(f"Hβ Absorption, fontsize=14")
plt.xlabel("Wavelength (Å)")
plt.ylabel("Absorption (1 - F / Fc)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Hβ_7000K_absorption_profile_comparison.png", dpi=300)
plt.show()

# 可选：单独绘制差值图
plt.figure(figsize=(10, 4))
plt.plot(wavelength, abs_diff, color='purple')
plt.fill_between(wavelength, 0, abs_diff, color='purple', alpha=0.3)
plt.xlabel("Wavelength (Å)")
plt.ylabel("Absolute Difference")
plt.title("Absolute Difference between Two Absorption Profiles")
plt.grid(True)
plt.tight_layout()
plt.savefig("Hβ_7000K_absorption_difference_only.png", dpi=300)
plt.show()

# 输出差值的一些统计量
print(f"最大绝对差值: {np.max(abs_diff):.4f}")
print(f"平均绝对差值: {np.mean(abs_diff):.4f}")


