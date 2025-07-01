"""
代码功能：利用玻尔兹曼分布得出的理论公式，计算吸收线强度和等效宽度，并生成结果
可调参数：logg,metallicity,T_list

"""
import numpy as np
from scipy.special import wofz
from astropy import constants as const
from astropy import units as u
from scipy.integrate import simpson 

# 物理常数
e = const.e.gauss.value  # 电子电荷 (esu)
me = const.m_e.cgs.value  # 电子质量 (g)
c = const.c.cgs.value    # 光速 (cm/s)
h = const.h.cgs.value    # 普朗克常数 (erg·s)
k = const.k_B.cgs.value   # 玻尔兹曼常数 (erg/K)
m_H = const.m_p.cgs.value  # 氢原子质量 (g)
Ry = 13.605693122994 * 1.60218e-12  # 里德伯能量 (erg)  

T_list = [5000, 6000, 7000, 8000]    # 温度列表 (K)
logg = 4.0                                  # logg值 (cgs)
  
metallicity = 0.0  # 金属丰度
# 巴尔末线系定义
balmer_series = {
    'Hβ': {'wavelength': 4861.3, 'n_l': 2, 'n_u': 4, 'f_lu': 0.1193, 'gl': 8},
    'Hγ': {'wavelength': 4340.5, 'n_l': 2, 'n_u': 5, 'f_lu': 0.0447, 'gl': 8},
    'Hδ': {'wavelength': 4101.7, 'n_l': 2, 'n_u': 6, 'f_lu': 0.0221, 'gl': 8},
    'Hε': {'wavelength': 3970.1, 'n_l': 2, 'n_u': 7, 'f_lu': 0.0128, 'gl': 8},
    'Hζ': {'wavelength': 3889.0, 'n_l': 2, 'n_u': 8, 'f_lu': 0.0082, 'gl': 8}
}

def partition_function(T):   # 计算配分函数
    n_max = min(20, int(1 + np.sqrt(Ry / (k*T))))
    Z = 0
    for n in range(1, n_max+1):
        g_n = 2 * n**2
        E_n = -Ry / n**2
        Z += g_n * np.exp(-E_n / (k * T))
    return Z

def calculate_P_e(T, N_total, tol=1e-4, max_iter=100):
    """迭代求解电子压强 P_e (cgs单位)"""
    #P_e = 100    #可固定参数看其影响
    P_e = N_total * k * T * 0.1  # 初始猜测 (dyn/cm^2)
    for _ in range(max_iter):
        N_HI = saha_equation(T, N_total, P_e)
        n_e = (N_total - N_HI)
        P_e_new = n_e * k * T
        if abs((P_e_new - P_e) / P_e_new) < tol:
            return P_e_new
            
        P_e = P_e_new
    
    raise RuntimeError(f"迭代未收敛，在{max_iter}次迭代后相对误差为{(P_e_new - P_e) / P_e_new}")



def calc_N_total(Teff, metallicity):
    # 金属丰度影响电子密度，进而影响总粒子数
    N_total = 10**(19.5 - 0.5*np.log10(Teff) - 0.4*logg + 0.1*metallicity)
    if N_total < 1e12 or N_total > 1e20:
        print(f"警告: 计算的粒子数密度可能不合理: {N_total:.2e} cm^-3")
    
    return N_total

def scale_height(T, logg, N_HI, N_total):  #计算标高
    g = 10**logg
    x = 1 - N_HI / N_total  # 电离度
    mu = 1 / (1 + x)        # 考虑电离的平均分子量
    return (k * T) / (mu * m_H * g)

# 萨哈方程计算中性氢比例
def saha_equation(T, N_total, P_e):
    """计算中性氢比例"""
    # 氢的电离能 (erg)
    chi_H = 13.6 * 1.60218e-12
    
    # 氢的配分函数 (中性氢Z0=2, 一阶电离氢Z1=1)
    Z0 = 2.0
    Z1 = 1.0
    
    # 萨哈方程
    n_e = P_e / (k * T)  # 电子数密度 (cm^-3)
    S = (2 * np.pi * me * k * T)**1.5 / h**3 * 2 * Z1 / Z0 * np.exp(-chi_H / (k * T))
    x = np.sqrt(S**2 + 4 * S * N_total / n_e) - S
    N_HI = N_total * x / (x + 1)
    return N_HI


# 玻尔兹曼分布
def boltzmann_column_density(gl, El, T, Z, N_HI, logg, N_total):
    """计算能级粒子柱密度 (cm^-2)"""
    N_l = (gl * np.exp(-El/(k*T)) / Z) * N_HI  # 原始数密度 (cm^-3)
    H = scale_height(T, logg, N_HI, N_total)   # 大气标高 (cm)
    return N_l * H                              # 柱密度 (cm^-2)

# 高斯型谱线轮廓 (多普勒展宽)
def gaussian_profile(wavelength, line_center, T):
    """返回高斯型谱线轮廓函数值"""
    delta_lambda_D = (line_center / c) * np.sqrt(2*k*T/m_H)   # 多普勒宽度 (Å)
    return (1/(delta_lambda_D * np.sqrt(np.pi))) * np.exp(-((wavelength - line_center)/delta_lambda_D)**2)

# Voigt型谱线轮廓 (多普勒+自然展宽)
def voigt_profile(wavelength, line_center, T, P_e):
    """返回Voigt型谱线轮廓函数值"""
    delta_lambda_D = (line_center / c) * np.sqrt(2*k*T/m_H)   # 多普勒宽度 (Å)
    gamma = 6.265e8 / (line_center**2)  # 自然展宽参数 (Å)
    sigma = delta_lambda_D / np.sqrt(2)
    n_e = P_e / (k * T)
    stark_width = 1e-8 * (n_e / 1e14)**0.4 * (T/10000)**0.1  # 示例修正   
    total_gamma = gamma + stark_width
    
    sigma = delta_lambda_D / np.sqrt(2)
    x = (wavelength - line_center) / sigma
    y = total_gamma / (2 * sigma)  # 使用总展宽参数
    z = x + 1j*y
    w = wofz(z)
    return w.real / (sigma * np.sqrt(2*np.pi))

# 黑体辐射谱
def blackbody_spectrum(wavelength, T):
    """返回黑体辐射强度 (erg/s/cm^2/Å/sr)"""
    wavelength_cm = wavelength * 1e-8  # 转换为cm
    return (2*h*c**2 / (wavelength_cm**5)) / (np.exp(h*c/(wavelength_cm*k*T)) - 1)

def get_5sigma_lambda_range(line_center, T):  #计算5σ区间
    # 多普勒展宽 (Å): Δλ_D = λ_0 * √(2kT/mc²)
    delta_lambda_D = line_center * np.sqrt(2 * k * T / (m_H * c**2))
    # 标准差 (Å): σ = Δλ_D / (2√(2ln2))
    sigma = delta_lambda_D / (2 * np.sqrt(2 * np.log(2)))
    # 5σ 范围
    lambda_min = line_center - 200 * sigma
    lambda_max = line_center + 200 * sigma
    return np.linspace(lambda_min, lambda_max, 2000)

# 计算吸收线强度
def calc_absorption_line(wavelength_range, line_center, f_lu, gl, El, T, profile_type='voigt'):
    """计算吸收线强度"""
    
    N_total = calc_N_total(T, metallicity)
    P_e = calculate_P_e(T, N_total)
    N_HI = saha_equation(T, N_total, P_e)
    Z = partition_function(T)
    N_col = boltzmann_column_density(gl, El, T, Z, N_HI, logg, N_total)
    
    # 吸收系数中的常数部分
    const_part = (np.pi * e**2) / (me * c) * f_lu   
    delta_lambda_D = (line_center * 1e-8) * np.sqrt(2 * k * T / m_H) / c  # 多普勒宽度 (cm)
    delta_nu_D = (c / (line_center * 1e-8)**2) * delta_lambda_D
    nu_mn = c / (wavelength_range * 1e-8)  # 线心频率 (Hz)
    stimulation_correction = 1 - np.exp(-h * nu_mn / (k * T))
    # 选择谱线轮廓
    if profile_type == 'gaussian':
        phi = gaussian_profile(wavelength_range, line_center, T)  # Å^-1
    else:
        phi = voigt_profile(wavelength_range, line_center, T, P_e)  # Å^-1
    
    # 连续谱强度
    I_cont = blackbody_spectrum(wavelength_range, T)  # erg/s/cm^2/Å/sr
    
    # 吸收线强度
    #tau = const_part * N_col / delta_nu_D * stimulation_correction   #量级估算(可选)
    tau = const_part * N_col * phi * 1e-8   # 转换为Å^-1 (因为phi是Å^-1而const_part是cm^-1)
    I_abs = I_cont * np.exp(-tau)  # erg/s/cm^2/Å/sr
    print(f"T={T}K: N_HI={N_HI:.2e} cm^-3, N_col={N_col:.2e} cm^-2, tau_max={tau.max():.2f}, ")   #tau_wing={tau[len(tau)//2]:.2e}
    #print(f"T={T}K: N_HI={N_HI:.2e} cm^-3, H={H:.2e} cm, N_l={N_l:.2e} cm^-3, N_col={N_col:.2e} cm^-2")
    print(f"P_e={P_e:.2e} dyn/cm^2, n_e={P_e/(k*T):.2e} cm^-3, x={1 - N_HI/N_total:.2f}")
    print(f"phi_max={phi.max():.2e} Å^-1, delta_nu_D={delta_nu_D:.2e} Hz")
    return I_abs, I_cont

# 计算等效宽度
def calc_equivalent_width(wavelengths, I_cont, I_abs):
    """计算等效宽度 (Å)，使用Simpson积分提高精度"""
    y = 1 - I_abs / I_cont
    EW = simpson(y, wavelengths)
    return EW

# 主计算程序
def main():
    results = []
    
    for line_name, line_data in balmer_series.items():
        line_center = line_data['wavelength']
        f_lu = line_data['f_lu']
        gl = line_data['gl']
        n_l = line_data['n_l']
        El = Ry * (1/n_l**2)    # 巴尔末线系基态能量 (erg)  
        
        line_results = []
        
        # 为每条谱线创建单独的结果文件
        with open(f'{line_name}_results.txt', 'w', encoding='utf-8') as line_file:
            line_file.write(f"# {line_name} 线计算结果 (波长: {line_center} Å)\n")
            line_file.write("# 温度(K) | 总粒子数密度(cm^-3) | 电子压强(dyn/cm^2) | 配分函数 | 波长范围(Å) | 等效宽度(Å)\n")
            
            for T in T_list:
                print(f"\n计算 {line_name} 线，温度 T = {T} K")
                
                N_total = calc_N_total(T, metallicity)
                P_e = calculate_P_e(T, N_total)
                Z = partition_function(T)
                
                wavelengths = get_5sigma_lambda_range(line_center, T)
                wavelength_range = (min(wavelengths), max(wavelengths))
                
                I_abs, I_cont = calc_absorption_line(wavelengths, line_center, f_lu, gl, El, T)
                EW = calc_equivalent_width(wavelengths, I_cont, I_abs)
                
                # 保存谱线轮廓数据
                np.savetxt(
                    f"{line_name}_{T}K_profile.txt",
                    np.column_stack((wavelengths, I_cont, I_abs)),
                    header="Wavelength(Å) I_cont(erg/s/cm2/Å/sr) I_abs(erg/s/cm2/Å/sr)",
                    fmt="%.4f %.4e %.4e",
                    encoding="utf-8"
                )

                line_results.append({
                    'Temperature': T,
                    'Equivalent_Width': EW,
                    'N_total': N_total,
                    'P_e': P_e,
                    'Z': Z,
                    'Wavelength_Range': wavelength_range
                })
                
                # 将当前温度的结果写入线结果文件
                line_file.write(f"{T} {N_total:.2e} {P_e:.2e} {Z:.4f} {wavelength_range[0]:.1f}-{wavelength_range[1]:.1f} {EW:.4f}\n")
                
                print(f"{line_name} T={T}K: 等效宽度 EW={EW:.4f} Å")
        
        results.append({
            'Line': line_name,
            'Wavelength': line_center,
            'Results': line_results
        })
    
    # 创建汇总结果文件，包含波长范围信息
    with open('balmer_line_summary.txt', 'w', encoding='utf-8') as summary_file:
        summary_file.write("# 氢巴尔末线系计算汇总结果\n")
        summary_file.write("# 线名 | 波长(Å) | 温度(K) | 波长范围(Å) | 等效宽度(Å)\n")
        
        for res in results:
            line_name = res['Line']
            wavelength = res['Wavelength']
            
            for temp_res in res['Results']:
                T = temp_res['Temperature']
                EW = temp_res['Equivalent_Width']
                wl_range = temp_res['Wavelength_Range']
                summary_file.write(f"{line_name} {wavelength:.1f} {T} {wl_range[0]:.1f}-{wl_range[1]:.1f} {EW:.4f}\n")
    
    print("\n计算完成，结果已保存到:")
    print("- 各谱线详细结果文件: {LineName}_results.txt")
    print("- 谱线轮廓数据: {LineName}_{Temperature}K_profile.txt")
    print("- 汇总结果文件: balmer_line_summary.txt")
    
if __name__ == "__main__":
    main()
