import numpy as np
import matplotlib.pyplot as plt
import glob

# 设置图形
plt.figure(figsize=(10, 6))
plt.title("Normalized Absorption Spectra")
plt.xlabel("Wavelength")
plt.ylabel("Normalized Intensity")

# 获取所有txt文件（根据实际情况修改路径）
file_list = glob.glob("H*profile.txt",)  # 假设所有txt文件在当前目录

# 处理每个文件
for file in file_list:
    # 加载数据
    data = np.loadtxt(file)
    wavelength = data[:, 0]
    continuum = data[:, 1]
    intensity = data[:, 2]
    
    # 计算归一化强度
    normalized = intensity / continuum
    
    # 绘制曲线（使用文件名作为标签）
    plt.plot(wavelength, normalized, '-', label=file.replace('.txt', ''))

# 添加图例和网格
# 添加图例和网格
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))  # 图例放在图外右侧
plt.grid(alpha=0.3)
plt.tight_layout(rect=[0, 0, 0.85, 1])  # 为图例留出空间

# 显示或保存图形
plt.show()
#plt.savefig("normalized_spectra.png")  