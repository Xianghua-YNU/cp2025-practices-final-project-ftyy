# 2_code/ -源代码目录

**目的:** 存放所有用于模拟、分析和可视化的代码。代码的质量是评分的重要组成部分。

**内容:**
- **main_simulation.py:** 模型模拟计算。负责设置参数、用玻尔兹曼模型计算吸收线强度和等效宽度。
- **absorption_profile.py:** 吸收线谱型可视化。负责读取main_simulation.py生成的吸收线强度的数据，生成谱型可视化结果。
- **numerical_methods.py:** 核心算法。实现具体的数值方法，例如龙格-库塔法、Crank-Nicolson、Metropolis算法等。这个模块应该具有通用性，不依赖于特定的模拟场景。
- **data_analysis.py:** 数据分析。存放用于处理 raw_data 的函数，例如计算平均值、误差、傅里叶变换、拟合曲线等。
- **visualization.py:** 可视化。存放所有绘图函数。这些函数应该接收数据作为输入，然后生成图表。这使得绘图逻辑与计算逻辑分离。
- **requirements.txt :** 项目依赖。通过 pip freeze > requirements.txt 命令生成，确保任何人都可以在新环境中复现你们的计算环境。

  - 配置环境：运行pip install -r requirements.txt。
  - 运行主程序：python main_simulation.py。


## 注意：
你也可以不按照预设的目录结构，按照自己的喜好组织代码，只要能够清晰地说明每个文件的功能即可。
