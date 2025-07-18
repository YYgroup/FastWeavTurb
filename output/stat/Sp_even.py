import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator, FuncFormatter
from matplotlib.font_manager import FontProperties

# 启用 LaTeX 渲染
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.rc('text.latex', preamble=r'\usepackage{newtxtext}\usepackage{newtxmath}')

# 定义数据文件和变量名
data_files = ['Sp.dat', 'Sp_DNS.dat']

# 读取数据
data_sets = []
for file in data_files:
    data = np.loadtxt(file)
    data_sets.append(data)

# 设置图形大小和布局
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

# 设置坐标轴范围和比例
ax.set_xlim(0.01, 2)
ax.set_ylim(1E-04, 100)
ax.set_xscale('log')
ax.set_yscale('log')

# 设置坐标轴标签
ax.set_xlabel(r'$r$', fontsize=46, fontfamily='Times New Roman', style='italic')
ax.set_ylabel(r'$S_p(r)$', fontsize=46, fontfamily='Times New Roman')
# 设置 x 轴和 y 轴标签与刻度线的间距
ax.xaxis.labelpad = -5  # 可以根据需要调整这个值
ax.yaxis.labelpad = 5  # 可以根据需要调整这个值

# 设置坐标轴刻度和标签
ax.xaxis.set_major_locator(LogLocator(base=10, numticks=10))
ax.xaxis.set_minor_locator(LogLocator(base=10, subs='all', numticks=10))
ax.yaxis.set_major_locator(LogLocator(base=10, numticks=10))
ax.yaxis.set_minor_locator(LogLocator(base=10, subs='all', numticks=10))
ax.yaxis.set_minor_formatter(NullFormatter())
def custom_y_formatter(x, pos):
    target_powers = [-6, -4, -2, 0, 2]
    power = np.log10(x)
    if power in target_powers:
        return r'$10^{{{}}}$'.format(int(power))
    return ''
# 创建自定义格式化器对象
y_formatter = FuncFormatter(custom_y_formatter)
# 将自定义格式化器应用到 y 轴主刻度上
ax.yaxis.set_major_formatter(y_formatter)
# 创建 FontProperties 对象指定字体家族和大小
tick_font = FontProperties(family='Times New Roman', size=46)
# 设置 x 轴刻度标签与刻度线的间距
ax.tick_params(axis='x', pad=15)  # 可根据需要调整这个值
# 设置 y 轴刻度标签与刻度线的间距
ax.tick_params(axis='y', pad=5)  # 可根据需要调整这个值
# 设置坐标框的线宽
for spine in ax.spines.values():
    spine.set_linewidth(2)  # 这里将线宽设置为 2，可以根据需要调整
# 设置主刻度线的粗细和长短
ax.tick_params(axis='both', which='major', width=2, length=8)
# 设置次刻度线的粗细和长短
ax.tick_params(axis='both', which='minor', width=1, length=4)
# 在上方和右方显示刻度线
ax.tick_params(axis='x', top=True, which='both')
ax.tick_params(axis='y', right=True, which='both')

# 设置刻度标签的字体和刻度线方向
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontproperties(tick_font)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontproperties(tick_font)
ax.tick_params(axis='both', which='major', direction='in')
ax.tick_params(axis='both', which='minor', direction='in')

# 设置网格线
# ax.grid(True, which='major', linestyle='-', linewidth=0.5)
# ax.grid(True, which='minor', linestyle='--', linewidth=0.2)

# 定义线条和符号的颜色和样式
line_colors = ['#6A5ACD', '#cab2d6', '#d73027', '#8dd3c7']
line_thickness = [4, 1]
symbol_size = [10, 20]
symbol_line_thickness = [0, 3]
symbol_skipping = [30, 1]
zorders = [3, 2]

# 存储绘图对象和标签
handles = []
labels = []

# 绘制线条和符号
for i, data in enumerate(data_sets):
    x = data[:, 0]
    for j in range(2,7):
        y = data[:, j]
        up_DNS = 0.7782
        if i < 1:
            line_label = fr'$S_{{{j}}}$ WT2'
            if j % 2 == 1:
                y = -y
                line_label = fr'$-S_{{{j}}}$ WT2'
        else:
            line_label = fr'$S_{{{j}}}$ DNS2'
            y = y/up_DNS**j
            if j % 2 == 1:
                y = -y
                line_label = fr'$-S_{{{j}}}$ DNS2'

        line_label = ''
        # 筛选出 y 大于等于 0 的点
        valid_indices = y >= 1E-08
        x_valid = x[valid_indices]
        y_valid = y[valid_indices]


        if j % 2 == 0:
            linestyle = '-' if i == 0 else 'None'

            line, = ax.plot(x_valid, y_valid, color=line_colors[j//2-1], linewidth=line_thickness[i], label=line_label, marker='o',
                            markersize=symbol_size[i],
                            markevery=symbol_skipping[i], markerfacecolor='none', markeredgewidth=symbol_line_thickness[i],
                            zorder=zorders[i], linestyle=linestyle)
            handles.append(line)
            labels.append(line_label)


# 定义新的顺序，这里将 Gaussian 放到最上面
new_order = []
# 按照新顺序重新排列 handles 和 labels
new_handles = [handles[i] for i in new_order]
new_labels = [labels[i] for i in new_order]

# 设置图例
# 增大字体大小
font = FontProperties(family='Times New Roman', size=44)
# 调整标记大小
legend = ax.legend(new_handles, new_labels, prop=font, loc='lower center', bbox_to_anchor=(0.74, -0.04), frameon=False, markerscale=1,
                   handletextpad=0.25,  # 调整图例标记和文本之间的间距
                   labelspacing=-0.05)  # 调整图例项之间的垂直间距
# 添加文本注释
ax.text(0.1, 0.8, r'$r_L=4$', transform=ax.transAxes, fontsize=46, fontfamily='Times New Roman',color='black')

# 显示图形
plt.tight_layout()
# 先保存图形
plt.savefig('./Sp_even.pdf', format='pdf')
# 再显示图形
plt.show()