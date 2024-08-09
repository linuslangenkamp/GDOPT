import pandas as pd
import matplotlib.pyplot as plt

path = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData"
model = "dieselMotor"
it = 0
df = pd.read_csv(path + "/" + model + str(it) + ".csv", sep=",")
print(df.head())
specifCols = None # [col for col in df.columns if col.startswith('u')] 
interval = [0, 0.5]

plt.rcParams.update({
    'font.serif': ['Times New Roman'],
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'legend.frameon': True,
    'legend.loc': 'best',
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'figure.figsize': (12, 8),
    'axes.titlepad': 20,
    'axes.labelpad': 10
})

if specifCols is None:
    columns_to_plot = df.columns[1:]
else:
    columns_to_plot = specifCols

num_plots = len(columns_to_plot)
fig, axs = plt.subplots(num_plots, 1, figsize=(12, 8 * num_plots), sharex=True)

if num_plots == 1:
    axs = [axs]

for idx, column in enumerate(columns_to_plot):
    ax = axs[idx]
    ax.plot(df['time'], df[column], label=column, linewidth=2, linestyle='-', color='steelblue')
    ax.scatter(df['time'], df[column], color='red', s=30, edgecolor='black', alpha=0.8, zorder=5)
    ax.set_xlabel('time')
    ax.set_ylabel(column)
    ax.set_xlim(interval[0], interval[1])
    ax.legend(frameon=True, loc='best')
    ax.grid(True)
    ax.title.set_fontsize(16)

plt.tight_layout()
plt.show()
