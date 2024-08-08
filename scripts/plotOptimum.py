import pandas as pd
import matplotlib.pyplot as plt

path = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData"
model = "satellite"
it = 0
specifCol = "u1"
interval = [0, 100]

df = pd.read_csv(path + "/" + model + str(it) + ".csv", sep=",")
print(df.head())

plt.rcParams.update({
    'font.family': 'serif',
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

if specifCol is None:
    num_plots = len(df.columns) - 1
    fig, axs = plt.subplots(num_plots, 1, figsize=(12, 8 * num_plots), sharex=True)

    if num_plots == 1:
        axs = [axs]

    for idx, column in enumerate(df.columns[1:]):
        ax = axs[idx]
        ax.plot(df['time'], df[column], label=column, linewidth=2, linestyle='-', color='steelblue')
        ax.scatter(df['time'], df[column], color='red', s=30, edgecolor='black', alpha=0.8, zorder=5)
        ax.set_xlabel('Time')
        ax.set_ylabel(column)
        ax.set_xlim(interval[0], interval[1])
        ax.legend(frameon=True, loc='best')
        ax.grid(True)
        ax.title.set_fontsize(16)

    plt.tight_layout()
    plt.show()

else:
    plt.figure(figsize=(12, 8))
    plt.plot(df['time'], df[specifCol], label=specifCol, linewidth=2, linestyle='-', color='cornflowerblue')
    plt.scatter(df['time'], df[specifCol], color='red', s=30, edgecolor='black', alpha=0.8, zorder=5)
    plt.xlabel('Time')
    plt.ylabel(specifCol)
    plt.title(f'{specifCol} over Time', fontsize=16)
    plt.xlim(interval[0], interval[1])
    plt.legend(frameon=True, loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
