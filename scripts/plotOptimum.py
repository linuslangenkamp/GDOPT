import pandas as pd
import matplotlib.pyplot as plt

it = 3
specifCol = 'u0'

df = pd.read_csv("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData/batchReactorRefinement/BatchReactor" + str(it) + ".csv" , sep=",")
print(df.head())
plt.figure(figsize=(10, 6))

if specifCol == None:
    for column in df.columns[1:]:
        plt.plot(df['time'], df[column], label=column)
        plt.scatter(df['time'], df[column], color='red', s=10)
        plt.xlabel('Time')
        plt.ylabel(column)
        plt.legend()
        plt.grid(True)
        plt.show()
else:
    plt.plot(df['time'], df[specifCol], label=specifCol)
    plt.scatter(df['time'], df[specifCol], color='red', s=10)
    plt.xlabel('Time')
    plt.ylabel(specifCol)
    plt.xlim(0, 1)
    plt.legend()
    plt.grid(True)
    plt.show()
