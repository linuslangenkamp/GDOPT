import pandas as pd
import matplotlib.pyplot as plt

path = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData"
model = "analyticHypersensitive"
it = 14
specifCol = 'u0'
interval = [9985, 10000]

df = pd.read_csv(path + "/" + model + str(it) + ".csv" , sep=",")
print(df.head())
plt.figure(figsize=(10, 6))

if specifCol == None:
    for column in df.columns[1:]:
        plt.plot(df['time'], df[column], label=column)
        plt.scatter(df['time'], df[column], color='red', s=1)
        plt.scatter(df['time'], [0] * len(df['time']), color='red', s=1)
        plt.xlabel('Time')
        plt.ylabel(column)
        plt.xlim(interval[0], interval[1])
        plt.legend()
        plt.grid(True)
        plt.show() 
else:
    plt.plot(df['time'], df[specifCol], label=specifCol)
    plt.scatter(df['time'], df[specifCol], color='red', s=1)
    plt.scatter(df['time'], [0] * len(df['time']), color='red', s=1)
    plt.xlabel('Time')
    plt.ylabel(specifCol)
    plt.xlim(interval[0], interval[1])
    plt.legend()
    plt.grid(True)
    plt.show()
