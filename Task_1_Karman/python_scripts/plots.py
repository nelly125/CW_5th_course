import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('./results/result.csv', delimiter='\t', names=['F', 'R', 'G', 'res_T', 'H'])

data["x"] = data.index

fig = plt.figure(figsize=(15, 7))
plt.plot(data.F, linewidth=4, markersize=1.0, label = "F")
# plt.plot(data.R, linewidth=2, markersize=1.0, label = "R")
plt.plot(data.G, linewidth=4, markersize=1.0, label = "G")
# plt.plot(data.res_T, linewidth=2, markersize=1.0, label = "T")
plt.plot(data.H, linewidth=4, markersize=1.0, label = "H")

# data.to_csv("../results/result_" + str(data.x.max()) + ".txt")


plt.xlabel("zeta", fontsize=20)
plt.ylabel("Solution", fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
plt.grid()
plt.savefig("./results/plots/RK_4_plots_" + str(data.x.max()) + ".png")
# plt.show()
