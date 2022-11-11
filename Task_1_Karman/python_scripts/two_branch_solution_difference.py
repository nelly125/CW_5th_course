# Скрипт для построения эпюр скоростей в области неоднозначного решения

import matplotlib.pyplot as plt
import pandas as pd

first_sol = pd.read_csv("./results/first_solution/eps_-0.1575/" + 'result.csv', delimiter=' ',
                        names=['F', 'R', 'G', 'res_T', 'H'])
second_sol = pd.read_csv("./results/second_solution/eps_-0.1575/" + 'result.csv', delimiter=' ',
                         names=['F', 'R', 'G', 'res_T', 'H'])

dir_to_output="./results/plots_to_report/"

fig, ax = plt.subplots(1, 1, figsize=(20, 10))
ax.plot(first_sol.F, label="first solution branch", color="r")
ax.plot(second_sol.F, label="second solution branch", color="b")
ax.set_title("F", fontsize=16)
ax.set_xlabel("zeta", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
ax.grid()
plt.savefig(dir_to_output + "F"+"_difference.png")
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(20, 10))
ax.plot(first_sol.G, label="first solution branch", color="g")
ax.plot(second_sol.G, label="second solution branch", color="y")
ax.set_title("G", fontsize=16)
ax.set_xlabel("zeta", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
ax.grid()
plt.savefig(dir_to_output + "G"+"_difference.png")
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(20, 10))
ax.plot(first_sol.H, label="first solution branch", color="orange")
ax.plot(second_sol.H, label="second solution branch", color="pink")
ax.set_title("H", fontsize=16)
ax.set_xlabel("zeta", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)
ax.grid()
plt.savefig(dir_to_output + "H"+"_difference.png")
plt.close()
