# Скрипт для построения графика зависимости пристрелосных параметров от скорости вращения жидкости
# и для получения таблицы, содержащей значения пристрелосных параметров и нормальной компоненты скорости на бесконечности

import os
import re

import matplotlib.pyplot as plt
import pandas as pd

dir_to_output = "./results/plots_to_report/"

dirs = [f for f in os.listdir('./results/first_solution/') if re.match(r'eps*', f)]
i = 0
for dir in dirs:
    path = "./results/first_solution/" + dir
    data_i = pd.read_csv(path + "/result.csv", delimiter=' ', names=['ksi_infty', 'F', 'R', 'G', 'T', 'H']).tail(1)
    shoot_data_i = pd.read_csv(path + "/shooting_parameters.csv", delimiter='\t',
                               names=['ksi', 'alpha', 'beta']).tail(1)
    s = float(dir[dir.find("_") + 1:])
    data_i["s"] = s
    shoot_data_i["s"] = s
    # print(-s)
    if (s < -0.159):
        continue;
    data_i = pd.merge(data_i, shoot_data_i, on=["s"], how='inner', )
    if (i == 0):
        data = data_i
    else:
        pass
        data = pd.concat([data, data_i])
    i += 1

dirs = [f for f in os.listdir('./results/second_solution/') if re.match(r'eps*', f)]
i = 0
for dir in dirs:
    path = "./results/second_solution/" + dir
    data_second_i = pd.read_csv(path + "/result.csv", delimiter=' ', names=['ksi_infty', 'F', 'R', 'G', 'T', 'H']).tail(
        1)
    shoot_data_second_i = pd.read_csv(path + "/shooting_parameters.csv", delimiter='\t',
                                      names=['ksi', 'alpha', 'beta']).tail(1)
    s = float(dir[dir.find("_") + 1:])
    data_second_i["s"] = s
    shoot_data_second_i["s"] = s
    # print(-s)
    data_second_i = pd.merge(data_second_i, shoot_data_second_i, on=["s"], how='inner', )
    if (i == 0):
        data_second = data_second_i
    else:
        data_second = pd.concat([data_second, data_second_i])
    i += 1
# print(shoot_data_i)
# print(data_i)
# print(data_second)

# data.to_csv("./results/all_s_parameters.csv", sep='\t', index=False)
#
data_s = data[["s", "H", "alpha", "beta"]].copy()
data_s = data_s.sort_values(by="s", ascending=False)

data_second_s = data_second[["s", "H", "alpha", "beta"]].copy()
data_second_s = data_second_s.sort_values(by="s", ascending=False)
# print(data_s)
#
# data_s.to_csv("./results/all_s_parameters.csv", sep='\t', index=False)
#
# data_s=data_s.iloc[25:]

fig = plt.figure(figsize=(15, 7))
plt.plot(-data_s.s, data_s.alpha, label="f'(0) first solution branch")
plt.plot(-data_s.s, -data_s.beta, label="g'(0) first solution branch")
plt.plot(-data_second_s.s, data_second_s.alpha, label="f'(0) second solution branch")
plt.plot(-data_second_s.s, -data_second_s.beta, label="g'(0) second solution branch")
plt.xlabel("-s")
plt.grid()
plt.legend(fontsize=12)
plt.savefig(dir_to_output + "two_solutions.png")
plt.close()

data_s=data_s.iloc[25:]
fig = plt.figure(figsize=(15, 7))
plt.plot(-data_s.s, data_s.alpha, label="f'(0) first solution branch")
plt.plot(-data_s.s, -data_s.beta, label="g'(0) first solution branch")
plt.plot(-data_second_s.s, data_second_s.alpha, label="f'(0) second solution branch")
plt.plot(-data_second_s.s, -data_second_s.beta, label="g'(0) second solution branch")
plt.xlabel("-s")
plt.legend(fontsize=12)
plt.grid()
plt.savefig(dir_to_output + "two_solutions_mini.png")
plt.close()

path_first_sol = './results/first_solution/'
dirs_first_sol = [f for f in os.listdir(path_first_sol) if re.match(r'eps*', f)]

path_second_sol = './results/second_solution/'
dirs_second_sol = [f for f in os.listdir(path_second_sol) if re.match(r'eps*', f)]

mask = list(set(dirs_second_sol) & set(dirs_first_sol))
dirs_first_sol = [path_first_sol + d + "/" for d in mask]
dirs_second_sol = [path_second_sol + d + "/" for d in mask]

both_solutions = pd.DataFrame(columns=['s', 'alpha_1', 'alpha_2', 'beta_1', 'beta_2', 'H_1', 'H_2'])

i = 0
for dir in mask:
    first_sol = pd.read_csv(path_first_sol + dir + '/' + 'result.csv', delimiter=' ',
                            names=['F', 'R', 'G', 'res_T', 'H'])
    second_sol = pd.read_csv(path_second_sol + dir + '/' + 'result.csv', delimiter=' ',
                             names=['F', 'R', 'G', 'res_T', 'H'])
    shoot_data_1 = pd.read_csv(path_first_sol + dir + '/' + "/shooting_parameters.csv", delimiter='\t',
                               names=['ksi', 'alpha', 'beta']).tail(1)
    shoot_data_2 = pd.read_csv(path_second_sol + dir + '/' + "/shooting_parameters.csv", delimiter='\t',
                               names=['ksi', 'alpha', 'beta']).tail(1)
    alpha_1 = shoot_data_1.iloc[[-1]].alpha.values[0]
    alpha_2 = shoot_data_2.iloc[[-1]].alpha.values[0]
    beta_1 = shoot_data_1.iloc[[-1]].beta.values[0]
    beta_2 = shoot_data_2.iloc[[-1]].beta.values[0]
    H_1 = first_sol.iloc[[-1]].H.values[0]
    H_2 = second_sol.iloc[[-1]].H.values[0]
    speed = float(dir[dir.find("_") + 1:])
    # print(shoot_data_1)
    # print(shoot_data_1['alpha'][0])
    both_solutions=both_solutions.append({"s": speed, "alpha_1": alpha_1, "alpha_2": alpha_2, "beta_1": beta_1, "beta_2": beta_2, "H_1": H_1, "H_2": H_2}, ignore_index=True)

# print(both_solutions)
both_solutions=both_solutions.sort_values(by="s", ascending=False)
both_solutions.to_csv(dir_to_output + "two_branches.csv", index=False, sep="\t")