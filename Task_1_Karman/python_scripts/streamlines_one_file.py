# Скрипт для визуализации течения для одной линиии тока в тому случае, когда у нас появляется разделяющая поверхность

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fig = plt.figure(figsize=(15, 15))
ax = plt.axes(projection='3d')

dir = './results/first_solution/'
eps = -0.1605
eps_file = "eps_" + str(eps)
streamline_file = eps_file + "/streamlines_1.0000_0.csv"
rphiz_file = eps_file + "/rphiz_1.0000_0.csv"

dir_to_output="./results/plots_to_report/"

streamlines_data = pd.read_csv(dir + streamline_file, delimiter='\t', names=['x', 'y', 'z'])
mask = streamline_file[streamline_file.find('_') + 1:]
rphiz_data = pd.read_csv(dir + rphiz_file, delimiter='\t', names=['r', 'phi', 'zeta', 'F', 'G', 'H'])
rphiz_data["v"] = rphiz_data.r * rphiz_data.r * rphiz_data.F * rphiz_data.F + \
                  rphiz_data.r * rphiz_data.r * rphiz_data.G * rphiz_data.G + \
                  rphiz_data.G * rphiz_data.G
size = len(streamlines_data)

low_coeff = 0.
high_coeff = 0.99
streamlines_data_0 = streamlines_data.iloc[int((size + 1) * low_coeff):int((size + 1) * high_coeff)]
rphiz_data_0 = rphiz_data.iloc[int((size + 1) * low_coeff):int((size + 1) * high_coeff)]

min_x = streamlines_data_0.x.min()
max_x = streamlines_data_0.x.max()
min_y = streamlines_data_0.y.min()
max_y = streamlines_data_0.y.max()
x = np.linspace (min_x, max_x, 100)
y = np.linspace (min_y, max_y, 100)
xgrid, ygrid = np.meshgrid(x, y)
mask_H = rphiz_data[rphiz_data["H"] > 0].index[-1]
zeta_plane = rphiz_data["zeta"][mask_H]
z = xgrid * 0 + zeta_plane

ax.plot(streamlines_data_0.x, streamlines_data_0.y, streamlines_data_0.z, linewidth=1, alpha=0.3)
ax.scatter3D(streamlines_data_0.x, streamlines_data_0.y, streamlines_data_0.z, s=2, c=rphiz_data_0.v)
ax.plot_surface(xgrid, ygrid, z, alpha=0.2)

ax.set_xlabel("x", fontsize=17)
ax.set_ylabel("y", fontsize=17)
ax.set_zlabel("z", fontsize=17)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='z', labelsize=15)

plt.savefig(dir_to_output + "streamlines_" + eps_file + ".png")
# plt.show()
plt.close()

fig = plt.figure(figsize=(15, 15))
ax = plt.axes(projection='3d')


streamlines_data_1 = streamlines_data.iloc[int((size + 1) * low_coeff):mask_H]
rphiz_data_1 = rphiz_data.iloc[int((size + 1) * low_coeff):mask_H]
min_x = streamlines_data_1.x.min()
max_x = streamlines_data_1.x.max()
min_y = streamlines_data_1.y.min()
max_y = streamlines_data_1.y.max()
x = np.linspace (min_x, max_x, 100)
y = np.linspace (min_y, max_y, 100)
xgrid, ygrid = np.meshgrid(x, y)
mask_H = rphiz_data_1[rphiz_data_1["H"] > 0].index[-1]
zeta_plane = rphiz_data_1["zeta"][mask_H]
z = xgrid * 0 + zeta_plane

ax.plot(streamlines_data_1.x, streamlines_data_1.y, streamlines_data_1.z, linewidth=1, alpha=0.3)
ax.scatter3D(streamlines_data_1.x, streamlines_data_1.y, streamlines_data_1.z, s=2, c=rphiz_data_1.v)
ax.plot_surface(xgrid, ygrid, z, alpha=0.2)

ax.set_xlabel("x", fontsize=17)
ax.set_ylabel("y", fontsize=17)
ax.set_zlabel("z", fontsize=17)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='z', labelsize=15)

plt.savefig(dir_to_output + "streamlines_positive_" + eps_file + ".png")
# plt.show()
plt.close()


fig = plt.figure(figsize=(15, 15))
ax = plt.axes(projection='3d')

streamlines_data_1 = streamlines_data.iloc[mask_H:int((size + 1) * high_coeff)]
rphiz_data_1 = rphiz_data.iloc[mask_H:int((size + 1) * high_coeff)]
min_x = streamlines_data_1.x.min()
max_x = streamlines_data_1.x.max()
min_y = streamlines_data_1.y.min()
max_y = streamlines_data_1.y.max()
x = np.linspace (min_x, max_x, 100)
y = np.linspace (min_y, max_y, 100)
xgrid, ygrid = np.meshgrid(x, y)
mask_H = rphiz_data_1[rphiz_data_1["H"] > 0].index[-1]
zeta_plane = rphiz_data_1["zeta"][mask_H]
z = xgrid * 0 + zeta_plane

ax.plot(streamlines_data_1.x, streamlines_data_1.y, streamlines_data_1.z, linewidth=1, alpha=0.3)
ax.scatter3D(streamlines_data_1.x, streamlines_data_1.y, streamlines_data_1.z, s=2, c=rphiz_data_1.v)
ax.plot_surface(xgrid, ygrid, z, alpha=0.2)

ax.set_xlabel("x", fontsize=17)
ax.set_ylabel("y", fontsize=17)
ax.set_zlabel("z", fontsize=17)
ax.tick_params(axis='x', labelsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.tick_params(axis='z', labelsize=15)

plt.savefig(dir_to_output + "streamlines_negative_" + eps_file + ".png")
# plt.show()
plt.close()
