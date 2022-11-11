# Скрипт для построения линий тока во всех директориях

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fig = plt.figure(figsize=(15, 15))
ax = plt.axes(projection='3d')

path = './results/first_solution/'
# path = './results/second_solution/'

dirs = [f for f in os.listdir(path) if re.match(r'eps*', f)]
dirs = [path + d + "/" for d in dirs]
# dirs = [path + "eps_-0.1605/"]
for dir in dirs:
    streamline_files = [f for f in os.listdir(dir) if re.match(r'streamlines_*', f)]
    rphiz_files = [f for f in os.listdir(dir) if re.match(r'rphiz*', f)]
    # plots = [f for f in os.listdir(dir + "plots/") if f.rfind('.png') != -1]
    # for plot in plots:
    #     os.remove(dir + "plots/" + plot)

    fig1 = plt.figure(1, figsize=(15, 15))
    ax = plt.axes(projection='3d')
    N = len(streamline_files)
    low_coeff = 0.71
    high_coeff = 0.72
    xgrid = 0
    ygrid = 0
    z = 0

    for streamline_file in streamline_files:
        streamlines_data = pd.read_csv(dir + streamline_file, delimiter='\t', names=['x', 'y', 'z'])
        mask = streamline_file[streamline_file.find('_') + 1:]
        rphiz_file = [f for f in rphiz_files if f.find(mask) != -1][0]

        rphiz_data = pd.read_csv(dir + rphiz_file, delimiter='\t', names=['r', 'phi', 'zeta', 'F', 'G', 'H'])
        rphiz_data["v"] = rphiz_data.r * rphiz_data.r * rphiz_data.F * rphiz_data.F + \
                          rphiz_data.r * rphiz_data.r * rphiz_data.G * rphiz_data.G + \
                          rphiz_data.G * rphiz_data.G
        size = len(streamlines_data)

        mask_H = rphiz_data[rphiz_data["H"] > 0]

        if (len(mask_H) ):
            mask_H_index = mask_H.index[-1]
            zeta_plane = rphiz_data["zeta"][mask_H_index]
            z = xgrid * 0 + zeta_plane
        else:
            mask_H_index = int((size + 1) * high_coeff)


        # streamlines_data =streamlines_data.iloc[int((size + 1) * low_coeff):mask_H_index]
        # rphiz_data = rphiz_data.iloc[int((size + 1) * low_coeff):mask_H_index]

        # streamlines_data =streamlines_data.iloc[int((size + 1) * low_coeff):mask_H_index]
        # rphiz_data = rphiz_data.iloc[int((size + 1) * low_coeff):mask_H_index]


        streamlines_data =streamlines_data.iloc[int((size + 1) * low_coeff):int((size + 1) * high_coeff)]
        rphiz_data = rphiz_data.iloc[int((size + 1) * low_coeff):int((size + 1) * high_coeff)]

        min_x = streamlines_data.x.min()
        max_x = streamlines_data.x.max()
        min_y = streamlines_data.y.min()
        max_y = streamlines_data.y.max()
        x = np.linspace (min_x, max_x, 100)
        y = np.linspace (min_y, max_y, 100)
        xgrid, ygrid = np.meshgrid(x, y)
        # if (N == 1):
        # ax.scatter3D(streamlines_data.x, streamlines_data.y, streamlines_data.z, linewidth=3, c=rphiz_data.v)
        if (N == 1):
            ax.scatter3D(streamlines_data.x, streamlines_data.y, streamlines_data.z, s=2, linewidth=2, c=rphiz_data.v)
            ax.plot(streamlines_data.x, streamlines_data.y, streamlines_data.z, linewidth=1, alpha=0.3)
        else:
            # ax.scatter3D(streamlines_data.x, streamlines_data.y, streamlines_data.z, s=2, linewidth=2, c=rphiz_data.v)
            ax.plot(streamlines_data.x, streamlines_data.y, streamlines_data.z, linewidth=2, alpha=0.6, color="green")

    # if (len(mask_H)):
    #     ax.plot_surface(xgrid, ygrid, z, alpha=0.2)

    if (N):

        fig1.savefig(dir + "/plots/" + "plot_streamlines_3d_" + str(N) + ".png")
        plt.close()

    if (N < 5):

        fig = plt.figure(figsize=(15, 15))
        ax = plt.axes()
        ax.set_xlabel("r", fontsize=30)
        ax.set_ylabel("zeta", fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        for streamline_file in streamline_files:
            streamlines_data = pd.read_csv(dir + streamline_file, delimiter='\t', names=['x', 'y', 'z'])
            mask = streamline_file[streamline_file.find('_') + 1:]
            rphiz_file = [f for f in rphiz_files if f.find(mask) != -1][0]
            r0 = streamline_file[streamline_file.find('_') + 1:streamline_file.rfind('_')]

            rphiz_data = pd.read_csv(dir + rphiz_file, delimiter='\t', names=['r', 'phi', 'zeta', 'F', 'G', 'H'])
            rphiz_data["v"] = rphiz_data.r * rphiz_data.r * rphiz_data.F * rphiz_data.F + \
                              rphiz_data.r * rphiz_data.r * rphiz_data.G * rphiz_data.G + \
                              rphiz_data.G * rphiz_data.G
            size = len(streamlines_data)

            streamlines_data =streamlines_data.iloc[int((size + 1) * low_coeff):int((size + 1) * high_coeff)]
            rphiz_data = rphiz_data.iloc[int((size + 1) * low_coeff):int((size + 1) * high_coeff)]
            # rphiz_data.append(-rphiz_data.r, rphiz_data.zeta)

            plt.plot(rphiz_data.r, rphiz_data.zeta, linewidth=2, label=r0)
            plt.plot(-rphiz_data.r, rphiz_data.zeta, linewidth=2, label=r0)

        plt.legend()
        plt.savefig(dir + "/plots/" + "generating_lines" + str(N) + ".png")
        plt.close()

    for streamline_file in streamline_files:
        try:
            os.remove(dir + streamline_file)
        except:
            print("Error while deleting file ")

    for rphiz_file in rphiz_files:
        try:
            os.remove(dir + rphiz_file)
        except:
            print("Error while deleting file ")
