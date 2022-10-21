import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import re

fig = plt.figure(figsize=(15, 15))
ax = plt.axes(projection='3d')

# data_0 = pd.read_csv('./results/streamlines0.txt', delimiter='\t', names=['x', 'y', 'z'])
# data_1 = pd.read_csv('./results/streamlines1.txt', delimiter='\t', names=['x', 'y', 'z'])

files = [f for f in os.listdir('./results/') if re.match(r'streamlines*', f)]
# print(files)
for file in files:
    data = pd.read_csv('./results/' + file, delimiter='\t', names=['x', 'y', 'z'])
    data = data.iloc[:]
    ax.scatter3D(data.x, data.y, data.z, c=data.z)

plt.savefig("./results/plots/streamlines3d.png")


fig.show()