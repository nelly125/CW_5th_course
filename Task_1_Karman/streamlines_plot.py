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


# ax.scatter3D(data_0.x, data_0.y, data_0.z, c=data_0.z, cmap='Greens')
# ax.scatter3D(data_1.x, data_1.y, data_1.z, c=data_1.z, cmap='Oranges')
#
# plt.show()

# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# import numpy as np
#
# x = data_0.x.values
# y = data_0.y.values
# z = data_0.z.values
#
# u = np.zeros_like(x)
# v = np.zeros_like(y)
# w =  np.zeros_like(z)
#
#
#
# fig = go.Figure(go.Streamtube(x, y, z, u, v, w))
#
# fig.add_trace(go.Streamtube(x=x, y=y, z=z, u=u, v=v, w=w), 1, 1)
# fig.add_trace(go.Streamtube(x=x, y=y, z=z, u=w, v=v, w=u), 1, 2)
# fig.add_trace(go.Streamtube(x=x, y=y, z=z, u=u, v=w, w=v), 1, 3)
#
# fig.update_layout(scene_camera_eye=dict(x=2, y=2, z=2),
#                   scene2_camera_eye=dict(x=2, y=2, z=2),
#                   scene3_camera_eye=dict(x=2, y=2, z=2))
fig.show()