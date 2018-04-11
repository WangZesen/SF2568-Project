import mpl_toolkits.mplot3d as a3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy as sp
from random import randint 


ax = a3.Axes3D(plt.figure())

f = open("positions.txt", "r");
lines = f.readlines();
count = 0

for line in lines:
	color = "r"
	if count < 250:
		color = "b"
	count += 1
	x = float(line.split(",")[0])
	y = float(line.split(",")[1])
	z = float(line.split(",")[2])
	ax.scatter(x, y, z, color = color)

f.close()

f = open("output.txt", "r");
lines = f.readlines()
fixed_color = sp.rand(3);
for line in lines:
	triangle = []
	if len(line.split(",")) < 9:
		fixed_color = sp.rand(3)
		continue
	for i in range(3):
		triangle.append((float(line.split(",")[i * 3]), float(line.split(",")[i * 3 + 1]), float(line.split(",")[i * 3 + 2])))
		facet = a3.art3d.Poly3DCollection([triangle])
		facet.set_color(colors.rgb2hex(fixed_color))
		facet.set_edgecolor('k')
		facet.set_alpha(0.8)
		ax.add_collection3d(facet)

f.close()
'''
f = open("chs.txt", "r");
lines = f.readlines()
x = []
y = []
z = []
for line in lines:
	if len(line.split(",")) < 3:
		fixed_color = sp.rand(3)
		if x != []:
			ax.plot(x, y, z, linewidth = 5)
		x = []
		y = []
		z = []
		continue
	x.append(float(line.split(",")[0]))
	y.append(float(line.split(",")[1]))
	z.append(float(line.split(",")[2]))
	
if x != []:
	ax.plot(x, y, z, linewidth = 5)

x = [18.84, 62.71]
y = [15.37, 14.20]
z = [197.82, 6.95]
ax.plot(x, y, z, linewidth = 5);


fixed_color = sp.rand(3)
x = [18.84, 62.71, 5.72]
y = [15.37, 14.20, 29.92]
z = [197.82, 6.95, 55.78]
nodes = zip(x, y, z)
facet = a3.art3d.Poly3DCollection([nodes])
facet.set_color(colors.rgb2hex(fixed_color))
facet.set_edgecolor('k')
facet.set_alpha(0.8)
ax.add_collection3d(facet)


x = [62.71, 5.72, 26.04]
y = [14.20, 29.92, 78.88]
z = [6.95, 55.78, 16.34]
nodes = zip(x, y, z)
facet = a3.art3d.Poly3DCollection([nodes])
facet.set_color(colors.rgb2hex(fixed_color))
facet.set_edgecolor('k')
facet.set_alpha(0.8)
ax.add_collection3d(facet)

x = [58.64, 62.71, 26.04]
y = [32.58, 14.20, 78.88]
z = [4.77, 6.95, 16.34]
nodes = zip(x, y, z)
facet = a3.art3d.Poly3DCollection([nodes])
facet.set_color(colors.rgb2hex(fixed_color))
facet.set_edgecolor('k')
facet.set_alpha(0.8)
ax.add_collection3d(facet)

f.close()
'''

ax.set_xlabel("X axis")
ax.set_ylabel("Y axis")
ax.set_zlabel("Z axis")


plt.show()
