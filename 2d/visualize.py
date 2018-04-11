import matplotlib.pyplot as plt

f = open("positions.txt", "r")
lines = f.readlines()
points = []
x = []
y = []
for line in lines:
	x1 = float(line.split(", ")[0])
	y1 = float(line.split(", ")[1])
	points.append([x1, y1])
	x.append(x1)
	y.append(y1)
	plt.plot([x1], [y1], 'ro', markersize = 1)
f.close()
# plt.plot(x, y, 'r')

f = open("log.txt", "r")
lines = f.readlines()

index = 0
x = []
y = []

for i in range(2):
	index += 1
	a = int(lines[index].split(",")[0])
	b = int(lines[index].split(",")[1])
	for j in range(b - a + 1):
		index += 1;
		curIndex = int(lines[index])
		x.append(points[curIndex][0])
		y.append(points[curIndex][1])
	plt.plot(x, y)
	x = []
	y = []

# plt.plot([points[82][0]], [points[82][1]], 'ro')


plt.show()
