import matplotlib.pyplot as plt

path = 'C:/Users/denis/Documents/repos/flow_and_transport_2d/build/Release/outC.txt'

f = open(path, 'r')
t, C = [], []
for line in f.readlines():
    fields = line.split(' ')
    t.append(float(fields[0]))
    C.append(float(fields[1]))
f.close()

plt.plot(t, C)