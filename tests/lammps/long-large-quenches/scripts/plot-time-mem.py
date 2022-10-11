import  numpy as np
import matplotlib.pyplot as plt

import datetime

ini_time = 0
data     = []
with open("time-memory.log", "r") as f:

    for line in f.readlines():
        time, memory = line.split(",")

        new_time = datetime.datetime.strptime(time, "%Y-%m-%d %H:%M:%S") #2021-07-08 09:59:54

        timestamp = new_time.timestamp()
        if ini_time == 0:
            ini_time = timestamp

        delta_t = timestamp - ini_time

        data.append([int(delta_t), int(memory.strip())])

data = np.array(data)


fig, ax = plt.subplots(figsize=(10,10))

ax.plot(data[:,0] / 3600, data[:,1] / 1e3)

ax.grid()
ax.tick_params(axis="both", labelsize=16)

ax.set_xlabel("Time [h]", fontsize=16)
ax.set_ylabel("Memory Usage [GB]", fontsize=16)

fig.savefig("memory-usage.png")
