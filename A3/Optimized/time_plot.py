import matplotlib.pyplot as plt
import pandas as pd;

times = [0.001, 0.002, 0.003, 0.003, 0.005, 0.008, 0.030, 0.073, 0.141, 0.215, 0.625, 1.107, 2.497, \
    6.441, 12.728, 20.951]
N_s = [10, 30, 50, 70, 90, 150, 300, 500, 700, 900, 1500, 2000, 3000, 5000, 7000, 9000]

data = [[N, time] for N, time in zip(N_s, times)]
df = pd.DataFrame(data, columns = ['N', 'Time (s)'])
df.plot(x='N', y='Time (s)', kind='line', title='Time as a function of N, 100 time steps', \
     ylabel='Time (s)', legend=False)
plt.show()
