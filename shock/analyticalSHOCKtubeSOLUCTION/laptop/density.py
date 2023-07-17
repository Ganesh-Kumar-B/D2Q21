import matplotlib.pyplot as plt 
import numpy as np

#lets read whole file into numpy 2d array
data=np.loadtxt(fname='density.dat',delimiter = ' ')

#print(data)

#now let's split them into 2 different arrays x and y arrays
t=[]#time
y1=[]

for i in range(len(data)):
		for j in range(6):
	 		if j==0:
	 			t.append(data[i][j])
	 		if j==1:
	 			y1.append(data[i][j])




plt.plot(t,y1)
plt.title('density')



plt.show()
