import matplotlib.pyplot as plt 
import numpy as np

#lets read whole file into numpy 2d array
data=np.loadtxt(fname='properties.dat',delimiter = ' ')
#print(data)

#now let's split them into 2 different arrays x and y arrays
t=[]#time
y1=[]
y2=[]
y3=[]
y4=[]

for i in range(len(data)):
		for j in range(6):
	 		if j==0:
	 			t.append(data[i][j])
	 		if j==1:
	 			y1.append(data[i][j])
	 		if j==2:
	 			y2.append(data[i][j])
	 		if j==3:
	 			y3.append(data[i][j])
	 		if j==4:
	 			y4.append(data[i][j])



plt.subplot(2,2,1)
plt.plot(t,y1)
plt.title('velocity')


plt.subplot(2,2,2)
plt.plot(t,y2)
plt.title('sound speed')

plt.subplot(2,2,3)
plt.plot(t,y3)
plt.title('pressure')

plt.subplot(2,2,4)
plt.plot(t,y4)
plt.title('density')

# fig, axs = plt.subplots(2, 2)	
	 		
# axs[0, 0].plot(t, y1)
# axs[0, 0].set_title('y1')
# # axs[0, 0].set_xlim(0,5)
# # axs[0, 0].set_ylim(0,2)

# axs[0, 1].plot(t, y2)
# axs[0, 1].set_title('y2')
# # axs[0, 1].set_xlim(0,5)
# # axs[0, 1].set_ylim(0,2)

# axs[1, 0].plot(t, y3)
# axs[1, 0].set_title('y3')
# # axs[1, 0].set_xlim(0,10)
# # axs[0, 2].set_ylim(0,2)

# axs[1, 1].plot(t, y4)
# axs[1, 1].set_title('y4')
# # axs[1, 1].set_xlim(0,10)
# # axs[0, 3].set_ylim(0,2)




# axs[3, 1].plot(t, y8)
# axs[3, 1].set_title('y8')
# axs[3, 1].set_xlim(0,5)
# axs[1, 3].set_ylim(0,2)


#print(x)
#print(y)
#now let's plot
# plt.plot(x,y,'r',x,z,'b',x,p,'g')

# plt.xlabel('time')
# plt.ylabel('concentration')
# plt.title('q9')
# plt.legend(['y1','y2','B'])
plt.show()
