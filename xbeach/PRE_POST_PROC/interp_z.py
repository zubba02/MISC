import numpy as np  #import necessary packages
import matplotlib.pyplot as plt
import scipy.interpolate

x_new = np.genfromtxt('x_new.csv') #load our data
y_new = np.genfromtxt('y_new.csv')


depth = np.genfromtxt('bed.dep') #load their data as reference points
x_old = np.genfromtxt('x.grd')
y_old = np.genfromtxt('y.grd')


a = np.shape(x_old) #define the shape of their data
#print(a[0]) #should give (46,178)


x_all = [] #empty array to save new x values
for i in range(0, a[0], 1):
    for j in range(0, a[1], 1): #nested loop to go through 178 elements in each of the 46 indexes
        print('this is {}th element for x in index {}'.format (j, i))
        x_all.append(x_new[i][j])


y_all = [] #empty array to save new y values
for i in range(0, a[0], 1):
    for j in range(0, a[1], 1): #nested loop to go through 178 elements in each of the 46 indexes
        print('this is {}th element for y'.format (i))
        y_all.append(y_new[i][j])


z_all = [] #empty array to save new z values
for i in range(0, a[0], 1):
    for j in range(0, a[1], 1): #nested loop to go through 178 elements in each of the 46 indexes
        print('this is {}th element for z'.format (i))
        z_all.append(depth[i][j])


interpolator_mv = scipy.interpolate.NearestNDInterpolator((x_all, y_all), z_all)

new_vals = [] #empty array to save new values
for i in range(0, a[0], 1):
    for j in range(0, a[1], 1): #nested loop to go through 178 elements in each of the 46 indexes
        new_val = interpolator_mv(x_new[i][j],y_new[i][j]) #interpolates our data
        print (new_val)
        new_vals.append(new_val) #saves new z


new_vals = np.array(new_vals) #make into array
print (np.shape(new_vals))
new_vals = np.reshape(new_vals,(46,178)) #reshape to format being used (46, 178)

np.savetxt('new_depth.csv',new_vals, fmt='%1.7e', delimiter='  ')



'''x_min,x_max = datax.min(),datax.max()
y_min,y_max = datay.min(),datay.max()

x_range = np.linspace(x_min,x_min,1500)
y_range = np.linspace(y_min,y_max,1500)

gridxy = np.mgrid[x_min:x_max:1500j, y_min:y_max:1500j].T
gridxy = np.reshape(gridxy, (1500*1500, 2))

x_coord = []
y_coord = []

for i,xy in enumerate(gridxy):
    x_coord.append(xy[0])
    y_coord.append(xy[1])

len_x = np.shape(x_coord)'''

'''coords = np.column_stack([datax, datay])

bathy = []
abS = []

for i,xy in enumerate (coords):
    print (i,len_x,interpolator_mv(xy[0], xy[1]))
    bathy.append(interpolator_mv(xy[0], xy[1]))
    abS.append(np.abs(interpolator_mv(xy[0], xy[1])))

all_data_mesh = np.column_stack([abS, x_coord, y_coord])

all_data = np.column_stack([x_coord, y_coord, bathy])
np.savetxt('clippedbathy_1500_1500.csv',bathy, fmt='%f')
np.savetxt('clippedbathy_with_coord_1500_1500.csv',all_data, fmt='%f', delimiter=',')

np.savetxt('clippedbathy_with_coord_1500_1500_mesh.csv',all_data_mesh, fmt='%f', delimiter=',')'''

