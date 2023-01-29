import numpy as np 

#######################################

#Prepare the x grid 
x_vals = np.linspace(0,2.5036598e+03,178) 
y_vals = np.linspace(0,1.9600000e+03,46) 

x_new = np.tile(np.array(x_vals), (46,1)) 
np.shape(x_new) 
x_new 
np.savetxt('x_new.csv', x_new) 

#######################################

#Prepare the Y grid 
y_new = np.empty_like(x_new)
np.shape(y_new)

x_range = np.arange(0,np.shape(y_new)[0],1)
y_range = np.arange(0,np.shape(y_new)[1],1)


for i in (x_range):
    for j in (y_range):
        y_new[i,j] = y_vals[i]
        #print (i,j,y_new[i,j])


np.savetxt('y_new.csv', y_new) 

########################################

X_COORD = []
Y_COORD = []

# GET THE BATHYMETRY AS A SINGLE FILE

for i in (x_range):
    for j in (y_range):
        X_COORD.append(x_new[i,j])
        Y_COORD.append(y_new[i,j])
        #print (i,j,y_new[i,j])
