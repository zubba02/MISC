import xarray as xr
import cmocean
import matplotlib.pyplot as plt
import numpy as np


#read the ncdf file using xarray
ds = xr.open_dataset('xboutput.nc')


#get the summary of the file
#print("get the summary of the file")
#print (ds)


#get the variable as numpy array
zsnp = np.array(ds['v'])
#print(zsnp)


#print("get the shape of the numpy array")
print(np.shape(zsnp))
a = np.shape(zsnp)[0]
print ('The number of timesteps is {}'.format(a))


#get the snap shot at the third time output 
#pp = zsnp[3,:,:]
#print("get the snap shot at the third time output")
#print(pp)


#plot the image
#plt.imsave('test.png', pp, cmap='cmo.balance', aspect='auto', interpolation='quadric', vmin=0.0, vmax=0.4)
#plt.imsave('test.png', pp, cmap='cmo.balance', vmin=0.0, vmax=0.4)
#print("the plotted image will be saved as test.png")


#plt.imshow(pp, cmap='cmo.balance', aspect='auto', interpolation='quadric', vmin=0.0, vmax=0.4)
#plt.show()


#select = ['zsnp[{},:,:]'.format(i) for i in range (0,a)]
#print(select)

#print(np.shape(select))
#b = np.shape(select)
#for k in 

for k in range(0,a):
	plt.figure(figsize=(20,10)) 
	print (k)
	pp = zsnp[k,:,:]
	print(pp)
	plt.imshow(pp, cmap='cmo.balance', aspect='auto', interpolation='quadric', vmin=0.0, vmax=0.4)
	plt.savefig('TEST {}.png'.format(k), dpi=300)
	#plt.imsave( pp, cmap='cmo.balance', vmin=0.0, vmax=0.4, dpi=1000)




