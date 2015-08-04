from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


def draw_map(): 
	#This function creates a orthogonal map of the north pole and the surrounding 30' and returns the map
	#which can then be used to make contour or other types of plots with the following syntax
	#
	# longitude, latitude = my_map(X, Y)
	# my_map.contour(longitude, latitude, data)
	#
	# where X and Y is the arrays containing longitudes and latitudes respectively, in degrees.

	# First setting up orthogonal map of the northen hemisphere, then we create our map to be a the center part of that map,
	# llcrnrx keyword means latitude-lower-corner-r-x -> corner setting stuff.
	# and similar for the other keywords

	lat_0=90.; lon_0=-0.
	my_temp_map = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution= None)

	my_map = Basemap(projection='ortho',
					lon_0=lon_0,lat_0=lat_0,resolution='l',\
    				llcrnrx= - my_temp_map.urcrnrx/4 ,llcrnry= - my_temp_map.urcrnry/4. ,urcrnrx=my_temp_map.urcrnrx/4.,urcrnry=my_temp_map.urcrnry/4.)

	my_map.drawcoastlines()
	my_map.drawmeridians(np.arange(0, 360, 30), labels=[1,1,0,1])
	my_map.drawparallels(np.arange(-90, 90, 30))


	return my_map

if __name__ == '__main__':
	fig = plt.figure()

	draw_map()
	plt.show()