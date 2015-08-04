from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


def draw_map(data): 
	# make sure the value of resolution is a lowercase L,
	#  for 'low', not a numeral 1
	plt.figure()
	my_map = Basemap(projection='ortho', lat_0=90, lon_0=-100,
	              resolution='l', area_thresh=1000.0)
	 
	my_map.drawcoastlines()
	my_map.drawmeridians(np.arange(0, 360, 30))
	my_map.drawparallels(np.arange(-90, 90, 30))
	# my_map.drawcountries()
	# my_map.fillcontinents(color='coral')
	# my_map.drawmapboundary()

	plt.show()