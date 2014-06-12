import cPickle
import matplotlib.pyplot as plt

stats = cPickle.load(open('./pkl/stats_test.pkl', 'rb'))

areas = stats['plume_meanarea']

plt.hist(areas/(1000*1000), bins=50)
plt.xlabel('sub-cloud plume cross sectional area, averaged over lifetime and height (km^2)')
plt.ylabel('number of plumes')
plt.show()


