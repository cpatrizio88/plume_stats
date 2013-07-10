import numpy as np
import cPickle
import glob
import matplotlib.pyplot as plt

filenames = glob.glob('/home/cpatrizi/repos/cloudtracker/pkl/clusters*.pkl')
filenames.sort()

dt = 2 #timestep in min
clusters_list = []
maxid = 0

for fname in filenames:
    clusters = cPickle.load(open(fname, 'rb'))
    clusters_list.append(clusters)
    ids = clusters.keys()
    maxid_tmp = max(ids)
    if maxid_tmp > maxid:
        maxid = maxid_tmp

print 'max. cluster id: ', maxid
#array to keep track of lifetimes of clusters
lifetimes = np.zeros(maxid+1)

#clusters_list contains the clusters at each time step
#clusters is a dictionary, with cluster id as key,
#and a dictionary as value (containing information about the cluster)
for clusters in clusters_list:
    #iterate over each cluster at the current time
    #and increment the lifetime of the cluster if it has a core
    for id, cluster in clusters.iteritems():
        if len(cluster['core']):
           lifetimes[id] = lifetimes[id] + dt


#filter out noise from lifetimes array
noise_filter = np.where(lifetimes > dt)
lifetimes = lifetimes[noise_filter]

mean_lt = np.mean(lifetimes)
print "mean plume lifetime (min): ", mean_lt

bins = np.arange(21)
plt.hist(lifetimes, bins)
plt.show()

    
    
    







