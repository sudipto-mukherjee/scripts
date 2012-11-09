import sys,os
sys.path.append('User/tud51931/scripts/gfzhou')
import bootstrap
from cutoffs import Cutoff

#cutoff = Cutoff(3.0,6.0,0.1)
#cutoffs = cutoff.ext(2.0,3.0,0.5)
#cutoffs = cutoff.cutoffs
cutoffs = [1.5]
print 'cutoffs:',cutoffs
numstraps = 100
Numoftrajs = 'all'
for cutoff in cutoffs:
    Assignmentsfile = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster%0.1f/lagtime50/Assignments.Fixed.h5'%cutoff 
    for i in range(0,1):
        SavePath = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster%0.1f/lagtime50/bootstrap/bs%d"%(cutoff,i)
        bootstrap.bootstrap(Assignmentsfile,Numoftrajs, numstraps, SavePath)
       
