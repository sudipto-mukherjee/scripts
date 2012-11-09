import sys, os
import numpy as np
import scipy.io
from msmbuilder import arglib
from msmbuilder import Serializer
from msmbuilder import MSMLib

Assignments=Serializer.LoadData("%s"%(sys.argv[1]))
NumStates = max(Assignments.flatten()) + 1
LagTime = sys.argv[2]
Counts = MSMLib.GetCountMatrixFromAssignments(Assignments, NumStates, LagTime=LagTime, Slide=True)
scipy.io.mmwrite('%s'%(sys.argv[3]), Counts)
