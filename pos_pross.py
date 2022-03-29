# -*- coding: utf-8 -*-
import numpy as np
from odbAccess import openOdb
from abaqusConstants import *



def nodalAveraged(odbInstance,Frame,StressType,timestep):
    #element_nodal - é o ponto de integracao(setado no pre)

    # Get number of nodes
    Field = Frame[timestep].fieldOutputs['S'].getSubset(position = ELEMENT_NODAL).getScalarField(invariant = StressType)
    Values = Field.bulkDataBlocks[0].data
    #print(Values)
    #print(len(Values))
    NodeLabels = Field.bulkDataBlocks[0].nodeLabels
    #print(NodeLabels)
    #print(len(NodeLabels))
    
    #caso tenha varias parts
    for i in range(len(Field.bulkDataBlocks)-1):
        Values = np.vstack((Values,Field.bulkDataBlocks[i+1].data))
        
        NodeLabels = np.hstack((NodeLabels,Field.bulkDataBlocks[i+1].nodeLabels))
    
    # Nodes are shared across multiple elements.  Get unique node labels.
   
    NodeLabels_unique, unq_idx = np.unique(NodeLabels, return_inverse=True)
    
    #print(NumNodes)

        # Calculate nodal averaged stresses at timestep
    Values_Averaged=np.zeros((NodeLabels_unique.size,Values.shape[1]))
    unq_counts = np.bincount(unq_idx)
    for i in xrange(0,Values.shape[1]):
       ValuesTemp = [item[i] for item in Values]
       unq_sum = np.bincount(unq_idx, weights=ValuesTemp)
       Values_Averaged[:,i] = unq_sum / unq_counts
    
    return NodeLabels_unique, Values_Averaged

#==============================================================================
# RUN THE PROGRAM
#==============================================================================
filename = 'Carrgamento_min_max'

#
# LOAD ABAQUS SOLUTION DATA
#-------------------------------------------------------------------
odb = openOdb(filename+'.odb',readOnly=True)

# Get Instance
allInstances = (odb.rootAssembly.instances.keys())
odbInstance = odb.rootAssembly.instances[allInstances[-1]]
#
# PROCESS RESULTS
#-------------------------------------------------------------------

# Retrieve nodal averaged stresses at steady-state solution
# para carga minima
timestep = -1
Frame = odb.steps['load_min'].frames

nodeNum, pressure = nodalAveraged(odbInstance,Frame,PRESS,timestep)
nodeNum, vonMises = nodalAveraged(odbInstance,Frame,MISES,timestep)

# Create a signed von Mises stress
vonMisesSigned_1 = np.sign(-1.*pressure)*vonMises

# para carga maxima
Frame = odb.steps['load_max'].frames
timestep = -1

nodeNum, pressure = nodalAveraged(odbInstance,Frame,PRESS,timestep)
nodeNum, vonMises = nodalAveraged(odbInstance,Frame,MISES,timestep)

# Create a signed von Mises stress
vonMisesSigned_2 = np.sign(-1.*pressure)*vonMises

#calculando Sigma_a e Sigma_m
sigma_a = np.zeros((len(vonMisesSigned_1),1))
sigma_m = np.zeros((len(vonMisesSigned_1),1))
Safety_factor = np.zeros((len(vonMisesSigned_1),1))
S_ult = 380 # tensão ultimate
Se = 0.5*S_ult # tensão de fadiga

for i in range(len(vonMisesSigned_1)):

    if vonMisesSigned_1[i]>vonMisesSigned_2[i]:
        # vonMissesSigned_1 é max
        sigma_a[i] = 0.5*(vonMisesSigned_1[i] - vonMisesSigned_2[i])
        sigma_m[i] = 0.5*(vonMisesSigned_1[i] + vonMisesSigned_2[i])
        

    else:
        # vonMissesSigned_2 é max
        sigma_a[i] = 0.5*(vonMisesSigned_2[i] - vonMisesSigned_1[i])
        sigma_m[i] = 0.5*(vonMisesSigned_2[i] + vonMisesSigned_1[i])

    Safety_factor[i] = 1/(sigma_a[i]/Se+sigma_m[i]/S_ult)
    print(float(sigma_a[i]), float(sigma_m[i]), i, float(Safety_factor[i]), float(vonMisesSigned_1[i]), float(vonMisesSigned_2[i]))

"""
sf = Safety_factor.flatten()
nodeNum = nodeNum.flatten()

# Get nodal coordinates
nodeList = Frame[0].fieldOutputs['S'].values[0].instance.nodes
nodeCoord = np.zeros((len(nodeList),4))

for item in range(len(nodeList)):
    nodeCoord[item,0] = nodeList[item].label
    nodeCoord[item,1] = nodeList[item].coordinates[0]
    nodeCoord[item,2] = nodeList[item].coordinates[1]
    nodeCoord[item,3] = nodeList[item].coordinates[2]

odb.close()

#print(Safety_factor)
odb = openOdb(filename+'.odb',readOnly=False)

# Get Instance
allInstances = (odb.rootAssembly.instances.keys())
odbInstance = odb.rootAssembly.instances[allInstances[-1]]

# save Odb


sf =sf.transpose()
nodeNum =nodeNum.transpose()

vMNodes  = np.ascontiguousarray(nodeNum, dtype=np.int32)

sf    = np.ascontiguousarray(np.reshape(sf,(-1,1)), dtype=np.float32)

print(len(vMNodes))
print(len(sf))

newFieldOutput = odb.steps['load_min'].frames[-1].FieldOutput(name = 'SF', description = 'Safety factor', type = SCALAR)
#print(vMNodes)
#print(sf.tolist())
newFieldOutput.addData(position = NODAL, instance = odbInstance, labels = vMNodes, data = sf.tolist())



"""
odb.save()
odb.close()