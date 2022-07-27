# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 09:01:10 2022

@author: admin
"""
import networkx as nx
import numpy as np
#Common_Neighbors------------------------------
def Common_Neighbors( matrix ):
    sim = np.dot( matrix, matrix )
    return sim
#Resource Allocation---------------------------
def RA( matrix ):
    add_row = matrix.sum( axis = 1 )          
    add_row_matrix = np.tile( add_row, (( matrix.shape )[0]) )
    sim = adjacent_matrix / add_row_matrix    
    sim = matrix * sim
    return sim 
#Jaccard---------------------------------------
def Jaccard( matrix ):
    sim = matrix * matrix
    sim1 = sim.copy()                        
    sim1[ np.nonzero( sim1 ) ] = 1            
    deg_row = matrix.sum( axis = 0 )         
    deg_row_matrix = np.tile( deg_row, ( (matrix.shape)[0], 1) )
    deg_row_matrix = np.multiply( deg_row_matrix, sim1 )
    deg_row_matrix = np.triu( deg_row_matrix ) + np.triu( (deg_row_matrix.T) )
    sim = sim / ( np.multiply( deg_row_matrix, sim1 ) - sim )
    return sim

G = nx.read_adjlist( r'D:\MN.txt' )       
adjacent_matrix = nx.to_numpy_matrix( G )   
numNode = G.number_of_nodes() 

#node similarity
CN_sim = np.triu(Common_Neighbors( adjacent_matrix ), k=1) 
Jac_sim = np.triu(Jaccard( adjacent_matrix ), k=1) 
RA_sim = np.triu(RA( adjacent_matrix ), k=1)

#Mean values of the node similarity
RA_sum = np.sum(np.sum(RA_sim, axis=0), axis=0)
ave_node_sim = RA_sum / ((numNode*(numNode-1))/2) #平均节点相似性计算
print(ave_node_sim)
