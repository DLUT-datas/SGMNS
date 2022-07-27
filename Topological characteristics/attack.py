# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:33:36 2022

@author: admin
"""
import random
import networkx as nx
def attack(G):
    output_ratio=[]
    numNode = G.number_of_nodes() 
    numDelete = 49 
    maxComG = len(max(nx.connected_components(G), key = len)) 
    for i in range(0, numDelete):
        nodesDegree = G.degree() 
        nodesRank = sorted(G.nodes(), key=lambda n: nodesDegree[n]) 
        nodes = nodesRank[int(-0.02*numNode)::] 
        G.remove_nodes_from(nodes) 
        maxComI = len(max(nx.connected_components(G), key = len)) 
        ratio =  maxComI / maxComG 
        print(i, ratio)
        output_ratio.append([i,ratio])
    return output_ratio
def random_attack(G, py):
    output_aveRatio = []
    numNode = G.number_of_nodes() 
    numDelete = 49
    numRandom = 100 
    maxComG = len(max(nx.connected_components(G), key = len))
    for i in range(0, numDelete): 
        ratio= 0
        for j in range(0, numRandom):
            G = nx.read_adjlist(py)
            random.seed(j)
            nodes = random.sample(G.nodes, int((i+1)*0.02*numNode)) 
            G.remove_nodes_from(nodes)
            maxComJ = len(max(nx.connected_components(G), key = len)) 
            ratio = ratio + maxComJ / maxComG 
        aveRatio = ratio / numRandom 
        print(i, aveRatio)
        output_aveRatio.append([i,aveRatio])
    return output_aveRatio

py = r'D:\MN.txt'  
G = nx.read_adjlist(py)    
output_aveRatio = random_attack(G, py)
output_ratio = attack(G)


