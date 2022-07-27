# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:33:36 2022

@author: admin
"""
import networkx as nx
import pandas as pd
#degree of nodes
G = nx.read_adjlist(r'D:\MN.txt')     
degree_node = G.degree()
column1=['Node','Degree'] 
degree_node = pd.DataFrame(columns=column1,data = degree_node)
degree_node.to_csv(r'D:\MN_degree_node.csv')
#degree distribution of network
degree_distribution = nx.degree_histogram(G) 
column2=['Frequency'] 
degree_distribution = pd.DataFrame(columns=column2,data = degree_distribution)
degree_distribution.to_csv(r'D:\MN_degree_distribution.csv')

#network density
network_density = nx.density(G)  
