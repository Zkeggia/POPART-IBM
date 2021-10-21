#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 14:08:50 2021

@author: fra
"""

import graph_tool.all as gt
import igraph


#filename='/home/fra/Projects/POPART-IBM/PARAMS_COMMUNITY1_ACCEPTED/Output/Partnership_network_2001_CL01_Za_B_V1.2_patch0_Rand10_Run1_PCseed0_0.csv'


giant_size=[]
average_degree=[]
average_var=[]
bicomponentsize=[]
average_giant_degree=[]

for index in range(1,500):
    filename = '/home/fra/Projects/POPART-IBM/PARAMS_COMMUNITY1_ACCEPTED/Output/Partnership_network_2019_Za_B_V1.2_Rand10_Run%d_PCseed0_0.csv'%index
    #graph= gt.load_graph_from_csv(filename,eprop_types=['int','float'], eprop_names=['patch','duration'], 
    #                            csv_options={'delimiter':' '}, )
    
    
    #_,part=gt.is_bipartite(graph,partition=True)
    
    
    #gt.graph_draw(graph, vertex_fill_color= part, output_size=(1920,1080), output='network.pdf')
    
    import networkx as nx
    import numpy as np
    import scipy.sparse as sparse
    import pandas as pd
    
    df = pd.read_csv(filename,sep=',')
    
    df=df[df.id_patch==0]
    
    Graphtype = nx.Graph()
    G = nx.from_pandas_edgelist(df, source='n_id', target='partner_id', create_using=Graphtype)

    
    
    
    #M = nx.adjacency_matrix(G)
    #M = sparse.csr_matrix.todense(M)
    
    
    s = sorted(G.degree, key=lambda x: x[1], reverse=True)
    degrees = sorted([x[1] for x in s])
    
    print('Full network, average and variance')
    average_degree.append(np.mean(degrees))
    average_var.append(np.var(degrees))

    
    
    
    
    
    #Get connected component
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G0 = G.subgraph(Gcc[0])
    #M0= nx.adjacency_matrix(G0)
    
    #M0 = sparse.csr_matrix.todense(M0)
    Gc = max(nx.biconnected_components(G), key=len)
    
    s = sorted(G0.degree, key=lambda x: x[1], reverse=True)
    degrees = sorted([x[1] for x in s])
    print('Giant component')

    
    #print('average shortest path length')
    #print(nx.average_shortest_path_length(G0))
    
    import matplotlib.pyplot as plt
    
    import numpy as np
    def plot_degree_dist(G,label):
        plt.yscale('log')
        degree_hist = nx.degree_histogram(G) 
        degree_hist = np.array(degree_hist, dtype=float)
        degree_prob = degree_hist
        plt.scatter(np.arange(degree_prob.shape[0]),degree_prob, label=label)
        plt.xlabel('k')
        plt.ylabel('p(k)')
        plt.title('Number of partners distribution')
        plt.show()
    #plot_degree_dist(G, 'full')
    #plot_degree_dist(G0,'giant') 
       
    #plt.legend()
    #plt.show()
    
    #plt.figure()
    nid=df.n_id.to_list() 
    npart=df.partner_id.to_list()
    
    males = np.unique(df.n_id.to_list())
    females = np.unique(df.partner_id.to_list())
    '''
    bottom_nodes, top_nodes = nx.bipartite.sets(G, top_nodes=males)
    
    
    Matrix=nx.algorithms.bipartite.matrix.biadjacency_matrix(G,bottom_nodes,top_nodes)
    
    bottom,top = nx.bipartite.sets(G0)
    
    
    
    giant=np.zeros(len(df.n_id))
    for index,i in enumerate(df.n_id):
        if i in bottom:
            giant[index]=1
            
    df['giant']=giant
    
    
    #Matrix1=nx.algorithms.bipartite.matrix.biadjacency_matrix(G0,bottom,top)
    
    #plt.spy(Matrix, markersize=0.4, color='blue', label='full')
    #plt.title('full')
    #plt.figure()
    #plt.spy(Matrix1, markersize=0.4, color= 'red', label='giant')
    #plt.title('giant')
    
    
    
    
    #Remove males aged less than 30
    condition1 = df.DoB_1<1991-25
    condition2 = df.DoB_1 > 1991-45
    rowstoremove = df[ condition1 & condition2]
    G_aux = G.copy()
    G_aux.remove_nodes_from(rowstoremove.n_id)
    Gcc2 = sorted(nx.connected_components(G_aux), key=len, reverse=True)
    G1 = G_aux.subgraph(Gcc2[0])
    
    
    #plt.figure()
    #plot_degree_dist(G1,'giant') 
    
    s = sorted(G1.degree, key=lambda x: x[1], reverse=True)
    degrees = sorted([x[1] for x in s])
    
    #print('Giant no males 25-45, average and var')
    #print(np.mean(degrees),np.var(degrees))
    
    
    
    #Remove high risk individuals
    low_risk_males = df['n_id'].value_counts()[df['n_id'].value_counts()<3].index.values
    rowstoremove=df[~df.n_id.isin(low_risk_males)]
    G_aux = G.copy()
    G_aux.remove_nodes_from(rowstoremove.n_id)
    Gcc2 = sorted(nx.connected_components(G_aux), key=len, reverse=True)
    G2 = G_aux.subgraph(Gcc2[0])    
    
    #plt.figure()
    #plot_degree_dist(G1,'giant') 
    
    s = sorted(G2.degree, key=lambda x: x[1], reverse=True)
    degrees = sorted([x[1] for x in s])
    
    #print('Giant no high risk males, average and var')
    #print(np.mean(degrees),np.var(degrees))
    '''
    
    
    
    
    
    
    
    
    #print('Giant sizes\n')
    #print(G.size())
    #print(G0.size())
    giant_size.append(G0.size()/G.size())
    bicomponentsize.append(len(Gc)/G0.size())
    
    #print(G1.size()/G0.size())
    #print(G2.size()/G0.size())
    #print(len(Gc)/G0.size())
    
        
    
    #bottom_nodes, top_nodes = nx.bipartite.sets(G1)
    
    #Matrix2=nx.algorithms.bipartite.matrix.biadjacency_matrix(G,bottom_nodes,top_nodes)
#plt.spy(Matrix2, markersize=0.4, color='green', label='young removed')
#plt.legend()
#plt.title('Giant component no 25-45')




'''
stuff = np.loadtxt('/home/fra/Projects/POPART-IBM/PARAMS_COMMUNITY1_ACCEPTED/param_processed_patch0_partnerships.csv', delimiter=' ',skiprows=1)
    
plt.hist(stuff[:,16],bins=30)
'''