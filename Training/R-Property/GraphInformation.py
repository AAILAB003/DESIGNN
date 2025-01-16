import torch
import torch.nn as nn
import torch.nn.functional 
import numpy as np
from torch import Tensor
from typing import Optional
from torch.utils.data import TensorDataset, DataLoader
import torch.optim as optim
import torch_geometric
import torch_geometric.nn as pyg_nn
import torch_geometric.data as datas
from torch_geometric.nn import MessagePassing
from torch.nn import Linear,Parameter
from torch_geometric.utils import add_self_loops, degree, to_dense_adj,dense_to_sparse
from torch_geometric.nn import global_add_pool, global_mean_pool, global_max_pool
from torch_geometric.data import Batch, Data
import torch_geometric.graphgym.register as register
from torch_geometric.graphgym.config import cfg
import statistics
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import math
from sklearn.metrics import mean_squared_error
   
def GraphInformation():
    
    #----------------------------------------------------------------#
    #  Reading of atomic-features and pair indices & the conversions #
    #----------------------------------------------------------------#

    f=open("atomic-features.txt")
    lines=f.read().splitlines()
    
    f1=open("nAtoms.txt")
    lines1=f1.read().splitlines()
    lines1=[int(i) for i in lines1]

    nAtoms_tensor=torch.tensor(lines1) 

    atomic_feature=[]
    s=0
    for n in range(len(lines1)):
        nAtoms=lines1[n]
        my_x_mol=[]
        
        for i,line in enumerate(lines[s:s+nAtoms]):
            line=line.split()
            tlist=[]
            for j in line:
                tlist.append(float(j))
            my_x_mol.append(tlist)
        s=s+nAtoms
        atomic_feature.append(my_x_mol)
    
    f2=open('Pairs.txt')
    lines2=f2.read().splitlines()
    
    f3=open('nPairs-list.txt')
    lines3=f3.read().splitlines()
    lines3=[int(i) for i in lines3]
    
    pair_indices=[]
    s=0
    for n in range(len(lines3)):
        npairs=lines3[n]
        pair_mol=[]
        for i,line in enumerate(lines2[s:s+npairs]):
            line=line.split()
            tlist=[]
            for j in line:
                tlist.append(int(j))
            pair_mol.append(tlist)
        s=s+npairs
        pair_indices.append(pair_mol)
    
    
    f4=open("r-value")
    lines4=f4.read().split()
    
    system_property=[]
    for n in range(len(lines3)):
        my_prop_mol=[]
        new_prop_mol=[]
        my_prop_mol.append(float(lines4[n]))
        new_prop_mol.append(my_prop_mol)
        system_property.append(new_prop_mol)
    
    # Prepare the face specific node tensor
    
    f5=open("face-index.txt")
    lines5=f5.read().split()
    
    face_index=[]
    s=0
    for n in range(len(lines1)):
        nAtoms=lines1[n]
        f_mol=[]
        for i,line in enumerate(lines5[s:s+nAtoms]):
            line=line.split()
            tlist=[]
            for j in line:
                tlist.append(int(j))
            f_mol.append(tlist)
        face_index.append(f_mol)
        s=s+nAtoms
    
    
    #-----------------------------------------------------------------#
    # PADDING THE DATA                                              #
    #---------------------------------------------------------------#
    target=atomic_feature
    max_length = max(len(row) for row in target)
    max_cols=max([len(cols) for batch in target for cols in batch])
    max_rows=max([len(batch) for batch in target])
    max_rows=max_rows+1
    padded=[batch + [[0] * (max_cols)] * (max_rows - len(batch)) for batch in target]
    padded=torch.tensor([row + [0]* (max_cols - len(row)) for batch in padded for row in batch])
    padded=padded.view(-1, max_rows, max_cols)
    tensor_x=padded

    #---Scaling the input---#                            
    
    x=tensor_x
    full_list=[]
    mfull_list=[]
    for i in range(x.size(0)):
        for j in range(x.size(1)):
            mxv=max(x[i][j][0::])
            mnv=min(x[i][j][0::])
            full_list.append(float(mxv))
            mfull_list.append(float(mnv))

    mx=max(full_list)
    mn=min(mfull_list)

    for i in range(x.size(0)):
        for j in range(x.size(1)):
            for k in range(x.size(2)):
                x[i][j][k]=(x[i][j][k]-mn)/(mx-mn)
    tensor_x=x

    #------End of Scaling the data-------#
    
    target=pair_indices
    max_length = max(len(row) for row in target)
    max_cols=max([len(cols) for batch in target for cols in batch])
    max_rows=max([len(batch) for batch in target])
    padded=[batch + [[0] * (max_cols)] * (max_rows - len(batch)) for batch in target]
    padded=torch.tensor([row + [0]* (max_cols - len(row)) for batch in padded for row in batch])
    padded=padded.view(-1, max_rows, max_cols)
    tensor_y=padded
    
    target = system_property
    max_length = max(len(row) for row in target)
    max_cols=max([len(cols) for batch in target for cols in batch])
    max_rows=max([len(batch) for batch in target])
    padded=[batch + [[0] * (max_cols)] * (max_rows - len(batch)) for batch in target]
    padded=torch.tensor([row + [0]* (max_cols - len(row)) for batch in padded for row in batch])
    padded=padded.view(-1, max_rows, max_cols)
    tensor_property=padded
    
    target = face_index
    max_length = max(len(row) for row in target)
    max_cols=max([len(cols) for batch in target for cols in batch])
    max_rows=max([len(batch) for batch in target])
    max_rows=max_rows+1
    padded=[batch + [[0] * (max_cols)] * (max_rows - len(batch)) for batch in target]
    padded=torch.tensor([row + [0]* (max_cols - len(row)) for batch in padded for row in batch])
    padded=padded.view(-1, max_rows, max_cols)
    tensor_face_index=padded

    return tensor_x,tensor_y,tensor_property,tensor_face_index,nAtoms_tensor
