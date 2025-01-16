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
from torch_geometric.nn import global_add_pool, global_mean_pool, global_max_pool
import torch_geometric.graphgym.register as register
from PredictionGraphInformation import GraphInformation
from TNETWORKDESIGNN import RMSELoss, designn

#-------------------------------------------------------------------------------------------------#
#                               PREPARATION OF THE DATASETS (PYTORCH)                             #
#-------------------------------------------------------------------------------------------------#
fpred=open('predicted_t','w')
in_channels=5
hidden_channels=128
num_message_passing_steps=2
modelmpnn = designn(in_channels,hidden_channels,hidden_channels,num_message_passing_steps)
modelmpnn.load_state_dict(torch.load('mpnn-DESIGNN_model-T.pth'))
nbatch=4

tensor_x,tensor_y,tensor_face_index,nAtoms_tensor=GraphInformation()
dataset = TensorDataset(tensor_x,tensor_y,tensor_face_index,nAtoms_tensor)
dataloader = DataLoader(dataset)
structure_loader = torch.utils.data.DataLoader(dataset,batch_size=nbatch,shuffle=False)

#------------------------------------------------------------------------------------------------#
#                           DESIGNN MODEL TRAINING AND VALIDATION                                #
#------------------------------------------------------------------------------------------------#

for ii,data in enumerate(structure_loader,0):
    inputs,labels,findex,nats= data
    outputs = modelmpnn(data)
    print(outputs,file=fpred)





