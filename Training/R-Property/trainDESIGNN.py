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
from GraphInformation import GraphInformation
from RNETWORK import RMSELoss, designn
from datetime import datetime
from torch.utils.data import SubsetRandomSampler
from torch.utils.data.dataset import random_split

#---------------------------------Parameters which are modifiable----------------------------------#

timefile=open('Total-Time.txt','w')
fl=open('Batch-Loss.txt','w')
now = datetime.now()
print("Start Time =", now,file=timefile)
fp=open('Prediction.txt','w')
in_channels=5
hidden_channels=64
num_message_passing_steps=2
lr=0.00025
training_iterations=250
nbatch=4
#weight_decay=0.0001
weight_decay=0.00001

criterion = RMSELoss
mlpmodel=designn(in_channels,hidden_channels,hidden_channels,num_message_passing_steps)
optimizer=optim.Adam(mlpmodel.parameters(), lr=lr,weight_decay=weight_decay)

#-------------------------------------------------------------------------------------------------#
#                               PREPARATION OF THE DATASETS (PYTORCH)                             #
#-------------------------------------------------------------------------------------------------#

tensor_x,tensor_y,tensor_property,tensor_face_index,nAtoms_tensor=GraphInformation()
dataset = TensorDataset(tensor_x,tensor_y,tensor_property,tensor_face_index,nAtoms_tensor)
dataloader = DataLoader(dataset)

testing_split=0.8
shuffle_dataset=True
random_seed=42
data_size=len(dataset)
indices=list(range(len(dataset)))

train_size=int(testing_split*data_size)
val_size=data_size-train_size

if shuffle_dataset:
    np.random.seed(random_seed)
    np.random.shuffle(indices)
train_indices,val_indices=random_split(indices,[train_size,val_size])

train_sampler=SubsetRandomSampler(train_indices)
val_sampler=SubsetRandomSampler(val_indices)

training_loader = torch.utils.data.DataLoader(dataset,batch_size=nbatch,sampler=train_sampler)
validation_loader = torch.utils.data.DataLoader(dataset,batch_size=nbatch,sampler=val_sampler)


#------------------------------------------------------------------------------------------------#
#                           DESIGNN MODEL TRAINING AND VALIDATION                                #
#------------------------------------------------------------------------------------------------#
min_val_loss=5.0
for epoch in range(training_iterations):

    running_loss=0.0
    running_batchloss=0.0
    numbatch=0.0
    for ii,data in enumerate(training_loader,0):

        inputs,labels,rval,findex,nats= data

        # resizing the target value list 
        rval=torch.reshape(rval,(rval.size(0)*rval.size(1),rval.size(2)))
        target=rval

        outputs = mlpmodel(data)

#        loss=criterion(outputs,target)
        loss=torch.nn.functional.mse_loss(outputs,target)
        loss.backward(retain_graph=True)

        optimizer.step()
        optimizer.zero_grad()
        running_loss+= loss.item()
        batchloss=loss.item()/len(data)
        running_batchloss+=batchloss
        numbatch=numbatch+1
    tloss=running_batchloss/numbatch
    print(epoch,"  ------------  RUNNING BATCH-LOSS --------> ",running_batchloss/numbatch,file=fl )

    #--------------------------------------------------------------------------------------------#
    #                   Validation set for the testing of the weights                            #
    #--------------------------------------------------------------------------------------------#

    vloss=0.0
    vlossrunning=0.0
    vbatch=0

    for ii,data in enumerate(validation_loader,0):
        inpv,labv,rvalv,fcindx,nats = data
        # resizing the target value list
        rvalv=torch.reshape(rvalv,(rvalv.size(0)*rvalv.size(1),rvalv.size(2)))
        prop_val=rvalv

        mpnnget=mlpmodel(data)

        lossv=torch.nn.functional.mse_loss(mpnnget,prop_val)
        vloss +=lossv.item()
        batch_vloss=lossv.item()/(len(data))
        vlossrunning+=batch_vloss
        vbatch+=1


        print(epoch,file=fp)
        for z in range(mpnnget.size(0)):
            print(mpnnget[z][0],prop_val[z][0],file=fp)

    f_loss=vlossrunning/vbatch
    if(f_loss < min_val_loss):
        torch.save(mlpmodel.state_dict(),'mpnn-DESIGNN_model-R.pth')
        min_val_loss=f_loss

now = datetime.now()
print("End Time =", now,file=timefile)