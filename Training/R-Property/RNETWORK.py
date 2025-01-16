import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree
from torch_geometric.nn import global_mean_pool, global_max_pool,global_add_pool
from GraphInformation import GraphInformation
import numpy as np



#---------------------------------------------------------------------#
#                    The message passing layer                        #
#---------------------------------------------------------------------#

class mpnnlayer(MessagePassing):
    def __init__(self,in_channels,out_channels):
        super().__init__(aggr='add')
        torch.manual_seed(5)
        self.lin= nn.Sequential(
                nn.Linear(in_channels,out_channels),
                nn.Tanh())
        self.lin1g=nn.Sequential(
                nn.Linear(out_channels,out_channels),
                nn.Tanh(),
                nn.Linear(out_channels,in_channels)
                )
        self.reset_parameters()

    def reset_parameters(self):
        for layer in self.lin:
            if hasattr(layer, 'reset_parameters'):
               layer.reset_parameters()

    def forward(self,x,edge_index,gc):
        if gc>0:
            x=self.lin(x)
            nt=x.size(1)
            bnn=nn.BatchNorm1d(nt)
            x=self.lin1g(x)

        x=self.propagate(edge_index,x=x)
        return x

    def message(self,x_j,edge_index):
        return x_j

    def update(self,aggr_out):
        return aggr_out

# The DESIGNN class for complete functioning

class designn(nn.Module):
    def __init__(self,in_channels,hidden_channels,out_channels,num_message_passing_steps):
        super().__init__()
        torch.manual_seed(5)
        self.in_channels=in_channels
        self.hidden_channels=hidden_channels
        self.num_message_passing_steps=num_message_passing_steps
        self.out_channels=out_channels
        self.gnn_layers=nn.ModuleList(mpnnlayer(in_channels,hidden_channels) for _ in range(num_message_passing_steps))

        self.lin=nn.Sequential(
                nn.Linear(5,64),
                nn.Tanh()
                )
        self.lin1mlp=nn.Sequential(
                nn.Linear(64,32),
                nn.Tanh(),
                nn.Dropout(0.2),
                nn.Linear(32,16),
                nn.Tanh(),
                nn.Linear(16,1)
                )

    def forward(self,data):
        inputs,labels,rval,findex,nats= data
        x,edge_index=pygdata(inputs,labels)
        edge_index,_ = add_self_loops(edge_index,num_nodes=x.size(0))
        for gc,gnn in enumerate(self.gnn_layers):
            x=gnn(x,edge_index,gc)
       
        
        bs = inputs.size(0)*inputs.size(1)
        batch=np.empty(bs)
        bach=pooling_choices(inputs,batch,findex,nats)
        #----------------------> 2. Carry out node pooling for final MLP <-----------------------#
       
        outputs=x
        node_pooling=global_max_pool(outputs,bach)
        ctr=0
        node_doublepool=torch.empty((int(node_pooling.size(0)/3),5))
        for i in range(0,node_pooling.size(0),3):
             node_doublepool[ctr][0::]=(node_pooling[i][0::])
             ctr=ctr+1

        node_doublepool.clone().detach()
        x=node_doublepool
        x=self.lin(x)
        nt=x.size(1)
        bnn=nn.BatchNorm1d(nt)
        x=self.lin1mlp(x)

        return x

#------Specifically written RMSE loss function----#

def RMSELoss(yhat,y):
    return torch.sqrt(torch.mean((yhat-y)**2))

#-----Convert to the format relevant to pytorch geometric------------#

def pygdata(x,edge_index):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #       Converts the batch of edges to the [2,E] dimension       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #________________________________________________________________#
    #1. readjust the edge_indexes by adding the indexes based on the #
    #   number of atoms                                              #
    #----------------------------------------------------------------#
    num_nodes=x.size(1)
    p=0
    patomcount=0
    tEdges=edge_index.size(1)*edge_index.size(0)
    eg_inx=np.empty((2,edge_index.size(1)*edge_index.size(0)))

    for f in range(edge_index.size(0)):
        if f>0:
            patomcount=patomcount+num_nodes
        else:
            patomcount=0
        for j in range(edge_index.size(1)):
            eg_inx[0][p]=edge_index[f][j][0]+patomcount
            eg_inx[1][p]=edge_index[f][j][1]+patomcount
            if(eg_inx[0][p]==eg_inx[1][p]):
                eg_inx[0][p]=0
                eg_inx[1][p]=0
            p=p+1

    edge_indextorch=torch.from_numpy(eg_inx)
    edge_indextorch=edge_indextorch.int()
    edge_index= edge_indextorch
    edge_index=edge_index.type(torch.int64)
    x=torch.reshape(x,(x.size(0)*x.size(1),x.size(2)))
    return x,edge_index


def pooling_choices(inputs,batch,findex,nats):
          #----------------------------------------------------------------------------------------#
          #                   Carry out node pooling for the data in the batch                     #
          #----------------------------------------------------------------------------------------#
  
          #--------------------> 1.Preparation of the batch for pooling <--------------------------#
          #--------------------> 1(a). Individual structures batched <-----------------------------#
          sval=0
          for p in range(inputs.size(0)):
              evals=sval+inputs.size(1)
  
              for k in range(sval,evals):
                  batch[k]=p
              sval=evals
  
          #-------------------> 1(b). Individual structures batched <------------------------------#
          #-------------------> based on the relevance of the face nodes.<-------------------------#
  
          ctr=0
          sval=0
  
          for p in range(inputs.size(0)):
              tmp=sval
              evals=sval+inputs.size(1)
              evalspc=sval+nats[p]
              ct=0
              for k in range(0,inputs.size(1)):
                  if(findex[p][k][0]==1):
                      batch[tmp]=ctr
                      ct=ct+1
                  tmp=tmp+1
              ctr=ctr+1
              tmp=sval
              ct=0
              for k in range(0,nats[p]):
                  if(findex[p][k][0]==0):
                      batch[tmp]=ctr
                      ct=ct+1
                  tmp=tmp+1
  
              ctr=ctr+1
              tmp=evalspc
              ct=0
              for k in range(nats[p],inputs.size(1)):
                      batch[tmp]=ctr
                      ct=ct+1
                      tmp=tmp+1
              ctr=ctr+1
              sval=evals
  
          batch=torch.tensor(batch).int()
          batch=batch.type(torch.int64)
  
          return batch

