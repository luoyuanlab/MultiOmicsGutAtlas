import gc
import math
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from collections import OrderedDict
from typing import Dict, Mapping, Optional, Tuple, Any, Union
import time
import copy
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score, multilabel_confusion_matrix
import matplotlib.pyplot as plt
from collections import defaultdict
from torch.nn.parameter import Parameter



def glorot_init(input_dim, output_dim):
    init_range = np.sqrt(6.0/(input_dim + output_dim))
    initial = torch.rand(input_dim, output_dim)*2*init_range - init_range
    return nn.Parameter(initial)



class GraphConvSparse(nn.Module):
    def __init__(self, input_dim, output_dim, adj, activation=F.relu, dropout=0, norm=nn.LayerNorm, bias=True, **kwargs):
        super(GraphConvSparse, self).__init__(**kwargs)
        
        self.adj = adj
        self.activation = activation
        self.dropout = dropout
        self.norm = norm
        
        #self.weight = glorot_init(input_dim, output_dim) 
        #self.weight = Parameter(torch.FloatTensor(input_dim, output_dim))
        #self.reset_parameters()
        
        self.weight = Parameter(torch.FloatTensor(input_dim, output_dim))
        if bias:
            self.bias = Parameter(torch.FloatTensor(output_dim))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    #def reset_parameters(self):
    #    torch.nn.init.xavier_uniform_(self.weight)
        
    def reset_parameters(self):
        stdv = 1. / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, inputs):
        x = inputs
        x = torch.mm(x, self.weight)
        x = torch.mm(self.adj, x)
        
        if self.bias is not None:
            return x + self.bias
        else:
            return x
        x = self.norm(x)
        
        outputs = self.activation(x)
        return F.dropout(x, p=self.dropout, training=self.training)

    

class CosMxGAE(nn.Module):
    def __init__(self, adj, input_dim, output_dim, hidden_dim = None, dropout=0):
        super(CosMxGAE, self).__init__()

        self.input_dim = input_dim
        self.output_dim = output_dim
        if hidden_dim == None:
            hidden_dim = (input_dim + output_dim) // 2
        self.hidden_dim = hidden_dim

        ## Ecoder for RNA VAE
        self.encoder = nn.Sequential(
            GraphConvSparse(input_dim, hidden_dim, adj, activation=F.leaky_relu, dropout=dropout, norm=lambda x:x),
            #GraphConvSparse(input_RNA_dim, hidden1_RNA_dim, adj_RNA, activation=F.leaky_relu, dropout=dropout, norm=nn.LayerNorm),
            #nn.BatchNorm1d(hidden1_RNA_dim),
            nn.LayerNorm(hidden_dim),
            GraphConvSparse(hidden_dim, output_dim, adj, activation=F.leaky_relu, dropout=dropout, norm=lambda x:x),
            #GraphConvSparse(hidden1_RNA_dim, hidden2_dim, adj_RNA, activation=F.leaky_relu, dropout=dropout, norm=nn.LayerNorm),
            #nn.BatchNorm1d(hidden2_dim),
            #nn.LayerNorm(output_dim)
        )
        #self.encoder = nn.DataParallel(self.encoder)

    def to(self, device):
        return super(CosMxGAE, self).to(device)
        
    def forward(self, X_RNA):
        
        X_out =  self.encoder(X_RNA)
        
        return X_out

