#!/bin/env python

import sys
from math import sin
import torch
import numpy
import subgrid

#load the neural network
nn = subgrid.load_paper_net('cpu')


def my_testNN(x, y):
    print(x[:10])
    x = x.reshape((1, 2, 25, 25), order='F')
    x = x.astype(numpy.float32)
    x = torch.tensor(x)
    with torch.no_grad():
        out = nn(x)
    out = out.numpy().astype(numpy.float64)
    out = out.flatten(order='F')
    print(out[:10])
    y[:] = out[:]
    sys.stdout.flush()



def my_test(x) :
    print(' From Python test: {}'.format(x.size))
    for i in range(x.size) :
        x[i] += 1.0

    a = torch.from_numpy(x)
    print(a)

    b = torch.from_numpy(x.astype(numpy.float32))
    print(b)

    if torch.cuda.is_available() :
        b_dev = b.cuda()
        print(b_dev)

    my_test2()

    sys.stdout.flush()
    
def my_test2() :
    print(' **** From my_test2 ****')
    return

def MOM6_testNN(u,v): 
#    print(u.shape,v.shape)
    x = np.array([np.squeeze(u),np.squeeze(v)])
    x = x[np.newaxis,:]
    x = x.astype(numpy.float32)
#    print(x.shape)
    x = torch.tensor(x)
    with torch.no_grad():
        out = nn(x)
    out = out.numpy().astype(numpy.float64)
    # at this point, python shape is 1,4,ni,nj
#    print(out.shape)
    dim = np.shape(out)
    out = out.flatten(order='F')
    out = out.reshape(dim[0],dim[1],dim[2],dim[3], order='F')
    sys.stdout.flush()
    return out

#if __name__ == '__main__':
import numpy as np
import torch
x = np.arange(1, 1251).astype(numpy.float32)
x = x / 100
print(x[:10])
x = x.reshape((1, 2, 25, 25), order='F')
x = torch.tensor(x)
with torch.no_grad():
    out = nn(x)
out = out.numpy()
out = out.flatten(order='F')
print("BEGINNING OF PYTHON")
print(out[:10])
print("END OF PYTHON")

