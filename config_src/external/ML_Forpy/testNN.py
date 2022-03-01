#!/bin/env python

import torch
import numpy as np
import subgrid
import time

# GPU setup
args_no_cuda = True
use_cuda = not args_no_cuda and torch.cuda.is_available()
device = torch.device("cuda" if use_cuda else "cpu")
print('device for inference is',device)

#load the neural network
print(torch.cuda.is_available())
nn = subgrid.load_paper_net('cpu')
if use_cuda:
    nn = nn.cuda()

def MOM6_testNN(u,v): 
   # start_time = time.time()
#    print(u.shape,v.shape)
   #normaliza the input by 10
   u = u*10.0
   v = v*10.0
   x = np.array([np.squeeze(u),np.squeeze(v)])
   x = x.astype(np.float32)
   x = x.transpose((3,0,1,2)) # new the shape is (4,ni,nj,nk)
#    print(x.shape)
   x = torch.tensor(x)
   if use_cuda:
       x = x.to(device)
   with torch.no_grad():
       start_time = time.time()
       out = nn(x)
       end_time = time.time()
   if use_cuda:
       out = out.to('cpu')
   out = out.numpy().astype(np.float64)
   # At this point, python out shape is (nk,4,ni,nj)
   # Comment-out is tranferring arraies into F order
   """
   print(out.shape)
   dim = np.shape(out)
   out = out.flatten(order='F')
   out = out.reshape(dim[0],dim[1],dim[2],dim[3], order='F')
   """
   # convert out to (ni,nj,nk)
   out = out.transpose((1,2,3,0)) # new the shape is (4,ni,nj,nk)
   dim = np.shape(out)
#    print(dim)
   Sxy = np.zeros((2,dim[1],dim[2],dim[3])) # the shape is (2,ni,nj,nk)
   epsilon_x = np.random.normal(0, 1, size=(dim[1],dim[2]))
   epsilon_x = np.dstack([epsilon_x]*dim[3])
   epsilon_y = np.random.normal(0, 1, size=(dim[1],dim[2]))
   epsilon_y = np.dstack([epsilon_y]*dim[3])
   scaling = 1e-7
   Sxy[0,:,:,:] = (out[0,:,:,:] + epsilon_x/out[2,:,:,:])*scaling
   Sxy[1,:,:,:] = (out[1,:,:,:] + epsilon_y/out[3,:,:,:])*scaling
   """
   np.savetxt('Sx_mean.txt',out[0,:,:,0])
   np.savetxt('Sx_std.txt',out[2,:,:,0])
   np.savetxt('WH_u.txt',u[:,:,1])
   np.savetxt('Sx.txt',Sxy[0,:,:,0])
   """
   # end_time = time.time()
   print("--- %s seconds for CNN ---" % (end_time - start_time))
   # print(nn)
   print(Sxy.shape)
   return Sxy 

# if __name__ == '__main__':
#   start_time = time.time()
#   x = np.arange(1, 1251).astype(np.float32)
#   x = x / 100
#   print(x[:10])
#   x = x.reshape((1, 2, 25, 25), order='F')
#   x = torch.tensor(x)
#   if use_cuda:
#       x = x.to(device)
#   with torch.no_grad():
#       out = nn(x)
#   if use_cuda:
#       out = out.to('cpu')
#   out = out.numpy()
#   out = out.flatten(order='F')
#   print("BEGINNING OF PYTHON")
#   print(out[:10])
#   print("END OF PYTHON")
#   end_time = time.time()
#   print("time elapse with", device, "is", end_time-start_time, "s")