#!/usr/bin/python
# -*- coding: utf-8 -*-
from evtk.hl import gridToVTK
import numpy as np

def meshgrid_nd(*arrs):
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []    
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = np.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j) 
        ans.append(arr2)

    return tuple(ans)

def opg_op2vtk(data, fname):
    x = data['dens'][:]
    y = data['temp'][:]
    z = 0.5*(data['groups'][1:] + data['groups'][:-1])
    #X, Y, Z = meshgrid_nd(x, y,z)
    pointData = {key: data[key][:] for key in ['opp_mg', 'opr_mg', 'emp_mg']}

    gridToVTK(fname, x, y, z, pointData = pointData)
