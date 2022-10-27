#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 15:12:12 2022

@author: avilledieu
"""

import gurobipy as gp
from gurobipy import GRB
import numpy as np

# task = neighbor, wroker = sight_cell

def ilp_choose_subface(induced_cross_A,induced_cross_B,edge_length_dif_a,edge_length_dif_b):
    
    (W_a,T) = induced_cross_A.shape
    (W_b,T) = induced_cross_B.shape
    
    
    # Create model
    m = gp.Model("subfacechoice")
    
    
    # Create variables
    # c_i = 1   :   worker i is chosen
    cell_a=m.addVars(W_a, vtype=GRB.BINARY, name="ca")
    cell_b=m.addVars(W_b, vtype=GRB.BINARY, name="cb")
    
    # e_ij = 1   :   task j is assigned to worker i
    edge_a=m.addVars(W_a, T, vtype=GRB.BINARY, name="ea")
    edge_b=m.addVars(W_b, T, vtype=GRB.BINARY, name="eb")
    
    print(induced_cross_B.shape)
    print(edge_length_dif_b.shape)
    
    # Set objective
    obj = gp.LinExpr(0)
    for j in range(T):
        for i in range(W_a):
            obj.addTerms(induced_cross_A[i][j], edge_a[i,j])
            obj.addTerms(0.01*edge_length_dif_a[i][j], edge_a[i,j])
        for i in range(W_b):
            obj.addTerms(induced_cross_B[i][j], edge_b[i,j])
            obj.addTerms(0.01*edge_length_dif_b[i][j], edge_a[i,j])
    m.setObjective(obj, GRB.MINIMIZE)
    
    
    # Create constraints
    # have to select one subface per face
    m.addConstr(gp.quicksum(cell_a) <= 1)
    m.addConstr(gp.quicksum(cell_b) <= 1)
        
    # a task can only be done by one worker
    for j in range(T):
        assignment_amount = gp.LinExpr(0)
        
        for i in range(W_a):
            assignment_amount.addTerms(1,edge_a[i,j])
        for i in range(W_b):
            assignment_amount.addTerms(1,edge_b[i,j])
        
        m.addConstr(assignment_amount >= 1)
    
    # a task cannot be assigned to a unchosen worker
    for j in range(T):
        for i in range(W_a):
            m.addConstr(cell_a[i] >= edge_a[i,j])
        for i in range(W_b):
            m.addConstr(cell_b[i] >= edge_b[i,j])
    
    
    # Optimize model
    m.optimize()
    
    index_A = 0
    index_B = 0
    subface_A=-1
    subface_B=-1
    
    for v in m.getVars():
        if (v.x) == 1:
            print('%s %g' % (v.varName, v.x))
        if v.varName[0]=="c":
            if v.varName[1] == "a":
                if (v.x) == 1:
                    subface_A=index_A
                index_A += 1
                        
            if v.varName[1] == "b":
                if (v.x) == 1:
                        subface_B=index_B
                index_B += 1
        
    
    return(subface_A,subface_B)



#m = the number of tasks = the number of target vertex neighbors, m = 4

# n_a = 7 subfaces in face A
induced_edge_crossings_a = np.array(
    [
        [0, 7, 9, 3],
        [0, 1, 3, 4],
        [0, 3, 2, 2],
        [0, 4, 3, 1],
        [0, 1, 4, 1],
        [0, 1, 1, 3],
        [0, 2, 1, 4]
    ]
)

# n_b = 5 subfaces in face B
induced_edge_crossings_b = np.array(
    [
        [1, 0, 4, 1],
        [2, 0, 5, 9],
        [4, 0, 3, 2],
        [2, 0, 8, 7]
    ]
)

edge_length_dif_a = np.array(
    [
        [4.94866324, 1.92283011, 3.17342655, 3.96672705],
        [0.3594641,  0.30426609, 2.51659739, 4.37369457],
        [0.64312183, 3.93699876, 1.44410342, 1.49997662],
        [2.58012235, 1.33352778, 1.63271866, 1.67656283],
        [3.49056989, 0.95530932, 0.5160346,  4.76360044],
        [3.53438154, 1.17276964, 4.84642644, 2.17530969],
        [0.25036604, 1.09324651, 3.30273229, 1.09304908]
    ]
)

edge_length_dif_b = np.array(
    [
        [1.34734909, 3.65925183, 3.00975595, 1.69575616],
        [1.30095049, 0.11525443, 4.05927489, 2.31190416],
        [4.97803758, 0.35838847, 3.25700744, 2.64590638],
        [3.28782913, 2.24382805, 0.01995905, 1.52653519],
    ]
)


print(ilp_choose_subface(induced_edge_crossings_a,induced_edge_crossings_b,edge_length_dif_a,edge_length_dif_b))