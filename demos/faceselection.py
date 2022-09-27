#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:21:35 2022

@author: avilledieu
"""

import gurobipy as gp
from gurobipy import GRB
import numpy as np

# task = neighbor, wroker = cell

def ilp_choose_face(visibility_matrix):
    
    (W,T) = visibility_matrix.shape
    
    
    # Create model
    m = gp.Model("facechoice")
    
    
    # Create variables
    # c_i = 1   :   worker i is chosen
    cell=m.addVars(W, vtype=GRB.BINARY, name="c")
    
    # e_ij = 1   :   task j is assigned to worker i
    edge=m.addVars(W, T, vtype=GRB.BINARY, name="e")
    
    
    # Set objective
    obj = gp.quicksum(edge)
    m.setObjective(obj, GRB.MAXIMIZE)
    
    
    # Create constraints
    # can only select two workers
    m.addConstr(gp.quicksum(cell) <= 2)
        
    # a task can only be done by one worker
    for j in range(T):
        assignment_amount = gp.LinExpr(0)
        
        for i in range(W):
            assignment_amount.addTerms(1,edge[i,j])
        
        m.addConstr(assignment_amount <= 1)
    
    # a task cannot be assigned to a unchosen worker
    for i in range(W):
        for j in range(T):
            m.addConstr(cell[i] >= edge[i,j])
    
    # if a worker cant do the task, it is not assigned to him
    for i in range(W):
        for j in range(T):
            m.addConstr(edge[i,j] <= visibility_matrix[i][j])
    
    
    # Optimize model
    m.optimize()
    
    faces=[]
    index = 0
    
    for v in m.getVars():
        
        if v.varName[0]=="c":
            
            if (v.x) == 1:
                faces+=[index]
        index += 1
        
    
    return(faces)



# Create an n * m matrix, of n=7 workers (i.e. sight cells) and m=4 tasks (i.e. target neighbor vertices)
cost_matrx = np.array(
    [
        [0, 0, 1, 1],
        [1, 0, 0, 0],
        [0, 1, 1, 1],
        [1, 1, 0, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 0],
        [1, 0, 1, 0]
    ]
)

print(ilp_choose_face(cost_matrx))

