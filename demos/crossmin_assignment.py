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

def ilp_choose_subface(induced_cross_A,induced_cross_B):
    
    (W_a, T) = induced_cross_A.shape
    (W_b, T) = induced_cross_B.shape
    
    
    # Create model
    m = gp.Model("subfacechoice")

    # Create variables
    # c_i = 1   :   worker i is chosen
    cell_a = m.addVars(W_a, vtype=GRB.BINARY, name="ca")
    cell_b = m.addVars(W_b, vtype=GRB.BINARY, name="cb")
    
    # e_ij = 1   :   task j is assigned to worker i
    edge_a = m.addVars(W_a, T, vtype=GRB.BINARY, name="ea")
    edge_b = m.addVars(W_b, T, vtype=GRB.BINARY, name="eb")
    
    
    # Set objective
    obj = gp.LinExpr(0)
    for j in range(T):
        for i in range(W_a):
            obj.addTerms(induced_cross_A[i][j], edge_a[i, j])
        for i in range(W_b):
            obj.addTerms(induced_cross_B[i][j], edge_b[i, j])
    m.setObjective(obj, GRB.MINIMIZE)
    
    
    # Create constraints
    # have to select one subface per face
    m.addConstr(gp.quicksum(cell_a) <= 1)
    m.addConstr(gp.quicksum(cell_b) <= 1)
        
    # a task can only be done by one worker
    for j in range(T):
        assignment_amount = gp.LinExpr(0)
        
        for i in range(W_a):
            assignment_amount.addTerms(1, edge_a[i, j])
        for i in range(W_b):
            assignment_amount.addTerms(1, edge_b[i, j])
        
        m.addConstr(assignment_amount == 1)
    
    # a task cannot be assigned to a unchosen worker
    for j in range(T):
        for i in range(W_a):
            m.addConstr(cell_a[i] >= edge_a[i, j])
        for i in range(W_b):
            m.addConstr(cell_b[i] >= edge_b[i, j])
    
    
    # Optimize model
    m.optimize()
    
    index_A = 0
    index_B = 0
    subface_A=-1
    subface_B=-1
    
    for v in m.getVars():
        if (v.x) == 1:
            print('%s %g' % (v.varName, v.x))
        if v.varName[0] == "c":
            if v.varName[1] == "a":
                if v.x == 1:
                    subface_A = index_A
                index_A += 1
                        
            if v.varName[1] == "b":
                if v.x == 1:
                        subface_B = index_B
                index_B += 1

    # Get Assignments
    assignment_a = [None] * induced_cross_A.shape[1]
    assignment_b = [None] * induced_cross_A.shape[1]
    for v in m.getVars():
        if v.varName[0] != "e":
            continue
        if v.varName[1] == "a" and v.varName[3] == str(subface_A):
            assignment_a[int(v.varName[5])] = int(v.x)
        elif v.varName[1] == "b" and v.varName[3] == str(subface_B):
            assignment_b[int(v.varName[5])] = int(v.x)

    print(assignment_a)
    print(assignment_b)
    return subface_A, subface_B



induced_edge_crossings_a = np.array(
    [
        [0, 0, 0, 1],
        [0, 0, 0, 3],
        [0, 0, 0, 2]
    ]
)

induced_edge_crossings_b = np.array(
    [
        [0, 3, 1, 0],
        [0, 3, 3, 0],
        [0, 2, 2, 0],
        [0, 3, 1, 0],
        [0, 2, 1, 0]
    ]
)

print(ilp_choose_subface(induced_edge_crossings_a, induced_edge_crossings_b))