# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 08:46:14 2021

@author: M302212
"""
import numpy as np
from Adenventure.Synthesis.SynthInfra import Feed, Product, Flow, Waste

class Individual():
    def __init__(self, feed):
        self.incidence_matrix = np.ones((1,), dtype=bool)
        self.mandatory = [True, True]
        self.operations = [
            Feed(feed), Product(feed)
        ]
    
    @property
    def input_operations(self):
        inputs = []
        for operation in self.operations:
            inputs.extend([operation]*len(operation.inputs))
        return inputs
    
    @property
    def output_operations(self):
        outputs = []
        for operation in self.operations:
            outputs.extend([operation]*len(operation.outputs))
        return outputs
            
    def inputs(self, operation):
        outputs = []
        for i, op in enumerate(self.input_operations):
            if op is operation:
                absolute_pos = np.where(self.incidence_matrix[:,i])[0][0]
                relative_pos = -1
                for j in range(absolute_pos+1):
                    if self.output_operations[absolute_pos] is self.output_operations[j]:
                        relative_pos += 1
                outputs.append((self.output_operations[absolute_pos], relative_pos))
        return outputs
    
    def outputs(self, operation):
        inputs = []
        if self.incidence_matrix.shape == (1,):
            return self.input_operations[0]
        for i, op in enumerate(self.output_operations):
            if op is operation:
                output_stream = np.where(self.incidence_matrix[i,:])[0][0]
                inputs.append(self.input_operations[output_stream])
        return inputs
    
    def _fill_waste(self):
        waste = []
        for i, row in enumerate(self.incidence_matrix):
            if not np.any(row):
                waste.append(i)
        
        if len(waste) == 0:
            return None
        newx = self.incidence_matrix.shape[0]
        newy = self.incidence_matrix.shape[1] + len(waste)
        incidence_matrix = np.zeros((newx, newy), dtype=bool)
        incidence_matrix[
            :self.incidence_matrix.shape[0],
            :self.incidence_matrix.shape[1]
        ] = self.incidence_matrix
        for i, pos in enumerate(waste):
            incidence_matrix[pos, -i-1] = True
            self.operations.append(Waste(Flow(0)))
            self.mandatory.append(False)
        self.incidence_matrix = incidence_matrix
            
            
    
    def insert_unit_operation(self, unit_operation, i, j, product_yield, product_purity):
        feed = Flow(0)
        unit_operation = unit_operation(feed, product_yield, product_purity)
        
        dimx = self.incidence_matrix.shape[0]
        if self.incidence_matrix.shape != (1,):
            dimy = self.incidence_matrix.shape[1]
        else:
            dimy = 1
        
        newx = dimx + len(unit_operation.outputs)
        newy = dimy + len(unit_operation.inputs)

        incidence_matrix = np.zeros((newx, newy), dtype=bool)
        incidence_matrix[
            :dimx,
            :dimy
        ] = self.incidence_matrix
        incidence_matrix[i, j] = False
        # Assumption one feed - distillation test case
        incidence_matrix[i, -1] = True
        output = np.random.randint(0, len(unit_operation.outputs))+1
        incidence_matrix[-output, j] = True
        self.incidence_matrix = incidence_matrix
        self.operations.append(unit_operation)
        self.mandatory.append(False)
        self._fill_waste()
        
    def mutate_unit_operation(self):
        if np.all(self.mandatory) == False:
            found = False
            while not found:
                index = np.random.choice(len(self.operations))
                operation = self.operations[index]
                if (not operation.is_endpoint) and (not operation.is_startingpoint):
                    found = True
                    product_yield = operation.product_yield
                    product_purity = operation.product_purity
                    self.operations[index] = type(operation)(Flow(0), product_yield, product_purity)
    
    def pick_random_connection(self):
        if self.incidence_matrix.shape == (1,):
            return 0, 0
        options = np.where(self.incidence_matrix)
        sel = np.random.choice(len(options[0]))
        return options[0][sel], options[1][sel]
    
    def operation_position(self, operation):
        for i, op in enumerate(self.operations):
            if op is operation:
                return i
            
    def trace_operation(self, operation):
        inputs = self.inputs(operation)
        flows = []
        for this_input, pos in inputs:
            if not self.updated[self.operation_position(this_input)]:
                self.trace_operation(this_input)
                print('Something is wrong! Not if there is extraction.')
            flows.append(this_input.outputs[pos])
        operation.change_input(np.array(flows))
        operation.calculate_output()
        self.updated[self.operation_position(operation)] = True
        for output in self.outputs(operation):
            self.trace_operation(output)
        
        
    def calculate_material_flow(self):
        self.updated = np.zeros(len(self.operations), dtype=bool)
        # Assumption - first UO is main feed
        self.updated[0] = True
        #print(1)
        self.trace_operation(self.outputs(self.operations[0])[0])
        
    def cost(self):
        total = 0
        for operation in self.operations:
            total += operation.calculate_cost()
        return total
    
    def flow_at_pos(self, i, j):
        pos = -1
        for k, operation in enumerate(self.output_operations):
            if k > i:
                break
            if operation is self.output_operations[i]:
                pos += 1
        return self.output_operations[i].outputs[pos]
    
    
        
        
        