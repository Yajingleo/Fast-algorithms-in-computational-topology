# Fast-algorithms-in-computational-topology

In this repo, we give a python script for computing the generalized link Floer complexes for two-bridge links L(p,q), where p is an
even integer. 

All the methods are defined in the class ComputingAs_TwoBridgeLink. The user just need to instantialize an object by 
A=ComputingAs_TwoBridgeLink()
Then, the user can use the methods in A to compute the chain complexes for As. 

A chain complex is stored using a dictionary, which records the differentials. 
