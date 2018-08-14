# Fast-algorithms-in-computational-topology

Heegaard Floer homology is a [Topological Quantum Field Theory](https://en.wikipedia.org/wiki/Topological_quantum_field_theory) for 3-manifolds. In [this paper](https://arxiv.org/abs/1011.1317), Manolescu and Ozsvath give a way to compute Heegaard Floer homology combinatorically using surgeries on links. Generalized link Floer complexes is defined by Manolescu in this paper to compute Heegaard Floer homology of surgeries on links.  

In this repo, we give a python script for computing the generalized link Floer complexes for two-bridge links L(p,q), where p is an even integer. See [this paper.] (https://arxiv.org/abs/1402.5727) 



All the methods are defined in the class ComputingAs_TwoBridgeLink. The user just need to instantialize an object by 

*A = ComputingAs_TwoBridgeLink()*

Then, the user can use the methods in *A* to compute the chain complexes for *A*s. 

A chain complex is stored using a dictionary, which records the differentials. 
