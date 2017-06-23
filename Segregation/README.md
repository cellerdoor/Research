# Pattern formation in low-copy plasmids


Codes in this folder do a deterministic simulation of multiple plasmids on a substrate and show how the ParA/B interdynamics leads to the patterns observed in the cell. 
The variables of the model are:  

N - no of plasmids  
c - ratio of removal to attraction lenght scales  
L - Length of the substrate  
$K_{o}$ - rate at which the parA is renewed  
$\gamma$ - rate at which the parA is removed when plasmid is nearby  
$\delta$ - random rate at which the parA is removed   

# Steady state ParA  
  
In the absence of any ParB, assuming that the ParA substrate has reached steady state we would have   
$\frac{dA}{dt}$ = $\k_{o}$A - $\delta$A
