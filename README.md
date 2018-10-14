# Set-Cover-Problem-C-

Goal:
Get clear idea of,
- Exact Algorithm using Bitmask Dynamic Programming Technique.
- Problem formulation using Integer Programming.
- Approximation using Linear Programming.
- Deterministic Rounding Technique.


Input Format : 
First line contains T, number of test cases.
Following lines contain test cases in following format
N, number of elements in the parent set( let X ), that needed to be covered
M, number of subset ( let S0 S1  ....  SN-1 ) . [  S0 U S1 U .... U SN-1 = X condition should hold ]
Then M lines contains information of M subsets in format
Wi Ki e0 e1 e2 ..... ek-1
Wi is cost of the subset
Kiis size of the subset, ei0 ei1 ei2 ..... eik-1 are elements of that subset Si

Output Format :
Ids of the subset you are taking to cover the parent set. And the minimum total cost you have are paying to cover the parent set. 

Implemented in two different ways,
- Using Bitmask DP Approach 
- Using LP Approximation Approach 
