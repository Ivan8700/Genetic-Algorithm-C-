# Genetic Algorithm in C++
## An Heuristic to solve "Load scheduling with even sums" NP problem
**<ins>Problem :</ins>** Given N positive integers (also known as "jobs") and M machines, find a schedule for these jobs on the machines such that the makepsan will be minimal.  
**<ins>Algorithm :</ins>** Genetic/Evolutionary algorithm.  
**<ins>initialization :</ins>** 
Every solution for this problem can be seen as an array of size 'N' (as the number of jobs) where in every index we have an integer number in [0,M-1].  
index 'i' of the array will mean 'job in index i' and the number 'j' inside the array will mean 'the job in index i is scheduled in machine number j'.  
Set "Population" size - how many Genes will be in every generation.  
Random an array of size 'Population" which will contain different schedules of the input on the M machines.  
The algorithm will works recursively as follows :  
1) Evaluate each gene in the array based on an evaluation function
2) set for each gene the relative probabilities that fits it by the evaluation function  
3) set for each gene the interval of probabilities (first gene is [0,his probability], second gene is [previos upper probability, previos + this gene prob] and so on)
4) random 2 numbers in (0,1), find the genes associated with the random numbers, make a "switching procedure" which partially switches the genes but place them in the "next generation" (do not effect the current genes)  
5) repeat step 4 (population-mutation)/2 times  
6) random a number in (0,1), find the gene associated with the random number, make a "mutation procedure" which randomaly changes some indecies in the gene and place it in the "next generation"
7) repeat step 6 "mutation" times
8) check how many different genes in the "next generation" we are left with in the end of the procedure.
9) if too much, repeat from step 1. if small enought number - exit and print the best solution.  

Overall the point of the algorithm is - the better the solution is, the higher evaluation and probability it will have, the higher amount of times it will be chosen to participate in the next generation.  
For large inputs, the algorithm will never get an optimal solution, but a close approximation will be found in a pretty good run-time.
