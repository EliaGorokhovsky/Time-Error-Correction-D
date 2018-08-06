# Gaussianity
This code is for a project to correct misreported observation times in data assimilation. 
This branch is for an experimental change to the code that will allow the inferential likelihood algorithm to assume normality in time error. This could have several benefits:
* Computations should be faster
* The algorithm may be able to correct for time error internally
* More work can be done towards making inference entirely analytical
However, it adds an additional constraint to the problem.