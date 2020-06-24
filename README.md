# Feature-selection-using-SCAD-regression

This project is focused on optimization problems in SCAD regression "http://www.personal.psu.edu/ril4/research/penlike.pdf"
There are two main issues about SCAD regression. First one is hyper-parameters tuning("Lambda" and "a" value), 
and the second one is to optimize non-concave penalized likelihood function. So, this project includes two optimization problems and simulation results.

The main issuse of SCAD regression is its difficulty to optimize non-concave likelihood function. 
The authors approximated this likelihood function into a quadratic form and used Newton-Raphson method.
So, I followed their way for this with my own codes. And for hyper-parameters tuning I tried the other way, named "Sequential model-based optimization".
My way is one of the most popular way for hyper-parameter tuning these days.

I found that SCAD can be useful to select features. They proved it has oracle properties in the article, so, it can be a reliable way for feature selection.
And during the iterative calculation for beta coefficients, it is already very close to the final value such as within 3 or 4 iterations.


I attached my own codes and summary file for this work.
