# An Averaging Library

Unnecessarily many ways to get an average of a bunch of numbers

## todo

* math
  * double-check the correctness of the derivations in the PDFs
* decorators to make coding easier
  * single elem list -> just return the elem
  * zero elem list -> create new exception
  * list of all the same elem -> just return the elem
* docs
  * ordering of means via generalized mean, lehmer mean, stolarsky mean (??? but doesn't apply to multivariate case)
  * use a table with the dim in a column
* code
  * rewrite logarithmic mean to just calculate the constants (weights) at each step since that's probably faster
  * toader mean
  * all the two-variable means, but in a separate package?
