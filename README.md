# An Averaging Library

Unnecessarily many ways to get an average of a bunch of numbers

## todo

* decorators to make coding easier
  * single elem list -> just return the elem
  * zero elem list -> create new exception
  * list of all the same elem -> just return the elem
  * infinity -> return infinity
  * handle iterables and *args by converting to list
  * handle negatives correctly (for means that don't allow negatives)
  * handle complex numbers correctly
* docs
  * ordering of means via generalized mean, lehmer mean, stolarsky mean (??? but doesn't apply to multivariate case)
  * use a table with the dim in a column
* code
  * rewrite logarithmic mean to just calculate the constants (weights) at each step since that's probably faster
  * toader mean
  * all the two-variable means, but in a separate package?
  * automatically apply divided differences to 2-element means?
  * refer to muste.7z for the code for logarithmic means
