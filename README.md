# Codegebra
A math software that can solve equations, find derivatives, and more.
## Motivation
I created this so that I could learn how systems like Matlab and Mathematica solve equations. However, I continued
to add more features after reaching my initial goal.
## Commands
Note that commands are case-sensitive.
 - ```SOLVE```: solve linear, quadratic, and exponential equations. Exponential equations must be in simplest for this to work
 - ```DIFF```: Differentiate a function
 - ```SCALE```: Scale a vector
 - ```DOTPR```: Find the dot product of 2 vectors
 - ```ADDVEC```: Find the sum of 2 vectors
 - ```EIGVAL```: Find the eigenvalues of a matrix
 - ```TRANSP```: Transpose a matrix
 - ```FFT```: Run a Discrete Fourier Transform on a signal
 - ```IFFT```: Run an Inverse Fourier Transform on a signal
 - ```INTEGRATE```: Integrate a function
 - ```SUM```: The sum of a vector
 - ```AVG```: The average of a vector
 - ```WORD```: Define a word(similar to Wolfram|Alpha)
 - ```EXP```: Calculate e^x where x is an expression, vector, matrix, or number
 - ```SQRT```: Calculate the square root of an expression, vector, matrix, or number
 - ```SIN```: Calculate the sine of an expression, vector, matrix, or number
 - ```COS```: Calculate the cosine of an expression, vector, matrix, or number
 - ```TAN```: Calculate the tangent of an expression, vector, matrix, or number
 - ```ATAN```: Calculate the arc tangent of an expression, vector, matrix, or number
 - ```ASIN```: Calculate the arc sine of an expression, vector, matrix, or number
 - ```ACOS```: Calculate the arc cosine of an expression, vector, matrix, or number
 - ```TANH```: Calculate the hyperbolic tangent of an expression, vector, matrix, or number
 - ```SINH```: Calculate the hyperbolic sine of an expression, vector, matrix, or number
 - ```COSH```: Calculate the hyperbolic cosine of an expression, vector, matrix, or number
 - ```ATANH```: Calculate the inverse hyperbolic tangent of an expression, vector, matrix, or number
 - ```ASINH```: Calculate the inverse hyperbolic sine of an expression, vector, matrix, or number
 - ```ACOSH```: Calculate the inverse hyperbolic cosine of an expression, vector, matrix, or number
 - ```SEC```: Calculate the secant of an expression, vector, matrix, or number
 - ```SECH```: Calculate the hyperbolic secant of an expression, vector, matrix, or number
 - ```CSC```: Calculate the cosecant of an expression, vector, matrix, or number
 - ```CSCH```: Calculate the hyperbolic cosecant of an expression, vector, matrix, or number
 - ```LOG```: Calculate the natural logarithm of an expression, vector, matrix, or number
 - ```ABS```: Calculate the absolute value of an expression, vector, matrix, or number
 - ```DETR```: Calculate the determinant of a 2x2 or 3x3 matrix
 - ```CONJ```: Conjugate a  matrix
 - ```ROUND```: Round a number to the nearest whole number
 - ```EYE```: Create an identity matrix(a matrix with ones going diagonally and zero everywhere else)
 - ```EVAL```: Evaluate an expression
 - ```PROD```: Product of all elements in an array
 - ```ELEMENT```: Get info about an element in the periodic table(by atomic number)
 - ```MOVIE```: Get info about a movie(the data source I used had not been updated since 2013)
 - ```INV```: Compute the inverse of a square matrix
 - ```MXV```: Matrix-vector multiplication
 - ```ONES```: An NxM matrix of all ones
 - ```ZEROES```: An NxM matrix of all zeroes
 - ```LTRI```: Isolate the lower triangular portion of a square matrix and set the upper part to zeroes
 - ```UTRI```: Isolate the upper triangular portion of a square matrix and set the lower part to zeroes
 - ```FAC```: Factorial of a number
 - ```RAT```: Approximate the rational form of a number
 - ```DIST```: Find the distance between 2 points in a 2d or 3d plane
 - ```PATH```: Find the distance taken of a path that traverses user-specified points on a 2d or 3d plane
 - ```TRANS```: Translate a 2d or 3d shape
 - ```PTOL```: Distance from a point to a line
 - ```SYST```: Solve a system of linear equations
 - ```DILA```: Dilate a shape about point p by scale factor k
 - ```HILB```: An NxN Hilbert matrix
 - ```MAGIC```: Generate NxN magic matricies when N is odd
 - ```RCOND```: Reciprocal condition number of a matrix
 - ```FACTOR```: Factor a quadratic
## Matrices
Matrices are written in the following format:
```[1, 2, 3;4, 5, 6]```
The comma separates elements in the matrix, and the semicolon seperates rows. Matrix elements can be imaginary, complex, or real numbers. They can also be expressions, which will be simplified at runtime.
## Chaining Operations
Similar to a calculator, we provide the constants ANS and ANS2. ANS stores the result of all mathematical operations, and ANS2 stores secondary answers(such as when solving a quadratic or system of equations). Just make your input to a command ANS or ANS2 to include the result of the previous operation.
## Data credits
 - [Movie Data](https://github.com/reisanar/datasets/blob/master/HollywoodMovies.csv)
 - [Elements Data](https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee)
 - [English Dictionary Data](https://github.com/benjihillard/English-Dictionary-Database)