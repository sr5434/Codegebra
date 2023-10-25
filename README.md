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
 - ```EXP```: Calculate e^x where x is an expression or number
 - ```SQRT```: Calculate the square root of an expression or number
 - ```SIN```: Calculate the sine of an expression or number
 - ```COS```: Calculate the cosine of an expression or number
 - ```TAN```: Calculate the tangent of an expression or number
 - ```ATAN```: Calculate the arc tangent of an expression or number
 - ```LOG```: Calculate the natural logarithm of an expression or number
 - ```ABS```: Calculate the absolute value of an expression or number
 - ```DETR```: Calculate the determinant of a 2x2 or 3x3 matrix
 - ```CONJ```: Conjugate a  matrix
 - ```ROUND```: Round a number to the nearest whole number
 - ```EYE```: Create an identity matrix(a matrix with ones going diagonally and zero everywhere else)
 - ```EVAL```: Evaluate an expression
 - ```PROD```: Product of all elements in an array
 - ```ELEMENT```: Get info about an element in the periodic table(by atomic number)
 - ```MOVIE```: Get info about a movie(the data source I used had not been updated since 2013)
 - ```INV```: Compute the inverse of a square matrix
## Matrices
Matrices are written in the following format:
```[1, 2, 3;4, 5, 6]```
The comma separates elements in the matrix, and the semicolon seperates rows. Matrix elements can be imaginary, complex, or real numbers. They can also be expressions, which will be simplified at runtime.
## Data credits
 - [Movie Data](https://github.com/reisanar/datasets/blob/master/HollywoodMovies.csv)
 - [Elements Data](https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee)
 - [English Dictionary Data](https://github.com/benjihillard/English-Dictionary-Database)