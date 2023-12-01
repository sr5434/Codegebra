import math, cmath
import re
from scipy.linalg import lapack
from scipy.linalg import blas
import ast
import operator as op
import pandas as pd

# supported operators(for safely running eval)
operators = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul,
             ast.Div: op.truediv, ast.Pow: op.pow, ast.BitXor: op.xor,
             ast.USub: op.neg}

pattern = r'-?[0-9.]+\s*=\s*(?:[0-9.]+\s*\*\s*)?[0-9.]+\^\(x-[0-9.]+\)(?:\s*\+\s*[0-9.]+)?'


def simplifyTerms(terms):
    terms = [item for item in terms if item is not None and item is not False and item is not True]
    x_coeff = 0  # X coefficient(aka m in y=mx+b)
    x2_coeff = 0  # X^2 coefficient(aka A in ax^2+bx+c=0)
    # Quadratic and exponential flags
    is_quadratic = False
    is_expo = "^(x" in "".join(terms)
    const = 0  # Constant(aka b in y=mx+b)
    # Variables in the exponential equation
    y = 0
    a = 0
    b = 0
    h = 0
    k = 0
    for term in terms:
        if ("x" in term or "X" in term) and (term.find("x**2") == -1 and term.find("x^2") == -1) and not is_expo:
            # Take out the X to get the coefficient
            term = term.replace("x", "")
            term = term.replace("X", "")
            if term == "":
                # If it doesn't have a coefficient than set it to one
                x_coeff += 1
            elif term == "-":
                # take away 1 if the coefficient is just a negative sign
                x_coeff += -1
            else:
                x_coeff += float(term)  # If x has a coefficient then add it to the running tally
        elif "x**2" in term or "x^2" in term:
            is_quadratic = True
            # Take out the X^2 to get the coefficient
            term = term.replace("x**2", "")
            term = term.replace("x^2", "")
            if term == "":
                # If it doesn't have a coefficient than set it to one
                x2_coeff += 1
            elif term == "-":
                # take away 1 if the coefficient is just a negative sign
                x2_coeff += -1
            else:
                x2_coeff += float(term)  # If x has a coefficient then add it to the running tally
        elif is_expo:
            if len(terms) == 1:
                y = float(term)
            elif re.findall(r'\d+\^\(', term) != []:
                print(term)
                term = term.replace("^(", "")
                term = term.replace("x", "")
                if len(term.split("*")) == 2:
                    a = float(term.split("*")[0])
                    b = float(term.split("*")[1])
                else:
                    b = float(term)
            elif re.findall(r'\d\)', term):
                h = float(term.replace(")", "")) * -1
            else:
                k = float(term)
            # Take out the X^2 to get the coefficient
            term = term.replace("x**2", "")
            term = term.replace("x^2", "")
        else:
            const += float(term)  # Add to the running tally for the constant
    if x_coeff == 1:
        x_coeff = ""
    else:
        x_coeff = str(x_coeff)
    # If a is 0 in a quadratic, is it really a quadratic?
    if x2_coeff == 1:
        x2_coeff = "x^2"
    elif x2_coeff == 0:
        is_quadratic = False
        x2_coeff = None
    else:
        x2_coeff = str(x2_coeff) + "x^2"
    # Determine which return values should be set to null
    if is_expo:
        return [None, None, False, None, a, b, h, k]
    if x_coeff == "0" and const != "0":
        return [None, str(const), is_quadratic, x2_coeff]
    elif x_coeff != "0" and const == "0":
        return [None, None, is_quadratic, x2_coeff]
    elif x_coeff == "0" and const == "0":
        return [None, None, is_quadratic, x2_coeff]
    else:
        return [f"{x_coeff}x", str(const), is_quadratic, x2_coeff]

def simpStdForm(equation):
  equation = equation.replace(" ", "")  # Remove whitespace
  # Equation mode
  right, left = equation.split("=")  # Split string into left and right sides of the equal sign
  # Seperate terms on each side of the equation
  left = " ".join(left.split("+"))
  left = left.replace("-", " -")
  right = " ".join(right.split("+"))
  right = right.replace("-", " -")
  right = right.split()
  left = left.split()
  x_coeff = 0
  y_coeff = 0
  const = 0
  for term in right:
    if ("x" in term) or ("X" in term):
      # Take out the X to get the coefficient
      term = term.replace("x", "")
      term = term.replace("X", "")
      if term == "":
          # If it doesn't have a coefficient than set it to one
          x_coeff -= 1
      elif term == "-":
          # take away 1 if the coefficient is just a negative sign
          x_coeff += 1
      else:
          x_coeff -= float(term)  # If x has a coefficient then add it to the running tally
    elif ("y" in term) or ("Y" in term):
      # Take out the X to get the coefficient
      term = term.replace("y", "")
      term = term.replace("Y", "")
      if term == "":
          # If it doesn't have a coefficient than set it to one
          y_coeff -= 1
      elif term == "-":
          # take away 1 if the coefficient is just a negative sign
          y_coeff += 1
      else:
          y_coeff -= float(term)  # If y has a coefficient then add it to the running tally
    else:
      const -= float(term)  # Add to the running tally for the constant
  for term in left:
    if ("x" in term) or ("X" in term):
      # Take out the X to get the coefficient
      term = term.replace("x", "")
      term = term.replace("X", "")
      if term == "":
          # If it doesn't have a coefficient than set it to one
          x_coeff += 1
      elif term == "-":
          # take away 1 if the coefficient is just a negative sign
          x_coeff += -1
      else:
          x_coeff += float(term)  # If x has a coefficient then add it to the running tally
    elif ("y" in term) or ("Y" in term):
      # Take out the X to get the coefficient
      term = term.replace("y", "")
      term = term.replace("Y", "")
      if term == "":
          # If it doesn't have a coefficient than set it to one
          y_coeff += 1
      elif term == "-":
          # take away 1 if the coefficient is just a negative sign
          y_coeff += -1
      else:
          y_coeff += float(term)  # If y has a coefficient then add it to the running tally
    else:
      const += float(term)  # Add to the running tally for the constant
  return (x_coeff, y_coeff, const)

def parseEq(eq):
    eq = eq.replace(" ", "")  # Remove whitespace
    if "=" in eq:
        # Equation mode
        right, left = eq.split("=")  # Split string into left and right sides of the equal sign
        # Seperate terms on each side of the equation
        left = " ".join(left.split("+"))
        left = left.replace("-", " -")
        right = " ".join(right.split("+"))
        right = right.replace("-", " -")
        right = right.split()
        left = left.split()
        return simplifyTerms(right), simplifyTerms(left)
    else:
        # Expression mode
        eq = " ".join(eq.split("+"))
        eq = eq.replace("-", " -")
        eq = eq.split()
        return eq


def quadratic_formula(a, b, c):
    discriminant = b ** 2 - (4 * a * c)
    if discriminant >= 0:
        x_1 = ((-1 * b) + math.sqrt(discriminant)) / (2 * a)
        x_2 = ((-1 * b) - math.sqrt(discriminant)) / (2 * a)
    else:
        i_coeff = math.sqrt(abs(discriminant))  # The coefficient of i if the solution is complex
        if i_coeff >= 0:
            # If the i coefficient is positive
            if i_coeff / (2 * a) == 1:
                # If it is one
                x_1 = f"{(-1 * b) / (2 * a)}+i"
                x_2 = f"{(-1 * b) / (2 * a)}-i"
            else:
                x_1 = f"{(-1 * b) / (2 * a)}+{i_coeff / (2 * a)}i"
                x_2 = f"{(-1 * b) / (2 * a)}-{i_coeff / (2 * a)}i"
        else:
            # If the coefficient is negative
            x_1 = f"{(-1 * b) / (2 * a)}{i_coeff / (2 * a)}i"
            x_2 = f"{(-1 * b) / (2 * a)}+{(i_coeff / (2 * a)) * -1}i"
    return x_1, x_2


def exponential_solver(y, b, a=1, h=0, k=0):
    numerator = math.log((-1 * a) / (k - y)) - h * math.log(b) + (2 * math.pi)
    return f"{math.log((-1 * a) / (k - y)) - h * math.log(b) / math.log(b)}+{(2 * math.pi) / math.log(b)}i n"


def solve(equation):
    if "x**3" in equation or "x^3" in equation:
        print("Cannot solve polynomials with a degree above 2")
        return 0
    if ("x**2" in equation) or ("x^2" in equation):
        left, right = parseEq(equation)
        right_cpy = right.copy()  # Copy the equation to avoid an infinite loop
        # Move terms to left side of the equal sign
        if right[0] is not None:
            if "-" in right[0]:
                # Add to both sides
                right.append(right[0][1:])
                left.append(right[0][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(right[0]))
                left.append("-" + str(right[0]))
        if right[1] != None:
            if "-" in right[1]:
                # Add to both sides
                right.append(right[1][1:])
                left.append(right[1][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(right[1]))
                left.append("-" + str(right[1]))
        if right[3] != None:
            if "-" in right[3]:
                print(right[3])
                # Add to both sides
                right.append(right[3][1:] + "x^2")
                left.append(right[3][1:] + "x^2")
            else:
                # Subtract from both sides
                right.append("-" + str(right[3]) + "x^2")
                left.append("-" + str(right[3]) + "x^2")
        right, left = simplifyTerms(right), simplifyTerms(left)
        # Extract A, B, and C
        left[3] = left[3].replace("x^2", "")
        if left[3] == "":
            a = 1
        else:
            a = float(left[3].replace("x^2", ""))
        b = float(left[0].replace("x", ""))
        c = float(left[1])
        x_1, x_2 = quadratic_formula(a, b, c)  # Plug into the quadratic formula
        if x_1 == x_2:
            print(f"x = {x_1}")
        else:
            print(f"x = {x_1}")
            print(f"x = {x_2}")
    elif re.findall(pattern, equation) != []:
        print("Exponential detected")
        left, right = parseEq(equation)
        print(exponential_solver(float(left[1]), right[5], right[4], right[6], right[7]))
        print("n ∈ ℤ(ℤ is the set of integers)")
    else:
        left, right = parseEq(equation)
        # Goal: Get X on one side and the consts on the other side, then divide the coefficient by the constant
        if right[0] != None:
            # Isolate x on the left side of the equation
            if "-" in right[0]:
                # Add to both sides
                right.append(right[0][1:])
                left.append(right[0][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(right[0]))
                left.append("-" + str(right[0]))
            # Isolate the constants on the right side
            if "-" in left[1]:
                # Add to both sides
                right.append(left[1][1:])
                left.append(left[1][1:])
            else:
                # Subtract from both sides
                right.append("-" + str(left[1]))
                left.append("-" + str(left[1]))
        left, right = simplifyTerms(left), simplifyTerms(right)
        x_coeff = float(left[0].replace("x", ""))
        const = float(right[1])
        print(const / x_coeff)


def derivative(expression):
  og_expression = expression
  expression = parseEq(expression)
  new_expression = []
  for term in expression:
      if re.findall(r'(\d+(\*?))?x\^(\d+)', term):
          # Power rule
          term = term.split("x^")
          if term[0] == "":
              a = 1
          else:
              a = int(term[0])
          b = int(term[1])
          a = a * b
          b = b - 1
          if b == 1:
              b = "x"
          elif b == 0:
              b = ""
          else:
              b = f"x^{b}"
          if a == 1:
              a = ""
          term = f"{a}{b}"
          new_expression.append(term)
          # Special cases
      elif term == "log x" or term == "log(x)":
          new_expression.append("1/x")
      elif term == "sin x" or term == "sin(x)":
          new_expression.append("cos(x)")
      elif term == "tan x" or term == "tan(x)":
        new_expression.append("sec^2(x)")
      elif term == "sinh x" or term == "sinh(x)":
        new_expression.append("cosh(x)")
      elif term == "cosh x" or term == "cosh(x)":
        new_expression.append("sinh(x)")
      elif term == "tanh x" or term == "tanh(x)":
        new_expression.append("sech^2(x)")
      elif term == "-acos x" or term == "-acos(x)" or term == "acos x" or term == "acos(x)" or term == "asin(x)" or term == "asin x":
        new_expression.append("1/√(1-x^2)")
      elif term == "atan x" or term == "atan(x)":
        new_expression.append("1/(x^2-1)")
      elif term == "acosh x" or term == "acosh(x)":
        new_expression.append("1/√(x^2+1)*√(x^2-1)")
      elif term == "asinh x" or term == "asinh(x)":
        new_expression.append("1/√(x^2+1)")
      elif term == "atanh x" or term == "atanh(x)":
        new_expression.append("1/(1 - x^2)")
      elif term == "sec^2 x" or term == "sec^2(x)":
          new_expression.append("2*tan(x)*sec^2(x)")
      elif term == "sech^2 x" or term == "sech^2(x)":
          new_expression.append("-2*tanh(x)*sech^2(x)")
      elif term == "e^x" or term == "e**x":
          new_expression.append("e^x")
      elif term == "1/x":
          new_expression.append("-1/x^2")
      elif term == "cos x" or term == "cos(x)":
          new_expression.append("sin(x)")
      elif term == "-log x" or term == "-log(x)":
          new_expression.append("-1/x")
      elif term == "-sin x" or term == "-sin(x)":
          new_expression.append("-cos(x)")
      elif term == "-tan x" or term == "-tan(x)":
        new_expression.append("-sec^2(x)")
      elif term == "-sinh x" or term == "-sinh(x)":
        new_expression.append("-cosh(x)")
      elif term == "-cosh x" or term == "-cosh(x)":
        new_expression.append("-sinh(x)")
      elif term == "-tanh x" or term == "-tanh(x)":
        new_expression.append("-sech^2(x)")
      elif term == "-asin(x)" or term == "-asin x":
        new_expression.append("-1/√(1-x^2)")
      elif term == "-atan x" or term == "-atan(x)":
        new_expression.append("-1/(x^2+1)")
      elif term == "-acosh x" or term == "-acosh(x)":
        new_expression.append("-1/√(x^2+1)*√(x^2-1)")
      elif term == "-asinh x" or term == "-asinh(x)":
        new_expression.append("-1/√(x^2+1)")
      elif term == "-atanh x" or term == "-atanh(x)":
        new_expression.append("1/(x^2-1)")
      elif term == "-sec^2 x" or term == "-sec^2(x)":
          new_expression.append("-2*tan(x)*sec^2(x)")
      elif term == "-sech^2 x" or term == "-sech^2(x)":
          new_expression.append("2*tanh(x)*sech^2(x)")
      elif term == "-e^x" or term == "-e**x":
          new_expression.append("-e^x")
      elif term == "-1/x":
          new_expression.append("1/x^2")
      elif term == "-cos x" or term == "-cos(x)":
          new_expression.append("sin(x)")
      elif re.findall(r'\d+(\*?)x', term) != []:
          # Regular slope
          if "*" in term:
              new_expression.append(term[:-2])
          else:
              new_expression.append(term[:-1])
  new_expression_str = ""
  for term in new_expression:
      if new_expression_str == "":
          new_expression_str = term
      else:
          if term[0] == "-":
              new_expression_str = new_expression_str + term
          else:
              new_expression_str = new_expression_str + "+" + term
  print(f"d/dx ({og_expression}) = " + new_expression_str)


def integrate(expression):
  og_expression = expression
  expression = parseEq(expression)
  new_expression = []
  for term in expression:
      if re.findall(r'(\d+(\*?))?x\^(\d+)', term):
          # Inverse power rule
          term = term.split("x^")
          if term[0] == "":
              a = 1
          else:
              a = int(term[0])
          b = int(term[1])
          a = a / (b + 1)
          b = b + 1
          if b == 1:
              b = "x"
          elif b == 0:
              b = ""
          else:
              b = f"x^{b}"
          if a == 1:
              a = ""
          term = f"{a}{b}"
          new_expression.append(term)
      elif re.match(r'(\-?)\d+x', term):
          print(term)
          term = term.replace("x", "")
          if term == "":
              new_expression.append("x^2/2")
          else:
              new_expression.append(f"{int(term) / 2}x^2")
      # Cases that can't be solved by the inverse power rule
      elif term == "log x" or term == "log(x)":
          new_expression.append("x(log(x)-1)")
      elif term == "sin x" or term == "sin(x)":
          new_expression.append("-cos(x)")
      elif term == "e^x" or term == "e**x":
          new_expression.append("e^x")
      elif term == "1/x":
          new_expression.append("log(x)")
      elif term == "cos x" or term == "cos(x)":
          new_expression.append("sin(x)")
      elif term == "tan x" or term == "tan(x)":
        new_expression.append("-log(cos((x))")
      elif term == "sinh x" or term == "sinh(x)":
        new_expression.append("cosh(x)")
      elif term == "cosh x" or term == "cosh(x)":
        new_expression.append("sinh(x)")
      elif term == "tanh x" or term == "tanh(x)":
        new_expression.append("log(cosh(x))")
      elif term == "acos x" or term == "acos(x)":
        new_expression.append("x*acos(x)-√(1-x^2)")
      elif term == "asin x" or term == "asin(x)":
        new_expression.append("√(1-x^2)+x*asin(x)")
      elif term == "atan x" or term == "atan(x)":
        new_expression.append("x*atan(x)-0.5*log(x^2+1)")
      elif term == "acosh x" or term == "acosh(x)":
        new_expression.append("x*asinh(x)-√(x^2+1)*√(x^2-1)")
      elif term == "asinh x" or term == "asinh(x)":
        new_expression.append("x*asinh(x)-√(x^2+1)")
      elif term == "atanh x" or term == "atanh(x)":
        new_expression.append("0.5*log(1-x^2)+x*atanh(x)")
      elif term == "sec^2 x" or term == "sec^2(x)":
          new_expression.append("tan(x)")
      elif term == "sech^2 x" or term == "sech^2(x)":
          new_expression.append("tanh(x)")
      elif term == "-acos x" or term == "-acos(x)":
          new_expression.append("√(1-x^2)-x*acos(x)")
      elif term == "-log x" or term == "-log(x)":
          new_expression.append("x-x*log(x)")
      elif term == "-sin x" or term == "-sin(x)":
          new_expression.append("cos(x)")
      elif term == "-tan x" or term == "-tan(x)":
        new_expression.append("log(cos(x))")
      elif term == "-sinh x" or term == "-sinh(x)":
        new_expression.append("-cosh(x)")
      elif term == "-cosh x" or term == "-cosh(x)":
        new_expression.append("-sinh(x)")
      elif term == "-tanh x" or term == "-tanh(x)":
        new_expression.append("-log(cos(x))")
      elif term == "-asin(x)" or term == "-asin x":
        new_expression.append("-√(1-x^2)-x*asin(x)")
      elif term == "-atan x" or term == "-atan(x)":
        new_expression.append("0.5*log(x^2+1)-x*atan(x)")
      elif term == "-acosh x" or term == "-acosh(x)":
        new_expression.append("√(x^2+1)*√(x^2-1)-x*acosh(x)")
      elif term == "-asinh x" or term == "-asinh(x)":
        new_expression.append("√(x^2+1)-x*asinh(x)")
      elif term == "-atanh x" or term == "-atanh(x)":
        new_expression.append("-0.5*log(1-x^2)-x*atanh(x)")
      elif term == "-sec^2 x" or term == "-sec^2(x)":
          new_expression.append("-tan(x)")
      elif term == "-sech^2 x" or term == "-sech^2(x)":
          new_expression.append("-tanh(x)")
      elif term == "-e^x" or term == "-e**x":
          new_expression.append("-e^x")
      elif term == "-1/x":
          new_expression.append("-log(x)")
      elif term == "-cos x" or term == "-cos(x)":
          new_expression.append("-sin(x)")
      elif re.match(r'(\-?)\d+', term):
          # Reverse slope
          new_expression.append(term + "x")
  # Combine all the terms into a single string
  new_expression_str = ""
  for term in new_expression:
      if new_expression_str == "":
          new_expression_str = term
      else:
          if term[0] == "-":
              new_expression_str = new_expression_str + term
          else:
              new_expression_str = new_expression_str + "+" + term
  print(f"∫ ({og_expression}) dx = {new_expression_str} + C")


def scale_vector(scalar, vector):
    return list(map(lambda x: x * scalar, vector))


def dot_product(vector1, vector2):
    return sum([a * b for a, b in zip(vector1, vector2)])


def matrix_parse(matrix):
    matrix = matrix.split(";")
    matrixlist = []
    def eval_(node):
        if isinstance(node, ast.Constant):  # <number>
            return node.value
        elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
            return operators[type(node.op)](eval_(node.left), eval_(node.right))
        elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
            return operators[type(node.op)](eval_(node.operand))
        else:
            raise TypeError(node)
    for row in matrix:
        if row[0] == "[":
            row = row[1:]
        if row[len(row) - 1] == "]":
            row = row[:len(row) - 1]
        row = "[" + row + "]"
        row = list(map(lambda x: eval_(ast.parse(x, mode='eval').body), row.strip('][').split(', ')))
        #row = list(map(lambda x: complex(x) if "j" in x else int(x), row))
        matrixlist.append(row)
    return matrixlist


def verify_square(input_matrix):
    rows = len(input_matrix)
    first_row_len = len(input_matrix[0])
    if first_row_len != rows:
        return False
    else:
        for row in input_matrix:
            if len(row) != first_row_len:
                return False
            else:
                continue
        return True


def transpose(matrix):
    new_matrix = []
    for col in range(len(matrix[0])):
        new_matrix.append([])
    for i in range(len(new_matrix)):
        for row in matrix:
            new_matrix[i].append(row[i])
    return new_matrix


def fft(x):
    N = len(x)

    if N == 1:
        return [x[0]]

    transformed_x = [0] * N

    even = fft(x[:N:2])
    odd = fft(x[1:N:2])

    for k in range(N // 2):
        w = cmath.exp(-2j * math.pi * k / N)
        transformed_x[k] = even[k] + w * odd[k]
        transformed_x[k + N // 2] = even[k] - w * odd[k]
    return transformed_x


def ifft(x):
    N = len(x)

    if N == 1:
        return [x[0]]

    transformed_x = [0] * N

    even = ifft(x[:N:2])
    odd = ifft(x[1:N:2])

    for k in range(N // 2):
        w = cmath.exp(2j * math.pi * k / N)
        transformed_x[k] = (even[k] + w * odd[k]) / N
        transformed_x[k + N // 2] = (even[k] - w * odd[k]) / N
    return transformed_x


def pretty_print_matrix(matrix):
    if type(matrix) != list:
        matrix = matrix[0].tolist()
    for i in range(len(matrix)):
        if i == 0:
            print(f"[ {str(matrix[0]).replace('[', '').replace(']', '')}")
        elif i == len(matrix) - 1:
            print(f"  {str(matrix[i]).replace('[', '').replace(']', '')} ]")
        else:
            print(f"  {str(matrix[i]).replace('[', '').replace(']', '')}")


def conjugate(num):
    new_matrix = [[x.conjugate() if type(x) == complex else x for x in row] for row in matrix]
    return new_matrix


def eval_(node):
    if isinstance(node, ast.Constant):  # <number>
        return node.value
    elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
        return operators[type(node.op)](eval_(node.left), eval_(node.right))
    elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
        return operators[type(node.op)](eval_(node.operand))
    else:
        raise TypeError(node)

def matrix_multiplication(A, B):
    # Check if the matrices can be multiplied
    if len(A[0]) != len(B):
        raise ValueError("Matrix dimensions are not compatible for multiplication")

    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]

    return result

ascii_art = """  ______                   __                                __                          
 /      \                 /  |                              /  |                         
/$$$$$$  |  ______    ____$$ |  ______    ______    ______  $$ |____    ______   ______  
$$ |  $$/  /      \  /    $$ | /      \  /      \  /      \ $$      \  /      \ /      \ 
$$ |      /$$$$$$  |/$$$$$$$ |/$$$$$$  |/$$$$$$  |/$$$$$$  |$$$$$$$  |/$$$$$$  |$$$$$$  |
$$ |   __ $$ |  $$ |$$ |  $$ |$$    $$ |$$ |  $$ |$$    $$ |$$ |  $$ |$$ |  $$/ /    $$ |
$$ \__/  |$$ \__$$ |$$ \__$$ |$$$$$$$$/ $$ \__$$ |$$$$$$$$/ $$ |__$$ |$$ |     /$$$$$$$ |
$$    $$/ $$    $$/ $$    $$ |$$       |$$    $$ |$$       |$$    $$/ $$ |     $$    $$ |
 $$$$$$/   $$$$$$/   $$$$$$$/  $$$$$$$/  $$$$$$$ | $$$$$$$/ $$$$$$$/  $$/       $$$$$$$/ 
                                        /  \__$$ |                                       
                                        $$    $$/                                        
                                         $$$$$$/                                         """

print(ascii_art)
print("Created in 2023 by Samir")

while True:
    cmd = input("COMMAND>")
    if cmd == "SOLVE":
        equation = input("EQUATION TO BE SOLVED>")
        solve(equation)
    elif cmd == "DIFF":
        equation = input("EXPRESSION TO BE DIFFERENTIATED>")
        derivative(equation)
    elif cmd == "SCALE":
        vector = input("VECTOR>")
        vector = list(map(int, vector.strip('][').split(', ')))
        scalar = int(input("SCALAR>"))
        print(scale_vector(scalar, vector))
    elif cmd == "DOTPR":
        vectorA = input("VECTOR A>")
        vectorA = list(map(float, vectorA.strip('][').split(', ')))

        vectorB = input("VECTOR B>")
        vectorB = list(map(float, vectorB.strip('][').split(', ')))

        print(dot_product(vectorA, vectorB))
    elif cmd == "ADDVEC":
        vectorA = input("VECTOR A>")
        vectorA = list(map(float, vectorA.strip('][').split(', ')))

        vectorB = input("VECTOR B>")
        vectorB = list(map(float, vectorB.strip('][').split(', ')))

        print(blas.saxpy(vectorA, vectorB))
    elif cmd == "EIGVAL":
        matrix = input("Matrix>")
        matrix = matrix_parse(matrix)
        eigenvalues = lapack.sgeev(matrix)
        print("EIGENVALUES:")
        for eigenvalue in eigenvalues:
            print(eigenvalue)
    elif cmd == "TRANSP":
        matrix = input("Matrix>")
        matrix = matrix_parse(matrix)
        matrix = transpose(matrix)
        pretty_print_matrix(matrix)
    elif cmd == "FFT":
        signal = input("SIGNAL(AS VECTOR)>")
        signal = list(map(int, signal.strip('][').split(', ')))
        print(fft(signal))
    elif cmd == "IFFT":
        signal = input("SIGNAL(AS VECTOR)>")
        signal = list(map(complex, signal.strip('][').split(', ')))
        print(ifft(signal))
    elif cmd == "INTEGRATE":
        expression = input("EXPRESSION TO INTEGRATE>")
        integrate(expression)
    elif cmd == "SUM":
        vector = input("VECTOR>")
        vector = list(map(int, vector.strip('][').split(', ')))
        print(sum(vector))
    elif cmd == "AVG":
        vector = input("VECTOR>")
        vector = list(map(int, vector.strip('][').split(', ')))
        print(sum(vector) / len(vector))
    elif cmd == "WORD":
        word = input("WORD>")
        df = pd.read_csv('english Dictionary.csv')
        rows = df.loc[df['word'] == word]
        print(f"DEFINITION OF {word}")
        i = 1
        for row in rows["def"]:
            print(row)
            i += 1
            print(f"{i}. - {row}")
    elif cmd == "EXP":



        expression = input("EXPONENT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.exp(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.exp(i) for i in j] for j in expression])
        else:
            print(cmath.exp(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "SQRT":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.sqrt(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.sqrt(i) for i in j] for j in expression])
        else:
            print(cmath.sqrt(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "SIN":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.sin(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.sin(i) for i in j] for j in expression])
        else:
            print(cmath.sin(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "COS":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.cos(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.cos(i) for i in j] for j in expression])
        else:
            print(cmath.cos(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "TAN":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.tan(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.tan(i) for i in j] for j in expression])
        else:
            print(cmath.tan(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ATAN":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.atan(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.atan(i) for i in j] for j in expression])
        else:
            print(cmath.atan(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ASIN":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.asin(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.asin(i) for i in j] for j in expression])
        else:
            print(cmath.asin(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ACOS":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.acos(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.acos(i) for i in j] for j in expression])
        else:
            print(cmath.acos(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "LOG":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.log(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.log(i) for i in j] for j in expression])
        else:
            print(cmath.log(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ABS":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([abs(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[abs(i) for i in j] for j in expression])
        else:
            print(abs(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "DETR":
        matrix = input("MATRIX>")
        matrix = matrix_parse(matrix)
        if verify_square(matrix):
            if len(matrix) == 2:
                determinant = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]
                print(determinant)
            elif len(matrix) == 3:
                # aei+bfg+cdh-ceg-bdi-afh
                determinant = matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] + \
                              matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0] - \
                              matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1]
                print(determinant)
            else:
                print("ERROR: CAN ONLY HANDLE 2x2 or 3x3 MATRICES")
        else:
            print("ERROR: CAN ONLY HANDLE SQUARE MATRICES")
    elif cmd == "CONJ":
        matrix = input("MATRIX>")
        matrix = matrix_parse(matrix)
        pretty_print_matrix(conjugate(matrix))
    elif cmd == "ROUND":
        number = round(float(input("NUMBER>")))
        print(number)
    elif cmd == "EYE":
        size = input("LENGTH>")
        id_matrix = []
        for i in range(int(size)):
            row = []
            for q in range(int(size)):
                if i == q:
                    row.append(1)
                else:
                    row.append(0)
            id_matrix.append(row)
        pretty_print_matrix(id_matrix)
    elif cmd == "EVAL":



        expression = input("EXPRESSION>")
        print(eval_(ast.parse(expression, mode='eval').body))
    elif cmd == "PROD":
        vector = input("VECTOR>")
        vector = list(map(int, vector.strip('][').split(', ')))
        prod = vector[0]
        for element in vector[1:]:
            prod = prod * element
        print(prod)
    elif cmd == "ELEMENT":
        number = input("ATOMIC NUMBER>")
        df = pd.read_csv('elements.csv')
        row = df.loc[df['AtomicNumber'] == int(number)]
        print(row['Element'].tolist()[0])
        print(f"Atomic Mass: {row['AtomicMass'].tolist()[0]}")
        print(f"Neutrons: {row['NumberofNeutrons'].tolist()[0]}")
        print(f"Protons: {row['NumberofProtons'].tolist()[0]}")
        print(f"Electrons: {row['NumberofElectrons'].tolist()[0]}")
        print(f"Type: {row['Type'].tolist()[0]}")
        print(f"Radius: {row['AtomicRadius'].tolist()[0]} pm")
        print(f"Density: {row['Density'].tolist()[0]} g/cm3")
        if not math.isnan(row['MeltingPoint'].tolist()[0]):
            print(f"Melts at: {row['MeltingPoint'].tolist()[0]}°K")
        if not math.isnan(row['BoilingPoint'].tolist()[0]):
            print(f"Boils at: {row['BoilingPoint'].tolist()[0]}°K")
    elif cmd == "MOVIE":
        print("WARNING: DATABASE NOT UPDATED SINCE 2013")
        title = input("TITLE>")
        df = pd.read_csv('movies.csv')
        row = df.loc[df['Movie'] == title]
        print(f"Studio: {row['LeadStudio'].tolist()[0]}")
        print(f"Rotten Tomatoes: {row['RottenTomatoes'].tolist()[0]}%")
        print(f"Genre: {row['Genre'].tolist()[0]}")
        print(f"Year: {row['Year'].tolist()[0]}")
    elif cmd == "INV":
        matrix = input("Matrix>")
        matrix = matrix_parse(matrix)
        LU, PIV, _ = lapack.sgetrf(matrix)
        inv_a = lapack.sgetri(LU, PIV)
        pretty_print_matrix(inv_a)
    elif cmd == "SINH":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.sinh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.sinh(i) for i in j] for j in expression])
        else:
            print(cmath.sinh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "COSH":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.cosh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.cosh(i) for i in j] for j in expression])
        else:
            print(cmath.cosh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "TANH":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.tanh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.tanh(i) for i in j] for j in expression])
        else:
            print(cmath.tanh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ATANH":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.atanh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.atanh(i) for i in j] for j in expression])
        else:
            print(cmath.atanh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ASINH":



        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.asinh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.asinh(i) for i in j] for j in expression])
        else:
            print(cmath.asinh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "ACOSH":

        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([cmath.acosh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[cmath.acosh(i) for i in j] for j in expression])
        else:
            print(cmath.acosh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "SEC":
        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([1 / cmath.cos(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[1 / cmath.cos(i) for i in j] for j in expression])
        else:
            print(1 / cmath.cos(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "SECH":
        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([1 / cmath.cosh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[1 / cmath.cosh(i) for i in j] for j in expression])
        else:
            print(1 / cmath.cosh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "CSC":
        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([1 / cmath.sin(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[1 / cmath.sin(i) for i in j] for j in expression])
        else:
            print(1 / cmath.sin(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "CSCH":
        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([1 / cmath.sinh(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[1 / cmath.sinh(i) for i in j] for j in expression])
        else:
            print(1 / cmath.sinh(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "NORM":
        # Euclidean Normalization
        mode = input("1, 2, OR INFINITY NORM>")
        if mode == "1":
            vector = input("VECTOR>").strip('][').split(', ')
            vector = [abs(eval_(ast.parse(i, mode='eval').body)) for i in vector]
            print(sum(vector))
        elif mode == "INFINITY":
            vector = input("VECTOR>").strip('][').split(', ')
            vector = [abs(eval_(ast.parse(i, mode='eval').body)) for i in vector]
            print(max(vector))
        elif mode == "2":
            vector = input("VECTOR>").strip('][').split(', ')
            vector = [eval_(ast.parse(i, mode='eval').body) for i in vector]
            radicand = 0
            for j in vector:
                radicand += j ** 2
            print(math.sqrt(radicand))
    elif cmd == "MXV":
        matrix = input("MATRIX>")
        matrix = matrix_parse(matrix)
        vector = input("VECTOR>")
        vector = vector.strip('][').split(', ')
        vector = [eval_(ast.parse(i, mode='eval').body) for i in vector]
        nvec = blas.sgemv(1, matrix, vector)
        print(f"RESULT: {str(nvec)}")
    elif cmd == "ONES":
        n = int(input("N>"))
        m = int(input("M>"))
        ones = [[1 for i in range(m)] for j in range(n)]
        pretty_print_matrix(ones)
    elif cmd == "ZEROES":
        n = input("N>")
        m = input("M>")
        zeroes = [[1 for i in range(m)] for j in range(n)]
        pretty_print_matrix(zeroes)
    elif cmd == "UTRI":
        matrix = input("MATRIX>")
        matrix = matrix_parse(matrix)
        new_matrix = [[matrix[i][j] if j >= i else 0 for j in range(len(matrix[i]))] for i in range(len(matrix))]
        pretty_print_matrix(new_matrix)
    elif cmd == "LTRI":
        matrix = input("MATRIX>")
        matrix = matrix_parse(matrix)
        new_matrix = [[matrix[i][j] if j <= i else 0 for j in range(len(matrix[i]))] for i in range(len(matrix))]
        pretty_print_matrix(new_matrix)
    elif cmd == "FAC":
        expression = input("INPUT>")
        if "[" in expression and not ";" in expression:
            expression = expression.strip('][').split(', ')
            print([math.factorial(eval_(ast.parse(i, mode='eval').body)) for i in expression])
        elif "[" in expression and ";" in expression:
            expression = matrix_parse(expression)
            pretty_print_matrix([[math.factorial(i) for i in j] for j in expression])
        else:
            print(math.factorial(eval_(ast.parse(expression, mode='eval').body)))
    elif cmd == "RAT":
        num = input("DECIMAL TO BE RATIONALIZED>")
        num = num.replace(".", "")
        num = list(num)
        expr = f"{num[0]} + 1/"
        num = num[1:]
        for i in range(len(num)):
            if i != len(num)-1:
                expr += f"({num[i]} + 1/"
            else:
                expr += f"({num[i]}"
        expr += ")"*(len(num)-1)
        print(expr)
    elif cmd == "DIST":
        dims = input("DIMENSIONS(2 OR 3)>")
        if dims == "2":
            x1 = float(input("X1>"))
            y1 = float(input("Y1>"))
            x2 = float(input("X2>"))
            y2 = float(input("Y2>"))
            print(math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2))
        elif dims == "3":
            x1 = float(input("X1>"))
            y1 = float(input("Y1>"))
            z1 = float(input("Z1>"))
            x2 = float(input("X2>"))
            y2 = float(input("Y2>"))
            z2 = float(input("Z2>"))
            print(math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2))
        else:
            print("ERROR: DIMS NOT VALID")
    elif cmd == "PATH":
        dims = input("DIMENSIONS(2 OR 3)>")
        nodes = int(input("Number of nodes: "))
        print(nodes)
        if dims == "2":
            dist = 0
            x1 = float(input("X1>"))
            y1 = float(input("Y1>"))
            for i in range(nodes-1):
                x2 = float(input(f"X{i+2}>"))
                y2 = float(input(f"Y{i+2}>"))
                dist += math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                x1=x2
                y1=y2
            print("Dist:" + str(dist))
        elif dims == "3":
            dist = 0
            x1 = float(input("X1>"))
            y1 = float(input("Y1>"))
            z1 = float(input("Z1>"))
            for i in range(nodes):
                x2 = float(input(f"X{i+2}>"))
                y2 = float(input(f"Y{i+2}>"))
                z2 = float(input(f"Z{i+2}>"))
                dist += math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
                x1=x2
                y1=y2
            print("Dist:" + str(dist))
        else:
            print("ERROR: DIMS NOT VALID")
    elif cmd == "TRANS":
        dims = input("DIMENSIONS(2 OR 3)>")
        nodes = int(input("Number of points on the shape: "))
        if dims == "2":
          x = []
          y = []
          trans_x = float(input("ΔX>"))
          trans_y = float(input("ΔY>"))
          x.append(float(input("X1>")))
          y.append(float(input("Y1>")))
          for i in range(nodes-1):
            x.append(float(input(f"X{i+2}>")))
            y.append(float(input(f"Y{i+2}>")))
          print("TRANSLATED POINTS:")
          for i in range(len(x)):
            x[i] += trans_x
            y[i] += trans_y
            print(f"({x[i]}, {y[i]})")
        elif dims == "3":
          x = []
          y = []
          z = []
          trans_x = float(input("ΔX>"))
          trans_y = float(input("ΔY>"))
          trans_z = float(input("ΔZ>"))
          x.append(float(input("X1>")))
          y.append(float(input("Y1>")))
          z.append(float(input("Z1>")))
          for i in range(nodes-1):
            x.append(float(input(f"X{i+2}>")))
            y.append(float(input(f"Y{i+2}>")))
            z.append(float(input(f"Z{i+2}>")))
          print("TRANSLATED POINTS:")
          for i in range(len(x)):
            x[i] += trans_x
            y[i] += trans_y
            z[i] += trans_z
            print(f"({x[i]}, {y[i]}, {z[i]})")
        else:
          print("ERROR: DIMS NOT VALID")
    elif cmd == "PTOL":
      #https://math.stackexchange.com/questions/1013230/how-to-find-coordinates-of-reflected-point
      line = input("LINEAR EQUATION>")
      x_coeff, y_coeff, const = simpStdForm(line)
      x = float(input("X>"))
      y = float(input("Y>"))
      dist = abs((x_coeff*x)+(y_coeff*y)+const)/math.sqrt((x_coeff**2)+(y_coeff**2))
      print(f"Distance: {dist}")
    elif cmd == "SYST":
      #Use elimination to solve the equation
      eq1 = input("LINEAR EQUATION 1>")
      eq2 = input("LINEAR EQUATION 2>")
      x_coeff1, y_coeff1, const1 = simpStdForm(eq1)
      x_coeff2, y_coeff2, const2 = simpStdForm(eq2)
      y_coeff1_copy = y_coeff1
      x_coeff1_copy = x_coeff1
      const1_copy = const1
      y_coeff1 *= -1 * y_coeff2
      x_coeff1 *= -1 * y_coeff2
      const1 *= -1 * y_coeff2
      y_coeff2 *= y_coeff1_copy
      x_coeff2 *= y_coeff1_copy
      const2 *= y_coeff1_copy
      x_coeff = x_coeff1 + x_coeff2
      y_coeff = y_coeff1 + y_coeff2
      const = const1 + const2
      x = (const*-1)/x_coeff
      y = (-1 * const1_copy - x_coeff1_copy*x)/y_coeff1_copy
      print(f"x = {x}")
      print(f"y = {y}")
    elif cmd == "HELP":
        help_str = """Codegebra
A "computational intelligence system"(basically a fancy calculator that can also tell you data) that can solve equations, find derivatives, tell you about *some* movies, and more.
Motivation
I created this so that I could learn how systems like Matlab and Mathematica solve equations. However, I continued
to add more features after reaching my initial goal.
Commands
Note that commands are case-sensitive.
 - SOLVE: solve linear, quadratic, and exponential equations. Exponential equations must be in simplest for this to work
 - DIFF: Differentiate a function
 - SCALE: Scale a vector
 - DOTPR: Find the dot product of 2 vectors
 - ADDVEC: Find the sum of 2 vectors
 - EIGVAL: Find the eigenvalues of a matrix      
 - TRANSP: Transpose a matrix
 - FFT: Run a Discrete Fourier Transform on a signal
 - IFFT: Run an Inverse Fourier Transform on a signal
 - INTEGRATE: Integrate a function
 - SUM: The sum of a vector
 - AVG: The average of a vector
 - WORD: Define a word(similar to Wolfram|Alpha)
 - EXP: Calculate e^x where x is an expression, vector, matrix, or number
 - SQRT: Calculate the square root of an expression, vector, matrix, or number
 - SIN: Calculate the sine of an expression, vector, matrix, or number
 - COS: Calculate the cosine of an expression, vector, matrix, or number
 - TAN: Calculate the tangent of an expression, vector, matrix, or number
 - ATAN: Calculate the arc tangent of an expression, vector, matrix, or number
 - ASIN: Calculate the arc sine of an expression, vector, matrix, or number
 - ACOS: Calculate the arc cosine of an expression, vector, matrix, or number
 - TANH: Calculate the hyperbolic tangent of an expression, vector, matrix, or number
 - SINH: Calculate the hyperbolic sine of an expression, vector, matrix, or number
 - COSH: Calculate the hyperbolic cosine of an expression, vector, matrix, or number
 - ATANH: Calculate the inverse hyperbolic tangent of an expression, vector, matrix, or number
 - ASINH: Calculate the inverse hyperbolic sine of an expression, vector, matrix, or number
 - ACOSH: Calculate the inverse hyperbolic cosine of an expression, vector, matrix, or number
 - SEC: Calculate the secant of an expression, vector, matrix, or number
 - SECH: Calculate the hyperbolic secant of an expression, vector, matrix, or number
 - CSC: Calculate the cosecant of an expression, vector, matrix, or number
 - CSCH: Calculate the hyperbolic cosecant of an expression, vector, matrix, or number
 - LOG: Calculate the natural logarithm of an expression, vector, matrix, or number
 - ABS: Calculate the absolute value of an expression, vector, matrix, or number
 - DETR: Calculate the determinant of a 2x2 or 3x3 matrix
 - CONJ: Conjugate a  matrix
 - ROUND: Round a number to the nearest whole number
 - EYE: Create an identity matrix(a matrix with ones going diagonally and zero everywhere else)
 - EVAL: Evaluate an expression
 - PROD: Product of all elements in an array
 - ELEMENT: Get info about an element in the periodic table(by atomic number)
 - MOVIE: Get info about a movie(the data source I used had not been updated since 2013)
 - INV: Compute the inverse of a square matrix
 - MXV: Matrix-vector multiplication
 - ONES: An NxM matrix of all ones
 - ZEROES: An NxM matrix of all zeroes
 - LTRI: Isolate the lower triangular portion of a square matrix and set the upper part to zeroes
 - UTRI: Isolate the upper triangular portion of a square matrix and set the lower part to zeroes
 - FAC: Factorial of a number
 - RAT: Approximate the rational form of a number
 - DIST: Find the distance between 2 points in a 2d or 3d plane
 - PATH: Find the distance taken of a path that traverses user-specified points on a 2d or 3d plane
 - TRANS: Translate a 2d or 3d shape
 - PTOL: Distance from a point to a line
 - SYST: Solve a system of linear equations
Matrices
Matrices are written in the following format:
[1, 2, 3;4, 5, 6]
The comma separates elements in the matrix, and the semicolon seperates rows. Matrix elements can be imaginary, complex, or real numbers. They can also be expressions, which will be simplified at runtime.
Data credits
 - [Movie Data](https://github.com/reisanar/datasets/blob/master/HollywoodMovies.csv)
 - [Elements Data](https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee)
 - [English Dictionary Data](https://github.com/benjihillard/English-Dictionary-Database)"""
        print(help_str)
    else:
        print("ERROR: UNRECOGNIZED COMMAND")