# type: ignore
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2019 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
from math import sqrt

REQUIRE_DEPENDENCIES = False
if not REQUIRE_DEPENDENCIES:
    IS_PYPY = True
else:
    try:
        # The right way imports the platform module which costs to ms to load!
        # implementation = platform.python_implementation()
        IS_PYPY = "PyPy" in sys.version
    except AttributeError:
        IS_PYPY = False

#IS_PYPY = True # for testing

#if not IS_PYPY and not REQUIRE_DEPENDENCIES:
#    try:
#        import numpy as np
#    except ImportError:
#        np = None

__all__ = [
    "argsort1d",
    "array_as_tridiagonals",
    "det",
    "dot_product",
    "eye",
    "gelsd",
    "inv",
    "lu",
    "matrix_multiply",
    "matrix_vector_dot",
    "norm2",
    "null_space",
    "scalar_add_matrices",
    "scalar_divide_matrix",
    "scalar_multiply_matrix",
    "scalar_subtract_matrices",
    "shape",
    "solve",
    "solve_tridiagonal",
    "sort_paired_lists",
    "stack_vectors",
    "subset_matrix",
    "sum_matrix_cols",
    "sum_matrix_rows",
    "transpose",
]
primitive_containers = frozenset([list, tuple])

def transpose(matrix):
    """Convert a matrix into its transpose by switching rows and columns.

    Parameters
    ----------
    matrix : list[list[float]]
        Input matrix as a list of lists where each inner list represents a row.
        All rows must have the same length.

    Returns
    -------
    list[list[float]]
        The transposed matrix where element [i][j] in the input becomes [j][i]
        in the output.

    Raises
    ------
    ValueError
        If the input matrix has inconsistent row lengths.
    TypeError
        If the input is not a list of lists.

    Examples
    --------
    >>> transpose([[1, 2, 3], [4, 5, 6]])
    [[1, 4], [2, 5], [3, 6]]

    >>> transpose([[1, 2], [3, 4]])  # Square matrix
    [[1, 3], [2, 4]]

    >>> transpose([[1, 2, 3]])  # Single row matrix
    [[1], [2], [3]]

    Notes
    -----
    - Empty matrices are preserved as empty lists
    - The function creates a new matrix rather than modifying in place
    - For an MxN matrix, the result will be an NxM matrix
    """
    pass

def det(matrix):
    """Seems to work fine.

    >> from sympy import *
    >> from sympy.abc import *
    >> Matrix([[a, b], [c, d]]).det()
    a*d - b*c
    >> Matrix([[a, b, c], [d, e, f], [g, h, i]]).det()
    a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g

    A few terms can be slightly factored out of the 3x dim.

    >> Matrix([[a, b, c, d], [e, f, g, h], [i, j, k, l], [m, n, o, p]]).det()
    a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o - d*f*k*m - d*g*i*n + d*g*j*m

    72 mult vs ~48 in cse'd version'

    Commented out - takes a few seconds
    >> #Matrix([[a, b, c, d, e], [f, g, h, i, j], [k, l, m, n, o], [p, q, r, s, t], [u, v, w, x, y]]).det()

    260 multiplies with cse; 480 without it.
    """
    pass

# The inverse function below is generated via the following script
# import sympy as sp
# import re
# from sympy import Matrix, Symbol, simplify, zeros, cse

# def replace_power_with_multiplication(match):
#     """Replace x**n with x*x*...*x n times"""
#     var = match.group(1)
#     power = int(match.group(2))
#     if power <= 1:
#         return var
#     return '*'.join([var] * power)

# def generate_symbolic_matrix(n):
#     """Generate an nxn symbolic matrix with unique symbols"""
#     syms = [[Symbol(f'm_{i}{j}') for j in range(n)] for i in range(n)]
#     return Matrix(syms), syms

# def analyze_matrix(n):
#     """Generate symbolic expressions for determinant and inverse"""
#     M, syms = generate_symbolic_matrix(n)
#     det = M.det()
#     inv = M.inv()
#     return det, inv, syms

# def post_process_code(code_str):
#     """Apply optimizing transformations to the generated code"""
#     # Replace x**n patterns with x*x*x... (n times)
#     code_str = re.sub(r'([a-zA-Z_][a-zA-Z0-9_]*)\*\*(\d+)', replace_power_with_multiplication, code_str)
#     # Replace **0.5 with sqrt()
#     code_str = re.sub(r'\((.*?)\)\*\*0\.5', r'sqrt(\1)', code_str)
#     return code_str

# def generate_python_inv():
#     """Generate a single unified matrix inversion function with optimized 1x1, 2x2, and 3x3 cases"""
#     # Generate the specialized code for 2x2 and 3x3
#     size_specific_code = {}
#     for N in [2, 3, 4]:
#         det, inv, _ = analyze_matrix(N)
#         exprs = [det] + list(inv)
#         replacements, reduced = cse(exprs, optimizations='basic')
#         det_expr = reduced[0]
#         inv_exprs = reduced[1:]

#         # Build the size-specific code block
#         code = []

#         # Unpack matrix elements
#         unpack_rows = []
#         for i in range(N):
#             row_vars = [f"m_{i}{j}" for j in range(N)]
#             unpack_rows.append("(" + ", ".join(row_vars) + ")")
#         code.append(f"        {', '.join(unpack_rows)} = matrix")

#         # Common subexpressions
#         code.append("\n        # Common subexpressions")
#         for i, (temp, expr) in enumerate(replacements):
#             code.append(f"        x{i} = {expr}")

#         # Determinant check
#         code.append("\n        # Calculate determinant and check if we need to use LU decomposition")
#         code.append(f"        det = {det_expr}")
#         code.append("        if abs(det) <= 1e-7:")
#         code.append("            return inv_lu(matrix)")

#         # Return matrix
#         return_matrix = []
#         for i in range(N):
#             row = []
#             for j in range(N):
#                 idx = i * N + j
#                 row.append(str(inv_exprs[idx]))
#             return_matrix.append(f"            [{', '.join(row)}]")

#         code.append("\n        return [")
#         code.append(",\n".join(return_matrix))
#         code.append("        ]")

#         size_specific_code[N] = post_process_code("\n".join(code))

#     # Generate the complete function
#     complete_code = [
#         "def inv(matrix):",
#         "    size = len(matrix)",
#         "    if size == 1:",
#         "        return [[1.0/matrix[0][0]]]",
#         "    elif size == 2:",
#         size_specific_code[2],
#         "    elif size == 3:",
#         size_specific_code[3],
#         "    elif size == 4:",
#         size_specific_code[4],
#         "    else:",
#         "        return inv_lu(matrix)",
#         ""
#     ]

#     return "\n".join(complete_code)

# # Generate and print the complete function
# print(generate_python_inv())


def inv(matrix):
    pass

def shape(value):
    """Find and return the shape of an array, whether it is a numpy array or
    a list-of-lists or other combination of iterators.

    Parameters
    ----------
    value : various
        Input array, [-]

    Returns
    -------
    shape : tuple(int, dimension)
        Dimensions of array, [-]

    Notes
    -----
    It is assumed the shape is consistent - not something like [[1.1, 2.2], [2.4]]

    Examples
    --------
    >>> shape([])
    (0,)
    >>> shape([1.1, 2.2, 5.5])
    (3,)
    >>> shape([[1.1, 2.2, 5.5], [2.0, 1.1, 1.5]])
    (2, 3)
    >>> shape([[[1.1,], [2.0], [1.1]]])
    (1, 3, 1)
    >>> shape(['110-54-3'])
    (1,)
    """
    pass

def eye(N, dtype=float):
    """
    Return a 2-D array with ones on the diagonal and zeros elsewhere.

    Parameters
    ----------
    N : int
        Number of rows and columns in the output matrix.
    dtype : type, optional
        The type of the array elements. Defaults to float.

    Returns
    -------
    list[list]
        A N x N matrix with ones on the diagonal and zeros elsewhere.

    Examples
    --------
    >>> eye(3)
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

    >>> eye(2, dtype=int)
    [[1, 0], [0, 1]]

    Notes
    -----
    This function creates an identity matrix similar to numpy's eye function,
    but implemented in pure Python using nested lists.

    Raises
    ------
    ValueError
        If N is not a positive integer.
    TypeError
        If N is not an integer or dtype is not a valid type.
    """
    pass


def dot_product(a, b):
    """
    Compute the dot product (also known as scalar product or inner product) of two vectors.

    Calculates sum(a[i] * b[i]) for i in range(len(a)).

    Parameters
    ----------
    a : list[float]
        First vector
    b : list[float]
        Second vector of same length as a

    Returns
    -------
    float
        The dot product of vectors a and b

    Examples
    --------
    >>> dot_product([1, 2, 3], [4, 5, 6])
    32.0
    >>> dot_product([1, 0], [0, 1])
    0.0

    Notes
    -----

    Raises
    ------
    ValueError
        If vectors are not the same length
    TypeError
        If inputs are not valid vector types
    """
    pass

def matrix_vector_dot(matrix, vector):
    """
    Compute the product of a matrix and a vector.

    Parameters
    ----------
    matrix : list[list[float]]
        Input matrix represented as a list of lists.
    vector : list[float]
        Input vector represented as a list of floats.

    Returns
    -------
    list[float]
        The result of the matrix-vector multiplication as a vector.

    Raises
    ------
    ValueError
        If the number of columns in the matrix does not match the length of the vector.
    TypeError
        If inputs are not valid matrix and vector types.

    Examples
    --------
    >>> matrix_vector_dot([[1, 2, 3], [4, 5, 6]], [1, 0, 1])
    [4, 10]
    >>> matrix_vector_dot([[1.0, 2.0], [3.0, 4.0]], [0, 1])
    [2.0, 4.0]
    """
    pass

def matrix_multiply(A, B):
    r"""Multiply two matrices using pure Python.

    Computes the matrix product C = AÂ·B where A is an mxp matrix and B is a pxn matrix,
    resulting in an mxn matrix C.

    Parameters
    ----------
    A : list[list[float]]
        First matrix as list of lists, with shape (m, p)
    B : list[list[float]]
        Second matrix as list of lists, with shape (p, n)

    Returns
    -------
    list[list[float]]
        Resulting matrix C with shape (m, n)

    Examples
    --------
    >>> A = [[1, 2], [3, 4]]
    >>> B = [[5, 6], [7, 8]]
    >>> matrix_multiply(A, B)
    [[19.0, 22.0], [43.0, 50.0]]

    Notes
    -----
    Uses a straightforward three-loop implementation optimized for pure Python:
    C[i,j] = sum(A[i,k] * B[k,j] for k in range(p))

    The implementation avoids repeated len() calls and list accesses by caching
    frequently used values.

    Raises
    ------
    ValueError
        If matrices have incompatible dimensions for multiplication
        If input matrices are empty or irregular (rows of different lengths)
    TypeError
        If A or B contains non-numeric values or is not a list of lists.
    """
    pass

def sum_matrix_rows(matrix):
    """Sum a 2D matrix along rows, equivalent to numpy.sum(matrix, axis=1).

    Parameters
    ----------
    matrix : list[list[float]]
        Input matrix as a list of lists where each inner list is a row

    Returns
    -------
    list[float]
        List containing the sum of each row

    Examples
    --------
    >>> sum_matrix_rows([[1, 2, 3], [4, 5, 6]])
    [6.0, 15.0]
    >>> sum_matrix_rows([[1], [2]])
    [1.0, 2.0]

    Notes
    -----
    For a matrix with shape (m, n), returns a list of length m
    where each element is the sum of the corresponding row.

    Raises
    ------
    ValueError
        If matrix is empty or has irregular row lengths
    TypeError
        If matrix is not a list of lists of numbers
    """
    pass

def sum_matrix_cols(matrix):
    """Sum a 2D matrix along columns, equivalent to numpy.sum(matrix, axis=0).

    Parameters
    ----------
    matrix : list[list[float]]
        Input matrix as a list of lists where each inner list is a row

    Returns
    -------
    list[float]
        List containing the sum of each column

    Examples
    --------
    >>> sum_matrix_cols([[1, 2, 3], [4, 5, 6]])
    [5.0, 7.0, 9.0]
    >>> sum_matrix_cols([[1], [2]])
    [3.0]

    Notes
    -----
    For a matrix with shape (m, n), returns a list of length n
    where each element is the sum of the corresponding column.

    Raises
    ------
    ValueError
        If matrix is empty or has irregular row lengths
    TypeError
        If matrix is not a list of lists of numbers
    """
    pass

def scalar_add_matrices(A, B):
    """Add two matrices element-wise.

    Computes the element-wise sum of two matrices of the same dimensions.

    Parameters
    ----------
    A : list[list[float]]
        First matrix as a list of lists.
    B : list[list[float]]
        Second matrix as a list of lists.

    Returns
    -------
    list[list[float]]
        Resulting matrix after element-wise addition.

    Examples
    --------
    >>> A = [[1.0, 2.0], [3.0, 4.0]]
    >>> B = [[5.0, 6.0], [7.0, 8.0]]
    >>> scalar_add_matrices(A, B)
    [[6.0, 8.0], [10.0, 12.0]]

    Raises
    ------
    ValueError
        If matrices A and B have different shapes or if they are empty.
    TypeError
        If A or B contains non-numeric values or is not a list of lists.
    """
    pass


def scalar_subtract_matrices(A, B):
    """Subtract two matrices element-wise.

    Computes the element-wise difference of two matrices of the same dimensions.

    Parameters
    ----------
    A : list[list[float]]
        First matrix as a list of lists.
    B : list[list[float]]
        Second matrix as a list of lists.

    Returns
    -------
    list[list[float]]
        Resulting matrix after element-wise subtraction.

    Examples
    --------
    >>> A = [[5.0, 6.0], [7.0, 8.0]]
    >>> B = [[1.0, 2.0], [3.0, 4.0]]
    >>> scalar_subtract_matrices(A, B)
    [[4.0, 4.0], [4.0, 4.0]]

    Raises
    ------
    ValueError
        If matrices A and B have different shapes or if they are empty.
    TypeError
        If A or B contains non-numeric values or is not a list of lists.
    """
    pass


def scalar_multiply_matrix(scalar, matrix):
    """Multiply a matrix by a scalar.

    Multiplies each element of the matrix by the specified scalar.

    Parameters
    ----------
    scalar : float
        Scalar value to multiply each element by.
    matrix : list[list[float]]
        Input matrix as a list of lists.

    Returns
    -------
    list[list[float]]
        Resulting matrix after scalar multiplication.

    Examples
    --------
    >>> matrix = [[1, 2], [3, 4]]
    >>> scalar_multiply_matrix(2.0, matrix)
    [[2.0, 4.0], [6.0, 8.0]]

    Raises
    ------
    ValueError
        If the input matrix is empty.
    TypeError
        If the matrix contains non-numeric values or is not a list of lists.
    """
    pass


def scalar_divide_matrix(scalar, matrix):
    """Divide a matrix by a scalar.

    Divides each element of the matrix by the specified scalar.

    Parameters
    ----------
    scalar : float
        Scalar value to divide each element by (cannot be zero).
    matrix : list[list[float]]
        Input matrix as a list of lists.

    Returns
    -------
    list[list[float]]
        Resulting matrix after scalar division.

    Examples
    --------
    >>> matrix = [[2, 4], [6, 8]]
    >>> scalar_divide_matrix(2.0, matrix)
    [[1.0, 2.0], [3.0, 4.0]]

    Raises
    ------
    ValueError
        If the input matrix is empty or if the scalar is zero.
    TypeError
        If the matrix contains non-numeric values or is not a list of lists.
    ZeroDivisionError
        If scalar is zero.
    """
    pass

def stack_vectors(vectors):
    """Stack a list of vectors into a matrix, similar to numpy.stack.

    Parameters
    ----------
    vectors : list[list[float]]
        List of vectors to stack into rows of a matrix

    Returns
    -------
    list[list[float]]
        Matrix where each row is one of the input vectors

    Examples
    --------
    >>> stack_vectors([[1, 2], [3, 4]])
    [[1, 2], [3, 4]]
    """
    pass







def inplace_LU(A, ipivot):
    pass

def solve_from_lu(A, pivots, b):
    pass

def solve_LU_decomposition(A, b):
    pass

def inv_lu(a):
    pass

def lu(A):
    """
    Compute LU decomposition of a matrix with partial pivoting.
    Returns P, L, U such that PA = LU

    Parameters
    ----------
        A: list of lists representing square matrix

    Returns
    -------
        P: permutation matrix as list of lists
        L: lower triangular matrix with unit diagonal as list of lists
        U: upper triangular matrix as list of lists
    """
    pass


'''Script to generate solve function. Note that just like in inv the N = 4 case has too much numerical instability.
import sympy as sp
from sympy import Matrix, Symbol, simplify, solve_linear_system
import re

def generate_symbolic_system(n):
    """Generate an nxn symbolic matrix A and n-vector b"""
    A = Matrix([[Symbol(f'a_{i}{j}') for j in range(n)] for i in range(n)])
    b = Matrix([Symbol(f'b_{i}') for i in range(n)])
    return A, b

def generate_cramer_solution(n):
    """Generate symbolic solution using Cramer's rule for small matrices"""
    A, b = generate_symbolic_system(n)
    det_A = A.det()

    # Solve for each variable using Cramer's rule
    solutions = []
    for i in range(n):
        # Create matrix with i-th column replaced by b
        A_i = A.copy()
        A_i[:, i] = b
        det_i = A_i.det()
        # Store numerator only - we'll multiply by inv_det later
        solutions.append(det_i)

    return det_A, solutions

def generate_python_solve():
    """Generate a unified matrix solve function with optimized 1x1, 2x2, and 3x3 cases"""
    size_specific_code = {}

    # Special case for N=1
    size_specific_code[1] = """        # Direct solution for 1x1
        return [b[0]/matrix[0][0]]"""

    # Generate specialized code for sizes 2 and 3
    for N in [2, 3]:
        det, solutions = generate_cramer_solution(N)

        code = []

        # Unpack matrix elements
        unpack_rows = []
        for i in range(N):
            row_vars = [f"a_{i}{j}" for j in range(N)]
            unpack_rows.append("(" + ", ".join(row_vars) + ")")
        code.append(f"        {', '.join(unpack_rows)} = matrix")

        # Unpack b vector
        code.append(f"        {', '.join(f'b_{i}' for i in range(N))} = b")

        # Calculate determinant
        det_expr = str(det)
        code.append("\n        # Calculate determinant")
        code.append(f"        det = {det_expr}")

        # Check for singular matrix
        code.append("\n        # Check for singular matrix")
        code.append("        if abs(det) <= 1e-7:")
        code.append("            return solve_LU_decomposition(matrix, b)")

        # Calculate solution
        code.append("\n        # Calculate solution")
        code.append("        inv_det = 1.0/det")

        # Generate solution expressions (multiply by inv_det, don't divide by det)
        solution_lines = []
        for i, sol in enumerate(solutions):
            solution_lines.append(f"        x_{i} = ({sol}) * inv_det")
        code.append("\n".join(solution_lines))

        # Return solution
        code.append("\n        return [" + ", ".join(f"x_{i}" for i in range(N)) + "]")

        size_specific_code[N] = "\n".join(code)

    # Generate the complete function
    complete_code = [
        "def solve(matrix, b):",
        "    size = len(matrix)",
        "    if size == 1:",
        size_specific_code[1],
        "    elif size == 2:",
        size_specific_code[2],
        "    elif size == 3:",
        size_specific_code[3],
        "    else:",
        "        return solve_LU_decomposition(matrix, b)",
        ""
    ]

    return "\n".join(complete_code)

# Generate and print the optimized solve function
print(generate_python_solve())
'''






def solve(matrix, b):
    pass

def norm2(arr):
    pass

def array_as_tridiagonals(arr):
    """Extract the three diagonals from a tridiagonal matrix.

    A tridiagonal matrix is a matrix that has nonzero elements only on the
    main diagonal, the first diagonal below this (subdiagonal), and the first
    diagonal above this (superdiagonal).

    Parameters
    ----------
    arr : list[list[float]]
        Square matrix in tridiagonal form, where elements not on the three
        main diagonals are zero

    Returns
    -------
    tuple[list[float], list[float], list[float]]
        Three lists containing:
        a: subdiagonal elements (length n-1)
        b: main diagonal elements (length n)
        c: superdiagonal elements (length n-1)

    Examples
    --------
    >>> arr = [[2, 1, 0], [1, 2, 1], [0, 1, 2]]
    >>> a, b, c = array_as_tridiagonals(arr)
    >>> a  # subdiagonal
    [1, 1]
    >>> b  # main diagonal
    [2, 2, 2]
    >>> c  # superdiagonal
    [1, 1]

    Notes
    -----
    For a matrix of size nxn, returns:
    - a[i] contains elements at position (i+1,i) for i=0..n-2
    - b[i] contains elements at position (i,i) for i=0..n-1
    - c[i] contains elements at position (i,i+1) for i=0..n-2

    No validation is performed to ensure the input matrix is actually tridiagonal.
    Elements outside the three diagonals are ignored.
    """
    pass


def tridiagonals_as_array(a, b, c, zero=0.0):
    r"""Construct a square matrix from three diagonals.

    Creates a tridiagonal matrix using the provided sub-, main, and super-diagonal
    elements. All other elements are set to zero.

    Parameters
    ----------
    a : list[float]
        Subdiagonal elements (length n-1)
    b : list[float]
        Main diagonal elements (length n)
    c : list[float]
        Superdiagonal elements (length n-1)
    zero : float, optional
        Value to use for non-diagonal elements. Defaults to 0.0

    Returns
    -------
    list[list[float]]
        Square matrix of size nxn where n is the length of b

    Examples
    --------
    >>> a = [1, 1]  # subdiagonal
    >>> b = [2, 2, 2]  # main diagonal
    >>> c = [1, 1]  # superdiagonal
    >>> tridiagonals_as_array(a, b, c)
    [[2, 1, 0.0], [1, 2, 1], [0.0, 1, 2]]

    Notes
    -----
    For output matrix M of size nxn:
    - a[i] becomes M[i+1][i] for i=0..n-2
    - b[i] becomes M[i][i] for i=0..n-1
    - c[i] becomes M[i][i+1] for i=0..n-2

    No validation is performed on input lengths. For correct results:
    - len(b) should be n
    - len(a) and len(c) should be n-1

    The function is the inverse of array_as_tridiagonals() when zero=0.0
    """
    pass

def solve_tridiagonal(a, b, c, d):
    """Solve a tridiagonal system of equations using the Thomas algorithm.

    Solves the equation system Ax = d where A is a tridiagonal matrix composed of
    diagonals a, b, and c. This is an efficient O(n) method also known as the
    tridiagonal matrix algorithm (TDMA).

    The system of equations has the form:
    b[0]x[0] + c[0]x[1] = d[0]
    a[i]x[i-1] + b[i]x[i] + c[i]x[i+1] = d[i], for i=1..n-2
    a[n-1]x[n-2] + b[n-1]x[n-1] = d[n-1]

    Parameters
    ----------
    a : list[float]
        Lower diagonal (subdiagonal) elements a[i] at (i+1,i), length n-1, [-]
    b : list[float]
        Main diagonal elements b[i] at (i,i), length n, [-]
    c : list[float]
        Upper diagonal (superdiagonal) elements c[i] at (i,i+1), length n-1, [-]
    d : list[float]
        Right-hand side vector, length n, [-]

    Returns
    -------
    x : list[float]
        Solution vector, length n, [-]

    Examples
    --------
    >>> # Solve the system:
    >>> # [9 -1  0] [x0]   [1]
    >>> # [-1 2 -1] [x1] = [0]
    >>> # [0 -1  2] [x2]   [1]
    >>> a = [-1, -1]  # lower diagonal
    >>> b = [9, 2, 2]  # main diagonal
    >>> c = [-1, -1]  # upper diagonal
    >>> d = [1, 0, 1]  # right hand side
    >>> solve_tridiagonal(a, b, c, d)
    [0.16, 0.44, 0.72]

    Notes
    -----
    The algorithm modifies the input arrays b and d in-place to save memory,
    but makes copies first to preserve the originals.


    The algorithm fails if any diagonal element becomes zero during elimination.

    This implementation uses the Thomas algorithm, which is a specialized form
    of Gaussian elimination that exploits the tridiagonal structure for O(n)
    efficiency.

    No validation is performed on input lengths. For correct results:
    - len(b) should be n
    - len(a), len(c) should be n-1
    - len(d) should be n
    where n is the size of the system.

    References
    ----------
    .. [1] "Tridiagonal matrix algorithm", Wikipedia,
           https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    """
    pass



def subset_matrix(whole, subset):
    pass

def argsort1d(arr):
    """
    Returns the indices that would sort a 1D list.

    Parameters
    ----------
    arr : list
        Input array [-]

    Returns
    -------
    indices : list[int]
        List of indices that sort the input array [-]

    Notes
    -----
    This function uses the built-in sorted function with a custom key to get the indices.
    Note this does not match numpy's sorting for nan and inf values.

    Examples
    --------
    >>> arr = [3, 1, 2]
    >>> argsort1d(arr)
    [1, 2, 0]
    """
    pass

def sort_paired_lists(list1, list2):
    """
    Sort two lists based on the values in the first list while maintaining
    the relationship between corresponding elements.

    Parameters
    ----------
    list1 : list
        First list that determines the sorting order
    list2 : list
        Second list that will be sorted according to list1's ordering

    Returns
    -------
    tuple
        A tuple containing (sorted_list1, sorted_list2)

    Raises
    ------
    ValueError
        If the lists have different lengths
    TypeError
        If either input is not a list

    Examples
    --------
    >>> temps = [300, 100, 200]
    >>> props = ['hot', 'cold', 'warm']
    >>> sort_paired_lists(temps, props)
    ([100, 200, 300], ['cold', 'warm', 'hot'])

    Notes
    -----
    This function maintains the one-to-one relationship between elements
    in both lists while sorting them based on list1's values.
    """
    pass

def svd(matrix):
    """Compute the singular value decomposition of a matrix.

    This function wraps numpy.linalg.svd but maintains pure Python input/output
    interfaces.

    Parameters
    ----------
    matrix : list[list[float]]
        Input matrix A to decompose

    Returns
    -------
    tuple[list[list[float]], list[float], list[list[float]]]
        Returns (U, s, Vt) where:
        - U is the left singular vectors as a matrix
        - s is the singular values as a 1D array
        - Vt is the transpose of the right singular vectors as a matrix

    Notes
    -----

    Examples
    --------
    >>> A = [[1, 2], [3, 4]]
    >>> U, s, Vt = svd(A)
    """
    pass


def gelsd(a, b, rcond=None):
    """Solve a linear least-squares problem using SVD (Singular Value Decomposition).
    This is a simplified implementation that uses numpy's SVD internally.

    The function solves the equation arg min(|b - Ax|) for x, where A is
    an M x N matrix and b is a length M vector.

    Parameters
    ----------
    a : list[list[float]]
        Input matrix A of shape (M, N)
    b : list[float]
        Input vector b of length M
    rcond : float, optional
        Cutoff ratio for small singular values. Singular values smaller
        than rcond * largest_singular_value are considered zero.
        Default: max(M,N) * eps where eps is the machine precision

    Returns
    -------
    x : list[float]
        Solution vector of length N
    residuals : float
        Sum of squared residuals of the solution. Only computed for overdetermined
        systems (M > N)
    rank : int
        Effective rank of matrix A
    s : list[float]
        Singular values of A in descending order

    Notes
    -----
    The implementation uses numpy.linalg.svd for the core computation but
    maintains a pure Python interface for input and output.
    """
    pass

def null_space(a, rcond=None):
    """
    Construct an orthonormal basis for the null space of A using SVD.

    Parameters
    ----------
    a : list[list[float]]
        Input matrix A of shape (M, N)
    rcond : float, optional
        Relative condition number. Singular values ``s`` smaller than
        ``rcond * max(s)`` are considered zero.
        Default: floating point eps * max(M,N).

    Returns
    -------
    Z : list[list[float]]
        Orthonormal basis for the null space of A.
        K = dimension of effective null space, as determined by rcond
    """
    pass
