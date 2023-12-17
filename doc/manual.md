# Table of Contents
1. [Machine Precision](#Machine-Precision)
- [Single-Precision Machine Epsilon Value (smaceps)](#Single-Precision-Machine-Epsilon-Value)
- [Double-Precision Machine Epsilon Value (dmaceps)](#Double-Precision-Machine-Epsilon-Value)
2. [Norm Lengths & Distances](#Norm-Lengths-&-Distances)
- [L1 Norm Length (l1norm)](#L1-Norm-Length)
- [L2 Norm Length (l2norm)](#L2-Norm-Length)
- [L-Infinity Norm Length (linfnorm)](#L-Infinity-Norm-Length)
- [L1 Norm Distance (l1distance)](#L1-Norm-Distance)
- [L2 Norm Length (l2distance)](#L2-Norm-Distance)
- [L-Infinity Norm Length (linfdistance)](#L-Infinity-Norm-Distance)]
3. [Derivative Approximations](#Derivative-Approximations)
- [Forward Difference Quotient (forwarddf)](#Forward-Difference-Quotient)
- [Backward Difference Quotient (backwarddf)](#Backward-Difference-Quotient)
- [Central Difference Quotient (centraldf)](#Central-Difference-Quotient)
4. [Linear Systems of Equations](#Linear-Systems-of-Equations)
- [Gauss-Jordan Elimination (gaussjordan)](#Gauss-Jordan-Elimination)
- [LU Factorization (lufactorize)](#LU-Factorization)
- [Backward Substitution (backsub)](#Backward-Substitution)
- [Forward Substitution (forwardsub)](#Forward-Substitution)
- [Jacobi Iteration (jacobi)](#Jacobi-Iteration)
- [Gauss-Seidel Method (gaussseidel)](#Gauss-Seidel-Method)
5. [Statistics](#Statistics) 
- [Linear Regression (linreg)](#Linear-Regression)
6. [Root Finding]
- [Fixed Point Iteration (fixedpntiter)](#Fixed-Point-Iteration)
- [Bisection Method (bisectroot)](#Bisection-Method)
- [Newton Method (newtonroot)](#Newton-Method)
- [Secant Method (secantroot)](#Secant-Method)
- [Hybrid Bisection Secant Method (bisectsecantroot)](#Hybrid-Bisection-Secant-Method)
7. [Eigenvalues]
- [Power Method (powermethod)](#Power-Method)
- [Inverse Power Method (inversepowermethod)](#Inverse-Power-Method)
- [Shifted Inverse Power Method (shiftedinvpowermethod)](#Shifted-Inverse-Power-Method)
8. [Miscellaneous]
- [Matrix-Vector multiplication (matvec)](#Matrix-Vector-Multiplication)
- [Dot Product (dotproduct)](#Dot-Product)

# Machine Precision

<hr>

## Single-Precision Machine Epsilon Value

**Routine Name:** smaceps

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Calculates the machine-precision value of a single-precision floating-point number.

**Input:** No input.

**Output:** Returns the machine epsilon value of a single-precision float.

**Usage/Example:**
```
	printf("%E", smaceps());
```

Output:
```

```

**Implementation/Code:** The following is the code for smaceps.c:
```
float smaceps()
{
	float eps = 1.0f;
	float prevEps;   // Stores the last eps such that 1.0f + eps != 1.0f
	while (1.0f + eps != 1.0f)
	{
		prevEps = eps;
		eps /= 2.0f;
	}
	return prevEps;
}
```

**Last Modified:** October/2023

<hr>

## Double-Precision Machine Epsilon Value

**Routine Name:** dmaceps

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Calculates the machine-precision value of a double-precision floating point number.

**Input:** No input.

**Output:** Returns the machine epsilon value of a double-precision floating point number.

**Usage/Example:**
```
	printf("%E", dmaceps());
```

Output:
```

```

**Implementation/Code:** The following is the code for smaceps.c:
```
float dmaceps()
{
	double eps = 1.0;
	double prevEps;   // Stores the last eps such that 1.0 + eps != 1.0
	while (1.0 + eps != 1.0)
	{
		prevEps = eps;
		eps /= 2.0;
	}
	return prevEps;
}
```

**Last Modified:** October/2023

# Norm Lengths & Distances

## L1 Norm Length

**Routine Name:** l1norm

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Computes the l1-norm length of a vector.

**Input:** The vector (v) to compute the l1-norm length of, along with its size n.

**Output:** The l1-norm length of v.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
#include <math.h>

double l1norm(double v[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(v[i]);   // sum |v_i|
	}
	return sum;
}
```

**Last Modified:** October/2023

<hr>

## L2 Norm Length

**Routine Name:** l2norm

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Computes the l2-norm (Euclidean) length of a vector.

**Input:** The vector (v) to calculate the length of, along with its size n.

**Output:** The Euclidean length of the vector v.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
#include <math.h>

double l2norm(double v[], int n)
{
	double sum = 0.0;   // norm^2
	for (int i = 0; i < n; i++)
	{
		sum += v[i] * v[i]);
	}
	return sqrt(sum);
}
```

**Last Modified:** October/2023

<hr>

## L-Infinity Norm Length

**Routine Name:** linfnorm

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Computes the infinity-norm length of a vector.

**Input:** The vector (v) to calculate the length of, along with its size n.

**Output:** The infinity-norm length of the vector v.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
#include <math.h>

double linfnorm(double v[], int n)
{
	double maxAbs = 0.0;
	for (int i = 0; i < n; i++)
	{
		double abs = fabs(v[i]);
		if (x > norm)
		{
			maxAbs = abs;
		}
	}
	return norm;
}
```

**Last Modified:** October/2023

<hr>

## L1 Norm Distance

**Routine Name:** l1distance

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Calculates the distance between two vectors in the L1 norm.

**Input:** The vectors v1 and v2 to calculate the distance between, along with their size n.

**Output:** The l1 norm distance between the vectors.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
#include <math.h>

double l1distance(double v1[], double v2[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(v1[i] - v2[i]);
	}
	return sum;
}
```

**Last Modified:** October/2023

<hr>

## L2 Norm Distance

**Routine Name:** l2distance

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Calculates the distance between two vectors in the L2 norm.

**Input:** The vectors v1 and v2 to calculate the distance between, along with their size n.

**Output:** The l2 norm distance between the vectors.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
#include <math.h>

double l2distance(double v1[], double v2[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		double x = v1[i] - v2[i];
		sum += x * x;
	}
	return sqrt(sum);
```

**Last Modified:** October/2023

<hr>

## L-Infinity Norm Distance

**Routine Name:** linfdistance

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC).

**Description/Purpose:** Calculates the distance between two vectors in the L-infinity norm.

**Input:** The vectors v1 and v2 to calculate the distance between, along with their size n.

**Output:** The l-infinity norm distance between the vectors.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
#include <math.h>

double l2distance(double v1[], double v2[], int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		double x = v1[i] - v2[i];
		sum += x * x;
	}
	return sqrt(sum);
```

**Last Modified:** October/2023

# Derivative Approximations

## Forward Difference Quotient
**Routine Name:** forwarddf

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Approximates f'(x) utilizing the forward difference quotient with increment h.

**Input:** The function f, the value of x at which to approximate f'(x), and the increment h.

**Output:** An approximation of f'(x).

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
double forwarddf(double (*f)(double), double x, double h)
{
	return (f(x + h) - f(x)) / h;
}
```

**Last Modified:** October/2023

## Backward Difference Quotient
**Routine Name:** backwarddf

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Approximates f'(x) utilizing the backward difference quotient with increment h.

**Input:** The function f, the value of x at which to approximate f'(x), and the increment h.

**Output:** An approximation of f'(x).

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
double backwardDf(double (*f)(double), double x, double h)
{
	return (f(x) - f(x - h)) / h;
}	
```

**Last Modified:** October/2023

## Central Difference Quotient
**Routine Name:** centraldf

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Approximates f'(x) utilizing the central difference quotient with increment h.

**Input:** The function f, the value of x at which to approximate f'(x), and the increment h.

**Output:** An approximation of f'(x).

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
double centraldf(double (*f)(double), double(x), double h)
{
	return ((f(x + h) - f(x - h)) / (2 * h));
}
```

**Last Modified:** October/2023

# Linear Systems of Equations

## Gauss-Jordan Elimination
**Routine Name:** gaussjordan

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Given a system of equations ax=b, reduces the matrix a into upper-triangular form while modifying the solution vector b (utilizing Gauss-Jordan elimination).

**Input:** A matrix a and a vector b that are related by a system of equations ax = b. The size n of a and b.

**Output:** No return value, although the matrix a will be modified into upper-triangular form (without leading zeroes) and the vector b will be modified such that ax=b still holds.

**Usage/Example:** Can be utilized in conjunction with the `backsubstitution()` routine to solve a system of linear equations. For example:
```

```

**Implementation/Code:** The following is the code for ...
```
void gaussjordan(int n, double a[][n], double b[])
{
    // Make a upper triangular (modifying b accordingly)
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = a[i][k] / a[k][k];
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
}
```

**Last Modified:** October/2023

<hr>

## LU Factorization
**Routine Name:** lufactorize

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** LU factorizes a matrix in place. Upper-right and diagonal elements are those of the U matrix, while the lower-left elements are those of the L matrix (which has diagonal values of 1).

**Input:** The matrix a to factorize, along with its size n.

**Output:** The matrix a will be modified into an in-place lu-factorized form (refer to Description/Purpose).

**Usage/Example:** Can be used in conjunction with forward substitution and backward substitution to solve a system of linear equations.
```

```

**Implementation/Code:** The following is the code for ...
```
void lufactorize(int n, double a[][n])
{
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = a[i][k] / a[k][k];
            for (int j = k + 1; j < n; j++)
            {
                a[i][j] -= factor * a[k][j];
            }
            a[i][k] = factor;
        }
    }
}
```

<hr>

## Backward Substitution
**Routine Name:** backsub

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Utilizes back-substitution to solve for x in ax=b, where a is upper triangular.

**Input:** Upper triangular matrix a, vector b, vector x (empty), and the size n of all matrices and vectors mentioned.

**Output:** Fills the contents of the vector x with the solution to ax = b.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
void backsub(int n, double a[][n], double b[], double x[])
{
    x[n - 1] = b[n - 1] / a[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}
```

**Last Modified:** October/2023

<hr>

## Forward Substitution
**Routine Name:** forwardsub

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Utilizes forward substitution to solve a system of equations of the form ly = b, where l is a lower triangular matrix with diagonal values of 1.

**Input:** Lower triangular matrix l (diagonal values of 1), vector b, vector y (empty), and the size n of these inputs.

**Output:** The input vector y will be filled with the solution to ly = b.

**Usage/Example:**

**Implementation/Code:** The following is the code for forwardsub():
```
void forwardsub(int n, double l[][n], double b[], double y[])
{
    y[0] = b[0];
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum += l[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
}
```

**Last Modified:** December/2023

## Jacobi Iteration

**Routine Name:** 

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** 

**Input:** 

**Output:** 

**Usage/Example:**

**Implementation/Code:** The following is the code for forwardsub():
```

```

**Last Modified:** December/2023

## Gauss-Seidel Method

**Routine Name:** 

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** 

**Input:** 

**Output:** 

**Usage/Example:**

**Implementation/Code:** The following is the code for forwardsub():
```

```

**Last Modified:** December/2023

# Statistics

<hr>

## Linear Regression
**Routine Name:** linreg

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Calculates the coefficients for the linear regression of a vector of y values against a vector of x values.

**Input:** The vector x of x values, the vector y of y values, the vector A to store the results of the linear regression, and the size n of the x and y vectors.

**Output:** The vector A will be filled with the coefficients of the linear regression y = ax + b, where A[0] = a and A[1] = b.

**Usage/Example:**

**Implementation/Code:** The following is the code for ...
```
void linreg(double x[], double y[], int n, double A[]) {
	double s1 = 0;   // sum x_i
	double s2 = 0;   // sum x_i^2
	
	for (int i = 0; i < n; i++)
	{
		s1 += x[i];
		s2 += x[i] * x[i];
	}
	
	// y = mx + b
	A[0] = 0;   // m
	A[1] = 1;   // b
	
	// (X^T X)^-1 X^T Y
	double det = n * s2 - s1 * s1;
	for (int i = 0; i < n; i++)
	{
		A[0] += y[i] * (n * x[i] - s1);
		A[1] += y[i] * (s2 - s1 * x[i]);
	}
	A[0] /= det;
	A[1] /= det;
}
```

# Root Finding

## Fixed Point Iteration

**Routine Name:** fixedpntiter

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Attempts to find a root of a function f(x) utilizing fixed point iteration.

**Input:** The function must be modified into the form g(x) = x + cf(x) for some nonzero constant c. The function g(x) is used as input, along with the initial guess x for the root and the tolerance of the error (setting this too low may result in not finding the root).

**Output:** If an approximation is found to the root within the given tolerance, the x-value of this root is given. Otherwise, the routine returns NaN and prints an error message to the console.

**Usage/Example:** To find a root of $f(x) = x^2 - 5x + 6$, modify it into the form $g(x) = x + f(x) = x^2 - 4x + 6$. This function is used as input, along with an initial guess for the root of x = 1. Finally, the tolerance is set to 1/1000, a reasonable value that should not be too low.
```
#include "fixedpntiter.c"

double g(double x)
{
	return x * x - 4 * x + 6;
}

int main()
{
    double x = 1.0;
    double tol = 1e-3;
    double root = fixedpntiter(g, x, tol);

    printf("Root at x = %f", root);

    return 0;
}
``` 

Output:
```
Root at x = 3.000000
```

The program found a root at x = 3, which we know is a correct root. For an example of when the root is not found by changing the initial guess to x = 5:
```
#include "fixedpntiter.c"

double g(double x)
{
    return x * x - 4 * x + 6;
}

int main()
{
    double x = 5.0;
    double tol = 1e-9;
    fixedpntiter(g, x, tol);

    return 0;
}
```

Output:
```
Failed to find root
```

**Implementation/Code:** The following is the code for fixedpntiter():
```
/*
    Attempts to approximate the root of a function f(x) utilizing fixed point iteration
    
    Input:
        g: function of the form g(x) = x + cf(x) for some nonzero c
        x: initial guess for root
        tol: tolerance for approximation of root
        
    Output:
        return: approximation of root if found; otherwise, NaN and an error message to the console
*/
#include <stdio.h>
#include <math.h>

double fixedpntiter(double (*g)(double), double x, double tol)
{
    const int maxiter = 15;
    double x1 = x;
    for (int i = 0; i < maxiter; i++)
    {
        double x2 = g(x1);
        if (fabs(x2 - x1) <= tol)
        {
            return x2;
        }
        x1 = x2;
    }
    printf("Failed to find root");
    return NAN;
}
```

**Last Modified:** December/2023

## Bisection Method

**Routine Name:** bisectroot

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Utilizes the bisection method to solve for a root of a function f(x) within a given interval $[a, b]$.

**Input:** The function f(x) to find the root of, the values a and b in the interval $[a, b]$ on which the root is being found ($f(a)$ and $f(b)$ must be nonzero and of opposite sign), and the tolerance for which to approximate the root.

**Output:** Returns the approximated root if the values a and b are valid. Otherwise, the root does not exist, and the routine will return NaN.

**Usage/Example:** To find a root of $f(x) = x^2 - 2$, you can utilize the interval $[1, 2]$, as $f(1) = -1$ and $f(2) = 2$ are of opposite signs. Utilize $1/1000$ as a reasonable tolerance. 

```
#include "bisectroot.c"

double f(double x)
{
    return x * x - 2.0;
}

int main()
{
    double a = 1.0;
    double b = 2.0;
    double tol = 1e-3;
    
    double root = bisectroot(f, a, b, tol);
    printf("Root found at x = %f", root);
    
    return 0;
}
```

Output:
```
Root found at x = 1.414551
```

If we utilize the interval $[-2, 2]$, although this interval contains two roots of $f(x)$, it is invalid for the bisection method, as $f(-2)$ and $f(2)$ are of the same sign. Therefore, it should return NaN as the root.

```
#include "bisectroot.c"

double f(double x)
{
    return x * x - 2.0;
}

int main()
{
    double a = -2.0;
    double b = 2.0;
    double tol = 1e-3;
    
    double root = bisectroot(f, a, b, tol);
    printf("Root found at x = %f", root);
    
    return 0;
}
```

Output:
```
Root found at x = nan
```

**Implementation/Code:** The following is the code for bisectroot():
```
/*
	Utilizes the bisection method to solve for a root of a function n a given interval
	
	Input:
		f: function to find root of
		a: left-point of interval to find root in
		b: right-point of interval to find root in
		tol: determines the precision of the root
		
	Output:
		Approximate root or NaN if no root was found within the given tolerance
*/
double bisectroot(double (*f)(double), double a, double b, double tol)
{
    double x1 = a;
    double x2 = b;
    
    // Check that both values are at opposite signs
    double f1 = f(x1);
    double f2 = f(x2);
    if (f1 * f2 >= 0)
    {
        return NAN; // No root guaranteed within the interval
    }
    
    // Iterate
    double xmid = (x1 + x2) / 2.0; // Midpoint
    double fmid = f(xmid);
    while (x2 - x1 > tol)
    {
        if (fmid == 0)
        {
            return fmid;
        }
        else if (f1 * fmid < 0)
        {
            x2 = xmid;
            f2 = fmid;
        }
        else
        {
            x1 = xmid;
            f1 = fmid;
        }
        xmid = (x1 + x2) / 2.0; // Midpoint
        fmid = f(xmid);
    }
    return xmid;
}
```

**Last Modified:** December/2023

## Newton Method

**Routine Name:** newtonroot

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Approximates the root of a function $f(x)$ utilizing the Newton method.

**Input:** The function $f(x)$ and its exact derivative $f'(x)$ (in code as df), along with an initial guess x and the tolerance for which to approximate the root.

**Output:** If a root is found within the given tolerance, the root is returned. If the root is not found, the return value is NaN and an error message is printed to the console.

**Usage/Example:** To find a root of $f(x) = x^3 - 2x$ by utilizing its exact derivative $f'(x) = 3x^2 - 2$ with an initial guess of $x = 1$ and reasonable tolerance of $1/1000$.

```
#include "newtonroot.c"

double f(double x)
{
    return x * x * x - 2 * x;   
}

double df(double x)
{
    return 3 * x * 2 - 2;
}

int main()
{
    double x = 1.0;
    double tol = 1e-3;
    
    double root = newtonroot(f, df, x, tol);
    printf("Root at x = %f", root);
    
    return 0;
}
```

Output:
```
Root at x = 1.413666
```

If we utilize a large enough initial guess $x = 1000$, the method will not converge within the allotted number of iterations and shoul print an error message to the console.

```
#include "newtonroot.c"

double f(double x)
{
    return x * x * x - 2 * x;   
}

double df(double x)
{
    return 3 * x * 2 - 2;
}

int main()
{
    double x = 1000.0;
    double tol = 1e-3;
    
    newtonroot(f, df, x, tol);
    
    return 0;
}
```

Output:
```
Root not found
```

**Implementation/Code:** The following is the code for newtonroot():
```
#include <stdio.h>
#include <math.h>

/*
	Approximate the root of a function using Newton's method
	
	Input:
		f: function to approximate the root of
		df: derivative of f
		x: starting point for x
		tolerance: determines precision of calculated root
		
	Output:
		Approximated root or NaN if a root within the given tolerance was not found
*/
double newtonroot(double (*f)(double), double (*df)(double), double x, double tolerance)
{
    const int maxiter = 15;
    double x1 = x;
    for (int i = 0; i < maxiter; i++)
    {
        double x2 = x1 - f(x1) / df(x1);
        if (fabs(x2 - x1) <= tolerance)
        {
            return x2;
        }
        x1 = x2;
    }
    // If root is not found
    printf("Root not found");
    return NAN;
}

```

**Last Modified:** December/2023

## Secant Method

**Routine Name:** secantroot

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Approximates a root of a function $f(x)$ utilizing the secant method.

**Input:** The function $f(x)$ to calculate the root of, the interval $[a, b]$ that defines the starting bounds for the method, and the tolerance for which to approximate the root.

**Output:** If the root value is found, it is returned. Otherwise, the return value is NaN and an error message is printed to the console.

**Usage/Example:** To approximate the root of $f(x) = x^2 - 2$ starting on the interval $[1, 2], using a tolerance of $1/1000$:

```
#include "secantroot.c"

double f(double x)
{
    return x * x - 2.0;
}

int main()
{
    double a = 1.0;
    double b = 2.0;
    double tol = 1e-3;
    
    double root = secantroot(f, a, b, tol);
    printf("Root at x = %f", root);
    
    return 0;
}
```

Output:
```
Root at x = 1.414211
```

If we utilize the interval $[1000, 2000]$ in the following modified code, the root will not be found. The error message will print, and the value of the root will be NaN:
```
double f(double x)
{
    return x * x - 2.0;
}

int main()
{
    double a = 1000.0;
    double b = 2000.0;
    double tol = 1e-3;
    
    double root = secantroot(f, a, b, tol);
    printf("\nRoot at x = %f", root);
    
    return 0;
}
```

Output:
```
Root not found
Root at x = nan
```


**Implementation/Code:** The following is the code for secantroot():
```
#include <stdio.h>
#include <math.h>

/*
	Approximates the root of a function within a starting interval utilizing the secant method
	
	Input:
		f: function to calculate the root of
		a: left starting bound on x
		b: right starting bound on x
		tolerance: determines precision of root calculated
	
	Output:
		Root value if calculated, and NaN otherwise if a root with the given tolerance is not found
*/
double secantroot(double (*f)(double), double a, double b, double tolerance)
{
    const int maxiter = 15;
    double x0 = a;
    double x1 = b;
    double f1 = f(x1);
    for (int i = 0; i < maxiter; i++)
    {
        double x2 = x1 - f1 * (x1 - x0) / (f1 - f(x0));
        x0 = x1;
        x1 = x2;
        f1 = f(x1);
        if (fabs(x1 - x0) < tolerance)
        {
            return x1;
        }
    }
    printf("Root not found");
    return NAN;
}
```

**Last Modified:** December/2023

## Hybrid Bisection Secant Method

**Routine Name:** bisectsecantroot

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Utilizes a combination of the bisection and secant methods to increase the efficiency at which the root is found. The secant method is utilized until the bisection method becomes viable, at which point it becomes the main method.

**Input:** The function $f(x)$ to find the root of, the interval $[a, b]$ on which to begin the secant method, and the tolerance for which the root is approximated.

**Output:** Returns the root if found.

**Usage/Example:** To find the root of $f(x) = x^2 - 2$, utilize the interval $[1, 2]$ and a tolerance of $1/1000$:

```
#include "bisectsecantroot.c"

#include <stdio.h>

double f(double x)
{
    return x * x - 2.0;
}

int main()
{
    double a = 1.0;
    double b = 2.0;
    double tol = 1e-3;
    
    double root = bisectsecantroot(f, a, b, tol);
    printf("Root at x = %f", root);
    
    return 0;
}
```

Output:
```
Root at x = 1.414214
```

**Implementation/Code:** The following is the code for bisectsecantroot():
```
#include <math.h>

/*
	Hybrid root-finding method that combines the guaranteed success of bisection with the rapid convergence
	of the secant method
	
	Input:
		f: function to find the root of
		a: left bound of root
		b: right bound of root
		tolerance: precision of root
*/
double bisectsecantroot(double (*f)(double), double a, double b, double tolerance)
{
    double x0 = a;
    double x1 = b;
    while (1)
    {
        // Secant method
        double x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        while (x0 < x2 && x2 < x1) // Continue while x2 in [x0, x1]
        {
            // Pick interval with opposite signs on ends
            // Maintain x0 < x2 < x1
            if (f(x2) == 0.0 || fabs(x2 - x1) <= tolerance)
            {
                return x2;
            }
            else if (f(x0) * f(x2) < 0.0)
            {
                x1 = x2;
            }
            else
            {
                x0 = x2;
            }
            
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        }
        
        // Single iteration of bisection
        x2 = (x0 + x1) / 2.0;
        if (f(x2) == 0.0 || fabs(x2 - x1) <= tolerance)
        {
            return x2;
        }
        else if (f(x0) * f(x2) < 0)
        {
            x1 = x2;
        }
        else
        {
            x0 = x2;
        }
    }
}
```

**Last Modified:** December/2023

# Eigenvalues

## Power Method

**Routine Name:** powermethod

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Calculates the dominant eigenvalue of a given matrix.

**Input:** The matrix A to calculate the dominant eigenvalue of, the initial guess x0 for the dominant eigenvector, the maximum number of iterations to apply the power method, and the tolerance to which to approximate the dominant eigenvalue.

**Output:** The dominant eigenvalue is returned. If it is not found, an error message is printed and the final approximation for the dominant eigenvalue is returned.

**Usage/Example:**

```
#include "powermethod.c"

#include <stdio.h>

int main()
{
    double A[3][3] = {
        {-11, 0, 4},
        {14, -5, 4},
        {-98, 6, 25}
    };
    double x0[3] = {1.0, 1.0, 1.0};
    double tol = 1e-6;
    double maxiter = 1e3;
    double eigenvalue = powermethod(3, A, x0, tol, maxiter);
    printf("Dominant eigenvalue: %f", eigenvalue);
    
    return 0;
}
```

Output:
```
Dominant eigenvalue: 5.000006
```

**Implementation/Code:** The following is the code for powermethod():
```
#include "l2norm.c"
#include "dotproduct.c"
#include "matvec.c"

#include <stdio.h>

/*
    Copies elements from vector o into vector c
*/
double copyVector(int n, double o[n], double c[n])
{
    for (int i = 0; i < n; i++)
    {
        c[i] = o[i];
    }
}

/*
    Normalizes vector
*/
void normalize(int n, double v[n])
{
    double magnitude = l2norm(n, v);
    for (int i = 0; i < n; i++)
    {
        v[i] /= magnitude;
    }
}

/*
	Computes the dominant eigenvalue of a matrix utilizing the power Method
	
	Input:
		n: size of matrices and vectors
		A: square matrix to calculate dominant eigenvalue of
		x0: initial guess of dominant eigenvector
		tol: tolerance to approximate dominant eigenvalue to
		maxiter: maximum number of power method iterations
*/
double powermethod(int n, double A[n][n], double x0[n], double tol, int maxiter)
{
    double x1[n]; // x_i
    double x2[n]; // x_(i+1)
    copyVector(n, x0, x1); // Preserve original guess array
    normalize(n, x1); // Ensure initial guess is normalized
    matvec(n, A, x1, x2); // x_(i+1) = Ax_i; do not normalize x2 yet
    
    double lam = dotproduct(n, x1, x2); // lambda_k; Since x_(k+1) = Ax_k and ||x_k|| = 1 => Rayleigh quotient: lambda = x1 * x2;
    double res[n]; // Residual vector: res_k = Ax_k - l_k x
    for (int i = 0; i < n; i++)
    {
        res[i] = x2[i] - lam * x1[i];
    }
    double error = l2norm(n, res);
    
    int iter = 0;
    while (error > tol && iter < maxiter)
    {
        // Compute x_(k+1) = norm(Ax_k)
        copyVector(n, x2, x1);
        normalize(n, x1);
        matvec(n, A, x1, x2); // x_(i+1) = normalize(Ax_i)
        
        // Calculate residual and error
        lam = dotproduct(n, x1, x2);
        for (int i = 0; i < n; i++)
        {
            res[i] = x2[i] - lam * x1[i];
        }
        error = l2norm(n, res);
    }
    if (error > tol)
    {
        printf("Dominant eigenvalue not found with given tolerance and iterations");
    }
    return lam;
}
```

**Last Modified:** December/2023

## Inverse Power Method

**Routine Name:** inversepowermethod

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Utilizes the inverse power method to calculate the smallest eigenvalue (by magnitude) of a given matrix.

**Input:** The matrix A to find the smallest eigenvalue of, the initial guess x0 at the eigenvalue's corresponding eigenvector, and the meximum number of iterations and tolerance for which to approximate the eigenvalue.

**Output:** The smallest eigenvalue is returned. If it is not found, an error message is printed and the final approximation for the eigenvalue is returned.

**Usage/Example:**

```
int main()
{
    double A[3][3] = {
        {-11, 0, 4},
        {14, -5, 4},
        {-98, 6, 25}
    };
    double x0[3] = {1.0, 1.0, 1.0};
    double tol = 1e-6;
    double max_iter = 1e3;
    double eigenvalue = inversepowermethod(3, A, x0, tol, max_iter);
    printf("Smallest eigenvalue: %f", eigenvalue);
    
    return 0;
}
```

Output:
```
Smallest eigenvalue: 0.999997
```

**Implementation/Code:** The following is the code for inversepowermethod():
```
#include "l2norm.c"
#include "dotproduct.c"
#include "lufactorize.c"
#include "forwardsub.c"
#include "backsub.c"

#include <stdio.h>

/*
    Copies elements from vector o into vector c
*/
double copyVector(int n, double o[n], double c[n])
{
    for (int i = 0; i < n; i++)
    {
        c[i] = o[i];
    }
}

/*
    Normalizes vector
*/
void normalize(int n, double v[n])
{
    double magnitude = l2norm(n, v);
    for (int i = 0; i < n; i++)
    {
        v[i] /= magnitude;
    }
}

/*
    Copies the contents of a matrix
*/
void copyMatrix(int n, double o[n][n], double c[n][n])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            c[i][j] = o[i][j];
        }
    }
}

/*
    Solves a system of equations LUx = b
*/
void solveLU(int n, double lu[n][n], double b[n], double x[n])
{
    double y[n];
    forwardsub(n, lu, b, y); // Solve Ly = b for y
    backsub(n, lu, y, x); // Solve Ux = b for x
}

/*
	Computes the smallest (by magnitude) eigenvalue of a matrix utilizing the inverse power method
	
	Input:
		n: size of matrices and vectors
		A: square matrix to calculate dominant eigenvalue of
		x0: initial guess of dominant eigenvector
		tol: tolerance to approximate dominant eigenvalue to
		maxiter: maximum number of inverse power method iterations
*/
double inversepowermethod(int n, double A[n][n], double x0[n], double tol, int maxiter)
{
    // LU decompose A
    double lu[n][n];
    copyMatrix(n, A, lu);
    lufactorize(n, lu);
    
    // Compute lambda^-1
    double x1[n]; // x_i
    double x2[n]; // x_(i+1)
    copyVector(n, x0, x1); // Preserve original guess array
    normalize(n, x1); // Ensure initial guess is normalized
    solveLU(n, lu, x1, x2); // Ax_(i+1) = x_i; do not normalize x2 yet
    
    double invlam = dotproduct(n, x1, x2); // lambda_k^-1
    double res[n]; // Residual vector: res_k = Ax_k - l_k x
    for (int i = 0; i < n; i++)
    {
        res[i] = x2[i] - invlam * x1[i];
    }
    double error = l2norm(n, res);
    
    int iter = 0;
    while (error > tol && iter < maxiter)
    {
        // Compute x_(i+1) = norm(A^-1 x_i)
        copyVector(n, x2, x1);
        normalize(n, x1);
        solveLU(n, lu, x1, x2);
        
        // Calculate residual and error
        invlam = dotproduct(n, x1, x2);
        for (int i = 0; i < n; i++)
        {
            res[i] = x2[i] - invlam * x1[i];
        }
        error = l2norm(n, res);
    }
    if (error > tol)
    {
        printf("Dominant eigenvalue not found with given tolerance and iterations");
    }
    return 1 / invlam; // Invert lambda
}
```

**Last Modified:** December/2023

## Shifted Inverse Power Method

**Routine Name:** shiftedinvpowermethod

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Utilizes the shifted inverse power method to calculate the eigenvalue closest to a shift value of a matrix.

**Input:** The matrix A to find an eigenvalue of, the initial guess x0 at the eigenvalue's corresponding eigenvector, and the meximum number of iterations and tolerance for which to approximate the eigenvalue.

**Output:** The eigenvalue closest to the shift value if the method converges. Otherwise, an error message is printed to the console.

**Usage/Example:** Take a matrix with eigenvalues 1, 3, and 7. To find the middle eigenvalue, use a shift of 4.

```
#include "shiftedinvpowermethod.c"

#include <stdio.h>

int main()
{
    double A[3][3] = {
        {13, -30, 22},
        {-24, 67, -48},
        {-36, 96, -69}
    };
    
    double x0[3] = {1, 0, 0};
    double tol = 1e-6;
    double max_iter = 1e3;
    
    double shift = 4.0;
    double eigenvalue = shiftedinvpowermethod(3, shift, A, x0, tol, max_iter);
    printf("Middle eigenvalue: %.2f", eigenvalue);
    
    return 0;
}
```

**Implementation/Code:** The following is the code for shiftedinvpowermethod():
```
#include "l2norm.c"
#include "dotproduct.c"
#include "lufactorize.c"
#include "forwardsub.c"
#include "backsub.c"

#include <stdio.h>

/*
    Copies elements from vector o into vector c
*/
double copyVector(int n, double o[n], double c[n])
{
    for (int i = 0; i < n; i++)
    {
        c[i] = o[i];
    }
}

/*
    Normalizes vector
*/
void normalize(int n, double v[n])
{
    double magnitude = l2norm(n, v);
    for (int i = 0; i < n; i++)
    {
        v[i] /= magnitude;
    }
}

/*
    Copies the contents of a matrix
*/
void copyMatrix(int n, double o[n][n], double c[n][n])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            c[i][j] = o[i][j];
        }
    }
}

/*
    Solves a system of equations LUx = b
*/
void solveLU(int n, double lu[n][n], double b[n], double x[n])
{
    double y[n];
    forwardsub(n, lu, b, y); // Solve Ly = b for y
    backsub(n, lu, y, x); // Solve Ux = b for x
}

/*
	Computes the eigenvalue of a matrix closest to a shift value utilizing the shifted inverse power method
	
	Input:
		n: size of matrices and vectors
		shift: value to shift matrix A by
		A: square matrix to calculate dominant eigenvalue of
		x0: initial guess of dominant eigenvector
		tol: tolerance to approximate dominant eigenvalue to
		maxiter: maximum number of shifted inverse power method iterations
*/
double shiftedinvpowermethod(int n, double shift, double A[n][n], double x0[n], double tol, int maxiter)
{
    // LU decompose A
    double lu[n][n];
    copyMatrix(n, A, lu);
    for (int i = 0; i < n; i++) // Shift matrix before factoring
    {
        lu[i][i] -= shift;
    }
    lufactorize(n, lu);
    
    // Compute lambda^-1
    double x1[n]; // x_i
    double x2[n]; // x_(i+1)
    copyVector(n, x0, x1); // Preserve original guess array
    normalize(n, x1); // Ensure initial guess is normalized
    solveLU(n, lu, x1, x2); // Ax_(i+1) = x_i; do not normalize x2 yet
    
    double invlam = dotproduct(n, x1, x2); // lambda_k^-1
    double res[n]; // Residual vector: res_k = Ax_k - l_k x
    for (int i = 0; i < n; i++)
    {
        res[i] = x2[i] - invlam * x1[i];
    }
    double error = l2norm(n, res);
    
    int iter = 0;
    while (error > tol && iter < maxiter)
    {
        // Compute x_(i+1) = norm(A^-1 x_i)
        copyVector(n, x2, x1);
        normalize(n, x1);
        solveLU(n, lu, x1, x2);
        
        // Calculate residual and error
        invlam = dotproduct(n, x1, x2);
        for (int i = 0; i < n; i++)
        {
            res[i] = x2[i] - invlam * x1[i];
        }
        error = l2norm(n, res);
    }
    if (error > tol)
    {
        printf("Dominant eigenvalue not found with given tolerance and iterations");
    }
    return 1 / invlam + shift; // Invert lambda, then remove shift
}
```

**Last Modified:** December/2023

# Miscellaneous

## Matrix-Vector Multiplication

**Routine Name:** matvec

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Multiples a matrix and a vector.

**Input:** A square matrix and vector of the same size n, along with an array in which to store the resulting vector.

**Output:** The result of the multiplication is stored in the p vector.

**Usage/Example:** We can multiply a matrix and a vector as follows:

```
#include "matvec.c"

#include <stdio.h>

int main()
{
    int n = 20;
    
    // Fill matrix and vector with 1s
    double m[n][n];
    double v[n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            m[i][j] = 1.0;
        }
        v[i] = 1.0;
    }
    
    // Create vector p to store result of m * v and print contents of p to console
    double p[n];
    matvec(n, m, v, p);
    for (int i = 0; i < n; i++)
    {
        printf("%.2f ", p[i]);
    }
    
    return 0;
}
```

Output:
```
20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00 20.00
```

The product vector is full of 20s, as expected.

**Implementation/Code:** The following is the code for matvec():
```
/*
    Multiply a matrix and a vector
    
    Input:
        n: size of matrix and vectors
        m: square matrix
        v: vector
        
    Output:
        p: stores the result m * v
*/
void matvec(int n, double m[n][n], double v[n], double p[n])
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += m[i][j] * v[j];
        }
        p[i] = sum;
    }
}
```

**Last Modified:** December/2023

## Dot Product

**Routine Name:** dotproduct

**Author:** Bryan Armenta

**Language:** C. The code can be compiled using the GNU C compiler (GCC). 

**Description/Purpose:** Calculates the dot product of two vectors.

**Input:** The vectors v1 and v2 to find the dot product of, along with their size n.

**Output:** The dot product v1 * v2.

**Usage/Example:** The dot product of two vectors of size 100 filled entirely with the value 1 should be 100.

```
#include "dotproduct.c"

#include <stdio.h>

int main()
{
    // Create two vectors of size 100 filled with 1s
    int n = 100;
    double v1[n], v2[n];
    for (int i = 0; i < n; i++)
    {
        v1[i] = 1.0;
        v2[i] = 1.0;
    }
    
    // Their dot product should be 100
    double dot = dotproduct(n, v1, v2);
    printf("Dot product: %f", dot);
    
    return 0;
}
```

Output:
```
Dot product: 100.000000
```

**Implementation/Code:** The following is the code for dotproduct():
```
/*
    Calculate the dot product of two vectors
    
    Input:
        n: size of vectors
        v1: vector to calculate dot product
        v2: vector to calculate dot product
        
    Output:
        return: dot product of vectors
*/
double dotproduct(int n, double v1[n], double v2[n])
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += v1[i] * v2[i];
    }
    return sum;
}
```

**Last Modified:** December/2023