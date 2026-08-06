## 2026-05-13 - Integer Exponentiation is Faster than Floating-Point
**Learning:** In Fortran, raising a real number to a real power (e.g., `x ** 2.0_wp`) forces the compiler to evaluate it using the exponential and logarithm functions (`exp(2.0 * log(x))`). Changing the exponent to an integer (`x ** 2`) allows the compiler to optimize the operation to simple multiplication (`x * x`), which is significantly faster and often more precise.
**Action:** Always look for and convert floating-point exponentiations to integer exponentiations when the exponent is mathematically an integer.
