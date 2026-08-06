## 2026-08-04 - Optimize Fortran Exponentiation
**Learning:** In Fortran, floating-point exponentiation (e.g., `x**2.0_wp`) is evaluated using expensive math library functions (like `exp(y * log(x))`). Integer exponentiation (`x**2`) uses simple multiplication and is much faster. Also, floating-point exponents can cause NaN errors with negative bases.
**Action:** Always prefer integer exponents for integer powers in Fortran code.
