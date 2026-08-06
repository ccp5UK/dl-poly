## 2024-05-24 - Fortran Exponentiation Performance
**Learning:** Using floating-point exponentiation (like `x**2.0_wp`) instead of integer exponentiation (`x**2`) causes expensive math library function calls like `exp(y * log(x))`.
**Action:** Always prefer integer exponentiation (e.g., `x**2` or `x**3`) over floating-point exponentiation in Fortran files.
