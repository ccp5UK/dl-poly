## 2024-05-14 - Fortran Exponentiation Optimization
**Learning:** In Fortran, floating-point exponentiation (e.g., `x**2.0_wp` or `x**3.0_wp`) is significantly slower than integer exponentiation (e.g., `x**2` or `x**3`) because it uses an expensive math library routine `exp(y * log(x))` instead of a simple multiplication loop or inline multiply instructions.
**Action:** Always prefer integer types for exponents in Fortran expressions (e.g., `**2`, `**3`) when the exponent is an integer, avoiding `.0_wp` or `.0` suffixes for powers.
