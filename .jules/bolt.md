## 2026-07-26 - Use Integer Exponentiation
**Learning:** In Fortran, integer exponentiation (e.g. `x**2`) should be preferred over floating point exponentiation (e.g. `x**2.0_wp`) to avoid expensive math library calls like `exp(y * log(x))`.
**Action:** Always check for `**2.0` or similar floating point exponents and convert them to integer exponents when the power is a whole number to improve performance.
