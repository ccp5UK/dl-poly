## 2024-05-24 - Avoid Floating-Point Exponentiation
**Learning:** Fortran performance pattern: using floating-point exponentiation (e.g., `x**2.0_wp`) instead of integer exponentiation (`x**2`) invokes expensive math library function calls like `exp(y * log(x))` instead of simple sequential multiplications.
**Action:** Always prefer integer exponentiation (`x**2`, `x**3`) over floating-point exponents when raising variables to integer powers.
