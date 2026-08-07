## 2026-08-07 - Optimize integer powers in Fortran
**Learning:** Using floating-point exponents (e.g., `**2.0_wp`) for integer powers triggers expensive math library function calls like `exp(y * log(x))` and can cause NaN errors for negative bases.
**Action:** Always prefer integer exponentiation (e.g., `**2`, `**3`) over floating-point exponentiation in Fortran.
