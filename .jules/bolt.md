## 2026-08-03 - Optimize Fortran exponentiation
**Learning:** In Fortran, raising to a floating-point power (e.g., `x**2.0_wp`) translates to expensive math library function calls like `exp(y * log(x))` and can cause NaN errors with negative bases.
**Action:** Always prefer integer exponentiation (e.g., `x**2`) over floating point exponentiation.
