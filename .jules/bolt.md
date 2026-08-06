## 2024-10-27 - Fortran Exponentiation Optimization
**Learning:** Floating-point exponentiation (e.g., `**2.0_wp`) in Fortran is a performance anti-pattern because it causes expensive math library calls like `exp(y * log(x))` and can cause NaN errors with negative bases.
**Action:** Always prefer integer exponentiation (e.g., `**2`) when the exponent is a whole number.
