## 2024-05-27 - Fortran Exponentiation Optimization
**Learning:** In Fortran, replacing floating point exponents (e.g., `**2.0_wp`) with integer exponents (e.g., `**2`) avoids computationally expensive `exp(y * log(x))` math library function calls at runtime. This provides a measurable speedup and avoids potential NaN errors for negative bases.
**Action:** Use integer exponents for small, whole-number powers.
