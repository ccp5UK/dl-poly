## 2024-05-13 - Fortran Exponentiation Optimization
**Learning:** Fortran's `**` operator behaves differently depending on whether the exponent is an integer or a real number. `x**2` (integer exponent) translates to `x * x`, while `x**2.0` (real exponent) translates to `exp(2.0 * log(x))`, which is significantly slower due to math library function calls. The same applies for `.0_wp` kinds.
**Action:** Replace `**2.0_wp` with `**2` and `**3.0_wp` with `**3` in exponentiations.
