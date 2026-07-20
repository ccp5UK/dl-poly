## 2026-07-20 - Fortran Integer Exponentiation Optimization
**Learning:** In this Fortran codebase, using integer exponentiation (e.g., `x**2`) is faster and more robust than floating-point exponentiation (e.g., `x**2.0_wp`), which triggers expensive math library function calls like `exp(y * log(x))` and can cause NaN errors with negative bases.
**Action:** Prefer integer literals for exponents when the power is a whole number to improve computational efficiency.
