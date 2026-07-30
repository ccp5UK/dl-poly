## 2026-07-30 - [Performance: Integer vs Floating-Point Exponentiation]
**Learning:** In this Fortran codebase, evaluating exponents as floating point values (e.g., `x**2.0_wp`) translates to expensive math library calls like `exp(y * log(x))` and can cause NaN errors for negative bases, whereas integer exponentiation (e.g., `x**2`) is much faster and avoids these robust issues.
**Action:** Always prefer integer exponents (e.g., `**2`) over floating-point exponents (e.g., `**2.0_wp`) for small integer powers.
