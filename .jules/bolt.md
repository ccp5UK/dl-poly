## 2026-08-06 - Optimize integer exponentiation in Fortran
**Learning:** In Fortran, expressions like `x ** 2.0_wp` evaluate exponentiation using expensive floating-point mathematical library function calls (e.g., `exp(2.0 * log(x))`). Using integer exponents like `x ** 2` evaluates to direct multiplication (e.g., `x * x`), which is significantly faster and avoids potential NaN issues with negative bases.
**Action:** Prefer integer exponentiation over floating-point exponentiation for small, whole-number powers in Fortran code.
