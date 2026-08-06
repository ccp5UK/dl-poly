## 2024-06-11 - Fortran Floating-Point Exponentiation
**Learning:** Fortran codes often use floating-point exponentiation (e.g., `x**2.0_wp` instead of `x**2`) which calls expensive math library functions like `exp()` and `log()` under the hood instead of simple multiplication. This is a common performance anti-pattern in scientific computing.
**Action:** Always search for floating-point exponents like `**2.0` or `**3.0` and replace them with integer exponents (`**2`, `**3`) to avoid unnecessary generic math library calls.
