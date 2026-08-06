## 2024-05-24 - Bolt Initial Journal
**Learning:** For Fortran performance optimization in this codebase, prefer integer exponentiation (e.g., `x**2` or `x**3`) over floating-point exponentiation (e.g., `x**2.0_wp`), as it avoids expensive math library function calls like `exp(y * log(x))`.
**Action:** Replace `**2.0_wp`, `**3.0_wp`, etc. with `**2`, `**3`, etc.
