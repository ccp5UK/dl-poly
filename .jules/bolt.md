## 2024-04-23 - Integer Exponentiation
**Learning:** Avoid using floating point exponentiation for small integer powers like `**2.0_wp`. Fortran converts this into expensive math library calls like `exp(y * log(x))`.
**Action:** Replace `**2.0_wp` and `**3.0_wp` with `**2` and `**3`.
