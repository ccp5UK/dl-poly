## $(date +%Y-%m-%d) - Use integer exponentiation in Fortran
**Learning:** In Fortran, raising a floating-point number to a floating-point power (like `x**2.0_wp`) can be significantly slower than raising it to an integer power (`x**2`), as the compiler may use an expensive power function instead of simple multiplication.
**Action:** Always prefer integer exponentiation when the exponent is an integer.
