# Coding Style

The programming style of DL\_POLY is intended to ensure the code is a
consistent, readable and maintainable as possible. The following rules apply
throughout the code. Contributors to the code are urged to read and apply this
style before submitting code for inclusion consideration.

# Units

 All routines employ DL_POLY internal units. Output contains, when relevant, the units.

    SIMULATION CONTROL PARAMETERS

    simulation temperature (K)          1.0000E+01

    simulation pressure (katms)         0.0000E+00

    Integration : Leapfrog Verlet
    Ensemble : NVT Nose-Hoover
    thermostat relaxation time (ps)     1.0000E-01

    selected number of timesteps         5000

# General Style

-   Use modern Fortran free form syntax.
-   Code should be written in Fortran2003/2008 dialect if possible.
-   Check the deprecated features from Fortran 2003/2008 and prior.
    Avoid using these features.
-   Fortran2008 should be treated with care, since not all compilers
    implement it.
-   Do not use **Common** blocks and **Block Data**, use **Modules**
    with public data if constant or a user defined type to group data which
    changes during runtime.
-   Do not use **go to** statements
-   Do not use **format** statements
-   File extension shall be .F90 for any new written code.
-   Indent code blocks by two space characters. Do not use
    tabs as they are not part of the Fortran standard.

```fortran
Subroutine init(T,k,traj)
  Class(current_type) :: T
  Type(trajectory_type),Intent(in)  :: traj
  Real(Kind=dp),Intent(in)          :: k(3)

  Allocate(T%jk(traj%nFrames,3))
  T%k = k
  T%n = traj%nFrames
End Subroutine init
```

-   Indent **Contains** statements to the same level as the sorrounding
    **Module**, or procedure

```fortran
Module my_module
  Implicit None

  Integer( Kind = wi ), Parameter :: pi = 3.14...

  ...

Contains

  Subroutine my_subroutine(a,b)

    ...

End Module my_module
```

-   Do not use more than one blank line to separate blocks of code.
-   Code lines shall not exceed 132 characters (the modern Fortran standard
    limit).
-   Do not write Fortran keywords in ALL capitals. Capitalize the first
    letter of the statement to make them stand out for the reader.

```fortran
    Program dl_poly
      Use kinds, Only : wi,wp
      Implicit None
    End Program dl_poly
```
-   While Fortran supports multiple statements on a line separated by *;*
    try to use them sparingly and wisely.
-   Variables, constants and program units should be named in English or using
    common physics notation. Do not use cryptic names.
-   Try to group variables names based on physical meaning or usage.

```fortran
  ! Verlet neighbour list data
  !> Update cells flag
  Logical, Public :: update = .true.
  !> Unconditional update flag
  Logical, Public :: unconditional_update = .false.


  !> Tracking points for Verlet neighbour list
  Real( Kind = wp ), Allocatable, Public :: xbg(:),ybg(:),zbg(:)

  !> Largest vdw cutoff, defines Verlet neighbour list radius
  Real( Kind = wp ), Public :: cutoff
  !> Padding around cutoff
  Real( Kind = wp ), Public :: padding
  !> Actual Verlet neighbour list cutoff (cutoff+padding)
  Real( Kind = wp ), Public :: cutoff_extended

  !> Linked cell list
  Integer( Kind = wi ), Allocatable, Public :: list(:,:)
  !> Maximum rank of linked cell list
  Integer( Kind = wi ), Public :: max_list
  !> Maximum number of cells per domain
  Integer( Kind = wi ), Public :: max_cell
```

-   Avoid naming program units, variables, constants using
    intrinsincs routines or keywords from fortran standard.
-   While Fortran is case insensitive, use for lower case letters for program
    units, variables and constants.
-   Where more than a word is used in a name use _ as separator (snake
    case notation).
-   Use the mordern syntax for logical expressions,
    -   == not .eq.
    -   /= not .ne.
    -   &gt; not .gt.
    -   &lt; not .lt.
    -   &gt;= not .ge.
    -   &lt;= not .le. .
-   Prefer positive testing in logical blocks.

```fortran
  If (isw == 0) Then
     ! do some fancy code
  Else
    ! do some not so fancy code
  End If
! to
  If (isw /= 0) Then
    ! do some not so fancy code
  Else
    ! do some fancy code
  End If
```
-   One line **If** shall be avoided.
-   Always use the optional separation space for **End** constructs, *e.g.*
    **End If** not **Endif** and **End Do** not **Enddo**
-   For program units use the full **End**: **End Subroutine name**.
-   Never use **Go To** statements
-   Use **Do** ... **End Do** instead of **Do** ... Label **Continue**
-   Do not use **Format**.

```fortran

  Write(*,20) reclen

20 Format('(',i0,'a)')

! replace with

  Write(*,'(a,i0,a)')"(",reclen,"a)"
```

-   Do not use interactive commands, like **Read** from keyboard or
    **Pause**.
-   Use **Implicit None** to avoid unwanted side effects. This is only nessecary once per module.
-   Avoid implicit casting of data: use 2.0_wp instead of 2.0 in an
    expression.

```fortran
 If (sw == 1) Then
   factor = 2.0_wp*engke*rf
 End If
! and not
If (sw == 1) factor = 2.0_wp*engke*rf
```

-   Floating point comparisons shall be done with care. Consider the following
    examples using a tolerance of **Epsilon(a)** or **Tiny(a)**,
    -   a > tolerance instead of a > 0.0
    -   a - b > tolerance instead of a > b
    -   abs(a - b) < tolerance instead of a == b .
-   Avoid use of magic numbers, declare constants at the top level of a module and use those instead.

```fortran
If (integrator == 1 ) Then
! instead use

Integer( Kind = wi ), Parameter :: VELOCITY_VERLET = 1

...

If (integrator == VELOCITY_VERLET ) Then
```

-   Any new feature shall be turnable on/off from CONTROL file.
-   Any new feature shall be documented in the manual and will cite relevant
    literature.

## Modular structure

-   All subroutines/functions shall be enclosed by a module.
-   Modules may contain the following:
  -  Type declarations
  -  Subroutines and functions
  -  Paramter definitions (using the **Parameter** attribute)
  -  Interfaces for overloaded procedures
  -  Declaration of **Public** data
-  Modules may NOT contain the following:
  -  Variables (_i.e._ specifications without the **Parameter** attribute)
-   By default everything in a module shall be made private, explicitly
    make public what is needed outside the module or type using the **Public** statement.
-   Data which is used only in the defining module should be declared
    private.
-   Each module should be in a separate file.
-   Module names shall match their file names.
-   When using a module with the **Use** statement, **Only** must also be used
-   While overloading operators may be tempting, it is best avoided if
    you prefer performance to aesthetical beauty.

```fortran
Use domains_module, Only : map
```
-   User derived types should be created to contain data relevant to the module, preferably one per module.
-   Types may match the module name but append them with **_type**.
-   Types shall provide init and **Final** procedures, optionally a summary
    procedure.
-   If provided a summary procedure shall produce valid YAML 1.2 output.

## Specification Statements

-   Align declaration statements clearly separating types and attributes
    from names.
-   Separate variable declaration block from code block of a subroutine
    by a blank line.
-   Use separate declaration blocks for routine arguments and local
    variables.

```fortran
!> Calculate the factorial of n
Function factorial(n) Result(res)
  !> The integer to calculate the factorial of
  Integer( Kind = wi ), Intent( In    ) :: n
  !> The factorial of n
  Integer( Kind = wi ) :: res

  !> Running total
  Integer( Kind = wi ) :: total
  !> Loop iterator
  Integer( Kind = wi ) :: i

  ...
```

-   Always use :: to separate types and attributes from names.
-   Use :: as an alignment guide if the list of names is too long.
-   Separate logical declaration blocks can be aligned differently to
    save screen real estate, *e.g.* parameters vs internal variables of a
    routine.

```fortran
Integer,           Intent( In    ) :: imcon,mxshak
Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
Real( Kind = wp ), Intent(   Out ) :: strcon(1:9)
Real( Kind = wp ), Intent(   Out ) :: vircon

Logical           :: safe
Integer           :: fail(1:2),i,j,k,icyc
Real( Kind = wp ) :: amti,amtj,dli,dlj, &
                     gamma,gammi,gammj,tstep2
```

## Comments

-   Comments shall be indented to align with the code
-   Comments shall be written in lower case except for proper nouns and
    standard abreviations.
-   Comments shall explain what a code does and why, not how it does
    it. Let the code explain how it is done.

## Ford and Doxygen

-  By conforming to the following style useful developer documentation may be created automatically using either FORD or Doxygen.
-  Comments attached to program units, variables and derived types may be automatically documented.
-  Documentation must precede the declaration of the unit, variable or derived type.
-  Comments to be documented must use the tag "!>" (This is default tag in FORD for a comment preceding the content. In Doxygen, "!>" is the only tag which can both start and continue a comment. So this seems to be the best compromise to make the source compatible with both systems. FORD does not like inline comments).
-  To insert a line break in a multiline comment use a blank comment line.
-  Comments use markdown syntax for formatting.
-  LaTeX style maths may be used to include equations in the documentation.
-  See the example structure for more comprehensive examples of documentation comments.

```fortran
!># Example
!>
!> Author - John Smith
!>
!> An example program
Program example
  Use kinds, Only : wi
  Implicit None

  !> Integer counter
  Integer( Kind = wi ) :: i

  ...

  !> Calculate the factorial of n
  Function fact(n)

  ...
```

## Procedures

-   A **Function** shall be pure (with no side-effects). If side-effect are
    needed use a **Subroutine**.
-   All arguments of a **Function/Subroutine** shall have an **Intent**
    attribute (with the exception of the **Class** in type bound procedures).
-   Avoid using **Recursive** procedures.
-   When you are passing an array argument to a **Subroutine/Function**,
    and the **Subroutine/Function** does not change the size of the
    array, you should pass it as an assumed shape array. Memory
    management of such an array is automatically handled by the
    **Subroutine/Function**, and you do not have to worry about having
    to **Allocate** or **Deallocate** your array. It also helps the
    compiler to optimise the code.

## Allocatable Data and Pointers

-   If possible **Allocate** and **Deallocate** for an array or pointer
    shall appear in the same scope.
-   For **Allocatable** members of a user defined type, allocation shall
    happen in the **init** and deallocation in the **final** subroutine.
-   In all cases, **Deallocate** in the opposite order you did **Allocate**.

```fortran

  Allocate(xxx(n),Stat = stat)
  Allocate(yyy(n),Stat = stat)
  Allocate(zzz(n),Stat = stat)

  Deallocate(zzz)
  Deallocate(yyy)
  Deallocate(xxx)
```

-   If using **Pointer**, define it before usage by pointing to **Null**
    or using **Nullify**.
-   Similarly, when a pointer is not used anymore nullify it using the
    same techniques as above.
