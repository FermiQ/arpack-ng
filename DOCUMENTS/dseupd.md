## Overview
The `dseupd` subroutine is used after `dsaupd` has successfully completed. Its purpose is to return the converged approximations to eigenvalues (Ritz values) of the problem A*z = lambda*B*z and, optionally, the corresponding approximate eigenvectors (Ritz vectors) or an orthonormal (Lanczos) basis for the associated approximate invariant subspace.

These quantities are derived from the Lanczos factorization computed by `dsaupd` for the linear operator OP, which is defined by the `MODE` selection in `dsaupd`. The approximate eigenvalues of the original problem are returned in ascending algebraic order.

## Key Components
- **`dseupd` (subroutine):** The main routine responsible for extracting Ritz values and Ritz vectors (or Lanczos basis) from the results of a converged `dsaupd` run. It transforms these quantities back to the original problem's eigensystem if a spectral transformation was used in `dsaupd`.

## Important Variables/Constants
- **`RVEC` (LOGICAL, INPUT):** Specifies whether Ritz vectors are computed.
    - `.FALSE.`: Compute Ritz values only.
    - `.TRUE.`: Compute Ritz vectors.
- **`HOWMNY` (Character*1, INPUT):** Specifies how many Ritz vectors are wanted.
    - `'A'`: Compute all NEV Ritz vectors.
    - `'S'`: Compute a selection of Ritz vectors (specified by `SELECT` array). (Note: The documentation mentions 'S' is not fully implemented in some older versions).
- **`SELECT` (Logical array, INPUT/WORKSPACE):** If `HOWMNY = 'S'`, specifies which Ritz vectors to compute. Used as workspace if `HOWMNY = 'A'`.
- **`D` (Double precision array, OUTPUT):** Contains the Ritz value approximations to the eigenvalues of A*z = lambda*B*z, in ascending order.
- **`Z` (Double precision array, OUTPUT):** If `RVEC = .TRUE.`, contains the B-orthonormal Ritz vectors corresponding to the Ritz values in `D`.
- **`LDZ` (Integer, INPUT):** Leading dimension of `Z`.
- **`SIGMA` (Double precision, INPUT):** The shift value used if `IPARAM(7)` (MODE in `dsaupd`) was 3, 4, or 5.
- **`BMAT`, `N`, `WHICH`, `NEV`, `TOL`, `RESID`, `NCV`, `V`, `LDV`, `IPARAM`, `IPNTR`, `WORKD`, `WORKL`, `LWORKL`, `INFO`:** These arguments *must* be the same as those passed to the final successful call of `dsaupd`. They should not be modified between the `dsaupd` call and the `dseupd` call.

    - **`WORKL` (Double precision array, OUTPUT/WORKSPACE):** Contains information from `dsaupd` and is further used by `dseupd`. `IPNTR(8:10)` point to locations within `WORKL` for Ritz values of the original system, error bounds, and the eigenvector matrix of the tridiagonal matrix T.
    - **`INFO` (Integer, OUTPUT):** Error flag. Specific negative values indicate errors similar to `dsaupd`, with additional ones like:
        - `-14`: `dsaupd` did not find any eigenvalues to sufficient accuracy.
        - `-15`: `HOWMNY` invalid.
        - `-17`: Discrepancy in the count of converged Ritz values between `dsaupd` and `dseupd`.

## Usage Examples
`dseupd` is called *after* `dsaupd` has converged (i.e., `INFO = 0` from `dsaupd`).

```fortran
c ... (previous dsaupd reverse communication loop) ...

if (INFO .eq. 0) then
c   dsaupd converged, now call dseupd
    RVEC = .TRUE.     c To compute Ritz vectors
    HOWMNY = 'A'      c Compute all NEV vectors

    call dseupd ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA,
& BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,
& IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO_EUPD )

    if (INFO_EUPD .eq. 0) then
c       Successfully retrieved D (eigenvalues) and Z (eigenvectors)
c       Process the results...
    else
c       Error in dseupd
    end if
else if (INFO .eq. 1) then
c   Max iterations reached in dsaupd.
c   Optionally, can still call dseupd to get whatever converged.
    RVEC = .TRUE.
    HOWMNY = 'A'
c   Note: IPARAM(5) from dsaupd contains the number of converged Ritz values (NCONV)
c   NEV_actual = IPARAM(5)
    call dseupd ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA,
& BMAT, N, WHICH, IPARAM(5), TOL, RESID, NCV, V, LDV,
& IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO_EUPD )
c    ...
else
c   Error in dsaupd
end if
```
For specific examples, refer to the `EXAMPLES/SYM` directory in the ARPACK source (e.g., `dsdrv1.f`, `dsdrv2.f`). These drivers show the typical sequence: call `dsaupd` until convergence, then call `dseupd` to retrieve results.

## Dependencies and Interactions
- **Internal ARPACK routines called:**
    - `dsesrt`: Sorts an array and applies permutation to a matrix.
    - `dsortr`: Sorting routine.
    - `ivout`: Utility for printing integers.
    - `dvout`: Utility for printing vectors.
- **LAPACK routines called:**
    - `dgeqr2`: Computes QR factorization.
    - `dlacpy`: Matrix copy.
    - `dlamch`: Determines machine constants.
    - `dorm2r`: Applies an orthogonal matrix in factored form.
    - `dsteqr`: Computes eigenvalues/eigenvectors of a tridiagonal matrix.
- **BLAS routines called:**
    - `dger` (Level 2): Rank one update.
    - `dcopy` (Level 1): Vector copy.
    - `dnrm2` (Level 1): Vector norm.
    - `dscal` (Level 1): Vector scale.
    - `dswap` (Level 1): Vector swap.
- **Interactions:**
    - `dseupd` is critically dependent on the successful completion of `dsaupd`.
    - It uses many of the same arguments passed to `dsaupd`, which must not be altered.
    - The `WORKL` array, populated by `dsaupd` and further utilized by `dseupd`, is essential for storing intermediate results like the tridiagonal matrix H, its eigenvectors, and Ritz estimates.
    - `IPNTR` provides pointers to various data structures within `WORKL`.
    - If spectral transformations (shift-invert, etc.) were used in `dsaupd` (indicated by `IPARAM(7)`), `dseupd` transforms the Ritz values (and vectors if computed) back to the original problem's eigenpairs using the `SIGMA` parameter.
    - The results `D` (eigenvalues) and `Z` (eigenvectors) are for the original problem A*z = lambda*B*z.
```
