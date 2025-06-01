## Overview
The `cneupd` subroutine is invoked after a successful convergence of `cnaupd`. Its purpose is to extract the converged approximations to eigenvalues (Ritz values) for the complex eigenvalue problem A*z = lambda*B*z. Optionally, it can also compute:
1. The corresponding approximate eigenvectors (Ritz vectors).
2. An orthonormal basis for the associated approximate invariant subspace (Schur vectors).
3. Both Ritz and Schur vectors.

These quantities are derived from the Arnoldi factorization computed by `cnaupd` for the complex linear operator OP (defined by the `MODE` selection in `cnaupd`). `cneupd` handles the transformation of these OP-related values and vectors back to the original problem's eigensystem, particularly if a shift-and-invert strategy was employed.

## Key Components
- **`cneupd` (subroutine):** The main routine for processing results from a converged `cnaupd` run. It computes Ritz values and, if requested, Ritz vectors or Schur vectors, including necessary transformations from the OP domain (used in `cnaupd`) to the original problem A*z = lambda*B*z.

## Important Variables/Constants
- **`RVEC` (LOGICAL, INPUT):**
    - `.FALSE.`: Compute Ritz values only.
    - `.TRUE.`: Compute Ritz vectors or Schur vectors.
- **`HOWMNY` (Character*1, INPUT):** Specifies the form of the basis for the invariant subspace.
    - `'A'`: Compute NEV Ritz vectors.
    - `'P'`: Compute NEV Schur vectors.
    - `'S'`: Compute selected Ritz vectors (via `SELECT` array). (Note: `'S'` might not be fully implemented in all versions).
- **`SELECT` (Logical array, INPUT/WORKSPACE):** If `HOWMNY = 'S'`, indicates which Ritz vectors to compute. Used as workspace otherwise.
- **`D` (Complex array, OUTPUT):** Contains the Ritz value approximations to the eigenvalues of A*z = lambda*B*z.
- **`Z` (Complex array, OUTPUT):** If `RVEC = .TRUE.` and `HOWMNY = 'A'`, this N by NEV array holds the Ritz vectors.
- **`LDZ` (Integer, INPUT):** Leading dimension of `Z`.
- **`SIGMA` (Complex, INPUT):** The shift value used in `cnaupd` if `IPARAM(7)` (MODE) was 3 (shift-and-invert). Not referenced for Mode 1 or 2.
- **`WORKEV` (Complex array, WORKSPACE):** Work array of dimension 2*NCV.
- **`BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO`:** These arguments *must* be identical to those from the last successful call to `cnaupd` and must not be modified.

    - **`V` (Complex array, INPUT/OUTPUT):** On input, contains Arnoldi basis from `cnaupd`. On output (if `RVEC = .TRUE.`), the first `NCONV=IPARAM(5)` columns contain approximate Schur vectors. If `Z` shares memory with `V` and Ritz vectors (`HOWMNY='A'`) are requested, `V` is overwritten by Ritz vectors.
    - **`WORKL` (Complex array, OUTPUT/WORKSPACE):** Contains data from `cnaupd` and is further used/populated by `cneupd`. `IPNTR(9:13)` point to locations for original system Ritz values, error bounds, the Schur matrix T, and the eigenvector matrix of H.
    - **`RWORK` (Real array, WORKSPACE):** Passed from `cnaupd`.
    - **`INFO` (Integer, OUTPUT):** Error flag. Relevant values include:
        - `1`: Schur form reordering by `ctrsen` failed.
        - `-8`: Error in LAPACK eigenvalue calculation (e.g., `clahqr`).
        - `-9`: Error in LAPACK eigenvector calculation (`ctrevc`).
        - `-14`: `cnaupd` found no eigenvalues.
        - `-15`: Discrepancy in converged Ritz value count between `cnaupd` and `cneupd`.

## Usage Examples
`cneupd` is called subsequent to `cnaupd` converging (typically `INFO = 0` from `cnaupd`).

```fortran
C ... (previous cnaupd reverse communication loop) ...

      IF (INFO .EQ. 0) THEN
C        cnaupd converged, call cneupd
         RVEC = .TRUE.      ! Request vectors
         HOWMNY = 'A'       ! Request Ritz vectors

         CALL CNEUPD ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV,
     & BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,
     & IPARAM, IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO_EUPD )

         IF (INFO_EUPD .EQ. 0) THEN
C           Successfully retrieved D (eigenvalues)
C           and Z (eigenvectors if RVEC=.TRUE., HOWMNY='A')
C           Process results...
         ELSE
C           Error in cneupd
         END IF
      ELSE
C        Error or other condition in cnaupd
      END IF
```
For comprehensive examples, consult the `EXAMPLES/COMPLEX` directory in the ARPACK source code, specifically `cndrv1.f`. The `ex-complex.doc` file, if present, also offers guidance.

## Dependencies and Interactions
- **Internal ARPACK routines called:**
    - `ivout`: Prints integers.
    - `cmout`: Prints complex matrices.
    - `cvout`: Prints complex vectors.
- **LAPACK routines called:**
    - `cgeqr2`: QR factorization.
    - `clacpy`: Matrix copy.
    - `clahqr`: Schur form of a Hessenberg matrix.
    - `claset`: Matrix initialization.
    - `ctrevc`: Eigenvectors of an upper triangular matrix.
    - `ctrsen`: Reorders Schur form.
    - `cunm2r`: Applies unitary matrix.
    - `slamch`: Machine constants (real-valued).
- **BLAS routines called:**
    - `ctrmm` (Level 3): Matrix-triangular matrix product.
    - `cgeru` (Level 2): Rank-one update (unsymmetric).
    - `ccopy` (Level 1): Vector copy.
    - `cscal` (Level 1): Vector scale by complex scalar.
    - `csscal` (Level 1): Vector scale by real scalar.
    - `scnrm2` (Level 1): Vector 2-norm (real result for complex vector).
    - `ccdotc` (Level 1): Dot product of complex vectors, conjugating the first.
- **Interactions:**
    - `cneupd` is critically dependent on a successful run of `cnaupd`.
    - It uses the same arguments as `cnaupd`, which must be passed unmodified.
    - `WORKL` (complex) and `RWORK` (real) arrays carry essential data from `cnaupd` (like the Hessenberg matrix H) and are used by `cneupd` to store intermediate and final results (Schur form T, eigenvectors of H, transformed Ritz values/vectors).
    - `IPNTR` maps to various data sections within `WORKL`.
    - If `RVEC = .TRUE.`, `V` is altered to hold Schur vectors or potentially Ritz vectors if `Z` and `V` share memory.
    - For shift-and-invert mode (`IPARAM(7)=3`), `SIGMA` is used to transform eigenvalues from the OP system (1/(lambda-sigma)) back to the original problem's lambda. Ritz vectors are generally invariant under this transformation but might be scaled or subject to a final purification step.
```
