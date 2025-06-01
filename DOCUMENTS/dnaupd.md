## Overview
The `dnaupd` subroutine is a reverse communication interface for the Implicitly Restarted Arnoldi iteration. It is designed to compute approximations to a few eigenpairs (Ritz values and Ritz vectors) of a general non-symmetric linear operator "OP" with respect to a semi-inner product defined by a symmetric positive semi-definite real matrix B. B can be the identity matrix.

If the operator OP is real and symmetric with respect to B, the subroutine `dsaupd` should be used instead.

`dnaupd` is typically called iteratively to solve one of the following problems:
- **Mode 1:** A*x = lambda*x (standard non-symmetric eigenvalue problem)
  OP = A, B = I.
- **Mode 2:** A*x = lambda*M*x, M symmetric positive definite (generalized non-symmetric eigenvalue problem)
  OP = inv[M]*A, B = M.
- **Mode 3:** A*x = lambda*M*x, M symmetric semi-definite (shift-and-invert mode in real arithmetic)
  OP = Real_Part{inv[A - sigma*M]*M}, B = M.
- **Mode 4:** A*x = lambda*M*x, M symmetric semi-definite (shift-and-invert mode in real arithmetic)
  OP = Imaginary_Part{inv[A - sigma*M]*M}, B = M.

The choice of operator OP (especially in modes 3 and 4) allows for targeting eigenvalues near a (complex) shift `sigma`.

## Key Components
- **`dnaupd` (subroutine):** The primary reverse communication routine that orchestrates the Implicitly Restarted Arnoldi iteration. It manages the computational flow, returning control to the calling program to perform matrix-vector products (OP*x, B*x) or to supply shifts.

## Important Variables/Constants
- **`IDO` (Integer, INPUT/OUTPUT):** Reverse communication flag.
    - `0`: First call.
    - `-1`: Compute Y = OP * X (initialization phase).
    - `1`: Compute Y = OP * X.
    - `2`: Compute Y = B * X.
    - `3`: Compute NP real and imaginary parts of shifts.
    - `99`: Done.
- **`BMAT` (Character*1, INPUT):** Specifies the type of matrix B.
    - `'I'`: Standard eigenvalue problem (B=I).
    - `'G'`: Generalized eigenvalue problem.
- **`N` (Integer, INPUT):** Dimension of the eigenproblem.
- **`WHICH` (Character*2, INPUT):** Specifies which Ritz values to compute (e.g., 'LM', 'SM', 'LR', 'SR', 'LI', 'SI').
- **`NEV` (Integer, INPUT):** Number of eigenvalues to compute (0 < NEV < N-1).
- **`TOL` (Double precision, INPUT/OUTPUT):** Stopping criterion (relative accuracy of Ritz values).
- **`RESID` (Double precision array, INPUT/OUTPUT):** Initial/final residual vector.
- **`NCV` (Integer, INPUT):** Number of Arnoldi vectors. Must satisfy 2 <= NCV-NEV and NCV <= N. Recommended NCV >= 2*NEV+1.
- **`V` (Double precision array, OUTPUT):** Contains the final set of Arnoldi basis vectors.
- **`LDV` (Integer, INPUT):** Leading dimension of V.
- **`IPARAM` (Integer array, INPUT/OUTPUT):** Parameters for the run.
    - `IPARAM(1)`: `ISHIFT` (shift selection method: 0 for user-supplied, 1 for exact shifts).
    - `IPARAM(3)`: `MXITER` (max Arnoldi iterations).
    - `IPARAM(7)`: `MODE` (problem type: 1, 2, 3, or 4).
    - `IPARAM(8)`: `NP` (number of shifts to provide when `IDO=3` and `ISHIFT=0`).
- **`IPNTR` (Integer array, OUTPUT):** Pointers to locations in WORKD and WORKL.
    - `IPNTR(1)`: Pointer for X in WORKD.
    - `IPNTR(2)`: Pointer for Y in WORKD.
    - `IPNTR(3)`: Pointer for B*X in WORKD (Modes 3, 4).
    - `IPNTR(5)`: Pointer to Hessenberg matrix H in WORKL.
    - `IPNTR(6)`: Pointer to real part of Ritz values (RITZR) in WORKL.
    - `IPNTR(7)`: Pointer to imaginary part of Ritz values (RITZI) in WORKL.
    - `IPNTR(8)`: Pointer to Ritz estimates in WORKL.
    - `IPNTR(14)`: Pointer to shifts in WORKL.
- **`WORKD` (Double precision array, REVERSE COMMUNICATION):** Distributed work array for Arnoldi iteration.
- **`WORKL` (Double precision array, OUTPUT/WORKSPACE):** Private work array. LWORKL must be at least 3*NCV**2 + 6*NCV.
- **`LWORKL` (Integer, INPUT):** Length of WORKL.
- **`INFO` (Integer, INPUT/OUTPUT):** Error/status flag.

## Usage Examples
`dnaupd` is used in a reverse communication loop. The user's code must handle `IDO` flags to perform computations or supply data.

```fortran
c Simplified structure of a dnaupd loop:
      IDO = 0
      INFO = 0  c Or non-zero if RESID is pre-set

10    CONTINUE
      CALL DNAUPD ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV,
     & IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )

      IF (IDO .EQ. -1 .OR. IDO .EQ. 1) THEN
C        Compute Y = OP * X
C        X is in WORKD(IPNTR(1)), result Y must be placed in WORKD(IPNTR(2))
C        If Mode 3 or 4, B*X is available in WORKD(IPNTR(3))
         GO TO 10
      ELSE IF (IDO .EQ. 2) THEN
C        Compute Y = B * X
C        X is in WORKD(IPNTR(1)), result Y must be placed in WORKD(IPNTR(2))
         GO TO 10
      ELSE IF (IDO .EQ. 3) THEN
C        Provide NP shifts if IPARAM(1) = 0 (ISHIFT=0)
C        Real parts: WORKL(IPNTR(14) ... IPNTR(14)+NP-1)
C        Imag parts: WORKL(IPNTR(14)+NP ... IPNTR(14)+2*NP-1)
         GO TO 10
      END IF

C     Check INFO for convergence or errors
      IF (INFO .LT. 0) THEN
C        Error occurred
      ELSE
C        INFO = 0: Normal exit, converged. Call dneupd.
C        INFO = 1: Max iterations reached.
C        INFO = 3: No shifts could be applied.
C        ... etc.
      END IF
```
For detailed examples, consult the `EXAMPLES/NONSYM` directory in the ARPACK distribution, particularly files like `dndrv1.f` (standard non-symmetric) or `dndrv2.f` (generalized non-symmetric, for Mode 2 type problems). The file `ex-nonsym.doc` (often included or generated from source comments) also provides context.

## Dependencies and Interactions
- **Internal ARPACK routines called:**
    - `dnaup2`: Implements the core Implicitly Restarted Arnoldi Iteration for non-symmetric problems.
    - `dstatn`: Initializes timing and statistics variables for non-symmetric routines.
    - `ivout`: Utility for printing integers.
    - `arscnd`: Utility for timing.
    - `dvout`: Utility for printing vectors.
- **LAPACK routines called:**
    - `dlamch`: Determines machine constants.
- **Interactions:**
    - `dnaupd` is the primary computational routine for non-symmetric eigenvalue problems.
    - It relies on reverse communication for matrix-vector products (OP*x, B*x) and optionally for shift provision.
    - Upon successful convergence (`INFO=0`), `dneupd` must be called to extract the computed Ritz values and, if desired, Ritz vectors or Schur vectors.
    - `WORKD` and `WORKL` arrays are critical for data exchange and workspace.
    - `IPARAM` and `IPNTR` control the process and map data within work arrays.
    - The `MODE` setting in `IPARAM(7)` significantly alters the definition of OP and the interpretation of results, especially for shift-and-invert strategies.
```
