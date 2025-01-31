# Evaluation of Modular Polynomials

This repository contains a C++ implementation of the algorithms in the research paper *Evaluation of Modular Polynomials from Supersingular Elliptic Curves* by Maria Corte‑Real Santos, Jonathan Komada Eriksen, Antonin Leroux, Michael Meyer, and Lorenz Panny.

## Requirements

All the code in this repository is written in C++, and uses
* [NTL](https://libntl.org/) (tested on version 11.5.1)
* [fplll](https://github.com/fplll/fplll) (tested on version 5.4.5)

## Functionality

The code can be used to compute classical modular polynomials, respectively Weber polynomials, of prime level $\ell$ evaluated at a $j$-invariant, respectively Weber invariant, defined over $\mathbb{F}_p$. To compile:

```bash
$ make -j
```

The general usage is

```bash
./main bigl j [p] [l] [j]  # p = characteristic, l = level, j = j-invariant
./main bigl w [p] [l] [w]  # p = characteristic, l = level, w = Weber invariant
./main bigc w [p] [l] [w]  # p = characteristic, l = level, w = Weber invariant
```
Alternatively, run `make -j debug` and call `./debug` with the same possible arguments to run the code single-threaded and with assertions enabled.

For example, to compute the evaluated Weber polynomial $\Phi^\mathfrak{f}_{11}(X, 2)$ over $\mathbb{F}_p$ for $p = 1073741827$, run
```bash
$ ./main bigl w 1073741827 5 2
#....
Done!
Coefficient of x^0 is: 4096
Coefficient of x^1 is: 64
Coefficient of x^2 is: 0
Coefficient of x^3 is: 1073741123
Coefficient of x^4 is: 0
Coefficient of x^5 is: 2816
Coefficient of x^6 is: 0
Coefficient of x^7 is: 1073736195
Coefficient of x^8 is: 0
Coefficient of x^9 is: 5632
Coefficient of x^10 is: 0
Coefficient of x^11 is: 1073739779
Coefficient of x^12 is: 1
```
for evaluation using the "ModularEvaluationBigLevel" algorithm, or use the `bigc` argument to select the "ModularEvaluationBigCharacteristic" algorithm.

### How to reproduce the experiments

To reproduce the experiments detailed in the accompanying paper run the following command in the terminal:

```bash
$ ./experiments.sh
```

WARNING: The experiments should only be run on a machine with sufficient computing power. In total, the runtime exceeds 1000 core days.

## Additional Comments

Various parts of our codebase may find utility outside of their application in this paper. As a few suggestions:
- `interpolation.cpp`: an optimised implementation of polynomial interpolation, which outperforms NTL’s built-in function for polynomials of large degree
- `id2iso.cpp`, `klpt.cpp`, and others: a low-level optimised implementation of the Deuring correspondence for generic primes.
- `endring.hpp`: an optimised implementation of the computation of the endomorphism ring of a supersingular elliptic curve defined over $\mathbb{F}_p$

## Citations

- The file `crt.cpp` is heavily inspired by the file `crt.c` in Sutherland’s software [classpoly](https://math.mit.edu/~drew/classpoly.html), which implements the algorithms described in the following research papers:
    - Andrew V. Sutherland, [*Computing Hilbert class polynomials with the Chinese Remainder Theorem*](https://arxiv.org/abs/0903.2785), [Mathematics of Computation](https://www.ams.org/journals/mcom/2011-80-273/S0025-5718-2010-02373-7/home.html) **80** (2011), 501–538.
    - Andreas Enge and Andrew V. Sutherland, [*Class invariants by the CRT method*](https://arxiv.org/abs/1001.3394), [Ninth Algorithmic Number Theory Symposium (ANTS IX)](https://link.springer.com/chapter/10.1007/978-3-642-14518-6_14), Lecture Notes in Computer Science **6197** (2010), 142–156.
- The algorithms in `interpolation.cpp` follow §10 in [Modern Computer Algebra](https://www.cambridge.org/core/books/modern-computer-algebra/DB3563D4013401734851CF683D2F03F0) by Joachim von zur Gathen and Jürgen Gerhard.
- The code in `choosetorsion.cpp` and `costmodel.cpp` mainly translate the SageMath code from [Deuring from the People](https://github.com/friends-of-quaternions/deuring) by Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni.
- The product tree functionality constructed in `utils.hpp` is inspired from the [SCALLOP code](https://github.com/isogeny-scallop/scallop) due to Luca De Feo, Tako Boris Fouotsa, Péter Kutas, Antonin Leroux, Simon‑Philipp Merz, Lorenz Panny, and Benjamin Wesolowski.
- The prime sieve in `utils.cpp` is inspired from the code for the [PTE-Sieve](https://github.com/microsoft/twin-smooth-integers/tree/main) by Craig Costello, Patrick Longa, Michael Meyer, and Michael Naehrig.
