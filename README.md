# Nonlinear Optimization Algorithms Multilang [![Build Status](https://travis-ci.org/rgolubtsov/nonlinear-optimization-algorithms-multilang.svg?branch=master)](https://travis-ci.org/rgolubtsov/nonlinear-optimization-algorithms-multilang)

**Note:** The Travis CI build status badge above is applicable to the following packages ONLY (for the moment):
* The (original) Hooke and Jeeves :small_blue_diamond: Traditional (K&amp;R) C implementation
* The Hooke and Jeeves :small_blue_diamond: ISO C (C99/11) implementation
* The Hooke and Jeeves :small_blue_diamond: ISO C++ (C++98/03/11/14) implementation
* The Hooke and Jeeves :small_blue_diamond: Objective-C 2.0 implementation
* The Hooke and Jeeves :small_blue_diamond: FORTRAN 77 (MIL-STD 1753) implementation
* The Hooke and Jeeves :small_blue_diamond: ISO Fortran 95 implementation
* The Hooke and Jeeves :small_blue_diamond: Perl 5 (5.10+) implementation
* The Hooke and Jeeves :small_blue_diamond: JavaScript (ECMA-262 5.1) implementation
* The Hooke and Jeeves :small_blue_diamond: Java (Java SE 7) implementation
* The Hooke and Jeeves :small_blue_diamond: Go (golang) implementation
* The Hooke and Jeeves :small_blue_diamond: Python 2 &amp; 3 implementation

---

* The (original) Nelder-Mead :small_blue_diamond: FORTRAN 77 (MIL-STD 1753) implementation
* The Nelder-Mead :small_blue_diamond: ISO C (C99/11) implementation

---

This project is aimed at implementing nonlinear programming algorithms as the (un-)constrained minimization problems with the focus on their numerical expression using various programming languages.

The idea behind is to collect some of nonlinear programming (NLP) algorithms/methods which are well known and developed, and when they have low or moderate level of complexity during their coding in programming languages in an analytical form. Furthermore, such NLP-algorithms will have to be applied to one or more test problems in numerical expression, and then the user will be able to get solutions for those problems and watch the digits.

All of these algorithms should be coded in a series of programming languages. Not every subsequent selected algorithm must be implemented in every particular language chosen for the previous one, but in some of them (the standard C is an exception) which are widely used nowadays. The main requirement here is the resulting output containing the solution for any kind of test problem must be the same for all implementations. I.e. output blocks made after computations produced by the C code, Java code, Fortran code and so on should be identically equal each other with respect to computed precision, precision expression, layout and spacing emitted to standard out (or user console).

*Consider this project as a somewhat __educational__ approach to the subject of implementing math algorithms in programming languages rather than it might be considered otherwise as bringing something important to scientific applications and investigations.*

Why NLP-algorithms, exactly? &ndash; The answer is why not? Indeed, NLP optimization methods and techniques are widely used everywhere. Their complexity ranges are quite different. They are suitable and fit the most optimization models more precisely and accurately than their linear optimization counterparts. &ndash; That's why.

There are two main conventional groups of NLP optimization methods: constrained and unconstrained. In this project it is considered to be implemented first two unconstrained minimization methods: (1) The algorithm of **Hooke and Jeeves** and (2) The **Nelder-Mead** algorithm.

## Overview

The project has the following directory structure and logical parts and items.

```
.
|-- nlp-unconstrained-cli              <== Unconstrained methods impl. container
|   |                                      (including test problems)
|   |-- hooke-jeeves                   <== Hooke-Jeeves algorithm container
|   |   |                                  (all implementations)
|   |   |-- c                          <== ISO C (C99/11) impl. container
|   |   |   |-- Makefile
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- cc                         <== ISO C++ (C++98/03/11/14) impl. container
|   |   |   |-- Makefile
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- f77                        <== FORTRAN 77 (MIL-STD 1753) impl. container
|   |   |   |-- Makefile
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- f95                        <== ISO Fortran 95 impl. container
|   |   |   |-- Makefile
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- go                         <== Go (golang) impl. container
|   |   |   |-- Makefile
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- java                       <== Java (Java SE 7) impl. container
|   |   |   |-- Makefile
|   |   |   |-- pom.xml
|   |   |   `-- src
|   |   |       `-- main
|   |   |           `-- java
|   |   |               `-- optimization
|   |   |                   `-- nonlinear
|   |   |                       `-- unconstrained
|   |   |                           `-- cli
|   |   |                               `-- < sources >
|   |   |
|   |   |-- js                         <== JavaScript (ECMA-262 5.1) impl. container
|   |   |   |-- Makefile
|   |   |   `-- src
|   |   |       `-- < sources >
|   |   |
|   |   |-- objc                       <== Objective-C 2.0 impl. container
|   |   |   |-- Makefile                   (using the GNUstep Base library)
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- __orig                     <== Traditional (K&R) C impl. container
|   |   |   |-- Makefile                   (This is the original impl.)
|   |   |   `-- src
|   |   |       |-- < sources >
|   |   |       `-- Makefile
|   |   |
|   |   |-- perl                       <== Perl 5 (5.10+) impl. container
|   |   |    |-- Makefile
|   |   |    `-- src
|   |   |        |-- NLPUCCLIHooke
|   |   |        |   `-- < sources >
|   |   |        `-- < sources >
|   |   |
|   |   `-- python                     <== Python 2 & 3 impl. container
|   |       |-- Makefile
|   |       `-- src
|   |           |-- < sources >
|   |           `-- nlpucclihooko
|   |               `-- < sources >
|   |
|   `-- nelder-mead                    <== Nelder-Mead algorithm container
|       |                                  (all implementations)
|       |-- c                          <== ISO C (C99/11) impl. container
|       |   |-- Makefile
|       |   `-- src
|       |       |-- < sources >
|       |       `-- Makefile
|       |
|       `-- __orig                     <== FORTRAN 77 (MIL-STD 1753) impl. container
|           |-- Makefile                   (This is the original impl.)
|           `-- src
|               |-- < sources >
|               `-- Makefile
|
|-- [ nlp-unconstrained-api ]          <== [TODO] Unconstrained methods impl. container
|                                                 (API only)
|-- [ nlp-constrained-cli ]            <== [TODO] Constrained methods impl. container
|                                                 (including test problems)
`-- [ nlp-constrained-api ]            <== [TODO] Constrained methods impl. container
                                                  (API only)
```

As shown above, the directories that should contain ~~stuff for the Nelder-Mead algorithm implementations~~ as well as other three ones (`nlp-unconstrained-api`, `nlp-constrained-cli`, and `nlp-constrained-api`) are not yet exist. But it's planned they have to be created and populated accordingly somewhen during development process.
