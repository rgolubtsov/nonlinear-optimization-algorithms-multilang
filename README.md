# Nonlinear Optimization Algorithms Multilang [![Build Status](https://travis-ci.org/rgolubtsov/nonlinear-optimization-algorithms-multilang.svg?branch=master)](https://travis-ci.org/rgolubtsov/nonlinear-optimization-algorithms-multilang)

**Note:** The Travis CI build status badge above is applicable to the following packages ONLY (for the moment):
* The (original) Hooke and Jeeves :small_blue_diamond: K&R C implementation
* The Hooke and Jeeves :small_blue_diamond: ANSI C (C89/90/95) implementation

---

This project is aimed at implementing nonlinear programming algorithms as the (un-)constrained minimization problems with the focus on their numerical expression using various programming languages.

The idea behind is to collect some of nonlinear programming (NLP) algorithms/methods which are well known and developed, and when they have low or moderate level of complexity during their coding in programming languages in an analytical form. Furthermore, such NLP-algorithms will have to be applied to one or more test problems in numerical expression, and then the user will be able to get solutions for those problems and watch the digits.

All of these algorithms should be coded in a series of programming languages. Not every subsequent selected algorithm must be implemented in every particular language chosen for the previous one, but in some of them (the standard C is an exception) which are widely used nowadays. The main requirement here is the resulting output containing the solution for any kind of test problem must be the same for all implementations. I.e. output blocks made after computations produced by the C code, Java code, Fortran code and so on should be identically equal each other with respect to computed precision, precision expression, layout and spacing emitted to standard out (or user console).

> Consider this project as a somewhat **educational** approach to the subject of implementing math algorithms in programming languages rather than it might be considered otherwise as bringing something important to scientific applications and investigations.

Why NLP-algorithms, exactly? &ndash; The answer is why not? :smile: NLP optimization methods and techniques are widely used everywhere. Their complexity ranges are quite different. They are suitable and fit the most optimization models more precisely and accurately than their linear optimization counterparts. &ndash; That's why.

There are two main conventional groups of NLP optimization methods: constrained and unconstrained. In this project it is considered to be implemented first two unconstrained minimization methods: :one: The algorithm of **Hooke and Jeeves** and :two: The **Nelder-Mead** algorithm.
