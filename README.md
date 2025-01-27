# ForBES

**ForBES** (standing for **For**ward-**B**ackward **E**nvelope **S**olver) is a MATLAB solver for
nonsmooth optimization problems.

It is generic in the sense that the user can customize the problem to solve in an easy and flexible way.
It is efficient since it features very efficient algorithms, suited for large scale applications.

For full documentation refer to the [ForBES webpage](http://kul-optec.github.io/ForBES/) .

## Installation

Simply clone the git repository, or click on this [link](https://github.com/kul-forbes/ForBES/archive/master.zip)
to download it as a zip archive and decompress the archive. Then move with the MATLAB command line to
the directory of ForBES, and execute the following command:

```
>> forbes_setup
```

This will compile all the necessary source files and install the directory into MATLAB's path.

## How to use it

ForBES consists mainly of one MATLAB routine, `forbes`. In order to use it one
must provide a description of the problem and (optionally) a set of options:

```
>> out = forbes(f, g, init, aff, constr, opt);
```

Full documentation, explaining how to specify these arguments, can be
found at the [ForBES webpage](http://kul-forbes.github.io/ForBES/).

Examples on how to use `forbes` can be found in the [demos folder](https://github.com/kul-forbes/ForBES/tree/master/demos).
Furthermore, you can access the help file of the solvers directly from MATLAB with

```
>> help forbes
```

## References

1. L. Stella, A. Themelis, P. Patrinos, “Forward-backward quasi-Newton methods for nonsmooth optimization problems,” [arXiv:1604.08096](http://arxiv.org/abs/1604.08096) (2016).

2. A. Themelis, L. Stella, P. Patrinos, “Forward-backward envelope for the sum of two nonconvex functions: Further properties and nonmonotone line-search algorithms,” [arXiv:1606.06256](http://arxiv.org/abs/1606.06256) (2016).

## Credits and license

ForBES is developed by [Lorenzo Stella](https://lostella.github.io) [`lorenzo.stella-at-imtlucca.it`] and Panos Patrinos [`panos.patrinos-at-esat.kuleuven.be`]
at [IMT Lucca](http://www.imtlucca.it) and [KU Leuven](http://www.kuleuven.be).
Any feedback, bug report or suggestion for future improvements is more than welcome.
We recommend using the [issue tracker](https://github.com/kul-forbes/ForBES/issues) to report bugs.

ForBES is free software and distributed under the [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0.en.html). You may copy, distribute and modify the software provided that modifications are described and licensed for free under LGPL.
