# Efgen
Efgen is a collection of **ef**ficient and **gen**eric implementations of chosen data structures and algorithms.
There is always room for optimization, but Efgen is not about squeezing that last drop of performance.
Rather, these mini-libraries try to provide as generic of an interface as possible first, and then internally
implement algorithms using various datatype-agnostic optimizations such as eliminating recursion etc.
For the most part the intended use case is competetive programming, where problems might normally require
you to modify an existing implementation of a data structure in order to use it in a solution. Instead,
hopefully an Efgen implementation is going to be generic enough that copying and pasting it into a solution
is all that needs to be done.

All the code in this repository is public domain. That said, if you end up using any of it,
I'd appreciate it if you left a comment in your code saying where the fragment comes from.
There are some examples showcasing how to use each mini-library in the [`examples`](/examples) directory.
No build script, you'll figure it out ;)

| Name | Description |
| --- | --- |
| `dft` | Discrete (Fast) Fourier Transform |
| `segtree` | Segment tree over any monoid and set of endomorphisms |
