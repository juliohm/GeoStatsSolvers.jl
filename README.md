# GeoStatsSolvers.jl

[![][build-img]][build-url] [![][codecov-img]][codecov-url]

Built-in solvers for the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

## Estimation

#### IDW

This solver provides a high-performance implementation of the inverse distance weighting scheme
introduced in the very early days of geostatistics (see [Shepard 1968](https://dl.acm.org/citation.cfm?id=810616)).
It is perhaps the simplest first attempt in the literature to perform estimation based on the
notion of proximity to data locations.

This implementation makes use of k-d trees from the [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl)
package, which leads to a fast estimation method for large or high-resolution spatial domains.
Although this method is recommended for fast assessment of a new field, it has poor statistical
properties (lacks covariance model) and should mainly be used for qualitative purposes.

#### LWR

This solver provides an implementation of locally weighted regression (a.k.a. LOESS) introduced by
[Cleveland 1979](http://www.stat.washington.edu/courses/stat527/s13/readings/Cleveland_JASA_1979.pdf).
It is the most natural generalization of inverse distance weighting in which one is allowed to use a
custom weight function instead of distance-based weights.

Like in the inverse distance weighting solver, this solver makes use of k-d trees from the
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) package for fast data
lookup. Locally weighted regression (LWR or LOESS) is a popular non-parametric method, however
it still has poor statistical properties compared to other estimation methods such as Kriging
that explicitly model spatial correlation.

In the current implementation, the estimation variance is computed assuming Gaussian residuals.

#### Kriging

This polyalgorithm solver provides an implementation of various forms of Kriging introduced by
[Matheron 1971](https://books.google.com.br/books/about/The_Theory_of_Regionalized_Variables_and.html).
Kriging is a popular method in various industries due to its good statistical properties and flexibility.
Unlike the previous solvers, Kriging relies on the specification of a variogram model, which can be
fit to geospatial data.

## Simulation

#### LUGS

This solver provides an implementation of direct Gaussian simulation (a.k.a. LU simulation)
as described in [Alabert 1987](https://link.springer.com/article/10.1007/BF00897191). In this
method, the full covariance matrix is built to include all locations of the simulation domain,
and samples from the multivariate Gaussian are drawn via LU factorization.

The method, which is widely implemented in many packages for Gaussian processes (e.g.
[GaussianProcesses.jl](https://github.com/STOR-i/GaussianProcesses.jl)),
is appropriate for relatively small simulation domains (e.g. 100x100 grids) where it is feasible
to factorize the full covariance. For larger domains (e.g. 3D grids), other methods are available
such as sequential Gaussian simulation, spectral methods, and FFT moving averages.

#### FFTGS

This solver provides an implementation of spectral Gaussian simulation (a.k.a. FFT simulation)
as described in [Gutjahr 1997](https://link.springer.com/article/10.1007/BF02769641).
In this method, the covariance function is perturbed in the frequency
domain after a fast Fourier transform. White noise is added to the phase
of the spectrum, and a realization is produced by an inverse Fourier transform.

The method is limited to simulations on regular grids, and care must be taken
to make sure that the correlation length is small enough compared to the grid
size. As a general rule of thumb, avoid correlation lengths greater than 1/3
of the grid. The method is extremely fast, and can be used to generate large
3D realizations.

#### SGS

This solver provides an implementation of sequential Gaussian simulation as described in
[Gomez-Hernandez & Journel 1993](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8).
The method traverses all locations of the spatial domain according to a path, approximates the
conditional distribution function within a neighborhood using Kriging, and assigns a value to
the center of the neighborhood by sampling from this distribution.

#### SPDEGS

This solver provides an implementation of the stochastic partial differential equations
(SPDE) approach to Gaussian simulation [Lindgren et al. 2011.](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x).
The method relies on a discretization of the Laplace-Beltrami operator on meshes and is
adequate for highly curved domains (e.g. surfaces).

## Learning

#### PointwiseLearn

A geostatistical learning solver that takes a classical (non-spatial)
statistical learning model from the literature and applies it to each
point of the domain. Although this solver is widely used in geospatial
problems, it has various issues as described in
[Hoffimann et al. 2021](https://arxiv.org/abs/2102.08791).

## Asking for help

If you have any questions, please [contact our community](https://juliaearth.github.io/GeoStats.jl/stable/about/community.html).

[build-img]: https://img.shields.io/github/actions/workflow/status/JuliaEarth/GeoStatsSolvers.jl/CI.yml?branch=main&style=flat-square
[build-url]: https://github.com/JuliaEarth/GeoStatsSolvers.jl/actions

[codecov-img]: https://img.shields.io/codecov/c/github/JuliaEarth/GeoStatsSolvers.jl?style=flat-square
[codecov-url]: https://codecov.io/gh/JuliaEarth/GeoStatsSolvers.jl
