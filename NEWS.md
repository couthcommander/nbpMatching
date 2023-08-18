# nbpmatching 1.5.4
* fix old Fortran warning
* update maintainer email

# nbpmatching 1.5.3
* additional check for singular matrix in `gendistance`; thanks Jonathan Chipman
* new `gendistance` argument `outRawDist`

# nbpmatching 1.5.2
* properly register Fortran calls

# nbpmatching 1.5.1

* Bug fix, handles distancematrix with all zeroes
* Rework how `nonbimatch` handles precision when input includes Inf distances
* Bug fix, `nonbimatch` passes precision to Fortran call
* `nonbimatch` creates class of the same name
* add `subsetMatches` method
* include nonbimatch-method for `assign.grp`, `qom`, `get.sets`
* `gendistance` produces Inf distance to phantoms on elements with forced matches
* documentation fix for undocumented S4 methods warning
* moving suggested package testthat outside of package

# nbpmatching 1.5.0

* Updated NAMESPACE to include default packages other than base
* Removed distance scaling in `gendistance`
* `gendistance` and `make.phantoms` use max distance Inf
* Removed integer constraint on `distancematrix`
* Update documentation to reference Mahalanobis distance calculation

# nbpmatching 1.4.5

* Improved package description
* Bug fix, specifying threshold with nonbimatch failed on R 3.2.0; thanks Philippe Sulger.
