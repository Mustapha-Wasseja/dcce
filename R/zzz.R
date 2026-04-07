# Package-level hooks
# S3 method registrations for our own generics are handled via roxygen2
# @export tags in the method files, not here. The .onLoad hook is used
# solely to register S3 methods on OPTIONAL (Suggests) packages' generics
# when those packages are available.

.onLoad <- function(libname, pkgname) {
  # Register marginaleffects methods if that package is installed.
  # This keeps marginaleffects as a soft dependency in Suggests while
  # still providing seamless integration when it is available.
  tryCatch(
    .register_marginaleffects_methods(),
    error = function(e) invisible()
  )
}
