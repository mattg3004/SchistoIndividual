#' Setup SchistoIndividual
#'
#' This function initializes Julia and the DifferentialEquations.jl package.
#' The first time will be long since it includes precompilation.
#'
#' @param schist_path Full path (no tildes etc.) to the file SchistoIndividual.jl.
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @examples
#'
#' \donttest{ 
#' ## diffeq_setup() is time-consuming and requires Julia.
#'
#' }
#'
#' @export
SchistoIndividual <- function (...){
  julia <- JuliaCall::julia_setup(...)
  JuliaCall::julia_install_package_if_needed("Distributions")
  JuliaCall::julia_install_package_if_needed("Random")
  JuliaCall::julia_install_package_if_needed("PoissonRandom")
  JuliaCall::julia_install_package_if_needed("JLD")
  JuliaCall::julia_install_package_if_needed("Schistoxpkg")

  JuliaCall::julia_library("Distributions")
  JuliaCall::julia_library("Random")
  JuliaCall::julia_library("PoissonRandom")
  JuliaCall::julia_library("JLD")
  JuliaCall::julia_library("Schistoxpkg")
}
