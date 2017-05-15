



warp_fct_til_optim <- function(params, observations, ode_system, ode_params, x0, zero = F,
                               time_points, penalisering, dim = 3, pen_fct = NULL, ...) {
  if (!zero) params <- c(0, params)
  if (class(ode_system) == "ODE") {
    j <- attr(ode_system, 'address')
  set_warp_params(j, matrix(params))
   yfit <- curve_solver(attr(ode_system, 'address'), ode_params, x0, time_points)$traj[, 1:dim]
  ud <- sum((t(observations) - yfit)^2) 
  print(ud)
  ud <- ud +
  if (!is.null(pen_fct)) penalisering*pen_fct(params, ...) else penalisering*sum(params^2)
  }
  ud
}


# #' Title
# #'
# ##' @param params 
# #'
# #' @return
# #' @export
# #'
# #' @examples
#Square.penalty <- function(params) {
#  sum(diff(c(0, params))^2)
#}