#' @rdname nbinom2
#' @export
bell <- function(link="log") {
    r <- list(family="bell",
              variance=function(mu){
  V <- mu*(1+LambertW::W(mu))
  return(V)
})
    return(make_family(r,link))
}

dev_bell <- function(y, mu){
  W_y <- LambertW::W(y)
  W_mu <- LambertW::W(mu)
  dev <- ifelse(
    y == 0,
    exp(1-W_mu ),
    exp(W_mu) - exp(W_y) + y *( log(W_y) - log(W_mu) )
  )
  return(dev)
}

##.noDispersionFamilies <- c("binomial", "poisson", "truncated_poisson",
 ##                          "bell")
