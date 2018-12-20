#' Model outputs for GBD
#'
#' @param mod simulation model output
#' @param fp model fixed parameters
#'
#' @return A data.frame
#'
#' @examples
#' library(eppasm)
#'
#' pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
#' fp <- prepare_directincid(pjnz)
#' mod <- simmod(fp, "R")
#'
#' gbdout <- get_gbd_outputs(mod, fp)
#' write.csv(gbdout, "gbd-output-run1.csv", row.names=FALSE)
#'
#' @export

get_gbd_outputs <- function(mod, fp) {

  mod <- mod_dimnames(mod, fp$ss)
  hp1 <- hivpop_singleage(mod, fp$ss)
  
  pop <- as.data.frame.table(apply(mod, c(1, 2, 4), sum), responseName = "pop")
  hiv_deaths <- as.data.frame.table(attr(mod, "hivdeaths"),
                                    responseName = "hiv_deaths")
  non_hiv_deaths <- as.data.frame.table(attr(mod, "natdeaths"),
                                        responseName = "non_hiv_deaths")
  new_hiv <- as.data.frame.table(attr(mod, "infections"),
                                 responseName = "new_hiv")
  pop_neg <- as.data.frame.table(mod[,,fp$ss$hivn.idx,])
  
  total_births <- as.data.frame.table(get_births(mod, fp), responseName = "total_births")
  total_births$sex <- "female"
  pregprev <- as.data.frame.table(get_pregprev(mod, fp, hp1), responseName = "pregprev")
  pregprev$sex <- "female"
  
  pop_art <- as.data.frame.table(colSums(hp1$artpop1,,2), responseName = "pop_art")

  hivpop_daly <- as.data.frame.table(get_daly_hivpop(hp1$hivpop), responseName = "value")
  hivpop_daly <- data.table::dcast(hivpop_daly, ... ~ cd4daly)

  v <- pop
  v <- merge(v, hiv_deaths, all.x=TRUE)
  v <- merge(v, non_hiv_deaths, all.x=TRUE)
  v <- merge(v, new_hiv, all.x=TRUE)
  v <- merge(v, pop_neg, all.x=TRUE)
  v <- merge(v, total_births, all.x=TRUE)
  v <- merge(v, pregprev, all.x=TRUE)
  v$hiv_births <- v$total_births * v$pregprev  # number of births to HIV positive women
  v$birth_prev <- NA                           # HIV prevalence among newborns -- not yet output
  v <- merge(v, pop_art, all.x=TRUE)
  v <- merge(v, hivpop_daly, all.x=TRUE)

  v
}


#' Births by single age
#' 
get_births <- function(mod, fp){

  py <- fp$ss$PROJ_YEARS
  fertpop <- apply(mod[fp$ss$p.fert.idx, fp$ss$f.idx, , ] +
                   mod[fp$ss$p.fert.idx, fp$ss$f.idx, , c(1, 1:(py-1))],
                   c(1, 3), sum) / 2
  fertpop * fp$asfr
}


#' HIV prevalence among pregant women by single age
#' 
get_pregprev <- function(mod, fp, hp1){

  ss <- fp$ss
  py <- fp$ss$PROJ_YEARS
  expand_idx <- rep(fp$ss$h.fert.idx, fp$ss$h.ag.span[fp$ss$h.fert.idx])
  hivn_w <- mod[ss$p.fert.idx, ss$f.idx, ss$hivn.idx,  ]
  hivp_w <- colSums(hp1$hivpop[ , ss$p.fert.idx, ss$f.idx, ] * fp$frr_cd4[ , expand_idx, ])
  art_w <- colSums(hp1$artpop[ , , ss$p.fert.idx, ss$f.idx, ] * fp$frr_art[ , , expand_idx, ],,2)
  denom_w <- hivn_w + hivp_w + art_w
  pregprev_a <- 1 - (hivn_w + hivn_w[ , c(1, 1:(py-1))]) / (denom_w + denom_w[ , c(1, 1:(py-1))])

  pregprev_a
}


get_daly_hivpop <- function(hivpop1){
  idx <- rep(c("pop_gt350", "pop200to350", "pop_lt200"), c(2, 2, 3))
  v <- apply(hivpop1, 2:4, fastmatch::ctapply, idx, sum)
  names(dimnames(v))[1] <- "cd4daly"
  v
}
