#===================================================================
# Description code: Set of functions to obtain confidence intervals
#      and differences in correlations
#===================================================================

# 1. boot_fun: function to obtain the confidence intervals 
# input: x = vector of correlations
#        conf = significance level
# output : ci = dataframe with mean and 99% confidence intervals

boot_fun <- function(x, conf = 0.99){
  set.seed(541983)
  bt <- boot::boot(x, \(x,i) mean(x[i]), R=1000)
  ci <- boot::boot.ci(bt, conf = conf, type="perc")
  data.frame(mean=ci$t0,
             lower=ci$percent[4],
             upper=ci$percent[5])
}

# 2. calcSD: function taken from  https://github.com/schrauf/AccuracyComparer/ 

# The function gives the confidence intervals and the statistical difference in predictive abilities
# for the models compared
# input: ci = dataframe with mean and 99% confidence intervals
# output : out= dataframe with mean, 99% confidence intervals and significance letters


calcSD <- function(ci, make_plot=TRUE) {
  # 1
  comps = ci$contrast
  labels <- do.call(rbind, strsplit(comps, " - "))
  lower.CL = ci$lower
  upper.CL = ci$upper
  # 3
  sdH1  = lower.CL > 0 | upper.CL < 0
  # 4
  H1 <- sdH1
  V <- unique(c(labels))
  E <- c(t(labels[!H1,]))
  isolates <- setdiff(V, unique(c(E)))
  sdGr <- igraph::graph(E, isolates = isolates, directed = F)
  if (make_plot) plot(sdGr)
  mc <- igraph::max_cliques(sdGr)
  out <- character(length(V))
  for(i in seq_along(mc)) {
    out <- paste0(out, ifelse(V %in% names(mc[[i]]), letters[i], ""))
  }
  return(data.frame(sd = out, modc = V))
}

# 3. ff: function to compute the difference in performance between methods

ff <- function(df, keys) {
  pc <- df |>
    (\(df) inner_join(df,df, by=c("rep")))() |>
    filter(modc.x < modc.y) |>
    mutate(diff = y.y-y.x, comp=paste(modc.y,modc.x,sep=" - ")) |>
    group_by(comp) |>
    summarise(boot_fun(diff))|>
    mutate(signif = (0 < lower) | (0 > upper))
  out <- calcSD(with(pc, data.frame(contrast=comp, lower, upper)),
                make_plot=FALSE)
  df |> 
    group_by(modc) |>
    summarise(boot_fun(y))|>
    left_join(out, by = c("modc"))
}
