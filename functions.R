# calculate PFS/TTF/death random sample ------

pfs <- function(arms = c("pomd", "dara", "car"), N = 1000) {

  arms <- match.arg(arms)

  conv <- 100

  switch(arms,
         pomd = {
          shape_pfs <- 1.77
          rate_pfs <- 1/0.03
          rate_ttf <- 1.19*rate_pfs
         },
         dara = {
           shape_pfs <- 1.77
           rate_pfs <- 0.95*(1/0.03)
           rate_ttf <- 1.13*rate_pfs
         },
         car = {
           shape_pfs <- 1.77
           rate_pfs <- 0.83*(1/0.03)
           rate_ttf <- 1.03*rate_pfs
         })

  t_to_prog <- conv * flexsurv::rgompertz(N,
                              shape = shape_pfs,
                              rate = rate_pfs)

  t_to_tf <- conv * flexsurv::rgompertz(N,
                              shape = shape_pfs,
                              rate = rate_ttf)

  is_tf <- ifelse((t_to_prog > t_to_tf), 1, 0)

  prob_death_pfs <- 0.022

  death_ind <- sample(c(1, 0),
                size = ceiling(t_to_prog),
                replace = TRUE,
                prob = c(prob_death_pfs, 1 - prob_death_pfs))

  death_month <- match(1, death_ind)

  is_death <- ifelse((death_month < t_to_prog), 1, 0)

  for (i in seq_len(N)) {

  }

}
