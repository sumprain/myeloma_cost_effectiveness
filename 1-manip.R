library(tidyverse)
library(flexsurv)

Rcpp::sourceCpp("sim_functions.cpp")

# load patients
set.seed(100)

n_pats <- 100
n_sims <- 500

sim_one <- function(n_sim, n_pats) {

  df <- map2_df(.x = c(0, 1, 2), .y = c("pomd", "dara", "carf"),
                .f = function(.x, .y) {
                  filter(mutate(bind_rows(time_lines(arms = .x, n_pats),
                  .id = "pat_id"), init_reg = .y), cum_month <= 36)
                })

  mutate(df, sim_no = rep(n_sim, nrow(df)))

}

sim_all <- function(n_sims, n_pats) {
  map_df(seq_len(n_sims), function(x, y) {
    sim_one(x, y)
  }, y = n_pats)
}

df_final <- sim_all(n_sims = n_sims, n_pats = n_pats)

os <- df_final %>% group_by(init_reg, sim_no, pat_id) %>%
  arrange(cum_month) %>%
  summarise(last_fup = max(cum_month),
            status = 1*(last(state) == 3)) %>%
  ungroup()

os_surv_plot <- os %>% group_by(sim_no, init_reg) %>% nest() %>%
  mutate(os_df = map(data, ~ survfit(Surv(last_fup, status == 1) ~ 1, data = .)))

median_surv <- os_surv_plot %>% mutate(median_surv = map_dbl(os_df, ~ approx(.$surv, .$time, xout = 0.5)$y)) %>% select(sim_no, init_reg, median_surv)

surv_time <- os_surv_plot %>% mutate(survs = map(os_df, ~ approx(.$time, .$surv, xout = c(3, 6, 12, 18, 24, 36))$y), ts = map(os_df, ~ c(3,6,12,18,24,36))) %>% select(sim_no, init_reg, survs, ts) %>% unnest()

state_pfs <- df_final %>% group_by(sim_no, init_reg, pat_id) %>% summarise(status_pfs = 1*(any(state == 2) | any(state == 3))) %>% ungroup()

f_time_pfs <- function(t, s) {
  l <- (s == 2) | (s == 3)
  if (sum(1*l) == 0) {
    res <- t[l][1]
  } else {
    res <- t[length(t)]
  }
  res
}

t_pfs <- df_final %>% group_by(sim_no, init_reg, pat_id) %>%
  summarise(time_pfs = f_time_pfs(cum_month, state)) %>% ungroup()

df_pfs <- left_join(state_pfs, t_pfs)

rm(t_pfs, state_pfs)

pfs_plot <- df_pfs %>% group_by(sim_no, init_reg) %>% nest() %>%
  mutate(pfs_df = map(data, ~ survfit(Surv(time_pfs, status_pfs == 1) ~ 1, data = .)))

median_pfs <- pfs_plot %>% mutate(median_pfs = map_dbl(pfs_df, ~ approx(.$surv, .$time, xout = 0.5)$y)) %>% select(sim_no, init_reg, median_pfs)

rm(f_time_pfs)

pfs_time <- pfs_plot %>% mutate(survs = map(pfs_df, ~ approx(.$time, .$surv, xout = c(3, 6, 12, 18, 24, 36))$y), ts = map(pfs_df, ~ c(3,6,12,18,24,36))) %>% select(sim_no, init_reg, survs, ts) %>% unnest()

df_pflm <- df_final %>% group_by(sim_no, init_reg, pat_id) %>% filter(state == 0 | state == 1) %>% summarise(pfm = last(cum_month)) %>% ungroup()

df_qalm <- df_final %>% group_by(sim_no, init_reg, pat_id) %>% summarise(qalm = sum(disc_utility)) %>% ungroup()

df_cost_reg <- df_final %>% group_by(sim_no, init_reg, pat_id) %>% summarise(tot_cost_reg = sum(disc_cost_regimen), tot_months = max(cum_month), cost_month = tot_cost_reg/tot_months) %>% ungroup()

df_cost_ae <- df_final %>% mutate(disc_cost_ae = disc_cost_anemia + disc_cost_npenia + disc_cost_tpenia) %>% group_by(sim_no, init_reg, pat_id) %>% summarise(tot_cost_ae = sum(disc_cost_ae), tot_months = max(cum_month), cost_month = tot_cost_ae/tot_months)

df_cost_ha <- df_final %>% group_by(sim_no, init_reg, pat_id) %>% summarise(tot_cost = sum(disc_cost_health_adm), tot_months = max(cum_month), cost_month = tot_cost/tot_months)

df_cost_qalm <- df_final %>% mutate(disc_cost = disc_cost_regimen + disc_cost_anemia + disc_cost_npenia + disc_cost_tpenia) %>% group_by(sim_no, init_reg, pat_id) %>% summarise(tot_cost = sum(disc_cost), tot_utility = sum(disc_utility), cost_qalm = tot_cost/tot_utility) %>% ungroup()

save(list = ls(), file = "final_before_slides.Rdata")

# create animated gif for timeline

library(gganimate)

samp <- df_final %>%
  filter(sim_no == 2, init_reg == "pomd", pat_id %in% c(1,2,3,4,9,8)) %>%
  select(pat_id, cum_month, state) %>% mutate(state = factor(state, c(0,1,2,3)))

p <- ggplot(data = samp, aes(x = cum_month, y = pat_id,
      colour = state, frame = cum_month, cumulative = TRUE)) +
  geom_point(size = 6) + theme_bw() + labs(x = "Months", y = "Patient Number", colour = "State") + scale_color_manual(values = c("black", "orange", "blue", "red"), labels = c("PFS T+", "PFS T-", "PPS", "Death"), drop = FALSE)

old_opts <- animation::ani.options(interval = 0.5)
gganimate(p, filename = "anim_time_line.gif")
animation::ani.options(old_opts)
