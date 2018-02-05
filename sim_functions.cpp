#include <Rcpp.h>
#include <map>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

enum State {
  PFS_TPOS = 0,
  PFS_TNEG = 1,
  PPS = 2,
  DEATH = 3
};

enum Regimen {
  POM_d = 0,
  DARA = 1,
  CARF = 2,
  LEN_d = 3,
  BORT = 4,
  THAL_d = 5,
  PanVD = 6,
  NO = 7
};

enum Util_state {
  PFS_U = 0,
  PPS_U = 1,
  AE_U = 2
};

class Rng {
private:
  int m_min, m_max;
public:
  Rng(int min, int max): m_min{min}, m_max{max} {
    if (min > max) {
      Rcpp::stop("min cannot be more than max");
    }
  }

  const bool contains(int num) const {
    if ((num >= m_min) & (num <= m_max)) {
      return true;
    } else {
      return false;
    }
  }
};


class Entity {
private:

  std::vector<std::vector<double> > m_prices;
  std::vector<Rng> m_strata;

  int stratum_idx(int stratum) const {

    int idx = -1;

    for (int i = 0; i < m_strata.size(); ++i) {
      if (m_strata[i].contains(stratum)) {
        idx = i;
        break;
      }
    }
    return idx;
  }

public:

  Entity(std::vector<double> bp_strata_unit,
         std::vector<double> bp_strata,
         std::vector<Rng> strata,
         double bp, double unit, int t_sim,
         double infl_ll, double infl_ul): m_strata{strata} {

    NumericVector infl_vec = Rcpp::runif(t_sim - 1, infl_ll, infl_ul);
    infl_vec.push_front(0);

    auto iter_s = strata.cbegin();
    auto iter_bp_unit = bp_strata_unit.cbegin();
    auto iter_bp = bp_strata.cbegin();

    while (iter_s != strata.cend()) {
      std::vector<double> price_vector;
      double price = ((*iter_bp_unit) * unit) + (*iter_bp) + bp;
      for (R_xlen_t j = 0; j < t_sim; ++j) {
        price = price * (1 + infl_vec[j]);
        price_vector.push_back(price);
      }

      m_prices.push_back(price_vector);

      ++iter_s;
      ++iter_bp_unit;
      ++iter_bp;
    }         // end while
  }

  double cur_price(int stratum, int pres_t) const {

    int idx_stratum = stratum_idx(stratum);
    if (idx_stratum == -1) {
      return 0;
    } else {
      return (m_prices[idx_stratum])[pres_t - 1];
    }
  }

  double discounted_price(int stratum, int pres_t, double r_discount) const {
    double pres_price = cur_price(stratum, pres_t);
    return (pres_price / std::pow((1 + r_discount), (pres_t - 1)));
  }

};

void chng_reg_ctr(const std::vector<int>& regimen,
                  const Regimen& present, int* regimen_ctr) {

  //Rcout << "entered chng_reg_ctr" << std::endl;
  //Rcout << "regimen_ctr: " << *regimen_ctr << std::endl;
  //Rcout << "prev regimen: " << regimen.back() << std::endl;
  //Rcout << "present regimen: " << present << std::endl;

  if (regimen.back() == int(present)) {
    ++(*regimen_ctr);
  } else {
    (*regimen_ctr) = 0;
  }
}

Regimen select_second_line (Regimen arms) {

  double reg_ind = R::runif(0, 1);

  if (arms == POM_d) {
    if (reg_ind <= 0.075) {
      return DARA;
    } else if (reg_ind <= 0.459) {
      return CARF;
    } else if (reg_ind <= 0.681) {
      return LEN_d;
    } else if (reg_ind <= 0.859) {
      return BORT;
    } else if (reg_ind <= 0.925) {
      return THAL_d;
    } else {
      return PanVD;
    }
  } else if (arms == DARA) {
    if (reg_ind <= 0.402) {
      return POM_d;
    } else if (reg_ind <= 0.65) {
      return CARF;
    } else if (reg_ind <= 0.794) {
      return LEN_d;
    } else if (reg_ind <= 0.909) {
      return BORT;
    } else if (reg_ind <= 0.952) {
      return THAL_d;
    } else {
      return PanVD;
    }
  } else if (arms == CARF) {
    if (reg_ind <= 0.502) {
      return POM_d;
    } else if (reg_ind <= 0.562) {
      return DARA;
    } else if (reg_ind <= 0.742) {
      return LEN_d;
    } else if (reg_ind <= 0.886) {
      return BORT;
    } else if (reg_ind <= 0.94) {
      return THAL_d;
    } else {
      return PanVD;
    }
  }

}

// [[Rcpp::export]]
List time_lines (int arms, int N) {

  Function rgomp ("rgompertz");
  const double CONV_TO_MNTHS_PPS = 1000.0;
  const double CONV_TO_MNTHS_PFS = 200.0;
  const double PROB_DEATH_PFS = 0.022;
  const double MEAN_BSA = 1.91;         // kg/sq m
  const double SD_BSA = 0.03;
  const double MEAN_WT = 79.0;            // kg
  const double SD_WT = 1.7;
  const double UTIL_PFS_LL = 0.70;
  const double UTIL_PFS_UL = 0.76;
  const double UTIL_DELTA_PPS_LL = -0.084;
  const double UTIL_DELTA_PPS_UL = -0.025;
  const double UTIL_DELTA_AE_LL = -0.088;
  const double UTIL_DELTA_AE_UL = -0.009;
  const int MAX_IDX = 100000;
  const double INFLATION_LL = 0.035/12.0;
  const double INFLATION_UL = 0.06/12.0;
  const double NO_UNIT = 1.0;
  const double DISCOUNT_RATE = 0.03/12.0;

  double prob_anemia = 0.0;
  double prob_fatigue = 0.0;
  double prob_npenia = 0.0;
  double prob_tpenia = 0.0;

  //NumericVector mnths;
  List out;

  double shape_pfs = 0.0;
  double rate_pfs = 0.0;
  double rate_ttf = 0.0;
  double shape_pps = 2.94;
  double rate_pps = 1.0/0.02;

  Regimen first_line = NO;

  switch(arms) {
    case 0:             // pomd
      shape_pfs = 1.77;
      rate_pfs = 1.0/0.03;
      rate_ttf = (1.0 / 1.19) * rate_pfs;
      first_line = POM_d;
      prob_anemia = (0.22 * 113)/(14.2 * 113);
      //prob_fatigue = 0.14 / 14.2;
      prob_npenia = 0.1 / 14.2;
      prob_tpenia = 0.19 / 14.2;
      break;
    case 1:  // dara
      shape_pfs = 1.77;
      rate_pfs = (1.0 / 0.95) * ( 1.0 / 0.03 );
      rate_ttf = (1.0 / 1.13) * rate_pfs;
      first_line = DARA;
      prob_anemia = 25.0 / (9.3 * 106.0);
      //prob_fatigue = 3.0 / (9.3 * 106.0);
      prob_npenia = 13.0 / (9.3 * 106.0);
      prob_tpenia = 20.0 / (9.3 * 106.0);
      break;
    case 2:             // car
      shape_pfs = 1.77;
      rate_pfs = (1.0 / 0.83) * ( 1.0 / 0.03 );
      rate_ttf = (1.0 / 1.03) * rate_pfs;
      first_line = CARF;
      prob_anemia = 63.0 / (266.0 * 15.6);
      //prob_fatigue = 20.0 / (266.0 * 15.6);
      prob_npenia = 18.0 / (266.0 * 15.6);
      prob_tpenia = 77.0 / (266.0 * 15.6);

  }

  NumericVector t_to_prog = rgomp(N,
                      _["shape"] = shape_pfs,
                      _["rate"] = rate_pfs);

  NumericVector t_to_ttf = rgomp(N,
                      _["shape"] = shape_pfs,
                      _["rate"] = rate_ttf);

  NumericVector t_from_pps_to_death = rgomp(N,
                      _["shape"] = shape_pps,
                      _["rate"] = rate_pps);

  NumericVector bsa = Rcpp::rnorm(N, MEAN_BSA, SD_BSA);
  NumericVector wt = Rcpp::rnorm(N, MEAN_WT, SD_WT);

  NumericVector util_pfs = Rcpp::runif(N, UTIL_PFS_LL, UTIL_PFS_UL);
  NumericVector util_delta_pps = Rcpp::runif(N,
                          UTIL_DELTA_PPS_LL, UTIL_DELTA_PPS_UL);
  NumericVector util_delta_ae = Rcpp::runif(N,
                          UTIL_DELTA_AE_LL, UTIL_DELTA_AE_UL);

  std::vector<int> cum_mnths, state, regimen, ae_anemia,    //ae_fatigue
                    ae_npenia, ae_tpenia;
  std::vector<double> cost_regimen, disc_cost_regimen,
                      cost_anemia, disc_cost_anemia,
                      //cost_fatigue, disc_cost_fatigue,
                      cost_npenia, disc_cost_npenia,
                      cost_tpenia, disc_cost_tpenia,
                      utility, disc_utility,
                      cost_health_adm, disc_cost_health_adm;

  Entity* pom_d_reg;
  Entity* dara_reg;
  Entity* carf_reg;
  Entity* len_d_reg;
  Entity* bort_reg;
  Entity* thal_d_reg;
  Entity* panvd_reg;
  Entity* util;
  Entity* health_adm_cost;
  Entity* anemia_cost;
  //Entity* fatigue_cost;
  Entity* npenia_cost;
  Entity* tpenia_cost;

  // unit = 1
  std::vector<double> bp_str_unit_pomd {0};
  std::vector<double> bp_str_pomd {20080};
  std::vector<Rng> str_pomd {Rng {0, MAX_IDX}};

  // unit = kg
  std::vector<double> bp_str_unit_dara {2400, 1200, 600};
  std::vector<double> bp_str_dara {4200, 2100, 1050};
  std::vector<Rng> str_dara {Rng {0, 2}, Rng {3, 6}, Rng {7, MAX_IDX}};

  // unit = sqm
  std::vector<double> bp_str_unit_carf {24675.3, 26973, 17982};
  std::vector<double> bp_str_carf {0, 0, 0};
  std::vector<Rng> str_carf {Rng {0, 1}, Rng {2, 12}, Rng {13, MAX_IDX}};

  // unit = 1
  std::vector<double> bp_str_unit_lend {0};
  std::vector<double> bp_str_lend {15788};
  std::vector<Rng> str_lend {Rng {0, MAX_IDX}};

  // unit = 1
  std::vector<double> bp_str_unit_bor {0};
  std::vector<double> bp_str_bor {45200};
  std::vector<Rng> str_bor {Rng {0, MAX_IDX}};

  // unit = 1
  std::vector<double> bp_str_unit_thald {0};
  std::vector<double> bp_str_thald {7568};
  std::vector<Rng> str_thald {Rng {0, MAX_IDX}};

  // unit = 1
  std::vector<double> bp_str_unit_panvd {0, 0};
  std::vector<double> bp_str_panvd {110280, 87640};
  std::vector<Rng> str_panvd {Rng {0, 8}, Rng {9, MAX_IDX}};

  const double cost_vte = 1459.2;
  const double cost_herpes = 360;

  // unit = 1
  std::vector<double> bp_str_unit_health_cost {0, 0, 0, 0};
  std::vector<double> bp_str_health_cost {9566.67, 4783.33, 12316.67, 0};
  std::vector<Rng> str_health_cost {Rng {0, 0}, Rng {1, 1},
                                  Rng {2, 2}, Rng {3, 3}};

  // unit = 1
  std::vector<double> bp_str_unit_anemia {0};
  std::vector<double> bp_str_anemia {2100};
  std::vector<Rng> str_anemia {Rng {0, MAX_IDX}};

  // unit = 1
  //std::vector<double> bp_str_unit_fatigue {0};
  //std::vector<double> bp_str_fatigue {500};
  //std::vector<Rng> str_fatigue {Rng {0, MAX_IDX}};

  // unit = 1
  std::vector<double> bp_str_unit_npenia {0};
  std::vector<double> bp_str_npenia {5000};
  std::vector<Rng> str_npenia {Rng {0, MAX_IDX}};

  // unit = 1
  std::vector<double> bp_str_unit_tpenia {0};
  std::vector<double> bp_str_tpenia {1800};
  std::vector<Rng> str_tpenia {Rng {0, MAX_IDX}};

  Rcout << "Entering patient loop" << std::endl;
  // simulate for each patient out of N patients
  for (R_xlen_t i = 0; i < N; ++i) {

    double util_pfs_1 = util_pfs[i];
    double util_pps_1 = util_pfs[i] + util_delta_pps[i];
    double util_ae_1 = util_pfs[i] + util_delta_ae[i];

    Rcout << "Util: " << util_pfs_1 << std::endl;
    Rcout << "Util pps: " << util_pps_1 << std::endl;
    Rcout << "Util ae: " << util_ae_1 << std::endl;
    // unit = 1
    std::vector<double> bp_str_unit_util {0, 0, 0};
    std::vector<double> bp_str_util {util_pfs_1, util_pps_1, util_ae_1};
    std::vector<Rng> str_util {Rng {0, 0}, Rng {1, 1}, Rng {2, 2}};

    int pd_pfs = std::ceil(CONV_TO_MNTHS_PFS * t_to_prog[i]);
    int pd_ttf = std::ceil(CONV_TO_MNTHS_PFS * t_to_ttf[i]);
    int pd_pps = std::ceil(CONV_TO_MNTHS_PPS * t_from_pps_to_death[i]);
    int tot_dur_for_sim = pd_pfs + pd_pps;

    util = new Entity {bp_str_unit_util, bp_str_util,
                       str_util, 0, NO_UNIT, tot_dur_for_sim,
                       0, 0};

    pom_d_reg = new Entity {bp_str_unit_pomd, bp_str_pomd,
                            str_pomd, cost_vte, NO_UNIT,
                            tot_dur_for_sim, INFLATION_LL,
                            INFLATION_UL};

    dara_reg = new Entity {bp_str_unit_dara, bp_str_dara,
                           str_dara, cost_herpes, wt[i],
                           tot_dur_for_sim, INFLATION_LL,
                           INFLATION_UL};

    carf_reg = new Entity {bp_str_unit_carf, bp_str_carf,
                           str_carf, cost_herpes, bsa[i],
                           tot_dur_for_sim, INFLATION_LL,
                           INFLATION_UL};

    len_d_reg = new Entity {bp_str_unit_lend, bp_str_lend,
                            str_lend, cost_vte, NO_UNIT,
                            tot_dur_for_sim, INFLATION_LL,
                            INFLATION_UL};

    bort_reg = new Entity {bp_str_unit_bor, bp_str_bor,
                           str_bor, cost_herpes, NO_UNIT,
                           tot_dur_for_sim, INFLATION_LL,
                           INFLATION_UL};

    thal_d_reg = new Entity {bp_str_unit_thald, bp_str_thald,
                             str_thald, cost_vte, NO_UNIT,
                             tot_dur_for_sim, INFLATION_LL,
                             INFLATION_UL};

    panvd_reg = new Entity {bp_str_unit_panvd, bp_str_panvd,
                            str_panvd, cost_herpes, NO_UNIT,
                            tot_dur_for_sim, INFLATION_LL,
                            INFLATION_UL};

    health_adm_cost = new Entity {bp_str_unit_health_cost,
                                  bp_str_health_cost,
                                  str_health_cost, 0, NO_UNIT,
                                  tot_dur_for_sim, INFLATION_LL,
                                  INFLATION_UL};

    anemia_cost = new Entity {bp_str_unit_anemia, bp_str_anemia,
                              str_anemia, 0, NO_UNIT, tot_dur_for_sim,
                              INFLATION_LL, INFLATION_UL};

    //fatigue_cost = new Entity {bp_str_unit_fatigue, bp_str_fatigue,
    //                           str_fatigue, 0, NO_UNIT, tot_dur_for_sim,
    //                           INFLATION_LL, INFLATION_UL};

    npenia_cost = new Entity {bp_str_unit_npenia, bp_str_npenia,
                              str_npenia, 0, NO_UNIT, tot_dur_for_sim,
                              INFLATION_LL, INFLATION_UL};

    tpenia_cost = new Entity {bp_str_unit_tpenia, bp_str_tpenia,
                              str_tpenia, 0, NO_UNIT, tot_dur_for_sim,
                              INFLATION_LL, INFLATION_UL};

    bool is_dead = false;
    bool has_progressed = false;
    Regimen second_line = NO;
    int regimen_ctr = 0;
    bool is_ae = false;

    Rcout << "entering time line loop" << std::endl;

    for (R_xlen_t time = 1; time <= tot_dur_for_sim; ++time) {

      // check for progression, select second line therapy
      if ((time > pd_pfs) & (!has_progressed)) {
        second_line = select_second_line(first_line);
        has_progressed = true;
      }

      if (!has_progressed) {
        if (time <= pd_ttf) {
          // if repeated entry of same regimen regimen_ctr will be +1 else set 0.
          if (time == 1) {
            regimen.push_back(first_line);
          } else {
            chng_reg_ctr(regimen, first_line, &regimen_ctr);
            regimen.push_back(first_line);
          }
          //Rcout << "reg ctr: " << regimen_ctr << std::endl;
        } else {
          chng_reg_ctr(regimen, NO, &regimen_ctr);
          regimen.push_back(NO);
        }
      } else {
        chng_reg_ctr(regimen, second_line, &regimen_ctr);
        regimen.push_back(second_line);
      }
      //Rcout << "reached regimen costing" << std::endl;
      switch (regimen.back()) {
        case int(POM_d):
          cost_regimen.push_back(pom_d_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(pom_d_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(DARA):
          cost_regimen.push_back(dara_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(dara_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(CARF):
          cost_regimen.push_back(carf_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(carf_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(LEN_d):
          cost_regimen.push_back(len_d_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(len_d_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(BORT):
          cost_regimen.push_back(bort_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(bort_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(THAL_d):
          cost_regimen.push_back(thal_d_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(thal_d_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(PanVD):
          cost_regimen.push_back(panvd_reg->cur_price(regimen_ctr, time));
          disc_cost_regimen.push_back(panvd_reg->discounted_price(regimen_ctr, time, DISCOUNT_RATE));
          break;
        case int(NO):
          cost_regimen.push_back(0);
          disc_cost_regimen.push_back(0);
        }

      // check for aes
      if ((regimen.back() == POM_d) |
          (regimen.back() == DARA) |
          (regimen.back() == CARF)) {

        if (R::runif(0, 1) <= prob_anemia) {
          is_ae = true;
          ae_anemia.push_back(1);
          cost_anemia.push_back(anemia_cost->cur_price(1, time));
          disc_cost_anemia.push_back(anemia_cost->discounted_price(1, time, DISCOUNT_RATE));
        } else {
          ae_anemia.push_back(0);
          cost_anemia.push_back(0);
          disc_cost_anemia.push_back(0);
        }

        // if (R::runif(0, 1) <= prob_fatigue) {
        //   is_ae = true;
        //   ae_fatigue.push_back(1);
        //   cost_fatigue.push_back(fatigue_cost->cur_price(1, time));
        //   disc_cost_fatigue.push_back(fatigue_cost->discounted_price(1, time, DISCOUNT_RATE));
        // } else {
        //   ae_fatigue.push_back(0);
        //   cost_fatigue.push_back(0);
        //   disc_cost_fatigue.push_back(0);
        // }

        if (R::runif(0, 1) <= prob_npenia) {
          is_ae = true;
          ae_npenia.push_back(1);
          cost_npenia.push_back(npenia_cost->cur_price(1, time));
          disc_cost_npenia.push_back(npenia_cost->discounted_price(1, time, DISCOUNT_RATE));
        } else {
          ae_npenia.push_back(0);
          cost_npenia.push_back(0);
          disc_cost_npenia.push_back(0);
        }

        if (R::runif(0, 1) <= prob_tpenia) {
          is_ae = true;
          ae_tpenia.push_back(1);
          cost_tpenia.push_back(tpenia_cost->cur_price(1, time));
          disc_cost_tpenia.push_back(tpenia_cost->discounted_price(1, time, DISCOUNT_RATE));
        } else {
          ae_tpenia.push_back(0);
          cost_tpenia.push_back(0);
          disc_cost_tpenia.push_back(0);
        }

      } else {
        ae_anemia.push_back(0);
        //ae_fatigue.push_back(0);
        ae_npenia.push_back(0);
        ae_tpenia.push_back(0);
        cost_anemia.push_back(0);
        disc_cost_anemia.push_back(0);
        //cost_fatigue.push_back(0);
        //disc_cost_fatigue.push_back(0);
        cost_npenia.push_back(0);
        disc_cost_npenia.push_back(0);
        cost_tpenia.push_back(0);
        disc_cost_tpenia.push_back(0);
      }

      // check for death
      if (time <= pd_pfs) {
        is_dead = (R::runif(0, 1) <= PROB_DEATH_PFS);
      } else {
        if (time == tot_dur_for_sim) {
          is_dead = true;
        }
      }

    // assign a state to the patient at the particular time
      if (is_dead) {
        if (time == 1) {
          cost_health_adm.push_back(health_adm_cost->cur_price(PFS_TPOS, time));
          disc_cost_health_adm.push_back(health_adm_cost->discounted_price(PFS_TPOS, time, DISCOUNT_RATE));
          state.push_back(DEATH);
        } else {
          cost_health_adm.push_back(health_adm_cost->cur_price(state.back(), time));
          disc_cost_health_adm.push_back(health_adm_cost->discounted_price(state.back(), time, DISCOUNT_RATE));
          state.push_back(DEATH);
        }
      } else {
        if (has_progressed) {
          state.push_back(PPS);
          cost_health_adm.push_back(health_adm_cost->cur_price(PPS, time));
          disc_cost_health_adm.push_back(health_adm_cost->discounted_price(PPS, time, DISCOUNT_RATE));
        } else {
          if (time <= pd_ttf) {
            state.push_back(PFS_TPOS);
            cost_health_adm.push_back(health_adm_cost->cur_price(PFS_TPOS, time));
            disc_cost_health_adm.push_back(health_adm_cost->discounted_price(PFS_TPOS, time, DISCOUNT_RATE));
          } else {
            state.push_back(PFS_TNEG);
            cost_health_adm.push_back(health_adm_cost->cur_price(PFS_TNEG, time));
            disc_cost_health_adm.push_back(health_adm_cost->discounted_price(PFS_TNEG, time, DISCOUNT_RATE));
          }
        }
      }

      // fill up utility vector
      //Rcout << "status before util: " << regimen.back() << std::endl;
      //Rcout << "is_dead: " << is_dead << std::endl;
      //Rcout << "is_ae: " << is_ae << std::endl;

      if (is_dead) {
        utility.push_back(0);
        disc_utility.push_back(0);
      } else {
        if (is_ae) {
          utility.push_back(util->cur_price(AE_U, time));
          disc_utility.push_back(util->discounted_price(AE_U, time, DISCOUNT_RATE));
        } else {
          if (state.back() == PPS) {
            utility.push_back(util->cur_price(PPS_U, time));
            disc_utility.push_back(util->discounted_price(PPS_U, time, DISCOUNT_RATE));
          } else if ((state.back() == PFS_TNEG) | (state.back() == PFS_TPOS)) {
            utility.push_back(util->cur_price(PFS_U, time));
            disc_utility.push_back(util->discounted_price(PFS_U, time, DISCOUNT_RATE));
          }
        }
      }

      cum_mnths.push_back(time);

      is_ae = false;

      if (is_dead) {
        Rcout << "tot time for simulation: " << tot_dur_for_sim << std::endl;
        Rcout << "time of death: " << time << std::endl;
        break;
      }
    } // out of time loop

    Rcout << "out of time loop" << std::endl;

    Rcout << "cum mnths: " << cum_mnths.size() << std::endl;
    Rcout << "state: " << state.size() << std::endl;
    Rcout << "regimen: " << regimen.size() << std::endl;
    Rcout << "ae_anemia: " << ae_anemia.size() << std::endl;
    Rcout << "ae_npenia: " << ae_npenia.size() << std::endl;
    Rcout << "ae_tpenia: " << ae_tpenia.size() << std::endl;
    Rcout << "cost reg: " << cost_regimen.size() << std::endl;
    Rcout << "disc cost reg: " << disc_cost_regimen.size() << std::endl;
    Rcout << "cost anemia: " << cost_anemia.size() << std::endl;
    Rcout << "disc anemia: " << disc_cost_anemia.size() << std::endl;
    Rcout << "cost npenia: " << cost_npenia.size() << std::endl;
    Rcout << "disc npenia: " << disc_cost_npenia.size() << std::endl;
    Rcout << "cost tpenia: " << cost_tpenia.size() << std::endl;
    Rcout << "disc cost tpenia: " << disc_cost_tpenia.size() << std::endl;
    Rcout << "cost health: " << cost_health_adm.size() << std::endl;
    Rcout << "disc health: " << disc_cost_health_adm.size() << std::endl;
    Rcout << "utility: " << utility.size() << std::endl;
    Rcout << "disc utility: " << disc_utility.size() << std::endl;

    out.push_back(DataFrame::create(_["cum_month"] = cum_mnths,
                        _["state"] = state,
                        _["regimen"] = regimen,
                        _["ae_anemia"] = ae_anemia,
                        //_["ae_fatigue"] = ae_fatigue,
                        _["ae_npenia"] = ae_npenia,
                        _["ae_tpenia"] = ae_tpenia,
                        _["cost_regimen"] = cost_regimen,
                        _["disc_cost_regimen"] = disc_cost_regimen,
                        _["cost_anemia"] = cost_anemia,
                        _["disc_cost_anemia"] = disc_cost_anemia,
                        //_["cost_fatigue"] = cost_fatigue,
                        //_["disc_cost_fatigue"] = disc_cost_fatigue,
                        _["cost_npenia"] = cost_npenia,
                        _["disc_cost_npenia"] = disc_cost_npenia,
                        _["cost_tpenia"] = cost_tpenia,
                        _["disc_cost_tpenia"] = disc_cost_tpenia,
                        _["cost_health_adm"] = cost_health_adm,
                        _["disc_cost_health_adm"] = disc_cost_health_adm,
                        _["utility"] = utility,
                        _["disc_utility"] = disc_utility));

    cum_mnths.clear();
    state.clear();
    regimen.clear();
    ae_anemia.clear();
    //ae_fatigue.clear();
    ae_npenia.clear();
    ae_tpenia.clear();
    cost_regimen.clear();
    disc_cost_regimen.clear();
    cost_anemia.clear();
    disc_cost_anemia.clear();
    //cost_fatigue.clear();
    //disc_cost_fatigue.clear();
    cost_npenia.clear();
    disc_cost_npenia.clear();
    cost_tpenia.clear();
    disc_cost_tpenia.clear();
    utility.clear();
    disc_utility.clear();
    cost_health_adm.clear();
    disc_cost_health_adm.clear();

    delete pom_d_reg;
    delete dara_reg;
    delete carf_reg;
    delete len_d_reg;
    delete bort_reg;
    delete thal_d_reg;
    delete panvd_reg;
    delete util;
    delete health_adm_cost;
    delete anemia_cost;
    //delete fatigue_cost;
    delete npenia_cost;
    delete tpenia_cost;
  } // out of patient loop

  Rcout << "out of patient loop" << std::endl;

  pom_d_reg = nullptr;
  dara_reg = nullptr;
  carf_reg = nullptr;
  len_d_reg = nullptr;
  bort_reg = nullptr;
  thal_d_reg = nullptr;
  panvd_reg = nullptr;
  util = nullptr;
  health_adm_cost = nullptr;
  anemia_cost = nullptr;
  //fatigue_cost = nullptr;
  npenia_cost = nullptr;
  tpenia_cost = nullptr;

  return out;

}
