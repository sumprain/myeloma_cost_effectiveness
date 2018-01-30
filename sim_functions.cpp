#include <Rcpp.h>
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
      return Regimen::POM_d;
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
  const double CONV_TO_MNTHS_PPS = 1000;
  const double CONV_TO_MNTHS_PFS = 200;
  const double PROB_DEATH_PFS = 0.022;
  const double MEAN_BSA = 1.91;         // kg/sq m
  const double SD_BSA = 0.03;
  const double MEAN_WT = 79;            // kg
  const double SD_WT = 1.7;
  const double UTIL_PFS_LL = 0.70;
  const double UTIL_PFS_UL = 0.76;
  const double UTIL_DELTA_PPS_LL = -0.084;
  const double UTIL_DELTA_PPS_UL = -0.025;
  const double UTIL_DELTA_AE_LL = -0.088;
  const double UTIL_DELTA_AE_UL = -0.009;

  double prob_anemia = 0;
  double prob_fatigue = 0;
  double prob_npenia = 0;
  double prob_tpenia = 0;

  //NumericVector mnths;
  List out(N);

  double shape_pfs = 0;
  double rate_pfs = 0;
  double rate_ttf = 0;
  double shape_pps = 2.94;
  double rate_pps = 1.0/0.02;

  Regimen first_line = NO;

  //Rcout << "POM" << Regimen::POM_d << std::endl;
  //Rcout << "CARF" << Regimen::CARF << std::endl;
  //Rcout << "NO" << Regimen::NO << std::endl;

  switch(arms) {
    case 0:             // pomd
      shape_pfs = 1.77;
      rate_pfs = 1.0/0.03;
      rate_ttf = (1.0 / 1.19) * rate_pfs;
      first_line = POM_d;
      prob_anemia = (0.22 * 113)/(14.2 * 113);
      prob_fatigue = 0.14 / 14.2;
      prob_npenia = 0.1 / 14.2;
      prob_tpenia = 0.19 / 14.2;
      break;
    case 1:  // dara
      //Rcout << "entered dara 1" << std::endl;
      shape_pfs = 1.77;
      rate_pfs = (1.0 / 0.95) * ( 1.0 / 0.03 );
      rate_ttf = (1.0 / 1.13) * rate_pfs;
      first_line = DARA;
      prob_anemia = 25.0 / (9.3 * 106);
      prob_fatigue = 3.0 / (9.3 * 106);
      prob_npenia = 13.0 / (9.3 * 106);
      prob_tpenia = 20.0 / (9.3 * 106);
      break;
    case 2:             // car
      //Rcout << "entered car 2" << std::endl;
      shape_pfs = 1.77;
      rate_pfs = (1.0 / 0.83) * ( 1.0 / 0.03 );
      rate_ttf = (1.0 / 1.03) * rate_pfs;
      first_line = CARF;
      prob_anemia = 63.0 / (266 * 15.6);
      prob_fatigue = 20.0 / (266 * 15.6);
      prob_npenia = 18.0 / (266 * 15.6);
      prob_tpenia = 77.0 / (266 * 15.6);

  }

  //Rcout << "first line " << first_line << std::endl;
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

  std::vector<int> cum_mnths;
  std::vector<int> state;
  std::vector<int> regimen;
  std::vector<int> ae_anemia;
  std::vector<int> ae_fatigue;
  std::vector<int> ae_npenia;
  std::vector<int> ae_tpenia;

  // simulate for each patient out of N patients
  for (R_xlen_t i = 0; i < N; ++i) {

    double util_pfs_one = util_pfs[i];
    double util_pps_one = util_pfs[i] + util_delta_pps[i];
    double util_ae_one = util_pfs[i] + util_delta_ae[i];

    int pd_pfs = std::ceil(CONV_TO_MNTHS_PFS * t_to_prog[i]);
    int pd_ttf = std::ceil(CONV_TO_MNTHS_PFS * t_to_ttf[i]);
    int pd_pps = std::ceil(CONV_TO_MNTHS_PPS * t_from_pps_to_death[i]);
    int tot_dur_for_sim = pd_pfs + pd_pps;

    bool is_dead = false;
    bool has_progressed = false;
    Regimen second_line = NO;

    for (R_xlen_t time = 1; time <= tot_dur_for_sim; ++time) {

      // check for progression, select second line therapy
      if ((time > pd_pfs) & (!has_progressed)) {
        second_line = select_second_line(first_line);
        has_progressed = true;
      }

      if (!has_progressed) {
        if (time <= pd_ttf) {
          regimen.push_back(first_line);
        } else {
          regimen.push_back(NO);
        }
      } else {
        regimen.push_back(second_line);
      }

      // check for aes
      if ((regimen.back() == POM_d) |
          (regimen.back() == DARA) |
          (regimen.back() == CARF)) {

        if (R::runif(0, 1) <= prob_anemia) {
          ae_anemia.push_back(1);
        } else {
          ae_anemia.push_back(0);
        }

        if (R::runif(0, 1) <= prob_fatigue) {
          ae_fatigue.push_back(1);
        } else {
          ae_fatigue.push_back(0);
        }

        if (R::runif(0, 1) <= prob_npenia) {
          ae_npenia.push_back(1);
        } else {
          ae_npenia.push_back(0);
        }

        if (R::runif(0, 1) <= prob_tpenia) {
          ae_tpenia.push_back(1);
        } else {
          ae_tpenia.push_back(0);
        }

      } else {
        ae_anemia.push_back(0);
        ae_fatigue.push_back(0);
        ae_npenia.push_back(0);
        ae_tpenia.push_back(0);
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
        state.push_back(DEATH);
      } else {
        if (has_progressed) {
          state.push_back(PPS);
        } else {
          if (time <= pd_ttf) {
            state.push_back(PFS_TPOS);
          } else {
            state.push_back(PFS_TNEG);
          }
        }
      }

      cum_mnths.push_back(time);

      if (is_dead) {
        break;
      }
    } // out of time loop

    out[i] = DataFrame::create(_["cum_month"] = cum_mnths,
                        _["state"] = state,
                        _["regimen"] = regimen,
                        _["ae_anemia"] = ae_anemia,
                        _["ae_fatigue"] = ae_fatigue,
                        _["ae_npenia"] = ae_npenia,
                        _["ae_tpenia"] = ae_tpenia);

    cum_mnths.clear();
    state.clear();
    regimen.clear();
    ae_anemia.clear();
    ae_fatigue.clear();
    ae_npenia.clear();
    ae_tpenia.clear();

  } // out of patient loop
  return out;

}

