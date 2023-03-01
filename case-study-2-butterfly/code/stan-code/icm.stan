// icm.stan: stan code for running the integrated community model for
//           the butterfly case study.
// Author: Jeffrey W. Doser and Wendy Leuenberger
data {

	// Constants 
	int<lower=1> n_counties;	// number of counties 
	int<lower=1> n_counties_long;	// number of counties * species
	int<lower=1> n_years;		// number of years
	int<lower=1> n_years_long;	// number of years * species
	int<lower=1> n_sites;		// number of unique summer survey locations
	int<lower=1> n_sites_long;	// number of unique summer survey locations * species
	int<lower=1> n_surveys;		// number of surveys (aka number of data points for each species)
        int<lower=1> n_species;         // number of species
        int<lower=1> n_surveys_long;    // number of species * number of surveys

	// Indices
	int<lower=1,upper=n_years> year_id[n_surveys];
	int<lower=1,upper=n_counties> county_id[n_surveys];
	int<lower=1,upper=n_sites> site_id[n_surveys];

	// Covariates
        vector[n_surveys] ia_ind;
        vector[n_surveys] il_ind;
        vector[n_surveys] mi_ind;
        vector[n_surveys] oh_ind;
        vector[n_surveys] year_cov;
	vector[n_surveys] week_cov;
	vector[n_surveys] week2_cov;
        vector[n_surveys] survey_cov;
        
	// Butterfly counts
	int<lower=0> y[n_surveys_long];
}

parameters{

        // Random effects on raw scale (mean 0 and single variance for all species)
        vector[n_counties_long] betaStarCountyRaw;
        vector[n_years_long] betaStarWeekRaw;
        vector[n_years_long] betaStarWeek2Raw;
        vector[n_sites_long] betaStarSiteRaw;
        vector[n_years_long] betaStarYearRaw;

	// Species-specific effects
        vector[n_species] betaInt;
        vector[n_species] betaIA;
        vector[n_species] betaIL;
        vector[n_species] betaMI;
        vector[n_species] betaOH;
        vector[n_species] betaYear;
        vector[n_species] betaWeek;
        vector[n_species] betaWeek2;
        vector[n_species] betaSurvey;

        // Community-level means
        real muBetaInt;
        real muBetaIA;
        real muBetaIL;
        real muBetaMI;
        real muBetaOH;
        real muBetaYear;
        real muBetaWeek;
        real muBetaWeek2;
        real muBetaSurvey;

        // Community-level variances
        real<lower=0> tauSqInt;
        real<lower=0> tauSqIA;
        real<lower=0> tauSqIL;
        real<lower=0> tauSqMI;
        real<lower=0> tauSqOH;
        real<lower=0> tauSqYear;
        real<lower=0> tauSqWeek;
        real<lower=0> tauSqWeek2;
        real<lower=0> tauSqSurvey;

	// Variances for random effects
	real<lower=0> sigma2_county;
	real<lower=0> sigma2_week;
	real<lower=0> sigma2_week2; 
	real<lower=0> sigma2_site;
	real<lower=0> sigma2_year;

	// Stuff for Gamma-Poisson mixture (NB1)
	vector<lower=0>[n_surveys_long] rho;	// gamma variates
	vector<lower=0>[n_species] r_count;	// shape/rate parameter for Gamma.
}

transformed parameters {

	// Expected count on log-scale. 
	vector<lower=0>[n_surveys_long] mu;
	
	// mu * rho 
	vector<lower=0>[n_surveys_long] mu_star;
        
        vector[n_counties_long] betaStarCounty;
        vector[n_years_long] betaStarWeek;
        vector[n_years_long] betaStarWeek2;
        vector[n_sites_long] betaStarSite;
        vector[n_years_long] betaStarYear;

        betaStarCounty = betaStarCountyRaw * sqrt(sigma2_county);
        betaStarWeek = betaStarWeekRaw * sqrt(sigma2_week);
        betaStarWeek2 = betaStarWeek2Raw * sqrt(sigma2_week2);
        betaStarSite = betaStarSiteRaw * sqrt(sigma2_site);
        betaStarYear = betaStarYearRaw * sqrt(sigma2_year);

        // Model
        for (i in 1:n_species) {
          for (j in 1:n_surveys) {
            mu[(i - 1) * n_surveys + j] = exp(betaInt[i] + betaIA[i] * ia_ind[j] + betaIL[i] * il_ind[j] +
                                              betaMI[i] * mi_ind[j] + betaOH[i] * oh_ind[j] + 
                                              betaYear[i] * year_cov[j] + 
                                              betaWeek[i] * week_cov[j] + 
                                              betaStarWeek[(i - 1) * n_years + year_id[j]] * week_cov[j] +
                                              betaWeek2[i] * week2_cov[j] + 
                                              betaStarWeek2[(i - 1) * n_years + year_id[j]] * week2_cov[j] + 
                                              betaSurvey[i] * survey_cov[j] + 
                                              betaStarCounty[(i - 1) * n_counties + county_id[j]] + 
                                              betaStarSite[(i - 1) * n_sites + site_id[j]] + 
                                              betaStarYear[(i - 1) * n_years + year_id[j]]);
            mu_star[(i - 1) * n_surveys + j] = mu[(i - 1) * n_surveys + j] * rho[(i - 1) * n_surveys + j];
          } // j (survey)
        } // i (species)
}

model {
	// IG priors on random effect variances and community-level variances
        target += inv_gamma_lpdf(sigma2_county | 0.1, 0.1);
        target += inv_gamma_lpdf(sigma2_week | 0.1, 0.1);
        target += inv_gamma_lpdf(sigma2_week2 | 0.1, 0.1);
        target += inv_gamma_lpdf(sigma2_site | 0.1, 0.1);
        target += inv_gamma_lpdf(sigma2_year | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqInt | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqIA | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqIL | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqMI | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqOH | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqYear | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqWeek | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqWeek2 | 0.1, 0.1);
        target += inv_gamma_lpdf(tauSqSurvey | 0.1, 0.1);

        // Community-level normal priors
        target += normal_lpdf(muBetaInt | 0, sqrt(100));
        target += normal_lpdf(muBetaIA | 0, sqrt(100));
        target += normal_lpdf(muBetaIL | 0, sqrt(100));
        target += normal_lpdf(muBetaMI | 0, sqrt(100));
        target += normal_lpdf(muBetaOH | 0, sqrt(100));
        target += normal_lpdf(muBetaYear | 0, sqrt(100));
        target += normal_lpdf(muBetaWeek | 0, sqrt(100));
        target += normal_lpdf(muBetaWeek2 | 0, sqrt(100));
        target += normal_lpdf(muBetaSurvey | 0, sqrt(100));

        // Species-level random effects        
        target += normal_lpdf(betaInt | muBetaInt, sqrt(tauSqInt));
        target += normal_lpdf(betaIA | muBetaIA, sqrt(tauSqIA));
        target += normal_lpdf(betaIL | muBetaIL, sqrt(tauSqIL));
        target += normal_lpdf(betaMI | muBetaMI, sqrt(tauSqMI));
        target += normal_lpdf(betaOH | muBetaOH, sqrt(tauSqOH));
        target += normal_lpdf(betaYear | muBetaYear, sqrt(tauSqYear));
        target += normal_lpdf(betaWeek | muBetaWeek, sqrt(tauSqWeek));
        target += normal_lpdf(betaWeek2 | muBetaWeek2, sqrt(tauSqWeek2));
        target += normal_lpdf(betaSurvey | muBetaSurvey, sqrt(tauSqSurvey));

	// Uniform prior on shape/scale for Gamma
	target += uniform_lpdf(r_count | 0, 20);
	
	// Random effects
	target += std_normal_lpdf(betaStarCountyRaw);
	target += std_normal_lpdf(betaStarWeekRaw);
	target += std_normal_lpdf(betaStarWeek2Raw);
	target += std_normal_lpdf(betaStarSiteRaw);
	target += std_normal_lpdf(betaStarYearRaw);

        for (i in 1:n_species) {
	  // Gamma random variable 
	  target += gamma_lpdf(rho[((i - 1) * n_surveys + 1):(i * n_surveys)] | r_count[i], r_count[i]);
	  // Poisson likelihood
	  target += poisson_lpmf(y[((i - 1) * n_surveys + 1):(i * n_surveys)] | mu_star[((i -1) * n_surveys + 1):(i * n_surveys)]);
        }
}

