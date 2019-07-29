
// Set up a dataset based on the covariance matrix. Because Stata has limited
// support for analysis of covariance data, we will generate a dataset with
// matching covariances

matrix S = (4.00, ///
    3.23, 4.41, ///
    2.66, 2.67, 3.61, ///
    1.81, 1.50, 1.56, 4.41, ///
    1.85, 1.62, 1.63, 2.63, 4.84, ///
    1.24, 1.16, 1.09, 1.55, 1.72, 5.29, ///
    3.04, 2.83, 2.35, 2.07, 1.96, 1.31, 3.61, ///
    2.81, 2.57, 2.19, 1.55, 1.74, 1.12, 2.67, 3.24, ///
    1.73, 1.64, 1.32, 0.89, 1.05, 1.41, 1.62, 1.62, 2.89, ///
    2.66, 2.49, 2.12, 1.40, 1.43, 0.95, 2.22, 2.01, 1.25, 3.24, ///
    2.99, 2.86, 2.42, 1.52, 1.50, 1.01, 2.34, 2.14, 1.31, 3.05, 4.84)


corr2data acsi1 acsi2 acsi3 cuex1 cuex2 cuex3 perq1 perq2 perq3 perv1 perv2, ///
	n(10417) cov(S) cstorage(lower) clear



////////// Scale score correlations //////////////

gen acsi = acsi1 + acsi2 + acsi3
gen cuex = cuex1 + cuex2 + cuex3
gen perq = perq1 + perq2 + perq3
gen perv = perv1 + perv2

correlate acsi cuex perq perv

// Scale score correlations can also be calculated using unit weights,
// which can be implemented by standardizing the indicators before taking the
// sum. We use the user-written center command

center acsi1 acsi2 acsi3 cuex1 cuex2 cuex3 perq1 perq2 perq3 perv1 perv2, ///
	standardize prefix(z_)


gen z_acsi = z_acsi1 + z_acsi2 + z_acsi3
gen z_cuex = z_cuex1 + z_cuex2 + z_cuex3
gen z_perq = z_perq1 + z_perq2 + z_perq3
gen z_perv = z_perv1 + z_perv2

correlate z_acsi z_cuex z_perq z_perv


////////// Disattenuated correlation using tau-equivalent reliability //////////

alpha acsi1 acsi2 acsi3
scalar rel1 = r(alpha)

alpha perq1 perq2 perq3
scalar rel2 = r(alpha)

correlate acsi perq

display r(rho)/sqrt(rel1*rel2)


////////// Disattenuated correlation using parallel reliability (HTMT) //////////

alpha acsi1 acsi2 acsi3, std
scalar rel1 = r(alpha)

alpha perq1 perq2 perq3, std
scalar rel2 = r(alpha)

correlate z_acsi z_perq

// Disattenuated correlation
display r(rho)/sqrt(rel1*rel2)


////////// Disattenuated correlation using congeneric reliability //////////

sem (ACSI -> acsi1 acsi2 acsi3) ///
	(CUEX -> cuex1 cuex2 cuex3) ///
	(PERQ -> perq1 perq2 perq3) ///
	(PERV -> perv1 perv2), ///
	variance(ACSI@1 CUEX@1 PERQ@1 PERV@1) // Alternative scaling

estat framework
matrix loadings = r(Gamma)
matrix errors = r(Psi)
matrix factorCorrelations = r(Phi)

// Congeneric reliabilities

scalar rel1 = ((loadings[1,1]+loadings[2,1]+loadings[3,1])^2)/ ///
	((loadings[1,1]+loadings[2,1]+loadings[3,1])^2 + ///
	errors[1,1] + errors[2,2] + errors[3,3])

scalar rel1 = ((loadings[7,3]+loadings[8,3]+loadings[9,3])^2)/ ///
	((loadings[7,3]+loadings[8,3]+loadings[9,3])^2 + ///
	errors[7,7] + errors[8,8] + errors[9,9])


correlate acsi perq

// Disattenuated correlation
display r(rho)/sqrt(rel1*rel2)


////////// Cross-loadings (obtained from exploratory factor analysis) //////////


// Two factor solution, varimax rotation
factor acsi1-perv2, factors(2)

// Pattern coefficients are printed out by default
rotate, orthogonal varimax

// Structure coefficients equal  pattern coefficients
estat structure


// Two factor solution, promax rotation
factor acsi1-perv2, factors(2)

// Pattern coefficients are printed out by default
rotate, oblique promax

// Structure coefficients equal  pattern coefficients
estat structure


////////// Factor correlation (obtained from EFA) //////////


// Four factor solution, promax rotation. Factor correlations are printed out
// by default

factor acsi1-perv2, factors(4)
rotate, oblique promax

estat common

//////////  Techniques that require SEM software ////////// 


////////// Factor correlation (point estimate and confidence interval) ////////// 

sem (ACSI -> acsi1 acsi2 acsi3) ///
	(CUEX -> cuex1 cuex2 cuex3) ///
	(PERQ -> perq1 perq2 perq3) ///
	(PERV -> perv1 perv2), ///
	variance(ACSI@1 CUEX@1 PERQ@1 PERV@1) // Alternative scaling

estimates store m1



////////// Techniques using model fit indices: compared to ////////// 
////////// nested models with fewer factors                //////////

sem (ACSIPERQ -> acsi1 acsi2 acsi3 perq1 perq2 perq3) ///
	(CUEX -> cuex1 cuex2 cuex3) ///
	(PERV -> perv1 perv2), ///
	variance(ACSIPERQ@1 CUEX@1 PERV@1) // Alternative scaling

estimates store m2

lrtest m1 m2

// Same analysis implemented with constrained correlations

sem (ACSI -> acsi1 acsi2 acsi3) ///
	(CUEX -> cuex1 cuex2 cuex3) ///
	(PERQ -> perq1 perq2 perq3) ///
	(PERV -> perv1 perv2), ///
	variance(ACSI@1 CUEX@1 PERQ@1 PERV@1) /// Alternative scaling
	covariance(ACSI*PERQ@1 ACSI*CUEX@a PERQ*CUEX@a ACSI*PERV@b PERQ*PERV@b)

estimates store m3

lrtest m1 m3

////////// Techniques using model fit indices: comparison //////////
////////// against model with correlation fixed at 1      //////////

sem (ACSI -> acsi1 acsi2 acsi3) ///
	(CUEX -> cuex1 cuex2 cuex3) ///
	(PERQ -> perq1 perq2 perq3) ///
	(PERV -> perv1 perv2), ///
	variance(ACSI@1 CUEX@1 PERQ@1 PERV@1) /// Alternative scaling
	covariance(ACSI*PERQ@1)

estimates store m4

lrtest m1 m4


////////// Techniques using model fit indices: comparison //////////
////////// against model with correlation fixed at cutoff //////////

sem (ACSI -> acsi1 acsi2 acsi3) ///
	(CUEX -> cuex1 cuex2 cuex3) ///
	(PERQ -> perq1 perq2 perq3) ///
	(PERV -> perv1 perv2), ///
	variance(ACSI@1 CUEX@1 PERQ@1 PERV@1) /// Alternative scaling
	covariance(ACSI*PERQ@.96)

estimates store m5

lrtest m1 m5


////////// AVE: compared with the square of factor correlation //////////



scalar ave1 = (loadings[1,1]^2 + loadings[2,1]^2 + loadings[3,1]^2) / ///
	(loadings[1,1]^2 + loadings[2,1]^2 + loadings[3,1]^2 + ///
	errors[1,1] + errors[2,2] + errors[3,3])

scalar ave2 = (loadings[4,2]^2 + loadings[5,2]^2 + loadings[6,2]^2) / ///
	(loadings[4,2]^2 + loadings[5,2]^2 + loadings[6,2]^2 + ///
	errors[4,4] + errors[5,5] + errors[6,6])

scalar ave3 = (loadings[7,3]^2 + loadings[8,3]^2 + loadings[9,3]^2) / ///
	(loadings[7,3]^2 + loadings[8,3]^2 + loadings[9,3]^2 + ///
	errors[7,7] + errors[8,8] + errors[9,9])

scalar ave4 = (loadings[10,4]^2 + loadings[11,4]^2) / ///
	(loadings[10,4]^2 + loadings[11,4]^2 + ///
	errors[10,10] + errors[11,11])

// Raise each factor correlation to second power using mata
mata : st_matrix("comparisonMatrix", st_matrix("factorCorrelations") :^2)

matrix comparisonMatrix[1,1] = ave1
matrix comparisonMatrix[2,2] = ave2
matrix comparisonMatrix[3,3] = ave3
matrix comparisonMatrix[4,4] = ave4

// Print out the comparison matrix
matrix list comparisonMatrix

////////// AVE: compared with the square of scale score correlation //////////

correlate acsi cuex perq perv
matrix scaleScoreCorrelations = r(C)

// Raise each factor correlation to second power using mata
mata : st_matrix("comparisonMatrix", st_matrix("scaleScoreCorrelations") :^2)

matrix comparisonMatrix[1,1] = ave1
matrix comparisonMatrix[2,2] = ave2
matrix comparisonMatrix[3,3] = ave3
matrix comparisonMatrix[4,4] = ave4

// Print out the comparison matrix
matrix list comparisonMatrix


////////// Structure coefficients (obtained from CFA) //////////

matrix structureCoefficients = loadings * factorCorrelations
matrix list structureCoefficients

