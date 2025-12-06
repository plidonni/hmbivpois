clear all

cd "D:\Dropbox\ch1\code Ch1" //set your folder


***********************
**** simulate data ****
***********************
*** load programs to simulate data
qui do "hmbivpois_simdata.do"

*** set parameters 
*local het = "HET"
local het = "HOM"
*local het = "XTRANS"

local n = 15000 // sample size
local T = 8 // number of period including the intial time

if "`het'" == "XTRANS" {
    local k = 3
    
    matrix tau =  ( ///
        -1   ,  -1.5 , -2.5 \ ///
        -0.5 ,   0.5 ,  0   \ ///
         0.5 ,  -0.5 ,  1.5 \ ///
        -1   ,  -1.5 , -3   \ ///
        -0.5 ,  -0.5 ,  1.5 \ ///
        -0.5 ,  -0.5 ,  0   )
    
    matrix phi = ( ///
        -0.5 ,  1.5 , -1.5 \ ///
        -0.5 , -1.5 ,  1.5 )

    matrix al1 = (1 , 1.5 , 2)
    matrix be1 = (0.5 , -0.5)
    scalar ga11 = -0.1
    scalar ga12 = -0.1

    matrix al2 = (2 , 1.5 , 1)
    matrix be2 = (0.5 , -0.5)
    scalar ga21 = -0.1
    scalar ga22 = -0.1

    scalar lambda = 2
}

if "`het'" == "HET" {
    local k = 3

    matrix tau = ( ///
       -0.5 , -0.5 , -0.5 \ ///
       -0.8 , -0.8 , -0.8 \ ///
       -2.2 , -2.2 , -2.2 \ ///
       -1.5 , -1.5 , -1.5 \ ///
       -2.2 , -2.2 , -2.2 \ ///
       -1.5 , -1.5 , -1.5 \ ///
       -3   , -3   , -3   \ ///
       -0.5 , -0.5 , -0.5 \ ///
       -0.8 , -0.8 , -0.8 \ ///
       -2.2 , -2.2 , -2.2 \ ///
       -1.5 , -1.5 , -1.5 \ ///
       -2.2 , -2.2 , -2.2 \ ///
       -1.5 , -1.5 , -1.5 \ ///
       -3   , -3   , -3 )

    matrix phi = ( ///
        -0.5 ,  1.5 , -1.5 \ ///
        -0.5 , -1.5 ,  1.5 )

    matrix al1 = (1 , 1.5 , 2)
    matrix be1 = (0.5 , -0.5)
    scalar ga11 = -0.1
    scalar ga12 = -0.1

    matrix al2 = (2 , 1.5 , 1)
    matrix be2 = (0.5 , -0.5)
    scalar ga21 = -0.1
    scalar ga22 = -0.1

    scalar lambda = 2
}

if "`het'" == "HOM" {
    local k = 3

    matrix tau = ( ///
        -2   , -1   , -2.5 \ ///
         0   ,  0   ,  0   \ ///
         0   ,  0   ,  0   \ ///
        -2.2 , -1.5 , -2.5 \ ///
         0   ,  0   ,  0   \ ///
         0   ,  0   ,  0   )

    matrix phi = ( ///
        -0.5 ,  1.5 , -1.5 \ ///
        -0.5 , -1.5 ,  1.5 )

    matrix al1 = (1 , 1.5 , 2)
    matrix be1 = (0.5 , -0.5)
    scalar ga11 = -0.1
    scalar ga12 = -0.1

    matrix al2 = (2 , 1.5 , 1)
    matrix be2 = (0.5 , -0.5)
    scalar ga21 = -0.1
    scalar ga22 = -0.1

    scalar lambda = 2
}


sim_bivmc , nobs(`n') ntime(`T') k(`k') type("`het'") save("data") //simulate data


********************
**** estimation ****
********************
clear all
use "data.dta",clear

qui do "hmbivpois.do" // run estimation programs


hmbivpois	(dep1 = cov11 cov12 Ldep1 Ldep2 ) /*
	*/		(dep2 = cov21 cov22 Ldep1 Ldep2 ) /*
	*/,  k(3) id(id) time(year) /*
	*/  initcov(Ldep1 Ldep2)

	
*use option XTRANS to simulate data
hmbivpois	(dep1 = cov11 cov12 Ldep1 Ldep2 ) /*
	*/		(dep2 = cov21 cov22 Ldep1 Ldep2 ) /*
	*/,  k(3) id(id) time(year) /*
	*/  initcov(Ldep1 Ldep2) xtrans(xtrans1 xtrans2)  tracep em(20) save(mymod)
	

hmbivpois	(dep1 = cov11 cov12 Ldep1 Ldep2 ) /*
	*/		(dep2 = cov21 cov22 Ldep1 Ldep2 ) /*
	*/,  k(3) id(id) time(year) /*
	*/  initcov(Ldep1 Ldep2) xtrans(xtrans1 xtrans2)  tracep em(20) load(mymod) eval
	

