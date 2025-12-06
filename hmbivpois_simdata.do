*******************************************
**** program to simulate random number ****
**** from a biv. poisson distribution  ****
*******************************************
capture mata mata drop rbivpoisson()
mata:
function rbivpoisson(y1,y2,mu1,mu2,c1,c2,lambda) {

	/** simulate Y2 ***/
	y2=rpoisson(1,1,mu2)
	pr2=poissonp(mu2,y2)
	
	/** simualte y1|y2 **/
	cp=0; r=0;dum=0
	while (cp<=0.99999) {
		assoc=1 + lambda*(exp(-r)-c1)*(exp(-y2)-c2);
		pr1=poissonp(mu1,r)
		p=(pr1*pr2*assoc)/pr2;
        cp=cp+p; pcp=cp-p;
		if (y1>=pcp & y1<=cp & dum==0) {
			y1=r
			dum=1
		}
		r=r+1		
	}
	return((y1,y2))
}
end

capture program drop sim_bivmc
program define sim_bivmc
    version 18.0
    
    syntax , [nobs(integer 1000) ntime(integer 1) k(integer 2) type(string) ///
         save(string) seed(integer 123456789)]

    if ("`save'" == "") local save "data.dta"
    if ("`seed'" == "") local seed = 123456789


    mata: N = `nobs'
    mata: T = `ntime'
    mata: K = `k'
   
	mata: tau    = vec(st_matrix("tau"))  
	mata: phi    = st_matrix("phi")
	mata: be1    = st_matrix("be1")
	mata: be2    = st_matrix("be2")
	mata: al1    = st_matrix("al1")
	mata: al2    = st_matrix("al2")
 
	mata: ga11   = st_numscalar("ga11")
	mata: ga12   = st_numscalar("ga12")
	mata: ga21   = st_numscalar("ga21")
	mata: ga22   = st_numscalar("ga22")
	mata: lambda = st_numscalar("lambda")
		
	mata: output=_sim_bivmc(N,T,K,tau,phi,be1,be2,al1,al2,ga11,ga12,ga21,ga22,lambda,"`type'",`seed')
	mata: id      = *output[1]
	mata: year    = *output[2]
	mata: Y1F     = *output[3]
	mata: LY1F    = *output[4]
	mata: Y2F     = *output[5]
	mata: LY2F    = *output[6]
	mata: x1F     = *output[7]
	mata: x2F     = *output[8]
	mata: xtransF = *output[9]
	
	********************************
	**** store dataset in stata ****
	********************************
	preserve
	drop _all
	mata: st_addobs(N*(T+1))
	qui foreach w in id year Y1F LY1F  Y2F LY2F {	
		mata: st_addvar("double", "`w'")
		mata: st_store(., "`w'",`w')	
	}
	rename Y1F dep1
	rename LY1F Ldep1
	rename Y2F dep2
	rename LY2F Ldep2
	qui mata: /*
	*/ 	for (i=1; i<=cols(x1F); i++) {; 				/*
	*/ 		st_addvar("double", "cov1"+strofreal(i)) ; 	/*
	*/ 		st_store(., "cov1"+strofreal(i),x1F[.,i]) ;	/*
	*/ 	}
	qui mata: /*
	*/ 	for (i=1; i<=cols(x2F); i++) {; 				/*
	*/		st_addvar("double", "cov2"+strofreal(i)); 	/*
	*/		st_store(., "cov2"+strofreal(i),x2F[.,i]); 	/*
	*/ 	}
	qui mata: /*
	*/ for (i=1; i<=cols(xtransF); i++) {;				/*
	*/ 	st_addvar("double", "xtrans"+strofreal(i));		/*
	*/ 	st_store(., "xtrans"+strofreal(i),xtransF[.,i]);/*
	*/ }
	if "`type'"=="HOM"{
		cap drop xtrans*
	}
	save "`save'.dta",replace
	restore


end

capture mata: mata drop _sim_bivmc()
mata 
function _sim_bivmc(N,T,K,tau,phi,be1,be2,al1,al2,ga11,ga12,ga21,ga22,lambda,type,seed) {
	
	printf("{res}Computation time can be substantial and increases with N, T, and K; please adjust your expectations accordingly.\n")
	
    /******************************/
    /**** 1. Generate covariate ***/
    /******************************/
	rseed(seed)
	idxt = ( (1..N)' # J(T,1,1) , J(N,1,1) # (1..T)' )   /* individual-time indicator */
	
	x1   = runiform(N*T,cols(be1))
	x2   = runiform(N*T,cols(be2))

	y10 = rpoisson(N,1,2)
	y20 = rpoisson(N,1,2)

	/* Covariate di transizione (come nel tuo codice) */
	tauT = rowshape(tau, K)'
	xtrans    = runiform(N*T,cols(tauT)-1)
	
	q = 1-exp(-1)
	
    /*************************************/
    /**** 2. Transition probabilities ****/
    /*************************************/
	printf("{res}Simulate transition probabilities:\n")
	if (type=="HET") {
		tauT=rowshape(tau,K)'
		PI=J(N*T*K,K,0); 
		idxtc=idxt#J(K,1,1); idxtc=(idxtc,J(rows(idxt),1,1)#(1..K)');
		xtrans=idxt[.,2]:==(2..T)#J(rows(idxt),1,1)
		next=10
		for (i=1; i<=N; i++) {		
			for (t=2; t<=T; t++) {
				posit=selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==t)):==2 )			 
				for (c=1; c<=K; c++) {
					positc=selectindex( ((idxtc[.,1]:==i)+(idxtc[.,2]:==t)+(idxtc[.,3]:==c)):==3 )
					Wxtrans=I(K)#J(1,cols(xtrans[posit,.]),1)					
					Wxtrans[c,.]=J(1,cols(Wxtrans),0)
					Wxtrans=Wxtrans[.,selectindex(colsum(Wxtrans):==1)]				
					Wxtrans=Wxtrans:*J(K,K-1,1)#(xtrans[posit,.])
					PI[positc,.]=(exp(Wxtrans*tauT[.,c]):/quadsum(exp(Wxtrans*tauT[.,c])))'
				}			
			}
			progress = 100*(i/N)
			if (i==N) {
					printf("{res}100%%\n")
					displayflush()
			}
			else if (i==1) {
					printf("{res}%2.0f%%...", 0)
					displayflush()
			}				
			else if (progress >= next) {
					printf("{res}%2.0f%%...", next)
					displayflush()
					next = next + 10
			}
		} 

	}	
	else {
		tauT=rowshape(tau,K)'
		PI=J(N*T*K,K,0); 
		idxtc=idxt#J(K,1,1); idxtc=(idxtc,J(rows(idxt),1,1)#(1..K)'); 
		next=10
		for (i=1; i<=N; i++) {	
			for (t=2; t<=T; t++) {
				posit=selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==t)):==2 )			 
				for (c=1; c<=K; c++) {
					positc=selectindex( ((idxtc[.,1]:==i)+(idxtc[.,2]:==t)+(idxtc[.,3]:==c)):==3 )
					Wxtrans=I(K)#J(1,1+cols(xtrans[posit,.]),1)
					Wxtrans[c,.]=J(1,cols(Wxtrans),0)
					Wxtrans=Wxtrans[.,selectindex(colsum(Wxtrans):==1)]		
					Wxtrans=Wxtrans:*J(K,K-1,1)#(1,xtrans[posit,.])
					PI[positc,.]=(exp(Wxtrans*tauT[.,c]):/quadsum(exp(Wxtrans*tauT[.,c])))'				
				}
			}
			progress = 100*(i/N)
			if (i==N) {
					printf("{res}100%%\n")
					displayflush()
			}
			else if (i==1) {
					printf("{res}%2.0f%%...", 0)
					displayflush()
			}				
			else if (progress >= next) {
					printf("{res}%2.0f%%...", next)
					displayflush()
					next = next + 10
			}
		}		
	}

	if (type=="HET") {
		printf("{res}Transition matrix:\n")
	}
	else {
		printf("{res}Average transition matrix:\n")
	}
	for (c=1; c<=K; c++) {
		if (c==1)	PIavg=mean(PI[selectindex(((idxtc[.,2]:>1)+(idxtc[.,3]:==c)):==2),.])
		else		PIavg=PIavg\mean(PI[selectindex(((idxtc[.,2]:>1)+(idxtc[.,3]:==c)):==2),.])
	}
	PIavg

	/******************************/
	/**** Predict Markov chain ****/
	/******************************/
	printf("{res}Predict latent states:\n")
	U = J(N*T,1,0); 
	next=10
	for (i=1; i<=N; i++) {	
		t=1
		pos=selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==t)):==2 )	
		de = phi[.,1]+phi[.,2]:*y10[i]+phi[.,3]:*y20[i];
		piv = exp((0 \ de)); piv = piv/quadsum(piv);
		u = K+1-sum(runiform(1,1):<quadrunningsum(piv'))		
		U[pos] = u
		
		for (t=2; t<=T; t++) {
			pos=selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==t)):==2 )
			PIt=PI[selectindex( ((idxtc[.,1]:==i)+(idxtc[.,2]:==t)):==2 ),.]	
			u= K+1-sum(runiform(1,1):<quadrunningsum(PIt[u,.]))
			U[pos] = u 
		}
		
		progress = 100*(i/N)
		if (i==N) {
				printf("{res}100%%\n")
				displayflush()
		}
		else if (i==1) {
				printf("{res}%2.0f%%...", 0)
				displayflush()
		}				
		else if (progress >= next) {
				printf("{res}%2.0f%%...", next)
				displayflush()
				next = next + 10
		}		
	}
	
	/**************************************/
	/**** Predict Y for each i,t and c ****/
	/**************************************/
	printf("{res}Predict Ys for each i,t and c:\n")
	Y1=runiform(N*T,1);Y2=runiform(N*T,1)
	for (i=1; i<=N*T; i++) {
		if (Y1[i]>=0.99999) Y1[i]=0.99991
		if (Y2[i]>=0.99999) Y2[i]=0.99991
	}
	next=10
	for (i=1; i<=N; i++) { 
		for (t=1; t<=T; t++) {
			if (t==1) {
				pos=selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==t)):==2 )
				
				lin1=x1[pos,.]*be1'+y10[i]*ga11+y20[i]*ga12
				lin2=x2[pos,.]*be2'+y10[i]*ga21+y20[i]*ga22
			}
			if (t>1) {
				pos =selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==t)):==2 )
				pos1=selectindex( ((idxt[.,1]:==i)+(idxt[.,2]:==(t-1))):==2 )
				
				lin1=x1[pos,.]*be1'+Y1[pos1]*ga11+Y2[pos1]*ga12
				lin2=x2[pos,.]*be2'+Y1[pos1]*ga21+Y2[pos1]*ga22				
			}
			c=U[pos]
			eta1=al1[c] + lin1; mu1=exp(eta1)
			eta2=al2[c] + lin2; mu2=exp(eta2)
			c1=exp(-q:*mu1);c2=exp(-q:*mu2);
			temp=rbivpoisson(Y1[pos],Y2[pos],mu1,mu2,c1,c2,lambda);		
			Y1[pos]=temp[1];Y2[pos]=temp[2]		
		}
		progress = 100*(i/N)
		if (i==N) {
				printf("{res}100%%\n")
				displayflush()
		}
		else if (i==1) {
				printf("{res}%2.0f%%...", 0)
				displayflush()
		}				
		else if (progress >= next) {
				printf("{res}%2.0f%%...", next)
				displayflush()
				next = next + 10
		}	
	}
		
	LY1=.; LY2=. /* from predicted Y1 and Y2 set the vectors of lags */
	for (i=1; i<=N; i++) {
		temp=(y10[i] \ Y1[selectindex( idxt[.,1]:==i )]);temp=temp[1..(rows(temp)-1)] 
		LY1=(LY1 \ temp)
		temp=(y20[i] \ Y2[selectindex( idxt[.,1]:==i )]);temp=temp[1..(rows(temp)-1)] 
		LY2=(LY2 \ temp)	
	}
	LY1=LY1[2..rows(LY1)];LY2=LY2[2..rows(LY2)]

	/*******************************/
	/**** Collect final dataset ****/
	/*******************************/
	id=.;year=.;
	Y1F=.;LY1F=.;x1F=J(1,cols(x1),.);
	Y2F=.;LY2F=.;x2F=J(1,cols(x2),.);
	xtransF=J(1,cols(xtrans),.);
	for (i=1; i<=N; i++) {
		pos=selectindex(idxt[.,1]:==i)
		id=(id \ J(T+1,1,1):*i);
		year=(year \ (0..T)');
		Y1F=(Y1F \ y10[i] \ Y1[pos])
		LY1F=(LY1F \ . \ LY1[pos] )
		x1F=(x1F \ J(1,cols(x1),.) \ x1[pos,.] )
		Y2F=(Y2F \ y20[i] \ Y2[pos])
		LY2F=(LY2F \ . \ LY2[pos] )
		x2F=(x2F \ J(1,cols(x2),.) \ x2[pos,.] )
		xtransF=(xtransF \ J(1,cols(xtrans),.) \ xtrans[pos,.] )
	}
	id=id[2..rows(id)]; year=year[2..rows(year)]
	Y1F=Y1F[2..rows(Y1F)]; LY1F=LY1F[2..rows(LY1F)];x1F=x1F[2..rows(x1F),.]
	Y2F=Y2F[2..rows(Y2F)]; LY2F=LY2F[2..rows(LY2F)];x2F=x2F[2..rows(x2F),.]
	xtransF=xtransF[2..rows(xtransF),.]
	
	return(&id,&year,&Y1F,&LY1F,&Y2F,&LY2F,&x1F,&x2F,&xtransF)
}
end
