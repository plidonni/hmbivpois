clear mata
*****************************************
capture mata: mata drop hmbivpois_index()
mata /*hmbivpois_index*/
function hmbivpois_index(idt,idtc,k) {
    
    /*preliminaries*/
	id=idt[.,1]; time=idt[.,2]
	n=rows(id); T=rows(uniqrows(time))
	
	for (t=1; t<=max(time); t++) {
	   if (t==1) {
		   post_T=&.
		   post0_T=&.
		   post_pi_T=&.
	   }
	   else {
		   post_T=post_T,&.
		   post0_T=post0_T,&.
		   post_pi_T=post_pi_T,&.	   
	   }
	}
	temp=(J(1,cols(idt),.)\idt); temp=temp[1..rows(idt),.]
	post=selectindex( idt[.,1]:!=temp[.,1] )
	(*post_T[1])=post
	for (t=2; t<=max(time); t++) {
		post0=selectindex(idt[.,2]:==(t-1))
		post=selectindex(idt[.,2]:==t)
		post_pi=selectindex(idtc[.,2]:==t)
		/*adjustmet for unbalance data*/
		if (rows(post0)!=rows(post)) {
		    post=post[selectindex(rowsum(id[post]:==((id[post0]')#J(rows(post),1,1))):==1)]
			post0=post0[selectindex(rowsum(id[post0]:==((id[post]')#J(rows(post0),1,1))):==1)]
			post_pi=post_pi[selectindex(rowsum(idtc[post_pi,1]:==((id[post0]')#J(rows(post_pi),1,1))):==1)]
		}
		else if (id[post]!=id[post0]){
		    post=post[selectindex(rowsum(id[post]:==((id[post0]')#J(rows(post),1,1))):==1)]
			post0=post0[selectindex(rowsum(id[post0]:==((id[post]')#J(rows(post0),1,1))):==1)]
			post_pi=post_pi[selectindex(rowsum(idtc[post_pi,1]:==((id[post0]')#J(rows(post_pi),1,1))):==1)]		    
		}
	    	
		(*post_T[t])=post
		(*post0_T[t])=post0
		(*post_pi_T[t])=post_pi
	}
	return(&post_T,&post0_T,&post_pi_T)
}
end

*****************************************
capture mata: mata drop bivpoisson_est_HM()
mata /*bivpoisson_est_HM*/
void bivpoisson_est_HM(typeoptimi,start_par,npar,all_id,Y,X,option,output,/*
					*/ tech,maxiter,tracep,traceg,traceh,start_par_name,/*
					*/th,V,LL){ 

	F  = optimize_init() /*F is the function we want to maximize*/
	optimize_init_evaluator(F, &bivpoisson_HM())
	optimize_init_evaluatortype(F, typeoptimi)
	optimize_init_params(F, start_par) /*initial parameters */
	optimize_init_argument(F, 1, npar) /*additional entries in void*/
	optimize_init_argument(F, 2, all_id) /*additional entries in void*/
	optimize_init_argument(F, 3, Y) /*additional entries in void*/
	optimize_init_argument(F, 4, X) /*additional entries in void*/
	optimize_init_argument(F, 5, option) /*additional entries in void*/
	optimize_init_argument(F, 6, output) /*additional entries in void*/
	optimize_init_technique(F,tech)
	optimize_init_conv_ptol(F, 1e-4) /* tol on parameters */
	optimize_init_conv_vtol(F, 1e-4) /*tol in lok-lik */
	optimize_init_conv_maxiter(F, maxiter)
	optimize_init_trace_params(F, tracep)
	optimize_init_trace_gradient(F, traceg)
	optimize_init_trace_Hessian(F, traceh)
	optimize_init_colstripe(F , start_par_name)
	th = optimize(F)
	V = optimize_result_V(F)
	LL=optimize_result_value(F)
}
end

*****************************************
capture mata: mata drop bivpoisson_der_fun_HM()
mata /*bivpoisson_der_fun_HM*/
void bivpoisson_der_fun_HM(real vector b,  /*
		*/   npar,all_id,Y,X,option,output, /*
		*/ LL) {

	LL1=.;S=.;H=.;
	bivpoisson_HM(0,b,  /*
		*/   npar,all_id,Y,X,option,output, /*
		*/ LL1, S, H) 
	LL=quadsum(LL1)
}
end

*****************************************
capture mata: mata drop bivpoisson_der_HM()
mata /*bivpoisson_der_HM*/
void bivpoisson_der_HM(start_par,npar,all_id,Y,X,option,output,/*
					*/d0,d1,d2){ 
	
	/*
	LL=.		    
	bivpoisson_der_fun_HM(start_par,npar,all_id,Y,X,option,output,LL)
	LL
	*/
		
	D = deriv_init() /* this is used to compute second derivative */
	deriv_init_evaluator(D, &bivpoisson_der_fun_HM())
	deriv_init_evaluatortype(D, "d")
	deriv_init_params(D, start_par) /* it must be a rowvector */
	deriv_init_argument(D, 1, npar) /*additional entries in void*/
	deriv_init_argument(D, 2, all_id) /*additional entries in void*/
	deriv_init_argument(D, 3, Y) /*additional entries in void*/
	deriv_init_argument(D, 4, X) /*additional entries in void*/
	deriv_init_argument(D, 5, option) /*additional entries in void*/
	deriv_init_argument(D, 6, output) /*additional entries in void*/
	deriv_init_search(D, "off")
	d0=deriv(D, 0)
	d1=deriv(D, 1)/*
	d2 =deriv(D, 2)*/

}
end

*****************************************
capture mata: mata drop est_bivpoisson()
mata /*est_bivpoisson*/
function est_bivpoisson(all_id,Y,X,V,paro,npar,trace,tracep,tolstepsize,steplength,maxiterf,usemata) {
    
	/* perform  Fisher's scoring to estimate lambda,al1,al2,be1,be2 in the complete bivariate poisson log-likelihood function
	--> input
	all_id		--> pointer matrix of individual index
	Y			--> matrix nt x 2 containing observed counts Y1 and Y2
	X 			--> pointer matrix 1 x 2 containing individual covariates for the first and second count
	V			--> posterior individual probabilities of belonging at a certain type in each time period 
	paro		--> vector of inital paramters j x 1
	npar		--> vector containing the number of paramters lambda,al1,al2,be1,be2 reported in b
	trace		--> "on" show steps of covergece, "off" otherwise
	tracep		--> "on"=show estimates at each iteration; "off"=otherwise
	tolstepsize	--> tolerance for asjusting increment of b at each iteration
	steplength	--> set the lenght of maximization
	maxiterf	--> max. number of iterations
	usemata		--> use mata optimizer to perform maximization

	--> output
	parest		--> 1 x h vector of estimated parameters containing estimates of lambda,al1,al2,be1,be2
	*/
	if (trace=="on")	printf("{res}Biv. Poisson:\n")
	if (usemata=="on") {
		F  = optimize_init() /*F is the function we want to maximize*/
		optimize_init_evaluator(F, &LLmarg_bivpoisson_LM())
		optimize_init_evaluatortype(F, "d2")
		optimize_init_params(F, paro') /*initial parameters */
		optimize_init_argument(F, 1, npar) /*additional entries in void*/
		optimize_init_argument(F, 2, V) /*additional entries in void*/
		optimize_init_argument(F, 3, X) /*additional entries in void*/
		optimize_init_argument(F, 4, Y) /*additional entries in void*/
		optimize_init_argument(F, 5, all_id) /*additional entries in void*/
		optimize_init_technique(F,"nr")
		/*optimize_init_which(F, "min" )*/
		/*optimize_init_nmsimplexdeltas(F,J(1,sum(npar),1))*/
		/*optimize_init_singularHmethod(F, "hybrid")*/
		/*optimize_init_conv_ignorenrtol(F, "on")*/
		optimize_init_conv_ptol(F, 1e-4) /* tol on parameters */
		optimize_init_conv_vtol(F, 1e-4) /*tol in lok-lik */
		/*optimize_init_conv_nrtol(F, 1e-2) */
		/*optimize_init_tracelevel(F,"none")*/
		optimize_init_conv_maxiter(F, maxiterf)
		optimize_init_trace_gradient(F, tracep) /*show gradient in each step */
		optimize_init_trace_params(F, tracep)	/* show parameters in each step */
		optimize_init_trace_value(F, trace) /* displays the function value at each step */		
		parest = optimize(F)
		parest=parest'
	}	
	else {
	    out=LLmarg_bivpoisson_LMg(2,paro',npar,V,X,Y,all_id)
		lkp=*out[1];lkop=0; g=*out[2]; H=*out[3]
		parest=paro;itp=0;parop=0;lk0=lkp;
		if (trace == "on")	printf("%s {txt}Iteration = {res}%f, {txt}f(p) = {res}%f \n", "", itp, lkp)
		do {
			itp=itp+1; lkop=lkp; parop=parest
			dpar=pinv(H)*g'
			mdpar=max(abs(dpar))
		
			pos=selectindex(abs(dpar):>=tolstepsize)
			parop=parest
			
			step_size=J(length(dpar),1,1)
			if (length(pos)>0) {
				for (i=1; i<=length(pos); i++) {
					step_size[pos[i]]= 1/steplength
				}
			}
			parest=parop-(step_size:*dpar)
			lk1=*LLmarg_bivpoisson_LMg(0,parest',npar,V,X,Y,all_id)[1]
			if (abs(lk1-lk0)< 0.1) { /* long step compared to lk0 */
				delta=J(length(dpar),1,-1)
				dpar=pinv(H)*g'
				parest=parop-dpar
			}				
						
			out=LLmarg_bivpoisson_LMg(2,parest',npar,V,X,Y,all_id)
			lkp=*out[1];g=*out[2];H=*out[3];
			if (trace == "on")	printf("%s {txt}Iteration = {res}%f, {txt}f(p) = {res}%f \n", "", itp, lkp)
			testp=((max(abs(g))>10^-5)+(abs(lkp-lkop)>10^-5)+(max(abs(parest:-parop))>10^-5))*(itp<maxiterf)
		} while (testp>=1) /*continua fino a quando test =1 */		
	}
	return(parest)
}
end

*****************************************
capture mata: mata drop LLmarg_bivpoisson_LM()
mata /*LLmarg_bivpoisson_LM*/
void LLmarg_bivpoisson_LM(real scalar todo, real vector b,  /*
		*/   npar,V,X,Y,all_id, /*
		*/   LL, S, H) {
			
		/* comute log-likelihood, score and hessian for the complete bivariate log-likelihood function
		--> input
		todo	--> 0=compute only log-lik.; 1=compute log-lik and score; 2=compute log-lik, score and hessian
		b		--> vector of inital paramters 1 x j
		npar	--> vector containing the number of paramters lambda,al1,al2,be1,be2 reported in b
		V		--> posterior individual probabilities of belonging at a certain type in each time period 
		X 		--> pointer matrix 1 x 2 containing individual covariates for the first and second count
		Y		--> matrix nt x 2 containing observed counts Y1 and Y2
		all_id	--> pointer matrix of individual index

		--> output
		LL		--> 1 x 1 log-likelihood
		S		--> score vector 1 x h
		H		--> hessian matri h x h
		*/
		
		if (todo==0) {
			CHECK=LLmarg_bivpoisson_LMg(0,b,npar,V,X,Y,all_id)
			LL=*CHECK[1];
		}
		if (todo==1) {
			CHECK=LLmarg_bivpoisson_LMg(1,b,npar,V,X,Y,all_id)
			LL=*CHECK[1];S=*CHECK[2];
		}
		if (todo==2) {
			CHECK=LLmarg_bivpoisson_LMg(2,b,npar,V,X,Y,all_id)
			LL=*CHECK[1];S=*CHECK[2];H=*CHECK[3];
		}
}
end

*****************************************
capture mata: mata drop LLmarg_bivpoisson_LMg()
mata /*LLmarg_bivpoisson_LMg*/
function LLmarg_bivpoisson_LMg(todo, real vector b,  /*
		*/   npar,V,X,Y,all_id) {
			
	/* comute log-likelihood, score and hessian for the complete bivariate log-likelihood function
	--> input
	todo	--> 0=compute only log-lik.; 1=compute log-lik and score; 2=compute log-lik, score and hessian
	b		--> vector of inital paramters 1 x j
	npar	--> vector containing the number of paramters lambda,al1,al2,be1,be2 reported in b
	V		--> posterior individual probabilities of belonging at a certain type in each time period 
	X 		--> pointer matrix 1 x 2 containing individual covariates for the first and second count
	Y		--> matrix nt x 2 containing observed counts Y1 and Y2
	all_id	--> pointer matrix of individual index

	--> output
	LL		--> 1 x 1 log-likelihood
	S		--> score vector 1 x h
	H		--> hessian matri h x h
	*/			
	
	/* collect parameters in b to be estimated */
	rnpar=runningsum(npar)
	if (min(npar)==0) {
		rnpar=J(1,length(npar),0)
		for (x=1; x<=length(npar); x++) { 
			if (npar[x]> 0) rnpar[x]=sum(npar[1..x])
			if (npar[x]==0) rnpar[x]=rnpar[x-1]
		}
	}
	lambda =b[1..rnpar[1]];
	al1    =b[(rnpar[1]+1)..rnpar[2]];
	al2    =b[(rnpar[2]+1)..rnpar[3]];
	if (npar[4]>0) be1    =b[(rnpar[3]+1)..rnpar[4]]'; 
	if (npar[5]>0) be2    =b[(rnpar[4]+1)..rnpar[5]]'; 
	
	/* preliminaries */
	idt=*all_id[1]; id=idt[.,1]; time=idt[.,2]
	n=rows(uniqrows(id)); k=cols(al1); T=rows(uniqrows(time));
	idtc=*all_id[2]
	U=V'; X1=*X[1]; X2=*X[2]
	Phi=J(rows(idt),k,.)

	/* *** check if data are balance *** */
	checkbal=1-hasmissing(idtw)	
	
	/* marginal distribution */
	for (c=1; c<=k; c++) {
		eta1=al1[c]; eta2=al2[c];
		if (cols(X1)> 0) eta1=eta1 :+ X1*be1; 
		if (cols(X2)> 0) eta2=eta2 :+ X2*be2; 
		mu1=exp(eta1); mu2=exp(eta2)		
		q=1-exp(-1)
		c1=exp(-q:*mu1);c2=exp(-q:*mu2);
		prob=poissonp(mu1,Y[.,1]):*poissonp(mu2,Y[.,2]):*(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
		Phi[.,c]=prob
	}
	
	/* compute log-lik */
	LL=quadsum(log(Phi)':*U)
	
	if (todo >= 1) {
	    /* compute score */
		SSlambda=J(1,1,0);SSal1=J(1,k,0);SSal2=J(1,k,0);SSbe1=J(1,npar[4],0);SSbe2=J(1,npar[5],0)
		SSlambda=J(1,1,0);SSal1=J(1,k,0);SSal2=J(1,k,0);SSbe1=J(1,npar[4],0);SSbe2=J(1,npar[5],0)
		for (c=1; c<=k; c++) {
			/* *** */				
			Utc=U[c,.]'
			eta1=al1[c]; eta2=al2[c]; 
			if (cols(X1)> 0) eta1=eta1 :+ X1*be1;			
			if (cols(X2)> 0) eta2=eta2 :+ X2*be2; 
			mu1=exp(eta1); mu2=exp(eta2)
			c1=exp(-q:*mu1);c2=exp(-q:*mu2);				
			psi1 = (lambda:*(exp(-Y[.,1]):-c1)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
			psi2 = (lambda:*(exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
			/* *** */
			temp=((exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
			temp=(temp):*Utc
			SSlambda=SSlambda+quadsum(temp)
			temp=Y[.,1]:-mu1:+q:*psi2:*mu1:*c1
			temp=(temp):*Utc
			SSal1[c]=SSal1[c]:+quadsum(temp)
			temp=Y[.,2]:-mu2:+q:*psi1:*mu2:*c2
			temp=(temp):*Utc
			SSal2[c]=SSal2[c]:+quadsum(temp)
			if (npar[4]>0) {
				temp=( Y[.,1]:-mu1:+q:*psi2:*mu1:*c1 ):*X1:*Utc
				SSbe1=SSbe1:+quadcolsum(temp)
			}
			if (npar[5]>0) {
				temp=( Y[.,2]:-mu2:+q:*psi1:*mu2:*c2 ):*X2:*Utc
				SSbe2=SSbe2:+quadcolsum(temp)
			}				
			
		}
		S=(SSlambda,SSal1,SSal2) /* score without covariates */
		/* S(be1) */
		if (npar[4]>0) {
			S=(S,SSbe1)
		}
		/* S(be2) */
		if (npar[5]>0) {
			S=(S,SSbe2)
		}
		
		if (todo >= 2) {/* compute hessian */		    
			hll=J(rows(idt),k,0); hlx1=J(rows(idt),k,0); hlx2=J(rows(idt),k,0);
			hx1x1=J(rows(idt),k,0); hx1x2=J(rows(idt),k,0);	hx2x2=J(rows(idt),k,0); 
			for (c=1; c<=k; c++) {
				eta1=al1[c]; eta2=al2[c]; 
				if (cols(X1)> 0) eta1=eta1 :+ X1*be1;			
				if (cols(X2)> 0) eta2=eta2 :+ X2*be2; 
				mu1=exp(eta1); mu2=exp(eta2)
				c1=exp(-q:*mu1);c2=exp(-q:*mu2);
				
				psi1 = (lambda:*(exp(-Y[.,1]):-c1)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
				psi2 = (lambda:*(exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
				/* *** */
				hll[.,c] = ((exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
			
				phi3 = ((exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
				phi4 = (lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2):^2):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2)):^2;
				hlx1[.,c] = ( phi3 :- phi4 ):*q:*mu1:*c1;
				
				phi3 = ((exp(-Y[.,1]):-c1)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
				phi4 = (lambda:*((exp(-Y[.,1]):-c1):^2):*(exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2)):^2;
				hlx2[.,c] = ( phi3 :- phi4 ):*q:*mu2:*c2;	
				
				hx1x1[.,c] = -mu1 :+ psi2:*q:*mu1:*c1 :- psi2:*q:^2:*mu1:^2:*c1 :- (psi2:*q:*mu1:*c1):^2;
			
				phi3 = (lambda):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2));
				phi4 = (lambda^2:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2)):/(1:+lambda:*(exp(-Y[.,1]):-c1):*(exp(-Y[.,2]):-c2)):^2;
				hx1x2[.,c] = (phi3:-phi4):*q^2:*mu1:*mu2:*c1:*c2;
				
				hx2x2[.,c] = -mu2 :+ psi1:*q:*mu2:*c2 :- psi1:*q^2:*mu2:^2:*c2 :- (psi1:*q:*mu2:*c2):^2;
				
			
			}
			H=J(sum(npar),sum(npar),0)
			/* H(lambda) */
			nd=1;
			for (x=1; x<=npar[nd]; x++) {
				/*H(lambda,lambda)*/
				pos = x
				hess = (-hll:^2):*U'
				H[pos,pos]=H[pos,pos]+quadsum(hess)
				/*H(lambda,al1)*/
				nc=2
				pos1 = (1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
				hess = hlx1:*U'
				hess = quadcolsum(hess)
				H[pos,pos1]= H[pos,pos1]:+hess
				H[pos1,pos]= H[pos1,pos]:+hess'				
				/*H(lambda,al2)*/
				nc=3
				pos1 = (1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
				hess = hlx2:*U'
				hess = quadcolsum(hess)
				H[pos,pos1]= H[pos,pos1]:+hess
				H[pos1,pos]= H[pos1,pos]:+hess'
				/*H(lambda,be1)*/
				if (npar[4]>0) {
					nc=4
					pos1 = (1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					/*
					for (c=1; c<=k; c++) {
						Utc=U[c,.]'
						hess = hlx1[.,c]:*X1:*Utc
						hess = quadcolsum(hess)
						H[pos,pos1]= H[pos,pos1]:+hess
						H[pos1,pos]= H[pos1,pos]:+hess'
					}
					*/
					H[pos,pos1]= quadrowsum(X1'*(hlx1:*U'))'
					H[pos1,pos]= H[pos,pos1]'
				}
				/*H(lambda,be2)*/
				if (npar[5]>0) {
					nc=5
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					/*
					for (c=1; c<=k; c++) {
						Utc=U[c,.]'
						hess = hlx2[.,c]:*X2:*Utc
						hess = quadcolsum(hess)
						H[pos,pos1]= H[pos,pos1]:+hess
						H[pos1,pos]= H[pos1,pos]:+hess'
					}
					*/
					H[pos,pos1]= quadrowsum(X2'*(hlx2:*U'))'					
					H[pos1,pos]= H[pos,pos1]'
				}
			}
			/* H(al1) */
			nd=2;
			for (x=1; x<=npar[nd]; x++) {
				/*H(al1,al1)*/
				pos = x+sum(npar[1..(nd-1)]);
				c=x;
				Utc=U[c,.]'
				hess = hx1x1[.,c]
				H[pos,pos]=H[pos,pos]+quadsum( hess:*Utc )				
				/*H(al1,al2)*/
				nc=3
				for (xx=1; xx<=npar[nc]; xx++) {
					if (x==xx) {
						pos1 = xx+sum(npar[1..(nc-1)]);						
						Utc=U[c,.]'
						hess = hx1x2[.,c];
						H[pos,pos1]= H[pos,pos1]+quadsum(hess:*Utc);
						H[pos1,pos]= H[pos1,pos]+quadsum(hess:*Utc);					
					}
				}
				/*H(al1,be1)*/
				nc=4
				if (npar[nc]>0) {
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					Utc=U[c,.]'
					hess = hx1x1[.,c]:*X1:*Utc
					hess=quadcolsum(hess)
					H[pos,pos1]= H[pos,pos1]:+hess
					H[pos1,pos]= H[pos1,pos]:+hess'
				}
				/*H(al1,be2)*/
				nc=5
				if (npar[nc]>0) {
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					Utc=U[c,.]'
					hess = hx1x2[.,c]:*X2:*Utc
					hess = quadcolsum(hess)
					H[pos,pos1]= H[pos,pos1]:+hess
					H[pos1,pos]= H[pos1,pos]:+hess'
				}
			}
			/* H(al2) */
			nd=3;
			for (x=1; x<=npar[nd]; x++) {
				pos = x+sum(npar[1..(nd-1)]);
				c=x
				Utc=U[c,.]'
				hess = hx2x2[.,c];
				H[pos,pos]=H[pos,pos]+quadsum( hess:*Utc )
				/*H(al2,be1)*/
				nc=4
				if (npar[nc]>0) {
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					Utc=U[c,.]'
					hess = hx1x2[.,c]:*X1:*Utc
					hess = quadcolsum(hess)
					H[pos,pos1]= H[pos,pos1]:+hess
					H[pos1,pos]= H[pos1,pos]:+hess'
				}
				/*H(al2,be2)*/
				nc=5
				if (npar[nc]>0) {
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					Utc=U[c,.]'
					hess = hx2x2[.,c]:*X2:*Utc
					hess = quadcolsum(hess)
					H[pos,pos1]= H[pos,pos1]:+hess
					H[pos1,pos]= H[pos1,pos]:+hess'
				}
			}			
			/* H(be1) */
			nd=4;
			if (npar[nd]>0) {
				for (x=1; x<=npar[nd]; x++) {
					pos = x+sum(npar[1..(nd-1)]);					
					/*H(be1,be1+)*/
					nc=4
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					H[pos,pos1]=quadrowsum((X1[.,x]:*X1)'*(hx1x1:*U'))'
					/*H(be1,be2)*/
					nc=5
					if (npar[nc]>0) {
						pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
						H[pos,pos1]= quadrowsum((X1[.,x]:*X2)'*(hx1x2:*U'))'
						H[pos1,pos]= H[pos,pos1]'
					}
					/*
					for (c=1; c<=k; c++) {
						Utc=U[c,.]'
						/*H(be1,be1+)*/
						nc=4
						pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
						hess = hx1x1[.,c]:*X1[.,x]:*X1:*Utc
						hess = quadcolsum(hess)
						H[pos,pos1]= H[pos,pos1]:+hess
						/*H(be1,be2)*/
						nc=5
						if (npar[nc]>0) {
							pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
							hess = hx1x2[.,c]:*X1[.,x]:*X2:*Utc
							hess = quadcolsum(hess)
							H[pos,pos1]= H[pos,pos1]:+hess
							H[pos1,pos]= H[pos1,pos]:+hess'
						}
					}
					*/
				}
			}
			/* H(be2) */
			nd=5;
			if (npar[nd]>0) {
				for (x=1; x<=npar[nd]; x++) {
					pos = x+sum(npar[1..(nd-1)]);
					/*H(be2,be2+)*/
					nc=5
					pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
					H[pos,pos1]= quadrowsum((X2[.,x]:*X2)'*(hx2x2:*U'))'
					/*
					for (c=1; c<=k; c++) {
						Utc=U[c,.]'
						/*H(be2,be2+)*/
						nc=5
						pos1=(1+sum(npar[1..(nc-1)]))..(npar[nc]+sum(npar[1..(nc-1)]))
						hess = hx2x2[.,c]:*X2[.,x]:*X2:*Utc
						hess = quadcolsum(hess)
						H[pos,pos1]= H[pos,pos1]:+hess
					}
					*/
				} /*nd*/
			}
			/*H=makesymmetric(H)*/
			return((&LL,&S,&H))
		} /* todo 2 */	
		return((&LL,&S))
	} /* todo 1*/
	return((&LL))
	/*LL=-LL;S=-S;H=-H*/
}
/* mata mosave LLmarg_bivpoisson_LMg(),dir(PERSONAL) replace*/
end


*****************************************
capture mata: mata drop bivpoisson_est_EM_HM()
mata /*bivpoisson_est_EM_HM*/
void bivpoisson_est_EM_HM(b,npar,all_id,Y,X,option,output,/*
					*/th,LL){ 

	/* collect parameters in b to be estimated */
	rnpar=runningsum(npar)
	if (min(npar)==0) {
		rnpar=J(1,length(npar),0)
		for (x=1; x<=length(npar); x++) { 
			if (npar[x]> 0) rnpar[x]=sum(npar[1..x])
			if (npar[x]==0) rnpar[x]=rnpar[x-1]
		}
	}
	
	al1 					=b[1..rnpar[1]]';
	if (npar[2]>0) be1		=b[(rnpar[1]+1)..rnpar[2]]';
	al2 					=b[(rnpar[2]+1)..rnpar[3]]';
	if (npar[4]>0) be2		=b[(rnpar[3]+1)..rnpar[4]]';
	if (npar[5]>0) lambda	=b[(rnpar[4]+1)..rnpar[5]]';
	if (npar[6]>0) be3    	=b[(rnpar[5]+1)..rnpar[6]]'; 
	if (npar[7]>0) be4    	=b[(rnpar[6]+1)..rnpar[7]]';
	
	if (npar[5]==0) 	lambda=0

	/*preliminaries*/
	idt=*all_id[1]; id=idt[.,1]; time=idt[.,2]
	n=rows(id); k=rows(al1); T=rows(uniqrows(time));
	idtc=*all_id[2]
	index=*all_id[3];post_T=*index[1];post0_T=*index[2];post_pi_T=*index[3];	
	X1=*X[1]; X2=*X[2];	WXI=*X[3];WXT=*X[4]
	restr=*option[1];maxiter=*option[2]
	disparamEM=*option[3]
	
	/*initial likelihood*/
	LL=.
	bivpoisson_HM(0, b,  /*
		*/   npar,(&(id,time),&idtc,&index),Y,X,option,output, /*
		*/ LL, ., .)
	lk=quadsum(LL)
	printf("{txt}Performing EM algorithm: ")
	printf(" \n")
	printf("{txt} EM step = {res}%f, ", 0)
	printf("{txt}log-lik = {res}%12.0g ", lk)
	printf(" \n")
	Piv=(*output[1])
	PI=(*output[2])
	mu1=(*output[3])
	mu2=(*output[4])
	BPo=(*output[5])
	pv=(*output[6])
	Qt=(*output[7])
	/*--- EM ---*/
	it=0; th=b'; lko=lk-10^10; lkv=lko; /*(lk-lko)/abs(lk)*/	 
	while ( ( ((lk-lko)>0.001) & (it<maxiter) ) ==1) {
	    it=it+1; lko=lk; tho=th
		timer_on(1)
		
		/*--- E-step ---*/
		/*backward recursion to compute U and V*/
		out=hmbivpois_rec_backward(all_id,pv,Qt,BPo,Piv,PI)
		U=*out[1]; V=*out[2]
		
		/*--- M-step ---*/
		/*update be3 for Piv*/
		if (k>1) {
			tmp=0\idt[.,1];tmp=tmp[1..rows(idt)]
			pos=selectindex(tmp:!=idt[.,1])
			out=est_multilogit(V[pos,.],WXI,Piv,be3,k,1,100,"off")			
			be3=*out[1]; Piv=*out[2]
		}
		/*update be4 for PI*/
		if (k>1 & restr==0) {
		    be4T=rowshape(be4,k)'
			for (c=1; c<=k; c++) {		
				temp=(J(k,3,-99)\idtc); temp=temp[1..rows(idtc),.]
				pos=selectindex( (idtc[.,1]:==temp[.,1])+(idtc[.,3]:==c):==2 )
				out=est_multilogit(U[pos,.],WXT,PI[pos,.],be4T[.,c],k,c,100,"off")
				be4T[.,c]=*out[1]; PI[pos,.]=*out[2]
			}			
			be4=vec(be4T)
		}
		/*update al1,be1,al2,be2,lambda for Phi*/
		npar1=npar[5],npar[1],npar[3],npar[2],npar[4]
		parest=(lambda\al1\al2)
		
		if (npar[2]>0)	parest=(parest\be1)
		if (npar[4]>0)	parest=(parest\be2)
		parest=est_bivpoisson(all_id,Y,X,V,parest,npar1,"off","off",0.1,10,50,"off")
		rnpar1=runningsum(npar1)
		if (min(npar1)==0) {
			rnpar1=J(1,length(npar1),0)
			for (x=1; x<=length(npar1); x++) { 
				if (npar1[x]> 0) rnpar1[x]=sum(npar1[1..x])
				if (npar1[x]==0) rnpar1[x]=rnpar1[x-1]
			}
		}
		lambda =parest[1..rnpar1[1]];
		al1    =parest[(rnpar1[1]+1)..rnpar1[2]];
		al2    =parest[(rnpar1[2]+1)..rnpar1[3]];
		if (npar1[4]>0) be1    =parest[(rnpar1[3]+1)..rnpar1[4]]; 
		if (npar1[5]>0) be2    =parest[(rnpar1[4]+1)..rnpar1[5]];
		
		/*compute log-likelihood*/
		th=al1;
		if (npar[2]>0) th=th\be1;
		th=th\al2
		if (npar[4]>0) th=th\be2;
		if (npar[5]>0) th=th\lambda;
		if (npar[6]>0) th=th\be3;
		if (npar[7]>0) th=th\be4;
		bivpoisson_HM(0, th',  /*
		*/   npar,(&(id,time),&idtc,&index),Y,X,option,output, /*
		*/ LL, ., .)
		lk=quadsum(LL)
		Piv=(*output[1])
		PI=(*output[2])
		mu1=(*output[3])
		mu2=(*output[4])
		BPo=(*output[5])
		pv=(*output[6])
		Qt=(*output[7])
		
		timer_off(1)
		eti=timer_value(1)
		timer_clear(1)
		
/*
		printf(" \n")
		printf("{txt}k = {res}%f, ", k)
		printf("{txt}step = {res}%f, ", it)
		printf("{txt}log-lik = {res}%12.0g, ", lk)
		printf("{txt}diff.log-lik = {res}%7.0g, ", (lk-lko)) /*(lk-lko)/abs(lk)*/
		printf("{txt}diff.par = {res}%7.0g, ",max(abs(th:-tho)) )
		printf("{txt}eti = {res}%4.0f, ", eti[1])
*/
		NT=rows(id);
		bic=(log(NT)*sum(npar) - 2*(lk))
		printf("%s {txt}EM-iteration = {res}%f, {txt}f(p) = {res}%f, {txt}diff. = {res}%f, {txt}diff.par. = {res}%f, {txt}Ex.Time = {res}%f, {txt}BIC = {res}%f \n", "", it, lk,abs(lk-lko),abs(sum(th:-tho)),eti[1]/60,bic)
				
		if (disparamEM[1,1]!=0) {			
			st_matrix("betaEM", th')
			st_matrixcolstripe("betaEM", disparamEM)			
			st_matrix_list("betaEM")
		}
		/*
		printf(" \n")
		printf("{txt}k = {res}%f, ", k)
		printf("{txt}step = {res}%f, ", it)
		printf("{txt}log-lik = {res}%12.0g, ", lk)
		printf("{txt}diff.log-lik = {res}%7.0g, ", (lk-lko)) /*(lk-lko)/abs(lk)*/
		printf("{txt}diff.par = {res}%7.0g, ",max(abs(th:-tho)) )
		printf("{txt}eti = {res}%4.0f, ", eti[1])
		*/		
	 }
}
end

*****************************************
capture mata mata drop prob_multilogit()
mata /*prob_multilogit*/
function prob_multilogit(X,be,k) {
    
	/* comute probs. of multinomial logit
	--> input
	X		--> design matrix
	be		--> vector of initial parameters k x 1

	--> output
	Prob	-> matrix of estimated probs. n x j

	*/
	
	/* *** preliminaries *** */
	np=cols(rowshape(be,k-1))	
	
	Eta=X*(be:*(I(k-1)#J(np,1,1)))
	Eta=(J(rows(Eta),1,0),Eta)
	Prob=exp(Eta):/quadrowsum(exp(Eta))

	/*
	if (der == 1) {
		lbe=rows(be)
		dProb=J(n*(rows(selectindex(rowsum(id[1,.]:==id0):==d))),lbe,0)
		for (i=1; i<=n; i++) {
			pos=selectindex(rowsum(id[i,.]:==id0):==d)
			Op=diag(Prob[i,.])-Prob[i,.]'*Prob[i,.]
			dProb[pos,.]=Op*X[pos,.]
		}
	}
	else dProb=.
	return((&Prob,&dProb))	
	*/
	return(&Prob)

}
end

*****************************************
capture mata mata drop est_multilogit()
mata /*est_multilogit*/
function est_multilogit(Y,X,P_init,be_init,k,base,maxiterm,trace,|der) {

	/* estimates multinomial logit
	--> input
	X		--> dep. var
	X		--> design matrix
	P		--> matrix of initial prob
	be_init		--> vector of initial parameters k x 1
	k		--> number of alternatives
	base	--> define base alternative
	maxiterm--> max number of iterations in multinomial logit estimation
	trace	--> "on" show steps of covergence
	der		--> "deriv" compute only score and hessian without estimation

	--> output
	be		--> estimated parameters
	PF		-> matrix of estimated probs. n x j
	sc		--> score vector
	F		--> hessin matrix

	*/
	/*preliminaries*/
	be=be_init
	np=cols(rowshape(be,k-1)) 
	nnp=rows(be) /*total number of parameters to be estimated*/
	ind=(1..(k-1))'#J(np,1,1)
	ny= quadrowsum(Y)
	notbase=selectindex((1..k):!=base)
	
	P=P_init
	YF=Y;PF=P
	Y=YF[.,notbase]
	P=PF[.,notbase]
	
	/*initial log-likelihood*/
	lk=quadsum(YF:*log(PF))
	
	/* iterate until convergence */
	it=0; lko=lk
	if (trace=="on")	printf("%s {txt}Iteration = {res}%f, {txt}f(p) = {res}%f \n", "", it, lk)

 	while ( (lk-lko)>10^-6 && it<maxiterm || it==0) {
	    it=it+1; lko=lk
		F=J(nnp,nnp,0)
		Dp=(Y:-(ny:*P))#J(1,np,1)
		sc=quadcolsum(Dp:*X)'
		for (b=1; b<=nnp; b++) {
			for (bb=1; bb<=nnp; bb++) {		
				if (ind[b]==ind[bb]) 	p=P[.,ind[b]]:*(1:-P[.,ind[bb]])
				else					p=P[.,ind[b]]:*(  -P[.,ind[bb]])
				temp=ny:*p:*X[.,b]:*X[.,bb]			
				F[b,bb]=quadcolsum(temp)
			}
		}
		if (der=="deriv")	break
		dbe=pinv(F)*sc	
		mdbe=max(abs(dbe))
		
		if(mdbe>0.5)	dbe = dbe/mdbe*0.5
		be0=be
		flag=1
		while (flag==1) {
		    be=be0:+dbe
			Eta=X*(be:*(I(k-1)#J(np,1,1)))
			if(max(abs(Eta))>100){
			    dbe = dbe/2
				flag = 1
	        } else {
			    flag = 0
	        }
		}

		/* compute probabilities */
		out=*prob_multilogit(X,be,k)[1]
		PF[.,base]=out[.,1];
		PF[.,notbase]=out[.,2..k];
		P=PF[.,notbase]		
		lk=quadsum(YF:*log(PF))
		if (trace=="on")	printf("%s {txt}Iteration = {res}%f, {txt}f(p) = {res}%f \n", "", it, lk)
	}
	return(&be,&PF,&sc,&F)	

}
end


*****************************************
capture mata: mata drop bivpoisson_HM()
mata /*bivpoisson_HM*/
void bivpoisson_HM(real scalar todo,real vector b,  /*
		*/   npar,all_id,Y,X,option,output, /*
		*/ LL, S, H) {
			
		/* comute log-likelihood, score and hessian for the complete bivariate log-likelihood function
		--> input
		todo	--> 0=compute only log-lik.; 1=compute log-lik and score; 2=compute log-lik, score and hessian
		b		--> vector of inital paramters 1 x j
		npar	--> vector containing the number of paramters lambda,al1,al2,be1,be2 reported in b
		V		--> posterior individual probabilities of belonging at a certain type in each time period 
		X 		--> pointer matrix 1 x 2 containing individual covariates for the first and second count
		Y		--> matrix nt x 2 containing observed counts Y1 and Y2

		--> output
		LL		--> 1 x 1 log-likelihood
		S		--> score vector 1 x h
		H		--> hessian matri h x h
		*/
	
	/* collect parameters in b to be estimated */
	rnpar=runningsum(npar)
	if (min(npar)==0) {
		rnpar=J(1,length(npar),0)
		for (x=1; x<=length(npar); x++) { 
			if (npar[x]> 0) rnpar[x]=sum(npar[1..x])
			if (npar[x]==0) rnpar[x]=rnpar[x-1]
		}
	}
	
	al1 					=b[1..rnpar[1]];
	if (npar[2]>0) be1		=b[(rnpar[1]+1)..rnpar[2]]';
	al2 					=b[(rnpar[2]+1)..rnpar[3]];
	if (npar[4]>0) be2		=b[(rnpar[3]+1)..rnpar[4]]';
	if (npar[5]>0) lambda	=b[(rnpar[4]+1)..rnpar[5]]';
	if (npar[6]>0) be3    	=b[(rnpar[5]+1)..rnpar[6]]'; 
	if (npar[7]>0) be4    	=b[(rnpar[6]+1)..rnpar[7]]';
	
	if (npar[5]==0) lambda=0

	/*preliminaries*/
	idt=*all_id[1]; id=idt[.,1]; time=idt[.,2]
	n=rows(id); k=cols(al1); T=rows(uniqrows(time));
	idtc=*all_id[2]
	index=*all_id[3];post_T=*index[1];post0_T=*index[2];post_pi_T=*index[3];	
	X1=*X[1]; X2=*X[2];	WXI=*X[3];WXT=*X[4]
	restr=*option[1]

	/*initial prob.*/
	if (k>1) {
		np=cols(rowshape(be3,k-1))	
		Eta=WXI*(be3:*(I(k-1)#J(np,1,1)))
		Eta=(J(rows(Eta),1,0),Eta)
		Piv=exp(Eta):/quadrowsum(exp(Eta))
	}
	else Piv=J(rows(uniqrows(id)),1,1)

	/*transition prob.*/
	if (k>1) {
		if (restr==1){
			temp=(J(k,cols(idtc),.)\idtc); temp=temp[1..rows(idtc),.]
			pos=selectindex( (idtc[.,1]:==temp[.,1]) )
			PI=J(rows(idtc),k,.)
			temp=rows(uniqrows(idtc[pos,1..2]))
			PI[pos,.]=J(temp,1,1)#I(k)
		}		
		else {
			be4T=rowshape(be4,k)'
			PI=J(rows(idtc),k,.)
			temp=(J(k,cols(idtc),.)\idtc); temp=temp[1..rows(idtc),.]
			for (c=1; c<=k; c++) {
				pos=selectindex( (idtc[.,1]:==temp[.,1])+(idtc[.,3]:==c):==2 )
				np=cols(rowshape(be4T[.,c],k-1))
				Eta=WXT*(be4T[.,c]:*(I(k-1)#J(np,1,1)))
				Eta=(J(rows(Eta),1,0),Eta)
				temp1=exp(Eta):/quadrowsum(exp(Eta))
				PI[pos,c]=temp1[.,1]
				PI[pos,selectindex((1..k):!=c)]=temp1[.,2..k]				
			}
		}
	}
	else	PI=J(rows(idtc),k,1)
	
	/*marginal prob.*/
	q=1-exp(-1)
	c1=exp(-Y[.,1]);c2=exp(-Y[.,2]);
	mu1=J(rows(idt),k,.);mu2=J(rows(idt),k,.);
	d1=J(rows(idt),k,.);d2=J(rows(idt),k,.);
	Po1=J(rows(idt),k,.);Po2=J(rows(idt),k,.);	
	BPo=J(rows(idt),k,.);
	for (c=1; c<=k; c++) {
		eta1=J(n,1,al1[c]); eta2=J(n,1,al2[c]);
		if (cols(X1)> 0) eta1=eta1 :+ X1*be1; 
		if (cols(X2)> 0) eta2=eta2 :+ X2*be2;		
		mu1[.,c]=exp(eta1); mu2[.,c]=exp(eta2)
		d1[.,c]=exp(-q:*mu1[.,c]);d2[.,c]=exp(-q:*mu2[.,c]);		
		Po1[.,c]=poissonp(mu1[.,c],Y[.,1]); Po2[.,c]=poissonp(mu2[.,c],Y[.,2])
		BPo[.,c]=Po1[.,c]:*Po2[.,c]:*(1:+lambda:*(c1 - d1[.,c]):*(c2 - d2[.,c])) 
	}

	/*compute log-likelihood*/
	Qt=J(n,k,0)
	/*initial period*/
	post=(*post_T[1])
	Qt[post,.]=BPo[post,.]:*Piv

	/*subsequent periods*/
	for (t=2; t<=T; t++) {		 
		post0=(*post0_T[t])
		post=(*post_T[t])
		post_pi=(*post_pi_T[t])
		
		tmp=PI[post_pi,.]
		for (c=1; c<=k; c++) {
			Qt[post,c]=quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c]
		}
	}
	if (min(Qt)==0)	Qt=Qt:+((Qt:==0):*1e-323)

	/*individual likelihood at the last period*/
	temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]		
	post=selectindex( idt[.,1]:!=temp[.,1] )
	pv=quadrowsum(Qt[post,.])
	LL=log(pv)

	/*store some results*/
	(*output[1])=Piv
	(*output[2])=PI
	(*output[3])=mu1
	(*output[4])=mu2
	(*output[5])=BPo
	(*output[6])=pv
	(*output[7])=Qt

	if (todo==1) {
	    S=J(rows(uniqrows(id)),sum(npar),0)
		
		/*gradient for al1*/
		dBPo=Po1:*Po2:*( ((Y[.,1] :- mu1):*(1:+lambda:*(c1 :- d1):*(c2 :- d2))) :+ /*
			*/ (lambda:*q:*mu1:*d1:*(c2:-d2)) ):*1
		for (x=1; x<=npar[1]; x++) {			
			Qt_der=J(rows(idt),k,0)
			/*initial period*/
			post=(*post_T[1])			
			Qt_der[post,x]=dBPo[post,x]:*Piv[.,x]
			/*subsequent periods*/
			for (t=2; t<=T; t++) {
				post0=(*post0_T[t])
				post=(*post_T[t])
				post_pi=(*post_pi_T[t])
				tmp=PI[post_pi,.]			
				for (c=1; c<=k; c++) {
					if (c==x) Qt_der[post,c]=quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c] :+ /*
								*/ quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post))):*dBPo[post,c]
					if (c!=x) Qt_der[post,c]= quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c]
				}				
			}
			/*last period*/
			temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
			pos3=selectindex( idt[.,1]:!=temp[.,1] )
			Qt_der=quadrowsum(Qt_der[pos3,.]):/pv
			S[.,(1..rnpar[2])[x]]=Qt_der
		}			
						
		/*gradient for be1*/
		for (x=1; x<=npar[2]; x++) {
			dBPo=Po1:*Po2:*( ((Y[.,1] :- mu1):*(1:+lambda:*(c1 :- d1):*(c2 :- d2))) :+ /*
			*/ (lambda:*q:*mu1:*d1:*(c2:-d2)) ):*X1[.,x]
			
			Qt_der=J(rows(idt),k,0)
			/*initial period*/
			post=(*post_T[1])			
			Qt_der[post,.]=dBPo[post,.]:*Piv
			/*subsequent periods*/
			for (t=2; t<=T; t++) {
				post0=(*post0_T[t])
				post=(*post_T[t])
				post_pi=(*post_pi_T[t])
				tmp=PI[post_pi,.]			
				for (c=1; c<=k; c++) {
					Qt_der[post,c]=quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c] :+ /*
								*/ quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post))):*dBPo[post,c]
				}			
			}
			/*last period*/
			temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
			pos3=selectindex( idt[.,1]:!=temp[.,1] )
			Qt_der=quadrowsum(Qt_der[pos3,.]):/pv
			S[.,((rnpar[1]+1)..rnpar[2])[x]]=Qt_der
		}

		/*gradient for al2*/
		dBPo=Po1:*Po2:*( ((Y[.,2] :- mu2):*(1:+lambda:*(c1 :- d1):*(c2 :- d2))) :+ /*
			*/ (lambda:*q:*mu2:*d2:*(c1:-d1)) ):*1
		for (x=1; x<=npar[3]; x++) {			
			Qt_der=J(rows(idt),k,0)
			/*initial period*/
			post=(*post_T[1])			
			Qt_der[post,x]=dBPo[post,x]:*Piv[.,x]
			/*subsequent periods*/
			for (t=2; t<=T; t++) {
				post0=(*post0_T[t])
				post=(*post_T[t])
				post_pi=(*post_pi_T[t])
				tmp=PI[post_pi,.]			
				for (c=1; c<=k; c++) {
					if (c==x) Qt_der[post,c]=quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c] :+ /*
								*/ quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post))):*dBPo[post,c]
					if (c!=x) Qt_der[post,c]= quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c]
				}				
			}
			/*last period*/
			temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
			pos3=selectindex( idt[.,1]:!=temp[.,1] )
			Qt_der=quadrowsum(Qt_der[pos3,.]):/pv
			S[.,((rnpar[2]+1)..rnpar[3])[x]]=Qt_der
		}			
			
		/*gradient for be2*/
		for (x=1; x<=npar[4]; x++) {
			dBPo=Po1:*Po2:*( ((Y[.,2] :- mu2):*(1:+lambda:*(c1 :- d1):*(c2 :- d2))) :+ /*
			*/ (lambda:*q:*mu2:*d2:*(c1:-d1)) ):*X2[.,x]			
			Qt_der=J(rows(idt),k,0)
			/*initial period*/
			post=(*post_T[1])			
			Qt_der[post,.]=dBPo[post,.]:*Piv
			/*subsequent periods*/
			for (t=2; t<=T; t++) {
				post0=(*post0_T[t])
				post=(*post_T[t])
				post_pi=(*post_pi_T[t])
				tmp=PI[post_pi,.]			
				for (c=1; c<=k; c++) {
					Qt_der[post,c]=quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c] :+ /*
								*/ quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post))):*dBPo[post,c]
				}			
			}
			/*last period*/
			temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
			pos3=selectindex( idt[.,1]:!=temp[.,1] )
			Qt_der=quadrowsum(Qt_der[pos3,.]):/pv
			S[.,((rnpar[3]+1)..rnpar[4])[x]]=Qt_der
		}		
				
		/*gradient for lambda*/
		dBPo=Po1:*Po2:*((c1 :- d1):*(c2 :- d2))		
		Qt_der=J(rows(idt),k,0)
		/*initial period*/
		post=(*post_T[1])			
		Qt_der[post,.]=dBPo[post,.]:*Piv		
		/*subsequent periods*/
		for (t=2; t<=T; t++) {
			post0=(*post0_T[t])
			post=(*post_T[t])
			post_pi=(*post_pi_T[t])
			tmp=PI[post_pi,.]			
			for (c=1; c<=k; c++) {
				Qt_der[post,c]=quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c] :+ /*
							*/ quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post))):*dBPo[post,c]
			}			
		}
		/*last period*/
		temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
		pos3=selectindex( idt[.,1]:!=temp[.,1] )
		Qt_der=quadrowsum(Qt_der[pos3,.]):/pv	
		S[.,((rnpar[4]+1)..rnpar[5])]=Qt_der	
		
		/*gradient for be4*/
		for (x=1; x<=npar[6]; x++) {
			Qt_der=J(rows(idt),k,0)
			/*initial period*/
			post=(*post_T[1])
			if (k==2)	c=max(rowsum((1..npar[6]):==x):*(2..k)')
			if (k>2)	c=max(rowsum((rowshape(1..npar[6],k-1):==x):*(2..k)'))
			dlPiv=J(rows(Piv),k,1):*(-Piv[.,c])
			dlPiv[.,c]=1:+dlPiv[.,c]
			Qt_der[post,.]=BPo[post,.]:*Piv:*dlPiv:*WXI[.,x]
			/*subsequent periods*/
			for (t=2; t<=T; t++) {
				post0=(*post0_T[t])
				post=(*post_T[t])
				post_pi=(*post_pi_T[t])
				tmp=PI[post_pi,.]
				
				for (c=1; c<=k; c++) {
					Qt_der[post,c]=quadrowsum(Qt_der[post0,.]:*rowshape(tmp[.,c]',rows(post))):*BPo[post,c]
				}
			}
			/*last period*/
			temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
			pos3=selectindex( idt[.,1]:!=temp[.,1] )
			Qt_der=quadrowsum(Qt_der[pos3,.]):/pv			
			S[.,((rnpar[5]+1)..rnpar[6])[x]]=Qt_der
		}		
		
		/*gradient for be5*/
		for (x=1; x<=npar[7]; x++) {
			cbar=selectindex(rowsum(rowshape(1..npar[7],k):==x ))
			cc=J(k,npar[7]/(k*(k-1)),0)
			cc[selectindex(((1..k):==cbar):==0),.]=rowshape(rowshape(1..npar[7],k)[cbar,.],k-1)
			cc=max((cc:==x):*(1..k)')		
			xx=selectindex(rowshape(1..npar[7],k)[cbar,.]:==x)			
			dQt=J(rows(idt),k,0)
			temp=(J(1,cols(idt),.)\idt); temp=temp[1..rows(idt),.]
			post=selectindex( idt[.,1]:==temp[.,1] )
			WXTb=J(rows(idt),cols(WXT),.)
			WXTb[post,.]=WXT
			/*subsequent periods*/
			for (t=2; t<=T; t++) {
				post0=(*post0_T[t])
				post=(*post_T[t])
				post_pi=(*post_pi_T[t])				
				tmp=PI[post_pi,.]				
				pi_cbar=tmp[selectindex(idtc[post_pi,3]:==cbar),.]				
				for (c=1; c<=k; c++) {
				    pi_c_cbar=rowshape(tmp[.,c]',rows(post))					
					dQt[post,c]=quadrowsum(dQt[post0,.]:*pi_c_cbar):*BPo[post,c]					
					if (c==cc)	dQt[post,c]=dQt[post,c]:+pi_cbar[.,c]:*(1:-pi_cbar[.,cc]):*Qt[post0,cbar]:*BPo[post,c]:*WXTb[post,xx]
					else		dQt[post,c]=dQt[post,c]:+pi_cbar[.,c]:*(  -pi_cbar[.,cc]):*Qt[post0,cbar]:*BPo[post,c]:*WXTb[post,xx]
				}				
			}			
			/*last period*/
			temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]	
			pos3=selectindex( idt[.,1]:!=temp[.,1] )			
			dQt=quadrowsum(dQt[pos3,.]):/pv			
			S[.,((rnpar[6]+1)..rnpar[7])[x]]=dQt			
		}		
	}/*end todo==1*/
		
}
end

*****************************************
capture mata: mata drop hmbivpois_rec_backward()
mata /*hmbivpois_rec_backward*/
function hmbivpois_rec_backward(all_id,pv,L,Phi,Piv,PI) {
	
	/*Baum–Welch forward recursion backward recursion
	--> input
	all_id	--> pointer matrix of individual index
	pv		--> manifest distribution
	L		--> recursive manifest distribution by i,t,k
	Phi		--> conditional response probabilities by i,t,k
	Piv		--> initilal prob.
	PI		--> transition prob.
	
	--> output
	U		--> posterior probabilities for t=1,...,T and i=1,...,n  and ubar,u=1,..., k
	V		--> posterior probabilities for t=1,...,T and i=1,...,n	
	*/
	
	/*preliminaries*/
	idt		=*all_id[1]
	idtc	=*all_id[2]
	idu=uniqrows(idt[.,1]);id=idt[.,1]
	n=rows(idu); T=max(idt[.,2]); k=max(idtc[.,3])
	
	idtw=J(n,T,.);
	for (i=1; i<=n; i++) {
		pos=selectindex(idt[.,1]:==i)
		temp=idt[pos,2];idtw[i,temp] = pos'	
	}
	
	idtcw=&.; 
	for(t=1; t<=T; t++) {
		tmp=J(n,k,.)
		tmp[selectindex(idtw[.,t]:!=.),.]=colshape((1::rows(idtc))[selectindex(idtc[.,2]:==t),.],k)
		idtcw=(idtcw,&(tmp:+0)); 
	}
	idtcw=idtcw[2..(T+1)]
	
	
	V=J(rows(idt),k,0);U=J(rows(PI),k,0)
	V1=V
	AA=J(n,k,1)
	M=J(n,k,1)
	MP=J(n,k,.)
	
	/*last period*/
	temp=(idt\J(1,cols(idt),-99)); temp=temp[2..rows(temp),.]
	pos0=selectindex( idt[.,1]:!=temp[.,1] )
	V1[pos0,.]=L[pos0,.]:/pv
	V[pos0,.]=V1[pos0,.]
	AA=Phi[pos0,.]:/pv
	
	temp=(idtc\J(k,cols(idtc),-99));temp=temp[(k+1)..rows(temp),.]	
	pos1=selectindex( idtc[.,1]:!=temp[.,1])
	temp=(idt\J(2,cols(idt),-99));temp=temp[3..rows(temp),.]
	pos2=selectindex(  idt[.,1]:!=temp[.,1] )
	pos2=colshape(pos2',2)[.,1]	
	U[pos1,.]=(vec(L[pos2,.]')#J(1,k,1)):*(AA#J(k,1,1))
	

	for (t=(T-1); t>1; t--) {/*backward recursion*/
		pos=.
		pos0=idtw[pos,t]
		pos1=idtw[pos,t+1]
		pos2=idtw[pos,t-1]
		pos3=(*idtcw[t])[pos,.]
		pos4=(*idtcw[t+1])[pos,.]
		pos=(pos0,pos1,pos2,pos3,pos4); pos=selectindex(rowmissing(pos):==0)	
		pos0=idtw[pos,t]
		pos1=idtw[pos,t+1]
		pos2=idtw[pos,t-1]
		pos3=(*idtcw[t])[pos,.]
		pos4=(*idtcw[t+1])[pos,.]
		
		MP[pos,.]=Phi[pos1,.]:*M[pos,.]	
		for (c=1; c<=k; c++) {
			M[pos,c]=rowsum(PI[pos4[.,c],.]:*MP[pos,.])
		}
		AA[pos,.]=Phi[pos0,.]:/pv[pos]:*M[pos,.]
		
		U[vec(pos3'),.]=(vec(L[pos2,.]')#J(1,k,1)):*(AA[pos,.]#J(k,1,1))
		V1[pos0,.]=L[pos0,.]:*M[pos,.]:/pv[pos]
		V[pos0,.]=V1[pos0,.]
	}

	/*first period*/
	temp=(J(1,cols(idt),-99)\idt); temp=temp[1..rows(idt),.]
	pos0=selectindex( idt[.,1]:!=temp[.,1])
	temp=(J(2,cols(idt),-99)\idt); temp=temp[1..rows(idt),.]
	pos1=selectindex( idt[.,1]:!=temp[.,1])
	pos1=colshape(pos1',2)[.,2]
	temp=(J(2*k,cols(idtc),-99)\idtc);temp=temp[1..rows(idtc),.]
	pos2=selectindex( idtc[.,1]:!=temp[.,1])
	pos2=colshape(pos2,k);pos2=pos2[((1..n):*2),.];

	MP=Phi[pos1,.]:*M
	for (c=1; c<=k; c++) {
		M[.,c]=rowsum(PI[pos2[.,c],.]:*MP)
	}
	V1[pos0,.]=L[pos0,.]:*M:/pv
	V[pos0,.]=V1[pos0,.]	

	U=U:*PI	
	
	return(&U,&V)
}
end

*****************************************
capture program drop hmbivpois
program define hmbivpois, eclass sortpreserve
	version 15
	local i = 1
	local paren "("
	local CMD "`0'"
	while "`paren'" == "(" {    
		gettoken e`i' 0:0, parse(" ,[") match(paren)
		local left "`0'"
		local junk: subinstr local e`i' ":" ":", count(local number)
		if "`number'" == "1" {
			gettoken dep`i'n  e`i':  e`i', parse(":")
			gettoken junk  e`i':  e`i', parse(":")
		}
		local  e`i' : subinstr local  e`i' "=" " "
		gettoken dep`i' 0:  e`i', parse(" ,[") 
		unab dep`i': `dep`i'' 

		/* collecting together -- for e.g. marking sample */
		local deps "`deps' `dep`i''"
		confirm variable `dep`i''
		if "`dep`i'n'" == "" {
			local dep`i'n "`dep`i''"
		}
		syntax [varlist(default=none)] [, nadd ]
		local ind`i' `varlist'
		local  0 "`left'"
				/* reset `paren' to empty when thru eqns */

		gettoken check : left, parse(" ,[") match(paren)
		local i = `i' + 1
	}
	local 0 "`left'"
	syntax [if] ,[ k(integer 1) id(varname) time(varname) INITcov(varlist) xtrans(varlist) /*
	*/ HETerogeneous HOMogeneous restr INDPOISson/*
	*/ maxiter(integer 3000) em(integer 10) emonly /*
	*/ NOCONStant tracep traceg traceh /*
	*/ random poisson STARTfrom(name) prev(string) /*
	*/ load(string) save(string) EVALuate /*
	*/ optimize_init_technique(string) optimize_init_evaluatortype(integer 1)] 
	 
	if "`optimize_init_technique'" == ""	local optimize_init_technique = "bfgs"
	
	marksample touse
	markout `touse'
	
	tempvar _set_sample_HM _`id'_HM _`time'_HM _`time'1_HM _info_sample_HM _`id'1_HM
	
	if "`noconstant'" != "" & "`initcov'" == "" {
	    dis as err "You cannot set noconstant option without adding covariates in the measurement equation "
		exit(99)		
	}
	if `k' == 1 & "`initcov'" != "" {
		dis as txt "Option " as input "initcov " as txt "ignored since " as input "k = `k'"
	}
	if `k' > 1 {
		if "`noconstant'" != "" & "`initcov'" != ""   {
			local constant = 0
		}
		else {
			local constant = 1		
		}
	}
	else local constant = 0
	foreach o in tracep traceg traceh {
		if "``o''" != ""	{
			local `o' = "on"
		}
		else {
			local `o' = "off"
		}
	}
	**************
	/*set sample*/
	**************
	qui gen `_set_sample_HM' = 1 if `touse'
	local to_be_dropped = "`dep1n' `dep2n'"
	if  wordcount("`ind1'") > 0 {
		local to_be_dropped = "`to_be_dropped' `ind1'"
	}
	if  wordcount("`ind2'") > 0 {
		local to_be_dropped = "`to_be_dropped' `ind2'"
	}
	if "`initcov'" != "" {
	    local to_be_dropped = "`to_be_dropped' `initcov'"
	}
	if "`xtrans'" != "" {
	    local to_be_dropped = "`to_be_dropped' `xtrans'"
	}
	qui for var `to_be_dropped' : replace `_set_sample_HM' = . if X == .
	qui gen `_`id'1_HM'= `id'*`_set_sample_HM'
	qui egen `_`id'_HM' = group(`_`id'1_HM')
	qui gen `_`time'1_HM' = `time'*`_set_sample_HM'
	qui egen `_`time'_HM' = group(`_`time'1_HM')
	qui xtset `_`id'_HM' `_`time'_HM'
	local T = (r(tmax)-r(tmin))+1
	if r(gaps) == 1 {
		local bal = r(balanced)
		dis as err "panel is `bal', but with a gap"
		exit(198)
	}
	qui bys `_`id'_HM' : gen `_info_sample_HM' = _N if `_set_sample_HM' != .
	qui su `_info_sample_HM' if `_set_sample_HM' != .
	local g_min = r(min) 
	local g_max = r(max) 
	local g_avg = r(mean) 
	qui bys `_`id'_HM' : replace `_info_sample_HM' = _n if `_set_sample_HM' != .
	qui su `_info_sample_HM' if `_set_sample_HM' != . & `_info_sample_HM' == 1
	local n = r(N)
	qui sort `_set_sample_HM' `_`id'_HM' `_`time'_HM'
	qui count if `_set_sample_HM' == 1
	local N = r(N)
	mata: id 		= st_data((1..`N')', "`_`id'_HM'")
	mata: time 		= st_data((1..`N')', "`_`time'_HM'")
	mata: id_o		= st_data((1..`N')', "`id'")
	mata: time_o		= st_data((1..`N')', "`time'")
	mata: idtc		= (id,time)#J(`k',1,1),J(`N',1,1)#(1..`k')'
	mata: Y1 		= st_data((1..`N')', "`dep1n'")
	mata: Y2 		= st_data((1..`N')', "`dep2n'")	
	if  wordcount("`ind1'") > 0 {
		mata: X1 = st_data((1..`N')', "`ind1'")
	}
	else {
		mata: X1 = J(0,0,.)
	}
	if  wordcount("`ind2'") > 0 {
		mata: X2 = st_data((1..`N')', "`ind2'")
	}
	else {
		mata: X2 = J(0,0,.)
	}
	if "`initcov'" != "" & `k' > 1 {
		mata: WXI	= st_data((1..`N')', "`initcov'")
		mata: WXI=WXI[selectindex((.\id)[1..rows(id)]:!=id),.]
		if `constant' == 1	mata: WXI	= J(rows(WXI),1,1),WXI
		mata: WXI=J(1,`k'-1,1)#WXI
	}
	else {
		mata: WXI = J(rows(selectindex((.\id)[1..rows(id)]:!=id)),`k'-1,1)
	}
	if `k' > 1 {
	    mata: WXT = J(rows(selectindex((.\id)[1..rows(id)]:==id)),`k'-1,1)
		if "`xtrans'" != ""  {
			mata: WXT	= st_data((1..`N')', "`xtrans'")
			mata: WXT=WXT[selectindex((.\id)[1..rows(id)]:==id),.]
			mata: WXT=J(rows(id)-rows(uniqrows(id)),1,1),WXT
			mata: WXT=J(1,`k'-1,1)#WXT
		}
		if "`heterogeneous'" != "" {
			mata: WXT=time[selectindex((.\id)[1..rows(id)]:==id)]:==(uniqrows(time)'#J(rows(id)-rows(uniqrows(id)),1,1))
			mata: WXT=WXT[.,2..cols(WXT)]
			mata: WXT=J(1,`k'-1,1)#WXT
		}
	}
	else {
		mata: WXT = J(rows(Y1),`k'-1,1)
	}
	/* indexing */
	mata: index=hmbivpois_index((id,time),idtc,`k') /*(&post_T,&post0_T,&post_pi_T)*/
	
	/* set starting parameters */
	if "`random'" != "" & "`poisson'" == ""  & "`load'" == "" & "`prev'" == ""{	    
		mata: hmbiv_be1=-J(cols(cov1),1,0);
		mata: hmbiv_be2=J(cols(cov2),1,0);
		mata: hmbiv_al1=runiform(`k',1); fmmbiv_al2=runiform(`k',1);
		mata: hmbiv_lambda=J(1,1,0)
		/*initial prob.*/
		local lI = wordcount("`initcov'")+`constant'
		mata: hmbiv_be3 = J((`k'-1)*(`lI'),1,0)
		/*transition prob.*/
		local lT = 1
		if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
		if "`heterogeneous'" != ""		local lT = (`T'-1)
		mata: hmbiv_be4 = J((`k'-1)*(`lT')*(`k'),1,0); tmp=J(`lT',1,0); tmp[1]=-log(10); hmbiv_be4=hmbiv_be4+J((`k'-1)*(`k'),1,1)#tmp
	}
	else if "`random'" == "" & "`poisson'" != "" & "`load'" == "" & "`prev'" == ""{
		dis ""
		dis as text "Univariate Poisson regression "
		dis ""		    
	    qui poisson `dep1n' `ind1' if `touse'
		local lcomp = e(ll)
		local al1 = _b[`dep1n':_cons]
		mata: hmbiv_al1 = J(`k',1,`al1'):+((0..(`k'-1)):/10)'
		if wordcount("`ind1'") > 0 {
			local nxlist = wordcount("`ind1'") 
			mata: hmbiv_be1=J(`nxlist',1,0)
			local i = 1
			foreach x in `ind1' {
				local be1 = _b[`dep1n':`x']
				mata: hmbiv_be1[`i'] = `be1'
				local ++i
			}
		}		
		qui poisson `dep2n' `ind2' if `touse'
		local lcomp = e(ll) + `lcomp'
		local al2 = _b[`dep2n':_cons]
		mata: hmbiv_al2 =  J(`k',1,`al2'):+((0..(`k'-1)):/10)'
		if wordcount("`ind2'") > 0 {
			local nxlist = wordcount("`ind2'") 
			mata: hmbiv_be2=J(`nxlist',1,0)
			local i = 1
			foreach x in `ind2' {
				local be2 = _b[`dep2n':`x']
				mata: hmbiv_be2[`i'] = `be2'
				local ++i
			}
		}
		mata: hmbiv_lambda = J(1,1,0)
		/*initial prob.*/
		local lI = wordcount("`initcov'")+`constant'
		mata: hmbiv_be3 = J((`k'-1)*(`lI'),1,0)
		/*transition prob.*/
		local lT = 1
		if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
		if "`heterogeneous'" != ""		local lT = (`T'-1)
		mata: hmbiv_be4 = J((`k'-1)*(`lT')*(`k'),1,0); tmp=J(`lT',1,0); tmp[1]=-log(10); hmbiv_be4=hmbiv_be4+J((`k'-1)*(`k'),1,1)#tmp		
		if "`restr'" != ""	mata: hmbiv_be4 = J(0,1,0)
		dis as text "Comparison:    log likelihood = " as res `lcomp'	
		dis ""
	}
	else if "`random'" == "" & "`poisson'" == "" & "`load'" != "" & "`prev'" == ""{
		mata: fh = fopen("`load'.bin", "r")
		mata: hmbiv_al1 = fgetmatrix(fh)'
		mata: hmbiv_be1 = fgetmatrix(fh)'
		mata: hmbiv_al2 = fgetmatrix(fh)'
		mata: hmbiv_be2 = fgetmatrix(fh)'
		mata: hmbiv_lambda = fgetmatrix(fh)
		mata: hmbiv_be3 = fgetmatrix(fh)'
		mata: hmbiv_be4 = fgetmatrix(fh)'
		mata: hmbiv_varcov = fgetmatrix(fh)	
		mata: hmbiv_output = fgetmatrix(fh)	
		mata: hmbiv_U = fgetmatrix(fh)	
		mata: hmbiv_PI = fgetmatrix(fh)	
		mata: hmbiv_we_t = fgetmatrix(fh)		
		mata: fclose(fh)
		mata: st_local("canev", strofreal(hmbiv_varcov[1,1]) )
		if "`evaluate'" != "" dis as text "Evaluate log-likelihood from " as res "`load'.bin"
		else dis as text "Load starting parameters from matrix in the file " as res "`load'.bin"
		if "`canev'" == "." & "`evaluate'" != "" {
			dis as text "Variance-covariance is not provided by " as res "`load'.bin"
			dis as text "Option evaluate ignored" 
		}
		if "`canev'" != "." & "`evaluate'" == "" {
			dis as text "Variance-covariance matrix provided, but ingored"
			mata: hmbiv_varcov = .
		}
		if "`restr'" != ""		mata: hmbiv_be4=J(0,0,.)
	}
	else if "`random'" == "" & "`poisson'" == "" & "`load'" == "" & "`prev'" != ""{
		mata: fh = fopen("`prev'.bin", "r")
		mata: hmbiv_al1 = fgetmatrix(fh)'
		mata: hmbiv_be1 = fgetmatrix(fh)'
		mata: hmbiv_al2 = fgetmatrix(fh)'
		mata: hmbiv_be2 = fgetmatrix(fh)'
		mata: hmbiv_lambda = fgetmatrix(fh)
		mata: hmbiv_be3 = fgetmatrix(fh)'
		mata: hmbiv_be4 = fgetmatrix(fh)'
		mata: hmbiv_varcov = fgetmatrix(fh)	
		mata: hmbiv_output = fgetmatrix(fh)	
		mata: hmbiv_U = fgetmatrix(fh)	
		mata: hmbiv_PI = fgetmatrix(fh)	
		mata: hmbiv_we_t = fgetmatrix(fh)		
		mata: fclose(fh)
		local k1=`k'-1
		dis ""
		dis as text "Load starting values from the file " as res "`prev'.bin" as text " with " as res "`k1' " as text "latent classes"
		dis ""
		mata: hmbiv_varcov = .		
		mata: hmbiv_al1 = (hmbiv_al1\(hmbiv_al1[`k1']+0.5))
		mata: hmbiv_al2 = (hmbiv_al2\(hmbiv_al2[`k1']+0.5))
		mata: hmbiv_be3 = (hmbiv_be3\mean(rowshape(hmbiv_be3,`k'-2))')
		mata: hmbiv_be4 = rowshape(hmbiv_be4,`k1')',mean(rowshape(hmbiv_be4,`k1'))'
		mata: hmbiv_be4 = hmbiv_be4\mean(hmbiv_be4)
		mata: hmbiv_be4 = vec(hmbiv_be4)
		capture mata: mata drop hmbiv_varcov
		capture mata: mata drop hmbiv_output 
		capture mata: mata drop hmbiv_U 
		capture mata: mata drop hmbiv_PI 
		capture mata: mata drop hmbiv_we_t 
		if "`restr'" != ""		mata: hmbiv_be4=J(0,0,.)
	}	
	else if "`startfrom'" != "" {
		dis ""
		dis as text "Load starting values from the following matrix " as input "`startfrom'"
		dis ""
		mata: start_par = st_matrix("`startfrom'")'
		mata: hmbiv_al1 = start_par[1..`k']; start_par = start_par[(`k'+1)..rows(start_par)]
		if wordcount("`ind1'") > 0 {
			local nxlist = wordcount("`ind1'")
			mata: hmbiv_be1 = start_par[1..`nxlist']; start_par = start_par[(`nxlist'+1)..rows(start_par)]
		}
		mata: hmbiv_al2 = start_par[1..`k']; start_par = start_par[(`k'+1)..rows(start_par)]
		if wordcount("`ind2'") > 0 {
			local nxlist = wordcount("`ind2'")
			mata: hmbiv_be2 = start_par[1..`nxlist']; start_par = start_par[(`nxlist'+1)..rows(start_par)]
		}
		mata: hmbiv_lambda = start_par[1]; start_par = start_par[(1+1)..rows(start_par)]
		/*initial prob.*/
		local lI = wordcount("`initcov'")+`constant'
		local nxlist = (`k'-1)*(`lI')
		mata: hmbiv_be3 = start_par[1..`nxlist']; start_par = start_par[(`nxlist'+1)..rows(start_par)]
		/*transition prob.*/
		local lT = 1
		if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
		if "`heterogeneous'" != ""		local lT = (`T'-1)
		local nxlist = (`k'-1)*(`lT')*(`k')
		mata: hmbiv_be4 = start_par[1..`nxlist']
	}
	else { /* starting paramaters from a pooled Famoye biv. Poisson*/
		dis ""
		dis as text "Famoye bivariate Poisson regression "
		dis ""	
		qui bivcnto (`dep1n' = `ind1') (`dep2n' = `ind2'), pfamoye
		local lcomp = e(ll)
		local al1 = _b[`dep1n':_cons]
		mata: hmbiv_al1 = J(`k',1,`al1'):+((0..(`k'-1)):/10)'
		local al2 = _b[`dep2n':_cons]
		mata: hmbiv_al2 =  J(`k',1,`al2'):+((0..(`k'-1)):/10)'
		if wordcount("`ind1'") > 0 {
			local nxlist = wordcount("`ind1'") 
			mata: hmbiv_be1=J(`nxlist',1,0)
			local i = 1
			foreach x in `ind1' {
				local be1 = _b[`dep1n':`x']
				mata: hmbiv_be1[`i'] = `be1'
				local ++i
			}
		}
		if wordcount("`ind2'") > 0 {
			local nxlist = wordcount("`ind2'") 
			mata: hmbiv_be2=J(`nxlist',1,0)
			local i = 1
			foreach x in `ind2' {
				local be2 = _b[`dep2n':`x']
				mata: hmbiv_be2[`i'] = `be2'
				local ++i
			}
		}
		local lambda = _b[/:lambda]
		mata: hmbiv_lambda = `lambda'
		/*initial prob.*/
		local lI = wordcount("`initcov'")+`constant'
		mata: hmbiv_be3 = J((`k'-1)*(`lI'),1,0)
		/*transition prob.*/
		local lT=1
		if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
		if "`heterogeneous'" != ""		local lT = (`T'-1)
		mata: hmbiv_be4 = J((`k'-1)*(`lT')*(`k'),1,0); tmp=J(`lT',1,0); tmp[1]=-log(10); hmbiv_be4=hmbiv_be4+J((`k'-1)*(`k'),1,1)#tmp		
		dis as text "Comparison:    log likelihood = " as res `lcomp'
		dis ""
	}
	if  wordcount("`ind1'") == 0 {
		mata: hmbiv_be1 = J(0,0,.)
	}
	if  wordcount("`ind2'") == 0 {
		mata: hmbiv_be2 = J(0,0,.)
	}
	if `k' == 1 {
		if "`random'" == "" & "`poisson'" == "" & "`startfrom'" == "" {
			*local maxiter = 0
			local optimize_init_technique = "nr"
			*mata: hmbiv_al1=hmbiv_al1+(runiform(1,1)/10);hmbiv_al2=hmbiv_al2+(runiform(1,1)/10)
		}
		mata: hmbiv_be3=J(0,0,.)
		mata: hmbiv_be4=J(0,0,.)
	}
	if "`restr'" != "" mata: hmbiv_be4=J(0,0,.)
	if "`indpoisson'" != "" mata: hmbiv_lambda=J(0,0,.)
	
	mata: npar=(rows(hmbiv_al1),rows(hmbiv_be1),rows(hmbiv_al2),rows(hmbiv_be2),rows(hmbiv_lambda),rows(hmbiv_be3),rows(hmbiv_be4))
	mata: start_par=(hmbiv_al1')
	if  wordcount("`ind1'") > 0 {
		mata: start_par=(start_par,hmbiv_be1')
	}	
	mata: start_par=(start_par,hmbiv_al2')
	if  wordcount("`ind2'") > 0 {
		mata: start_par=(start_par,hmbiv_be2')
	}
	if "`indpoisson'" == "" mata: start_par=(start_par,hmbiv_lambda)
		
	if `k' > 1	{
		mata: start_par=(start_par,hmbiv_be3')
		if "`restr'" == ""	mata: start_par=(start_par,hmbiv_be4')
	}
	
	*******************************************************************************
	/* set name of matrix to be displayed at each iteration if trace is required */
	*******************************************************************************
	local eqname = ""
	local name_be = ""
	foreach y of numlist 1/2 {
		foreach c of numlist 1/`k' {
		    local eqname = "`eqname' `dep`y'n'"
			local name_be = "`name_be' alpha`y'_`c'"
		}
		if  wordcount("`ind`y''") > 0 {
			foreach x in `ind`y'' {
				local eqname = "`eqname' `dep`y'n'"
				local name_be = "`name_be' `x'"
			}	
		}
	}
	if "`indpoisson'" == "" {
	    local eqname = "`eqname' lambda"
		local name_be = "`name_be' lambda"
	}
	local name_be3 = ""
	if `k'>1 & `constant' == 1 {
		foreach c of numlist 2/`k' {
			foreach x in alpha `initcov' {
			    local eqname = "`eqname' Init.cond."
				local name_be3 = "`name_be3' `x'_`c'"
			}
		}
	}
	else if `k'>1 & `constant' == 0 {
		foreach c of numlist 2/`k' {
			foreach x in `initcov' {
			    local eqname = "`eqname' Init.cond."
				local name_be3 = "`name_be3' `x'_`c'"
			}
		}
	}	
	else {
	    if `k'>1 {
			foreach c of numlist `k'/`k' {
				foreach x in alpha `initcov' {
					local eqname = "`eqname' Init.cond."
					local name_be3 = "`name_be3' `x'_`c'"
				}
			}
		}
	}
	
	local name_be4 = ""
	if `k'>1 & "`restr'" == "" {
		if "`xtrans'" != "" {	
	    	
			foreach c of numlist 1/`k' {
				foreach cc of numlist 1/`k' {
					if `c'!=`cc' {
						foreach x in alpha `xtrans' {
							local eqname = "`eqname' `cc'vs`c'"
							local name_be4 = "`name_be4' `x'"
						}
					}
				}
			}
		}
		else if "`heterogeneous'" != "" {
			foreach c of numlist 1/`k' {
				foreach cc of numlist 1/`k' {
					if `c'!=`cc' {
						foreach t of numlist 2/`T' {
							local eqname = "`eqname' `cc'vs`c'"
							local name_be4 = "`name_be4' time_`t'"
						}
					}
				}
			}
		}
		else {
			foreach c of numlist 1/`k' {
				foreach cc of numlist 1/`k' {
					if `c'!=`cc' {
						foreach x in alpha {
							local eqname = "`eqname' `cc'vs`c'"
							local name_be4 = "`name_be4' `x'`cc'_`c'"
						}
					}
				}
			}
		}
	}
	
	tempname start_par
	mata: st_matrix("`start_par'", start_par)
	mat colnames `start_par' = `name_be' `name_be3' `name_be4'
	mat coleq    `start_par' = `eqname'
	if "`startfrom'" != "" mat list `start_par'
	
	if "`restr'" != ""	local restr = 1
	else				local restr = 0	

	if "`evaluate'" == ""	mata: hmbiv_output=(&.,&.,&.,&.,&.,&.,&.)
	
	/*
	*************
	/* log-lik */
	*************
	local restr = 0	
	mata: LL=.;S=.;H=.;hmbiv_output=(&.,&.,&.,&.,&.,&.,&.)
	mata: bivpoisson_HM(1, start_par,  /*
		*/   npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT),(&`restr'),hmbiv_output, /*
		*/ LL, S, H)
	*mata: LL
	*mata: hmbiv_U=*hmbivpois_rec_backward((&(id,time),&idtc),(*hmbiv_output[6]),(*hmbiv_output[7]),(*hmbiv_output[5]),(*hmbiv_output[1]),(*hmbiv_output[2]))[2]
exit(55)
	****************
	/* derivative */
	****************	
	mata: d0=.;d1=.;d2=.;
	mata: bivpoisson_der_HM(start_par,npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT), /*
						*/  (&`restr'),hmbiv_output,d0,d1,d2)
	mata: d0,quadsum(LL)
	mata: d1',quadcolsum(S)'
	exit(1)
	*/

	****************
	/* estimation */
	****************
	if "`evaluate'" != ""{
		mata: th=(hmbiv_al1)
		mata: if (npar[2]>0) {; th=th\hmbiv_be1; };
		mata: th=th\hmbiv_al2
		mata: if (npar[4]>0) {; th=th\hmbiv_be1; };
		mata: th=th\hmbiv_lambda
		mata: if (npar[6]>0) {; th=th\hmbiv_be3; };
		mata: if (npar[7]>0) {; th=th\hmbiv_be4; };
		mata: th=th'
		mata: V=hmbiv_varcov;LL=quadsum(log((*hmbiv_output[6])))
	}
	else {
		if "`emonly'" != ""{
			mata: th=.;V=.;LL=.;
			/*EM-step*/
			if `em' < 1000 {
				local em = 1000
			}
			mata disparamEM=0
			if "`tracep'" != "" {
				mata:disparamEM=st_matrixcolstripe("`start_par'")
			}
			mata: bivpoisson_est_EM_HM(start_par,npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT),(&`restr',&`em',&disparamEM),hmbiv_output,th,LL)
			mata: start_par=th'
			qui mata: bivpoisson_est_HM("gf`optimize_init_evaluatortype'",start_par,npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT), /*
			*/  (&`restr'),hmbiv_output,"`optimize_init_technique'",0,"`tracep'","`traceg'","`traceh'",st_matrixcolstripe("`start_par'"),th,V,LL)
		}
		else {
			mata: th=.;V=.;LL=.;
			mata disparamEM=0
			if "`tracep'" != "" {
				mata:disparamEM=st_matrixcolstripe("`start_par'")
			}
			/*EM-step*/
			mata: bivpoisson_est_EM_HM(start_par,npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT),(&`restr',&`em',&disparamEM),hmbiv_output,th,LL)
			mata: start_par=th'
			mata: printf(" \n")
			mata: printf("Maximized Log-likelihood \n")
			/*Direct maximization*/
			mata: bivpoisson_est_HM("gf`optimize_init_evaluatortype'",start_par,npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT), /*
			*/  (&`restr'),hmbiv_output,"`optimize_init_technique'",`maxiter',"`tracep'","`traceg'","`traceh'",st_matrixcolstripe("`start_par'"),th,V,LL)
		}
		dis "Estimation completed"
	}

	*************************************
	/* collect some additional outputs */
	*************************************
	if "`evaluate'" == "" {
		mata: hmbiv_output=(&.,&.,&.,&.,&.,&.,&.)
		mata: bivpoisson_HM(0, th,  /*
			*/   npar,(&(id,time),&idtc,&index),(Y1,Y2),(&X1,&X2,&WXI,&WXT),(&`restr'), hmbiv_output, /*
			*/ .,.,.)
		mata: hmbiv_PI=J(`k',`k',0)
		foreach c of numlist 1/`k' {
			mata: hmbiv_PI[`c',.]=mean((*hmbiv_output[2])[selectindex(idtc[.,3]:==`c'),.])
		}
		mata: hmbiv_U=*hmbivpois_rec_backward((&(id,time),&idtc),(*hmbiv_output[6]),(*hmbiv_output[7]),(*hmbiv_output[5]),(*hmbiv_output[1]),(*hmbiv_output[2]))[2]
		mata: hmbiv_we_t = J(`T',`k',.)
		foreach t of numlist 1/`T' {	
			mata: hmbiv_we_t[`t',.]=mean(hmbiv_U[selectindex(time:==`t'),.])
		}
	}

	tempname b V ll BIC AIC npar we we_t PI
	
	mata: st_matrix("`b'", th)
	mat colnames `b' = `name_be' `name_be3' `name_be4'
	mata: st_matrix("`V'", V)
	matrix colnames `V' = `name_be' `name_be3' `name_be4'
	matrix rownames `V' = `name_be' `name_be3' `name_be4'
	matrix coleq `b' = `eqname'
	matrix coleq `V' = `eqname'
	matrix roweq `V' = `eqname'	
	ereturn post `b' `V', obs(`N')					/* estimated vector of paramaters and var/cov matrix */					
	ereturn scalar k = `k'							/* number of latent states */
	mata: st_numscalar("`ll'", LL)
	ereturn scalar ll = `ll'						/* maximized log-lik */	
	mata: st_numscalar("`BIC'", (log(`N')*sum(npar) - 2*(LL)))
	ereturn scalar BIC = `BIC'						/* BIC */	
	mata: st_numscalar("`AIC'", (- 2*(LL) + 2*sum(npar)))
	ereturn scalar AIC = `AIC'						/* AIC */	
	mata: st_numscalar("`npar'", sum(npar))
	ereturn scalar npar = `npar'					/* npar */		
	ereturn scalar g_min = `g_min'					/* min. length time */
	ereturn scalar g_avg = `g_avg'					/* average length time */
	ereturn scalar g_max = `g_max'					/* max. length time */	
		
	mata: st_matrix("`we'", mean(hmbiv_U))	
	ereturn matrix we = `we'						/* we */
	mata: st_matrix("`we_t'", hmbiv_we_t)	
	ereturn matrix we_t = `we_t'						/* we */	
	mata: st_matrix("`PI'", hmbiv_PI)
	ereturn matrix PI = `PI'						/* Average transition prob. */
	
	ereturn local cmd = "hmbivpois"
	ereturn local depvar = "`dep1n' `dep2n'"
	
	****************
	/*Save results*/
	****************
	if "`save'" != "" {
		mata: rnpar=runningsum(npar)
		mata: if (min(npar)==0) { ; /*
			*/ rnpar=J(1,length(npar),0) ; /*
			*/ for (x=1; x<=length(npar); x++) {  ; /*
				*/ if (npar[x]> 0) rnpar[x]=sum(npar[1..x]) ; /*
				*/ if (npar[x]==0) rnpar[x]=rnpar[x-1] ; /*
			*/ } ; /*
		*/ } ;
		mata: 								           hmbiv_al1=th[1..rnpar[1]]		
		mata: hmbiv_be1=J(0,0,.);    if (npar[2]>0) {; hmbiv_be1=th[(rnpar[1]+1)..rnpar[2]]; };
		mata: 								           hmbiv_al2=th[(rnpar[2]+1)..rnpar[3]]
		mata: hmbiv_be2=J(0,0,.);    if (npar[4]>0) {; hmbiv_be2=th[(rnpar[3]+1)..rnpar[4]]; };
		mata: hmbiv_lambda=J(0,0,.); if (npar[5]>0) {; hmbiv_lambda=th[(rnpar[4]+1)..rnpar[5]]; };
		mata: hmbiv_be3=J(0,0,.);    if (npar[6]>0) {; hmbiv_be3=th[(rnpar[5]+1)..rnpar[6]]; };
		mata: hmbiv_be4=J(0,0,.);    if (npar[7]>0) {; hmbiv_be4=th[(rnpar[6]+1)..rnpar[7]]; };
		
		mata: fh = fopen("`save'.bin", "rw")
		mata: fputmatrix(fh, hmbiv_al1)
		mata: fputmatrix(fh, hmbiv_be1)
		mata: fputmatrix(fh, hmbiv_al2)
		mata: fputmatrix(fh, hmbiv_be2)
		mata: fputmatrix(fh, hmbiv_lambda)
		mata: fputmatrix(fh, hmbiv_be3)
		mata: fputmatrix(fh, hmbiv_be4)
		mata: fputmatrix(fh, V)
		mata: fputmatrix(fh, hmbiv_output)
		mata: fputmatrix(fh, hmbiv_U)
		mata: fputmatrix(fh, hmbiv_PI)		
		mata: fputmatrix(fh, hmbiv_we_t)
		mata: fclose(fh)
	}
	
	*******************
	/*Display results*/
	*******************
	#delimit ;
	di _n in gr "Hidden Markov Poisson Famoye  " 
		_col(49) in gr "Number of obs" _col(67) "=" 
		_col(69) in ye %10.0fc `N' ;
	di in gr "bivariate Poisson regression"  in gr
		_col(49) "Number of groups" _col(67) "="
		_col(69) in ye %10.0fc `n' _n ;		
	di in gr "Group variable: " in ye abbrev("`id'",12) _col(49) as txt "Obs per group:" ;
	di as txt "Latent states  = " as res %-10.0fc `k'
		_col(63) as txt "min" _col(67) "="
		_col(69) as res %10.0fc `g_min' ;
	di as txt "      " 
		_col(63) as txt "avg" _col(67) "="
		_col(69) as res %10.1fc `g_avg' ;
	di as txt "      " 
		_col(63) as txt "max" _col(67) "="
		_col(69) as res %10.0fc `g_max' _n ;
	di as txt "likelihood = " as res %10.3f `ll' /*10.3*/
		_col(62) as txt "BIC" _col(67) "="
		_col(69) as res %10.3fc `BIC'  ;	
	#delimit cr

	ereturn display
	ereturn repost, esample(`touse')	
end






/*
hmbivpois 	(dep1 = cov11 cov12 Ldep1 Ldep2 ) /*
	*/		(dep2 = cov21 cov22 Ldep1 Ldep2 )  /*
	*/,  k(1) id(id) time(year)    /*
	*/   initcov(Ldep1 Ldep2)   /*
	*/ optimize_init_technique(bfgs 30 nr 1000) save(pp) load(pp) eval
	
*/




/*
hmbivpois 	(viopr   =   )  /*
	*/		(peacepr =  $year )  /*
	*/,  k(2) id(cname1) time(yearorig) optimize_init_technique(bfgs 10 nr 1000) /*
	*/ initcov(Lviopr  Lpeacepr) 
*/

/*
hmbivpois 	(viopr   =   )  /*
	*/		(peacepr =  )  /*
	*/,  k(2) id(cname1) time(yearorig) optimize_init_technique(bfgs 10 nr 1000) /*
	*/ initcov(Lviopr  Lpeacepr)  xtrans(loggdp1) indpois
*/
	
	
/*
hmbivpois 	(totae = $year $indcov Ltotae Ltotgp ) /*
	*/		(totgp = $year $indcov Ltotae Ltotgp )  /*
	*/,  k(7) id(id) time(year)  em(20) /*
	*/  save(results/tempem_k7) initcov(Ltotae Ltotgp) prev(results/tempem_k6)
*/