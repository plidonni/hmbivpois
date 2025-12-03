clear mata

*******************************************
capture mata: mata drop idt_index()
mata /*idt_index*/
function idt_index(idt,idtc,k) {
    
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

*******************************************
capture mata: mata drop setintercept()
mata /*setintercept*/
function setintercept(k) {
	X0=J(1,k,NULL)
	for (v=1; v<=k; v++) { 
		du=J(1,k,0); du[v]=1
		X0[v]=&du[.,.]
	}

	return(X0)
}
end

*******************************************
capture mata: mata drop prlogit()
mata /*prlogit*/
function prlogit(X,b,lev,|base) {	
	/*
	 comute probs. of multinomial logit
	--> input
	X		--> design matrix n x k 
	be		--> col vector of initial parameters k x 1
    lev     --> number of levels
	--> output
	Prob	-> matrix of estimated probs. n x lev
	*/
	
    /*preliminaries*/
    np=cols(rowshape(b,lev-1)); /*number of parameters in equqation*/
    /*compute probabilities*/
	Eta=X*(b:*(I(lev-1)#J(np,1,1)))
	Eta=(J(rows(Eta),1,0),Eta)
    Prob=exp(Eta):/quadrowsum(exp(Eta));
	if (args()==4) {
		temp=Prob
		Prob[.,selectindex((1..lev):==base)]=temp[.,1]
		Prob[.,selectindex((1..lev):!=base)]=temp[.,2..lev]
	}	
	return(Prob)
}
end

*******************************************
capture mata: mata drop poissonfam_d()
mata: /*poissonfam_d*/
function poissonfam_d(Y1,Y2,X1,X2,be1,be2,lambda,BPotype) {
	
	if (BPotype=="normal") {
		theta = (exp(2*lambda)-1)/(exp(2*lambda)+1)
		mu1=exp(X1*be1');mu2=exp(X2*be2')
		n=rows(Y1)
		
		f10 = poisson(mu1, Y1)
		f11 = J(n,1,0); pos=selectindex(Y1:!=0); if (length(pos)>0) f11[pos]= poisson(mu1[pos], Y1[pos]:-1)
		f20 = poisson(mu2, Y2)
		f21 = J(n,1,0); pos=selectindex(Y2:!=0); if (length(pos)>0) f21[pos]= poisson(mu2[pos], Y2[pos]:-1)
		
		/*Calculate the Normal Copula function*/
		t1 =  binormal(invnormal(f10),invnormal(f20),theta)
		t2 =  binormal(invnormal(f11),invnormal(f20),theta)
		t3 =  binormal(invnormal(f10),invnormal(f21),theta)
		t4 =  binormal(invnormal(f11),invnormal(f21),theta)
		
		/* Set up for reasonable bounds of probabilities*/
		minv = J(n,1,1e-7) /*1e-7 -- 1e-6*/
		maxv = 1.0:-minv		
		
		/* Impose Frechet-Hoeffding bounds*/
		t1 = rowmin((rowmax((t1, rowmax((f10+f20:-1,minv)))),rowmin((f10, f20,maxv))))
		t2 = rowmin((rowmax((t2, rowmax((f11+f20:-1,minv)))),rowmin((f11, f20,maxv))))
		t3 = rowmin((rowmax((t3, rowmax((f10+f21:-1,minv)))),rowmin((f10, f21,maxv))))
		t4 = rowmin((rowmax((t4, rowmax((f11+f21:-1,minv)))),rowmin((f11, f21,maxv))))
		
		Phi = rowmin((maxv,rowmax((minv,t1-t2-t3+t4))))		
	}
	if (BPotype=="frank") {
		theta=lambda
		if (theta==0) theta=1e-6
		eta = 1 - exp(-theta)
		mu1=exp(X1*be1');mu2=exp(X2*be2')
		n=rows(Y1)
	
		f10 = poisson(mu1, Y1)
		f11 = J(n,1,0); pos=selectindex(Y1:!=0); if (length(pos)>0) f11[pos]= poisson(mu1[pos], Y1[pos]:-1)
		f20 = poisson(mu2, Y2)
		f21 = J(n,1,0); pos=selectindex(Y2:!=0); if (length(pos)>0) f21[pos]= poisson(mu2[pos], Y2[pos]:-1)
			
		/*Calculate the Frank Copula function*/
		t1 = - (1/theta) :* log((eta :- (1:-exp(-theta:*f10)):*(1:-exp(-theta:*f20))):/eta)
		t2 = - (1/theta) :* log((eta :- (1:-exp(-theta:*f11)):*(1:-exp(-theta:*f20))):/eta)
		t3 = - (1/theta) :* log((eta :- (1:-exp(-theta:*f10)):*(1:-exp(-theta:*f21))):/eta)
		t4 = - (1/theta) :* log((eta :- (1:-exp(-theta:*f11)):*(1:-exp(-theta:*f21))):/eta)
		
		/* Set up for reasonable bounds of probabilities*/
		minv = J(n,1,1e-7) /*1e-7 -- 1e-6*/
		maxv = 1.0:-minv

		/* Impose Frechet-Hoeffding bounds*/
		t1 = rowmin((rowmax((t1, rowmax((f10+f20:-1,minv)))),rowmin((f10, f20,maxv))))
		t2 = rowmin((rowmax((t2, rowmax((f11+f20:-1,minv)))),rowmin((f11, f20,maxv))))
		t3 = rowmin((rowmax((t3, rowmax((f10+f21:-1,minv)))),rowmin((f10, f21,maxv))))
		t4 = rowmin((rowmax((t4, rowmax((f11+f21:-1,minv)))),rowmin((f11, f21,maxv))))
		
		Phi = rowmin((maxv,rowmax((minv,t1-t2-t3+t4))))	
	}
	if (BPotype=="pfamoye") {
		q=1-exp(-1)
		c1=exp(-Y1);c2=exp(-Y2);	
		mu1=exp(X1*be1');mu2=exp(X2*be2')
		d1=exp(-q:*mu1);d2=exp(-q:*mu2);
		Po1=poissonp(mu1,Y1);Po2=poissonp(mu2,Y2)
		Phi=Po1:*Po2:*(1:+lambda:*(c1 - d1):*(c2 - d2))		
		if (min(Phi)<=0) {
			pos = selectindex(Phi:<=0)
			Phi[pos] = J(length(pos),1,10^-323)
		}		
	}	
	return(Phi)
}
end

*******************************************
capture mata: mata drop BWforback()
mata: /*BWforback*/
function BWforback(Phi,piv,PIm) {
    	
	/*preliminaries*/
	Tv=cols(Phi); k2=rows(Phi)
	TTv=(1..Tv)'#J(k2,1,1)
	tol=0.1; pot=100
	tol=tol^pot
	Qf=J(k2,Tv,.)
	Pp1=J(k2,Tv,0)
	
	/*forward recursion: compute probabilties*/

	//first time occasion
	ff=rowmax((Phi[.,1],J(k2,1,tol)))
	q=piv':*ff
	mm=quadsum(q)
	q=q:/mm
	norm=log(mm)
	Qf[.,1]=q

	//following occasions
	for (t=2; t<=Tv; t++) {
		ff=rowmax((Phi[.,t],J(k2,1,tol)))
		q0=q
		PI=PIm[selectindex(TTv:==t),.]
		for (c=1; c<=k2; c++) {
			q[c]=ff[c]*sum(q0:*PI[.,c])
		}
		mm=quadsum(q)
		q=q:/mm
		norm=norm+log(mm)
		Qf[.,t]=q
	}
	//compute log-likelihood
	lk=log(quadsum(q))+norm

	/*backward recursion: compute probabilties*/
	//for the last time occasion
	qb=J(1,k2,1)
	Pp1[.,Tv] = Qf[.,Tv]:/quadsum(Qf[.,Tv])
	if (Tv==1) Pp2=J(k2,k2,0)
	else {
		ff=rowmax((Phi[.,Tv],J(k2,1,tol)))
		qf0=Qf[.,Tv-1]
		temp=J(k2,k2,0)
		PI=PIm[selectindex(TTv:==Tv),.]
		for (c=1; c<=k2; c++) {
			temp[c,.]=qf0[c]:*PI[c,.]:*ff'
		}
		temp=temp:/quadsum(temp)
		Pp2=temp	

		//iterate in backward order
		for (t=(Tv-1); t>=1; t--) {
			qb0 = qb
			qb =J(1,k2,0)
			PI=PIm[selectindex(TTv:==(t+1)),.]
			for (c=1; c<=k2; c++) {
				qb[c]=quadsum(qb0:*ff':*PI[c,.])
			}
			qb=qb:/quadsum(qb)
			//compute first order posterior probabilities
			Pp1[.,t] = Qf[.,t]:*qb'
			Pp1[.,t] = Pp1[.,t]:/quadsum(Pp1[.,t])
			
			//compute second order probabilities
			if (t>1){
				ff=rowmax((Phi[.,t],J(k2,1,tol)))
				qf0=Qf[.,t-1]
				temp=J(k2,k2,0)
				PI=PIm[selectindex(TTv:==t),.]
				for (c=1; c<=k2; c++) {
					temp[c,.]=qf0[c]:*PI[c,.]:*ff':*qb
				}
				temp=temp:/quadsum(temp)
				Pp2=(temp\Pp2)
			}
		}
		Pp2=(J(k2,k2,0)\Pp2)
	}
	return((&(Pp1'),&Pp2))
	//return((&lk,&(Pp1'),&Pp2))
}
end

/*
*******************************************
capture mata: mata drop BWforback()
mata: /*BWforback*/
function BWforback(Phi,piv,PIm) {
    	
	/*preliminaries*/
	Tv=cols(Phi); k2=rows(Phi)
	TTv=(1..Tv)'#J(k2,1,1)
	tol=0.1; pot=100
	tol=tol^pot
	Qf=J(k2,Tv,.);Qb=J(k2,Tv,.)
	
	/*forward recursion: compute probabilties*/

	//first time occasion
	ff=rowmax((Phi[.,1],J(k2,1,tol)))
	q=piv':*ff
	Qf[.,1]=q
	
	//following occasions
	for (t=2; t<=Tv; t++) {
		ff=rowmax((Phi[.,t],J(k2,1,tol)))
		q0=q
		PI=PIm[selectindex(TTv:==t),.]
		for (c=1; c<=k2; c++) {
			q[c]=ff[c]*sum(q0:*PI[.,c])
		}
		Qf[.,t]=q
	}

	/*backward recursion: compute probabilties*/	
	for (t=(Tv); t>=1; t--) {
		if (t==Tv) {
			Qb[.,t]=J(k2,1,1)
		}
		if (t<Tv) {
			qb0 = Qb[.,t+1]
			PI=PIm[selectindex(TTv:==(t+1)),.]
			//ff=rowmax((Phi[.,t+1],J(k2,1,tol)))
			ff=Phi[.,t]
			qb=J(k2,1,.)
			for (c=1; c<=k2; c++) {
				qb[c]=quadsum(qb0:*ff:*PI[c,.]')
			}
			Qb[.,t]=qb
		}
	}
	
	/*compute posterior probabilties*/
	Pp1=(Qf':*Qb'):/quadrowsum(Qf':*Qb')
		
	if (Tv>1) {
		for (t=(Tv); t>1; t--) {
			qf=Qf[.,t-1]
			qb=Qb[.,t]
			//ff=rowmax((Phi[.,t],J(k2,1,tol)))
			ff=Phi[.,t]
			PI=PIm[selectindex(TTv:==t),.]
			temp=J(k2,k2,.)
			for (c=1; c<=k2; c++) {
				temp[c,.]=qf[c]:*PI[c,.]:*qb':*ff'
			}
			temp=temp:/sum(temp)
			
			if (t==Tv) Pp2=temp
			if (t!=Tv) Pp2=temp\Pp2
		}
		Pp2=J(k2,k2,.)\Pp2
	}
	else Pp2=J(k2,k2,0)
		
	return((&Pp1,&Pp2))
}
end

*/

*******************************************
capture mata mata drop postprob()
mata: /*postprob*/
function postprob(all_id,output) {
    
	/*preliminaries*/
	id=(*all_id[1]);idtc=(*all_id[3]);	
	N=rows(uniqrows(id)); NT=rows(id);	
	Piv=*output[1]; Pi=*output[2]; Phi=*output[3]; Qt=*output[4]; pu=*output[5];
	k=max(idtc[.,3])
	
	wv=J(NT,k,.); wvvb=J(NT*k,k,.)
	for (i=1; i<=N; i++) {
		pos=selectindex(id:==i)
		posc=selectindex(idtc[.,1]:==i)			
		a=BWforback(Phi[pos,.]',Piv[i,.],Pi[posc,.])
		wv[pos,.]=(*a[1])
		wvvb[posc,.]=(*a[2])
	}
		
	return(&wv,&wvvb)
}
end

*******************************************
capture mata: mata drop lk_fraclogit()
mata /*lk_fraclogit*/
void lk_fraclogit(real scalar todo,real vector b,W,X,lev,base,/*
		*/ LL, S, H) {
	/*preliminaries*/
	npp=cols(b)/(lev-1)
	n=rows(X)
	/* compute log-lik */
	Pr=prlogit(X,b',lev,base)
	LL=quadsum(log(Pr):*W)
	
	if (todo>0) {
	    S=J(1,cols(X),.)
	    bv=J(lev,npp,.);  bvv=1..(cols(X))		
		bv[selectindex((1..lev):!=base),.]=rowshape(bvv,lev-1)
		for (x=1; x<=cols(X); x++) {
			 pv=selectindex(quadrowsum(bv:==x):==1)
			 temp=J(1,lev,-Pr[.,pv]); temp[.,pv]=1:+temp[.,pv]
			 S[x]=quadsum(temp:*W:*X[.,x])
		}
	}
}
end

*******************************************
capture mata: mata drop lk_Phi()
mata /*lk_Phi*/
void lk_Phi(real scalar todo,real vector b,all_id,npar,W,Y,X,option,/*
		*/ LL, S, H) {
	
	/* collect parameters in b to be estimated */
	rnpar=runningsum(npar)
	if (min(npar)==0) {
		rnpar=J(1,length(npar),0)
		for (x=1; x<=length(npar); x++) { 
			if (npar[x]> 0) rnpar[x]=sum(npar[1..x])
			if (npar[x]==0) rnpar[x]=rnpar[x-1]
		}
	}
	be1    = b[1..rnpar[1]];
	be2    = b[(rnpar[1]+1)..rnpar[2]];
	lambda=0; if (npar[3]!=0) lambda = b[(rnpar[2]+1)..rnpar[3]];
	
	/*preliminaries*/
	id=(*all_id[1]);idtc=(*all_id[3]);
	Y1=Y[.,1]; Y2=Y[.,2];
	X1=*X[1]; X2=*X[2]; X0=*X[5]
	k=cols(X0); 
	N=rows(uniqrows(id)); NT=rows(id);
	BPotype=*option[1]
		
	Phi=J(NT,k,.)		
	for (c=1; c<=k; c++) {
		Xtc=(*X0[c]);
		Xtc1=J(NT,1,Xtc); Xtc2=Xtc1
		if (cols(X1)>0) Xtc1=(Xtc1,X1)
		if (cols(X2)>0) Xtc2=(Xtc2,X2)
		Phi[.,c]=poissonfam_d(Y1,Y2,Xtc1,Xtc2,be1,be2,lambda,BPotype)					
	}
	L=log(Phi):*W
	
	/* compute log-lik */
	LL=quadsum(L)
	
	if (todo>0) {
		if (BPotype=="pfamoye") {
			/* compute score */
			SSlambda=J(1,1,0);SSbe1=J(1,npar[1],0);SSbe2=J(1,npar[2],0)
			for (c=1; c<=k; c++) {
				/* *** */
				Utc=W[.,c]
				Xtc=(*X0[c]);
				Xtc1=J(NT,1,Xtc); Xtc2=Xtc1
				if (cols(X1)>0) Xtc1=(Xtc1,X1)
				if (cols(X2)>0) Xtc2=(Xtc2,X2)
				q=1-exp(-1)			
				eta1=Xtc1*be1'; eta2=Xtc2*be2'
				mu1=exp(eta1); mu2=exp(eta2)
				c1=exp(-q:*mu1);c2=exp(-q:*mu2);
				
				psi1 = (lambda:*(exp(-Y1):-c1)):/(1:+lambda:*(exp(-Y1):-c1):*(exp(-Y2):-c2));
				psi2 = (lambda:*(exp(-Y2):-c2)):/(1:+lambda:*(exp(-Y1):-c1):*(exp(-Y2):-c2));
				/* *** */
				temp=((exp(-Y1):-c1):*(exp(-Y2):-c2)):/(1:+lambda:*(exp(-Y1):-c1):*(exp(-Y2):-c2));
				temp=(temp):*Utc
				SSlambda=SSlambda+quadsum(temp)			
				temp=( Y1:-mu1:+q:*psi2:*mu1:*c1 ):*Xtc1:*Utc
				SSbe1=SSbe1:+quadcolsum(temp)
				temp=( Y2:-mu2:+q:*psi1:*mu2:*c2 ):*Xtc2:*Utc
				SSbe2=SSbe2:+quadcolsum(temp)			
			}			
			S=(SSbe1,SSbe2,SSlambda) /* score without covariates */	
		}
	}

}
end

*******************************************
capture mata: mata drop hmbivpois_EM()
mata: /*hmbivpois_EM*/
function hmbivpois_EM(b,npar,all_id,Y,X,option,emiter,mlem,emtol,tracep,start_par_name) {	
	
	/*collect parameters in b to be estimated*/
	rnpar=runningsum(npar)
	if (min(npar)==0) {
		rnpar=J(1,length(npar),0)
		for (x=1; x<=length(npar); x++) { 
			if (npar[x]> 0) rnpar[x]=sum(npar[1..x])
			if (npar[x]==0) rnpar[x]=rnpar[x-1]
		}
	}
	be1    = b[1..rnpar[1]];
	be2    = b[(rnpar[1]+1)..rnpar[2]];
	lambda=0; if (npar[3]>0) lambda = b[(rnpar[2]+1)..rnpar[3]];
	if (npar[4]>0) be3 = b[(rnpar[3]+1)..rnpar[4]];
	if (npar[5]>0) be4 = b[(rnpar[4]+1)..rnpar[5]];
	
	/*preliminaries*/
	id=(*all_id[1]);time=(*all_id[2]);idtc=(*all_id[3]);
	Y1=Y[.,1]; Y2=Y[.,2]; 
	k=cols(*X[5])
	N=rows(uniqrows(id)); NT=rows(id); T=max(time)
	BPotype=*option[1]
	 
	/*compute initial log-likelihood*/
	hmbivpois_lp(0,b,npar,all_id,Y,X,option,output=.,lk=., S=.,H=.) 
	lk0=lk
	Piv=*output[1]; Pi=*output[2]; Phi=*output[3]; Qt=*output[4]; pu=*output[5];
	
	/*EM algorithm*/	
	b0=b; lk0=lk-1000
	it=0
	printf("%s {txt}EM-iteration = {res}%f, {txt}f(p) = {res}%f, {txt}diff. = {res}%f \n", "", it, lk,.)
	if (tracep=="on") {
	    st_matrix("betaEM", b)
		st_matrixcolstripe("betaEM", start_par_name)
		st_matrix_list("betaEM")
	}
	while ( ( (abs(lk-lk0)>emtol) & (it<emiter) ) ==1) {
		timer_clear(16);timer_on(16); it=it+1; lk0=lk; b0=b;
	
		//E-step
		printf("{txt}Performing E-step... ")
		displayflush()
		
		wv=J(NT,k,.); wvvb=J(NT*k,k,.)
		for (i=1; i<=N; i++) {
			pos=selectindex(id:==i)
			posc=selectindex(idtc[.,1]:==i)			
			a=BWforback(Phi[pos,.]',Piv[i,.],Pi[posc,.])
			wv[pos,.]=(*a[1])
			wvvb[posc,.]=(*a[2])
		}
		printf("{txt}completed. \n")
		
		//M-step
		if (k>1) {
			//1.Piu
			temp=(-99\id); temp=temp[1..(rows(temp)-1)]		
			pos=selectindex( id:!=temp[.,1] )
			F  = optimize_init() /*F is the function we want to maximize*/
			optimize_init_evaluator(F, &lk_fraclogit())
			optimize_init_evaluatortype(F, "d1")
			optimize_init_params(F, be3) /*initial parameters */
			optimize_init_argument(F, 1, wv[pos,.]) /*additional entries in void*/
			optimize_init_argument(F, 2, (*X[3]))/*additional entries in void*/
			optimize_init_argument(F, 3, k)/*additional entries in void*/
			optimize_init_argument(F, 4, 1)/*additional entries in void*/
			optimize_init_technique(F,"nr")
			optimize_init_conv_maxiter(F, 25)
			be3 = optimize(F)
		}
		
		//2:Pi
		if (k>1) {
			be4T=rowshape(be4',k)'
			for (c=1; c<=k; c++) {
				pos=(idtc[.,2]:>1)+(idtc[.,3]:==c);pos=selectindex(pos:==2)
				F  = optimize_init() /*F is the function we want to maximize*/
				optimize_init_evaluator(F, &lk_fraclogit())
				optimize_init_evaluatortype(F, "d1")
				optimize_init_params(F, be4T[.,c]') /*initial parameters */
				optimize_init_argument(F, 1, wvvb[pos,.]) /*additional entries in void*/
				optimize_init_argument(F, 2, (*X[4]))/*additional entries in void*/
				optimize_init_argument(F, 3, k)/*additional entries in void*/
				optimize_init_argument(F, 4, c)/*additional entries in void*/
				optimize_init_technique(F,"nr")	
				optimize_init_conv_maxiter(F, 25)
				be4T[.,c] = optimize(F)'
			}
			be4=vec(be4T)'
		}

		//3:Phi
		F  = optimize_init() /*F is the function we want to maximize*/
		optimize_init_evaluator(F, &lk_Phi())
		if (BPotype=="normal") optimize_init_evaluatortype(F, "d0")
		if (BPotype=="frank")  optimize_init_evaluatortype(F, "d0")		
		if (BPotype=="pfamoye") optimize_init_evaluatortype(F, "d1")
		if (npar[3]>0) optimize_init_params(F, (be1,be2,lambda)) /*initial parameters */
		if (npar[3]==0) optimize_init_params(F, (be1,be2)) /*initial parameters */
		optimize_init_argument(F, 1, all_id) /*additional entries in void*/
		optimize_init_argument(F, 2, npar)/*additional entries in void*/
		optimize_init_argument(F, 3, wv)/*additional entries in void*/
		optimize_init_argument(F, 4, Y)/*additional entries in void*/
		optimize_init_argument(F, 5, X)/*additional entries in void*/
		optimize_init_argument(F, 6, option)/*additional entries in void*/
		optimize_init_technique(F,"nr")
		optimize_init_conv_maxiter(F, mlem)
		optimize_init_conv_ptol(F,1e-4)
		optimize_init_conv_vtol(F,1e-4)
		optimize_init_conv_nrtol(F,1e-4)		/*optimize_init_tracelevel(F,"gradient")*/
		be_Phi = optimize(F)
		be1    = be_Phi[1..rnpar[1]];
		be2    = be_Phi[(rnpar[1]+1)..rnpar[2]];
		if (npar[3]>0) lambda = be_Phi[(rnpar[2]+1)..rnpar[3]];
					
		//recollect parameters
		b=(be1,be2)
		if (npar[3]>0) b = (b,lambda)
		if (npar[4]>0) b = (b,be3);
		if (npar[5]>0) b = (b,be4);
	
		//compute log-likelihood
		hmbivpois_lp(0,b,npar,all_id,Y,X,option,output=.,lk=., S=.,H=.) 
		Piv=*output[1]; Pi=*output[2]; Phi=*output[3]; Qt=*output[4]; pu=*output[5];
		
		//display some results
		timer_off(16); eti=timer_value(16);timer_clear(16)
		bic=(log(NT)*sum(npar) - 2*(lk))
		printf(" \n")
		printf("%s {txt}EM-iteration = {res}%f, {txt}f(p) = {res}%f, {txt}diff. = {res}%f, {txt}diff.par. = {res}%f, {txt}Ex.Time = {res}%f, {txt}BIC = {res}%f \n", "", it, lk,abs(lk-lk0),abs(sum(b:-b0)),eti[1]/60,bic)
		if (tracep=="on") {
			st_matrix("betaEM", b)
			st_matrixcolstripe("betaEM", start_par_name)
			st_matrix_list("betaEM")
		}
		if ((*option[3])!="") {			
			fh = fopen((*option[3])+".bin", "rw")
			fputmatrix(fh, k)
			fputmatrix(fh, be1)
			fputmatrix(fh, be2)
			if (npar[3]>0) fputmatrix(fh, lambda); else fputmatrix(fh, J(0,0,0))
			if (npar[4]>0) fputmatrix(fh, be3);    else fputmatrix(fh, J(0,0,0))
			if (npar[5]>0) fputmatrix(fh, be4);    else fputmatrix(fh, J(0,0,0))
			fputmatrix(fh, V=.)
			fputmatrix(fh, output)
			fputmatrix(fh, U=.)
			fputmatrix(fh, lk)
			fputmatrix(fh, b)
			fputmatrix(fh, start_par_name)
			fclose(fh)		    
		}
	}
	
	th=b
	return(th)
}
end

*****************************************
capture mata: mata drop hmbivpois_der_fun_lp()
mata /*hmbivpois_der_fun_lp*/
void hmbivpois_der_fun_lp(real vector b,  /*
		*/   npar,all_id,Y,X,option,output, /*
		*/ LL) {

	LL1=.;S=.;H=.;
	hmbivpois_lp(0,b,  /*
		*/   npar,all_id,Y,X,option,output, /*
		*/ LL1, S, H) 
	LL=quadsum(LL1)
}
end

*****************************************
capture mata: mata drop hmbivpois_der_lp()
mata /*hmbivpois_der_lp*/
void hmbivpois_der_lp(start_par,npar,all_id,Y,X,option,output,/*
					*/d0,d1,d2){ 
		
	D = deriv_init() /* this is used to compute second derivative */
	deriv_init_evaluator(D, &hmbivpois_der_fun_lp())
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
	d1=deriv(D, 1)
	d2=deriv(D, 2)

}
end

****************************************
capture mata: mata drop hmbivpois_lpo()
mata /*hmbivpois_lpo*/
void hmbivpois_lpo(typeoptimi,start_par,npar,all_id,Y,X,option,output,/*
					*/ tech,maxiter,tracep,traceg,traceh,start_par_name,/*
					*/th,V,LL){ 
	//traceg="on";typeoptimi="d1debug";tech="nr"
	F  = optimize_init() /*F is the function we want to maximize*/
	optimize_init_evaluator(F, &hmbivpois_lp())
	optimize_init_evaluatortype(F, typeoptimi)
	optimize_init_params(F, start_par) /*initial parameters */
	optimize_init_argument(F, 1, npar) /*additional entries in void*/
	optimize_init_argument(F, 2, all_id) /*additional entries in void*/
	optimize_init_argument(F, 3, Y) /*additional entries in void*/
	optimize_init_argument(F, 4, X) /*additional entries in void*/
	optimize_init_argument(F, 5, option) /*additional entries in void*/
	optimize_init_argument(F, 6, output) /*additional entries in void*/
	optimize_init_technique(F,tech)
	optimize_init_conv_ptol(F, 1e-6) /* tol on parameters */
	optimize_init_conv_vtol(F, 1e-6) /*tol in lok-lik */
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

*******************************************
capture mata: mata drop hmbivpois_lp()
mata /*hmbivpois_lp*/
void hmbivpois_lp(real scalar todo,real vector b,  /*
		*/   npar,all_id,Y,X,option,output, /*
		*/ LL, S, H) {
	
	/* collect parameters in b to be estimated */
	rnpar=runningsum(npar)
	if (min(npar)==0) {
		rnpar=J(1,length(npar),0)
		for (x=1; x<=length(npar); x++) { 
			if (npar[x]> 0) rnpar[x]=sum(npar[1..x])
			if (npar[x]==0) rnpar[x]=rnpar[x-1]
		}
	}
	be1    = b[1..rnpar[1]];
	be2    = b[(rnpar[1]+1)..rnpar[2]];
	lambda=0; if (npar[3]>0) lambda = b[(rnpar[2]+1)..rnpar[3]];
	if (npar[4]>0) be3    = b[(rnpar[3]+1)..rnpar[4]];
	if (npar[5]>0) be4    = b[(rnpar[4]+1)..rnpar[5]];

	/*preliminareis*/
	id=(*all_id[1]);time=(*all_id[2]);idtc=(*all_id[3]);
	index=*all_id[4];post_T=*index[1];post0_T=*index[2];post_pi_T=*index[3];
	Y1=Y[.,1]; Y2=Y[.,2];
	X1=*X[1]; X2=*X[2]; X3=*X[3]; X4=*X[4]; X0=*X[5];
	k=cols(X0);
	N=rows(uniqrows(id)); NT=rows(id); T=max(time)
	BPotype=*option[1]	
	
	/*Pr(U=u)*/
	if (k>1) Piv=prlogit(X3,be3',k)
	else Piv=J(N,1,1)

	/*Pr(u|ub)*/
	if (k>1) {
	    be4T=rowshape(be4',k)'
		Pi=J(rows(idtc),k,.)
		temp=(J(k,cols(idtc),.)\idtc); temp=temp[1..rows(idtc),.]
		for (c=1; c<=k; c++) {
			pos=selectindex( ((idtc[.,1]:==temp[.,1])+(idtc[.,3]:==c)):==2 )			
			Pi[pos,.]=prlogit(X4,be4T[.,c],k,c)
		}
	}
	else	{
	    Pi=J(rows(idtc),k,1)
		Pi[selectindex(time:==1)]=J(N,1,.)
	}
	
	/*Pr(y1,y2|u)*/
	Phi=J(NT,k,.)
	for (c=1; c<=k; c++) {
		Xtc=J(NT,1,*X0[c])
		Xtc1=Xtc; Xtc2=Xtc; 
		if (cols(X1)>0) Xtc1=(Xtc,X1)
		if (cols(X2)>0) Xtc2=(Xtc,X2)
		Phi[.,c]=poissonfam_d(Y1,Y2,Xtc1,Xtc2,be1,be2,lambda,BPotype)
	}
	
	/*compute log-likelihood*/		
	Qt=J(NT,k,.);/*forward recursion given v*/
	for (t=1; t<=T; t++) {
		if (t==1) { //first time occasion
			post=(*post_T[1])
			Qt[post,.]=Phi[post,.]:*Piv
		}
		else { //subsequent periods
			post0=(*post0_T[t])
			post=(*post_T[t])
			post_pi=(*post_pi_T[t])		
			//Pi=J(NT,1,1)#(1,2,3\4,5,6\7,8,9)			
			tmp=Pi[post_pi,.]		
			for (c=1; c<=k; c++) {
				Qt[post,c]=quadrowsum(Qt[post0,.]:*rowshape(tmp[.,c]',rows(post)):*Phi[post,c])
			}			
		}
	}
	if (min(Qt)==0)	Qt=Qt:+((Qt:==0):*1e-323)		
	
	/*last period*/
	temp=(id\-99); temp=temp[2..rows(temp)]		
	post=selectindex( id:!=temp[.,1] )
	pv=quadrowsum(Qt[post,.])		
		
	/*individual likelihood*/
	if (min(pv)==0)	pv=pv:+((pv:==0):*1e-323)
	LL=quadsum(log(pv))

	theta=lambda
	if (BPotype=="normal") theta=(exp(2*lambda)-1)/(exp(2*lambda)+1)
	
	output=(&Piv,&Pi,&Phi,&Qt,&pv,&theta)
	
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
	syntax [if] ,[ k(integer 1) id(varname) time(varname) /*
	*/ INITcov(varlist)  EMiter(integer 50) emonly/*
	*/ xtrans(varlist) HETerogeneous restr BPotype(name) /*
	*/ maxiter(integer 3000) mlem(integer 25) emtol(real .0001)/*
	*/ tracep traceg traceh /*
	*/ STARTfrom(name) prev(string) STARTPOISson/*
	*/ load(string) save(string) EVALuate /*
	*/ optimize_init_technique(string) optimize_init_evaluatortype(integer 0)] 
	 
	if "`optimize_init_technique'" == ""	local optimize_init_technique = "bfgs"
	
	marksample touse
	markout `touse'
	
	tempvar _set_sample_HM _`id'_HM _`time'_HM _`time'1_HM _info_sample_HM _`id'1_HM
	
	if `k' == 1 & "`initcov'" != "" {
		dis as txt "Option " as input "initcov " as txt "ignored since " as input "k = `k'"
	}
	*if `k' == 1 dis as txt "EM algorithm is not used with " as input "k = `k'"
	/*
	if "`emonly'" != "" {
		local emiter = 1000
		dis as txt "Total number of EM iteration is set = " as input "1000 "
	}
	*/
	if "`xtrans'" != "" & "`heterogeneous'" != "" {
		dis as err "You cannot set both xtrans and heterogeneous options at the same time"
	}	
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
	mata: time_o	= st_data((1..`N')', "`time'")
	mata: idtc		= (id,time)#J(`k',1,1),J(`N',1,1)#(1..`k')'
	mata: Y1 		= st_data((1..`N')', "`dep1n'")
	mata: Y2 		= st_data((1..`N')', "`dep2n'")
	mata: X0        = setintercept(`k')
	mata: X1=J(0,0,0);X2=J(0,0,0);
	if  wordcount("`ind1'") > 0 {
		mata: X1 = st_data((1..`N')', "`ind1'")
	}
	if  wordcount("`ind2'") > 0 {
		mata: X2 = st_data((1..`N')', "`ind2'")
	}
	if "`initcov'" != "" & `k' > 1 {
		mata: X3=st_data((1..`N')', "`initcov'")
		mata: X3=X3[selectindex((.\id)[1..rows(id)]:!=id),.];X3=(J(rows(X3),1,1),X3)
		mata: X3=J(1,`k'-1,1)#X3
	}
	else {
		mata: X3 = J(rows(selectindex((.\id)[1..rows(id)]:!=id)),`k'-1,1)
	}
	if `k' > 1 {
	    mata: X4 = J(rows(selectindex((.\id)[1..rows(id)]:==id)),`k'-1,1)
		if "`xtrans'" != ""  {
			mata: X4=st_data((1..`N')', "`xtrans'")
			mata: X4=X4[selectindex((.\id)[1..rows(id)]:==id),.]
			mata: X4=J(rows(id)-rows(uniqrows(id)),1,1),X4
			mata: X4=J(1,`k'-1,1)#X4
		}
		if "`heterogeneous'" != "" {
			mata: X4=time[selectindex((.\id)[1..rows(id)]:==id)]:==(uniqrows(time)'#J(rows(id)-rows(uniqrows(id)),1,1))
			mata: X4=X5[.,2..cols(X4)]
			mata: X5=J(1,`k'-1,1)#X4
		}
	}
	else {
		mata: X4 = J(rows(Y1),`k'-1,1)
	}
	
	/* set starting parameters */
	if "`startfrom'" == "" {
		if "`prev'" == "" & "`load'" != "" & "`startpoisson'" == "" {		
			dis ""
			if "`evaluate'" == "" dis as text "Load starting values from the file " as res "`load'.bin"
			else                  dis as text "Evaluate log-likelihood using file " as res "`load'.bin"
			dis ""
			mata: fh = fopen("`load'.bin", "r")		
			mata: k0 = fgetmatrix(fh)
			mata: be1 = fgetmatrix(fh)
			mata: be2 = fgetmatrix(fh)
			mata: lambda = fgetmatrix(fh)
			mata: be3 = fgetmatrix(fh)
			mata: be4 = fgetmatrix(fh)
			mata: V = fgetmatrix(fh)	
			mata: output = fgetmatrix(fh)	
			mata: U = fgetmatrix(fh)	
			mata: LL = fgetmatrix(fh)
			mata: start_par_k0 = fgetmatrix(fh)
			mata: start_par_name_k0 = fgetmatrix(fh)
			mata: fclose(fh)
		}
		else if "`prev'" != "" & "`load'" == "" & "`startpoisson'" == "" {
			mata: fh = fopen("`prev'.bin", "r")
			mata: k0 = fgetmatrix(fh)
			mata: be1 = fgetmatrix(fh)
			mata: be2 = fgetmatrix(fh)
			mata: lambda = fgetmatrix(fh)
			mata: be3 = fgetmatrix(fh)
			mata: be4 = fgetmatrix(fh)
			mata: V = fgetmatrix(fh)	
			mata: output = fgetmatrix(fh)	
			mata: U = fgetmatrix(fh)	
			mata: LL = fgetmatrix(fh)
			mata: start_par_k0 = fgetmatrix(fh)
			mata: start_par_name_k0 = fgetmatrix(fh)		
			mata: fclose(fh)
			
			mata: st_local("k0",strofreal(k0));
			if `k0'!=(`k'-1)  {
				dis as err "Number of classes k in " as res "`prev'.bin" as err " are " as res "`k0' " as err " and not " as res "`=`k'-1' "
				exit(99)
			}
		
			dis ""
			dis as text "Load starting values from the file " as res "`prev'.bin" as text " with " as res "`k0' " as text "latent classes"
			dis ""
				
			mata: be1 = be1[1..`k0'],be1[`k0']+0.5,be1[(`k0'+1)..cols(be1)]
			mata: be2 = be2[1..`k0'],be2[`k0']+0.5,be2[(`k0'+1)..cols(be2)]

			/*initial prob.*/
			local lI = wordcount("`initcov'")+1
			mata: be3 = J((`k'-1)*(`lI'),1,0); be3=be3'

			/*transition prob.*/					
			local lT=1
			if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
			if "`heterogeneous'" != ""		local lT = (`T'-1)
			*mata: be5 = J((`k2'-1)*(`lT')*(`k2'),1,0); tmp=J(`lT',1,0); tmp[1]=-log(10); be5=be5+J((`k2'-1)*(`k2'),1,1)#tmp; be5=be5'
			mata: be4 = J(1,(`k'-1)*(`lT')*(`k'),0);
				
		}
		else if "`prev'" == "" & "`load'" == "" & "`startpoisson'" != "" {
			dis ""
			dis as text "Starting parameters from:"
			dis as text "Univariate Poisson regression for Y1 and Y2 respectively"			
			qui poisson `dep1n' `ind1'
			local lcomp = e(ll)
			local b = _b[`dep1n':_cons]
			mata: be1=J(1,`k',`b'):+((0..(`k'-1)):/10) //dep1
			if wordcount("`ind1'") > 0 {
				local i = 1
				foreach x in `ind1' {
					local b = _b[`dep1n':`x']
					mata: be1 = (be1,`b')
					local ++i
				}
			}			
			qui poisson `dep2n'  `ind2' 
			local lcomp = e(ll) + `lcomp'
			local b = _b[`dep2n':_cons]		
			*mata: be2=J(1,`k',1/`k'):*`b';be2[1]=`b'-0.5 //dep2
			mata: be2=J(1,`k',`b'):+((0..(`k'-1)):/10) //dep1
			if wordcount("`ind2'") > 0 {
				local i = 1
				foreach x in `ind2' {
					local b = _b[`dep2n':`x']
					mata: be2 = (be2,`b')
					local ++i
				}
			}
			mata: lambda = 0
			/*initial prob.*/
			local lI = wordcount("`initcov'")+1
			mata: be3 = J((`k'-1)*(`lI'),1,0); be3=be3'
			/*transition prob.*/
			local lT=1
			if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
			if "`heterogeneous'" != ""		local lT = (`T'-1)
			*mata: be5 = J((`k2'-1)*(`lT')*(`k2'),1,0); tmp=J(`lT',1,0); tmp[1]=-log(10); be5=be5+J((`k2'-1)*(`k2'),1,1)#tmp; be5=be5'
			mata: be4 = J(1,(`k'-1)*(`lT')*(`k'),0);	
			dis as text "Comparison:    log likelihood = " as res `lcomp'
			dis ""			
		}
		else {	    
			dis ""
			dis as text "Starting parameters from:"
			if "`bpotype'"== "normal" {
				dis as text "Bivariate Poisson regression for Y1 and Y2 with parameterization : `bpotype'"		
				qui bivcnto (`dep1n' = `ind1' ) (`dep2n' = `ind2' ), dist1(poisson) dist2(poisson) copula(normal)			
			}
			else if "`bpotype'"== "frank" {
				dis as text "Bivariate Poisson regression for Y1 and Y2 with parameterization : `bpotype'"		
				qui bivcnto (`dep1n' = `ind1' ) (`dep2n' = `ind2' ), dist1(poisson) dist2(poisson) copula(frank)
			}
			else {
				dis as text "Bivariate Poisson regression for Y1 and Y2 with parameterization : pfamoye"		
				qui bivcnto (`dep1n' = `ind1' ) (`dep2n' = `ind2' ), pfamoye
			}
			local lcomp = e(ll)
			local b = _b[`dep1n':_cons]		
			*mata: be1=J(1,`k',1/`k'):*`b';be1[1]=`b'-0.5 //dep1
			mata: be1=J(1,`k',`b'):+((0..(`k'-1)):/10) //dep1
			if wordcount("`ind1'") > 0 {
				local i = 1
				foreach x in `ind1' {
					local b = _b[`dep1n':`x']
					mata: be1 = (be1,`b')
					local ++i
				}
			}
			local b = _b[`dep2n':_cons]		
			*mata: be2=J(1,`k',1/`k'):*`b';be2[1]=`b'-0.5 //dep2
			mata: be2=J(1,`k',`b'):+((0..(`k'-1)):/10) //dep1
			if wordcount("`ind2'") > 0 {
				local i = 1
				foreach x in `ind2' {
					local b = _b[`dep2n':`x']
					mata: be2 = (be2,`b')
					local ++i
				}
			}

			if "`bpotype'"== "normal" local lambda = _b[/:atanhtheta]  //atanhtheta
			else if "`bpotype'"== "frank"  local lambda = _b[/:theta] //theta
			else                      local lambda = _b[/:lambda] //lambda
			mata: lambda = `lambda'
			/*initial prob.*/
			local lI = wordcount("`initcov'")+1
			mata: be3 = J((`k'-1)*(`lI'),1,0); be3=be3'
			/*transition prob.*/
			local lT=1
			if "`xtrans'" != ""				local lT = (wordcount("`xtrans'")+1)
			if "`heterogeneous'" != ""		local lT = (`T'-1)
			*mata: be5 = J((`k2'-1)*(`lT')*(`k2'),1,0); tmp=J(`lT',1,0); tmp[1]=-log(10); be5=be5+J((`k2'-1)*(`k2'),1,1)#tmp; be5=be5'
			mata: be4 = J(1,(`k'-1)*(`lT')*(`k'),0);	
			dis as text "Comparison:    log likelihood = " as res `lcomp'
			dis ""
		}
		if `k' ==1 mata: be1[1]=0;be2[1]=0;lambda=0
		mata: start_par=(be1,be2,lambda)
		if `k' > 1 mata: start_par=(start_par,be3,be4)
		
		*mata: npar=(cols(be1),cols(be2),cols(lambda),cols(be3),cols(be4))
	}
	else {
	    dis ""
		dis as text "Load starting values from the matrix " as res "`startfrom'"
		dis ""
		mat list `startfrom'
		mata: start_par = st_matrix("`startfrom'")
	}
	mata: npar=(cols(*X0[1])+cols(X1),cols(*X0[1])+cols(X2),1,cols(X3),`k'*cols(X4))
	
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
	if "`bpotype'" == "normal" {
		local eqname = "`eqname' atanhtheta"
		local name_be = "`name_be' atanhtheta"
	}
	else if "`bpotype'" == "frank" {
		local eqname = "`eqname' theta"
		local name_be = "`name_be' theta"
	}
	else {
		local eqname = "`eqname' lambda"
		local name_be = "`name_be' lambda"
	}
	local name_be3 = ""
	if `k'>1 {
		foreach c of numlist 2/`k' {
			foreach x in alpha `initcov' {
			    local eqname = "`eqname' Init.cond."
				local name_be3 = "`name_be3' `x'_`c'"
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
							if "`x'"=="alpha" local name_be4 = "`name_be4' `x'`cc'_`c'"
							else              local name_be4 = "`name_be4' `x'        "
						}
					}
				}
			}
		}	
		else if "`heterogeneous'" != "" {
			foreach c of numlist 1/`' {
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
	
	if "`load'"!="" & "`evaluate'" == "" & "`startfrom'"==""{
		dis as text "Check loaded starting values.."
		local full1 = "`eqname'"
		local full2 = "`name_be' `name_be3' `name_be4'"
		local nx = wordcount("`full1'")
		mata: th0=J(1,sum(npar),0)
		local coefnotfound = ""
		foreach x of numlist 1/`nx' {
			local x1 = word("`full1'",`x')
			local x2 = word("`full2'",`x')
			mata: pos=selectindex(quadrowsum(start_par_name_k0:==("`x1'","`x2'")):==2)
			mata: st_local("pos",strofreal(length(pos)));
			if `pos' == 1 mata: th0[`x']=start_par_k0[pos]
			else local coefnotfound = "`coefnotfound' `x1'-`x2';"
		}
		if wordcount("`coefnotfound'")>0 dis as txt "Starting values has been set to 0 for these vars:" as res "`coefnotfound'"
		mata: start_par=th0
	}
	
	tempname start_par
	mata: st_matrix("`start_par'", start_par)
	mat colnames `start_par' = `name_be' `name_be3' `name_be4'
	mat coleq    `start_par' = `eqname'
	*mat list `start_par'	

	if "`restr'" != ""	local restr = 1
	else				local restr = 0	

	mata: index=idt_index((id,time),idtc,`k') /*(&post_T,&post0_T,&post_pi_T)*/

	/* set options */
	if "`bpotype'" != "" mata: option=(&"`bpotype'",&`restr',&"`save'")
	else                 mata: option=(&"pfamoye",&`restr',&"`save'")
	/*
	*************
	/* log-lik */
	*************
	local restr = 0	
	*mata LL=.;S=.;H=.
	
	mata: hmbivpois_lp(0,start_par,  /*
			*/   npar,(&id,&time,&idtc,&index),(Y1,Y2),(&X1,&X2,&X3,&X4,&X0),option,output=., /*
			*/ LL=., S=.,H=.) 
	mata: LL
	mata: S
	
	****************
	/* derivative */
	****************	
	mata: d0=.;d1=.;d2=.;
	mata: hmbivpois_der_lp(start_par,npar,(&id,&time,&idtc,&index),(Y1,Y2),(&X1,&X2,&X3,&X4,&X0), /*
						*/  option,output,d0,d1,d2)
	mata: d0,quadsum(LL)
	mata: d1'
	exit(1)
	*/

	****************
	/* estimation */
	****************
	if "`evaluate'" != ""{
		mata: th=(be1,be2,lambda)
		if  `k' > 1 mata: th=(th,be3,be4)
	}
	else {
	    mata: th=.;V=.;LL=.;		
		/*EM*/
		*& (`k'>1 )
		if `emiter'>0  mata: start_par=hmbivpois_EM(start_par,npar,(&id,&time,&idtc,&index),(Y1,Y2),(&X1,&X2,&X3,&X4,&X0),option,`emiter',`mlem',`emtol',"`tracep'",st_matrixcolstripe("`start_par'"))
		/*Direct maximization*/
		if "`emonly'" != "" {
			mata: th=start_par
			mata: hmbivpois_der_lp(th,npar,(&id,&time,&idtc,&index),(Y1,Y2),(&X1,&X2,&X3,&X4,&X0),option,output=.,LL,S=.,V)
			mata: V=invsym(-V)
		}
		else {
			mata: hmbivpois_lpo("gf`optimize_init_evaluatortype'",start_par,npar,(&id,&time,&idtc,&index),(Y1,Y2),(&X1,&X2,&X3,&X4,&X0), /*
				*/  option,output=.,"`optimize_init_technique'",`maxiter',"`tracep'","`traceg'","`traceh'",st_matrixcolstripe("`start_par'"),th,V,LL)
		}
		dis "Estimation completed"
	}

	*************************************
	/* collect some additional outputs */
	*************************************
	mata: Piv=*output[1]; Pi=*output[2]; Phi=*output[3]; Qt=*output[4]; pu=*output[5];
	if "`evaluate'" == "" {
		mata: pprob=postprob((&id,&time,&idtc,&index),output); 
		mata: U=(*pprob[1]);
	}
	mata: U_we_t = J(`T',`k',.)		
	foreach t of numlist 1/`T' {	
		mata: U_we_t[`t',.]=mean(U[selectindex(time:==`t'),.])			
	}
	mata: Pim=J(`k',`k',.)
	foreach c of numlist 1/`k' {
		mata: Pim[`c',.]=mean(Pi[selectindex(((idtc[.,2]:>1):+(idtc[.,3]:==`c')):==2),`.'])
	}

	tempname b V ll BIC AIC npar we we_t Pi theta
	
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
	mata: theta=(*output[6]); st_numscalar("`theta'", theta )
	ereturn scalar theta = `theta'					/* theta */		
	mata: st_matrix("`we'", mean(U)) 
	ereturn matrix we = `we'						/* we */
	mata: st_matrix("`we_t'", U_we_t)	
	ereturn matrix we_t = `we_t'						/* we */	
	mata: st_matrix("`Pi'", Pim)
	ereturn matrix Pi = `Pi'						/* Average transition prob. */	

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
		mata:                                       be1=th[1..rnpar[1]]		
		mata:                                       be2=th[(rnpar[1]+1)..rnpar[2]];
		mata: lambda=J(0,0,.); if (npar[3]>0) {; lambda=th[(rnpar[2]+1)..rnpar[3]]; };
		mata: be3=J(0,0,.);    if (npar[4]>0) {;    be3=th[(rnpar[3]+1)..rnpar[4]]; };	
		mata: be4=J(0,0,.);    if (npar[5]>0) {;    be4=th[(rnpar[4]+1)..rnpar[5]]; };
		
		mata: fh = fopen("`save'.bin", "rw")
		mata: fputmatrix(fh, `k')
		mata: fputmatrix(fh, be1)
		mata: fputmatrix(fh, be2)
		mata: fputmatrix(fh, lambda)
		mata: fputmatrix(fh, be3)
		mata: fputmatrix(fh, be4)
		mata: fputmatrix(fh, V)
		mata: fputmatrix(fh, output)
		mata: fputmatrix(fh, U)
		mata: fputmatrix(fh, LL)
		mata: fputmatrix(fh, th)
		mata: fputmatrix(fh, st_matrixcolstripe("`start_par'"))
		mata: fclose(fh)
	}
	
	local bpotypev = "pfamoye"
	if "`bpotype'" != "" local bpotypev = "`bpotype'"
	
	*******************
	/*Display results*/
	*******************
	#delimit ;
	di _n in gr "Hidden Markov   " 
		_col(49) in gr "Number of obs" _col(67) "=" 
		_col(69) in ye %10.0fc `N' ;
	di in gr "bivariate Poisson regression"  in gr
		_col(49) "Number of groups" _col(67) "="
		_col(69) in ye %10.0fc `n' _n ;		
	di in gr "Group variable: " in ye abbrev("`id'",12) _col(49) as txt "Obs per group:" ;
	di as txt "Latent states = " as res %-10.0fc `k'
		_col(63) as txt "min" _col(67) "="
		_col(69) as res %10.0fc `g_min' ;
	di as txt " " as res %-10.0fc 
		_col(63) as txt "avg" _col(67) "="
		_col(69) as res %10.1fc `g_avg' ;
	di as txt "Biv. Poisson function: `bpotypev'" 
		_col(63) as txt "max" _col(67) "="
		_col(69) as res %10.0fc `g_max' _n ;
	di as txt "likelihood = " as res %10.3f `ll' /*10.3*/
		_col(62) as txt "BIC" _col(67) "="
		_col(69) as res %10.3fc `BIC'  ;	
	#delimit cr

	ereturn display
	ereturn repost, esample(`touse')
	

	foreach x in output	Phi	Pi	Pim	Piv	Pm	U	U_we_t	Ub	k0	lambda pu rnpar fh k0_1 k0_2 Phi_u Piu Pm_u Qt_u U1 U2 U2_we_t fd index intercept start_par_k0	x {
		capture mata: mata drop `x'
	}


end


/*
hmbivpois	(dep1 = cov11 cov12 Ldep1 Ldep2 ) /*
	*/		(dep2 = cov21 cov22 Ldep1 Ldep2 )  /*
	*/,  k(3) id(id) time(year) /*
	*/  initcov(Ldep1 Ldep2)   /*
	*/ optimize_init_technique(bfgs) optimize_init_evaluatortype(0) em(20) tracep /*
	*/ save(tempem_k3) bpotype(normal) xtrans(xtrans1 xtrans2) 
	
	load(tempem_k3)
	
	load(tempem_k3) eval
	
	prev(tempem_k2) 
*/

	