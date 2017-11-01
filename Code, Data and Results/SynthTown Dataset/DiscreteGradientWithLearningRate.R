#Using Gaussian Gradient algorithm with not aggressive learning rates
load("inference_1.RData")


upd_forward<-function(v1,inc,out,pn,dec,len){
  v2=v1*pn - 0:(len-1)*v1*out
  v2[2:len]=v2[2:len] + v1[1:(len-1)]*inc
  v2[1:(len-1)]=v2[1:(len-1)] + 1:(len-1)*v1[2:len]*dec
  v2
}

transition_forward_fra<-function(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2){
  m.inc=sapply(1:length(locations),function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
  m.eq=sapply(1:length(locations),function(n) sum(la1[[n]]*lb2[[n]]))
  m.eq[m.eq==0]=1e-20
  m.dec=sapply(1:length(locations),function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]] ))
  
  fra.inc=m.inc/m.eq
  fra.dec=m.dec/m.eq
  pinc=sapply(1:length(locations),function(n) sum(  ratein[[n]] * fra.dec[locin[[n]] ] )) 
  pdec= sapply(1:length(locations),function(n)  sum(  rateout[[n]] * fra.inc[locout[[n]]  ] ) )
  
  #for each link, calculate the prob of transition at all other links 
  tran=lapply(1:length(locations), function(n)  rateout[[n]] * fra.dec[n] * fra.inc[locout[[n]]  ] )
  alltran=sum(unlist(tran))
  trother=numeric(length = length(locations))
  trother[]=alltran
  trother=trother-sapply(tran,sum) # transition at other links = all transition - transition from local link - transition to local link
  for(n in 1:length(locations)) trother[locout[[n]]]=trother[locout[[n]]]-tran[[n]]
  pn=1-pnull+trother
  
  la2_tilde = lapply(1:length(locations), function(n) upd_forward(la1[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1) )
  lg2=lapply(1:length(locations), function(n) la2_tilde[[n]]*lb2[[n]] )
  K=sapply(lg2, sum )
  
  la2_tilde=lapply(1:length(locations), function(n) la2_tilde[[n]]/K[n])
  lg2=lapply(1:length(locations), function(n) lg2[[n]]/K[n] )
  
  list(la2_tilde=la2_tilde, lg2=lg2)
}

transition_forward_int<-function(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2){
  m.inc=numeric(length = length(locations))
  m.eq=numeric(length = length(locations))
  m.dec=numeric(length = length(locations))
  
  m.inc[unobservable]=sapply(unobservable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
  m.eq[unobservable]=sapply(unobservable,function(n) sum(la1[[n]]*lb2[[n]]))
  m.dec[unobservable]=sapply(unobservable,function(n) sum( 1:max.person[n] *la1[[n]][2:(max.person[n]+1)]* lb2[[n]][1:max.person[n]] ))
  
  m.inc[observable]=sapply(observable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)]  * obs.p2[[as.character(n) ]][2:(max.person[n]+1)] ))
  m.eq[observable]=sapply(observable,function(n) sum(la1[[n]]*lb2[[n]]*  obs.p2[[as.character(n) ]]))
  m.dec[observable]=sapply(observable,function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]]* obs.p2[[as.character(n) ]][1:max.person[n]] ) )
  
  m.eq[m.eq==0]=1e-20
  
  fra.inc=m.inc/m.eq
  fra.dec=m.dec/m.eq
  pinc=sapply(1:length(locations),function(n) sum(  ratein[[n]] * fra.dec[locin[[n]] ] )) 
  pdec= sapply(1:length(locations),function(n)  sum(  rateout[[n]] * fra.inc[locout[[n]]  ] ) )
  
  #for each link, calculate the prob of transition at all other links 
  tran=lapply(1:length(locations), function(n)  rateout[[n]] * fra.dec[n] * fra.inc[locout[[n]]  ] )
  alltran=sum(unlist(tran))
  trother=numeric(length = length(locations))
  trother[]=alltran
  trother=trother-sapply(tran,sum) # transition at other links = all transition - transition from local link - transition to local link
  for(n in 1:length(locations)) trother[locout[[n]]]=trother[locout[[n]]]-tran[[n]]
  pn=1-pnull+trother
  
  
  la2_tilde = lapply(1:length(locations), function(n) upd_forward(la1[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1) )
  
  la2_tilde[observable]=lapply(observable, function(n) la2_tilde[[n]] * obs.p2[[as.character(n) ]] )
  
  lg2=lapply(1:length(locations), function(n) la2_tilde[[n]]*lb2[[n]] )
  K=sapply(lg2, sum )
  
  la2_tilde=lapply(1:length(locations), function(n) la2_tilde[[n]]/K[n])
  lg2=lapply(1:length(locations), function(n) lg2[[n]]/K[n] )
  
  list(la2_tilde=la2_tilde, lg2=lg2)
}


upd_backward<-function(v1,inc,out,pn,dec,len){
  v2=v1*pn - 0:(len-1)*v1*out
  v2[1:(len-1)]=v2[1:(len-1)]+v1[2:len]*inc
  v2[2:len]=v2[2:len]+1:(len-1)*v1[1:(len-1)]*dec
  v2
}

transition_backward_fra<-function(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2){
  m.inc=sapply(1:length(locations),function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
  m.eq=sapply(1:length(locations),function(n) sum(la1[[n]]*lb2[[n]]))
  m.eq[m.eq==0]=1e-20
  m.dec=sapply(1:length(locations),function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]] ))
  
  fra.inc=m.inc/m.eq
  fra.dec=m.dec/m.eq
  pinc=sapply(1:length(locations),function(n) sum(  ratein[[n]] * fra.dec[locin[[n]] ] )) 
  pdec= sapply(1:length(locations),function(n)  sum(  rateout[[n]] * fra.inc[locout[[n]]  ] ) )
  
  #for each link, calculate the prob of transition at all other links 
  tran=lapply(1:length(locations), function(n)  rateout[[n]] * fra.dec[n] * fra.inc[locout[[n]]  ] )
  alltran=sum(unlist(tran))
  trother=numeric(length = length(locations))
  trother[]=alltran
  trother=trother-sapply(tran,sum) # transition at other links = all transition - transition from local link - transition to local link
  for(n in 1:length(locations)) trother[locout[[n]]]=trother[locout[[n]]]-tran[[n]]
  pn=1-pnull+trother
  
  lb1_tilde = lapply(1:length(locations), function(n) upd_backward(lb2[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1) )
  
  lg1=lapply(1:length(locations), function(n) la1[[n]]*lb1_tilde[[n]] )
  K=sapply(lg1, sum )
  
  lb1_tilde=lapply(1:length(locations), function(n) lb1_tilde[[n]]/K[n])
  lg1=lapply(1:length(locations), function(n) lg1[[n]] / K[n] )
  
  list(lb1_tilde=lb1_tilde, lg1=lg1)
}

transition_backward_int<-function(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2){
  m.inc=numeric(length = length(locations))
  m.eq=numeric(length = length(locations))
  m.dec=numeric(length = length(locations))
  
  m.inc[unobservable]=sapply(unobservable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
  m.eq[unobservable]=sapply(unobservable,function(n) sum(la1[[n]]*lb2[[n]]))
  m.dec[unobservable]=sapply(unobservable,function(n) sum( 1:max.person[n] *la1[[n]][2:(max.person[n]+1)]* lb2[[n]][1:max.person[n]] ))
  
  m.inc[observable]=sapply(observable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)]  * obs.p2[[as.character(n) ]][2:(max.person[n]+1)] ))
  m.eq[observable]=sapply(observable,function(n) sum(la1[[n]]*lb2[[n]]*  obs.p2[[as.character(n) ]]))
  m.dec[observable]=sapply(observable,function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]]* obs.p2[[as.character(n) ]][1:max.person[n]] ) )
  
  m.eq[m.eq==0]=1e-20
  
  fra.inc=m.inc/m.eq
  fra.dec=m.dec/m.eq
  pinc=sapply(1:length(locations),function(n) sum(  ratein[[n]] * fra.dec[locin[[n]] ] )) 
  pdec= sapply(1:length(locations),function(n)  sum(  rateout[[n]] * fra.inc[locout[[n]]  ] ) )
  
  #for each link, calculate the prob of transition at all other links 
  tran=lapply(1:length(locations), function(n)  rateout[[n]] * fra.dec[n] * fra.inc[locout[[n]]  ] )
  alltran=sum(unlist(tran))
  trother=numeric(length = length(locations))
  trother[]=alltran
  trother=trother-sapply(tran,sum) # transition at other links = all transition - transition from local link - transition to local link
  for(n in 1:length(locations)) trother[locout[[n]]]=trother[locout[[n]]]-tran[[n]]
  pn=1-pnull+trother
  
  lb2[observable]=lapply(observable, function(n) lb2[[n]] * obs.p2[[as.character(n) ]] )
  lb1_tilde = lapply(1:length(locations), function(n) upd_backward(lb2[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1) )
  
  lg1=lapply(1:length(locations), function(n) la1[[n]]*lb1_tilde[[n]] )
  K=sapply(lg1, sum )
  
  lb1_tilde=lapply(1:length(locations), function(n) lb1_tilde[[n]]/K[n])
  lg1=lapply(1:length(locations), function(n) lg1[[n]] / K[n] )
  
  list(lb1_tilde=lb1_tilde,lg1=lg1)
}


forward2 = function(la, lb, lg, obs, rate_in_f, rate_out_f, max.person ){
  
  la_old=la
  lb_old=lb
  
  new.t = c()
  
  length.la = length(la)
  length.lb = length(lb)
  length.lg = length(lg)
  
  
  for(i in 1:1440 ){ 
    #print(i)
    
    ratein=rate_in_f(i)
    rateout=rate_out_f(i)
    locin=loc_in_f(i)
    locout=loc_out_f(i)
    
    la1=la[[i]]
    lb2=lb_old[[i+1]]
    
    m.eq=numeric(length = length(locations))
    m.eq[unobservable]=sapply(unobservable,function(n) sum(la1[[n]]*lb2[[n]]))
    obs.p2=lapply(observable, function(n) obs.matrix[ 1:(max.person[n]+1) ,observation[[as.character(n)]][i+1]+1 ] )
    names(obs.p2)=observable_nominal
    m.eq[observable]=sapply(observable,function(n) sum(la1[[n]]*lb2[[n]]*obs.p2[[as.character(n)]]))
    m.eq[m.eq==0]=1e-20
    m.eq.x=numeric(length = length(locations))
    m.eq.x[unobservable]=sapply(unobservable,function(n) sum(0:max.person[n]*la1[[n]]*lb2[[n]]))
    m.eq.x[observable]=sapply(observable,function(n) sum(0:max.person[n]*la1[[n]]*lb2[[n]]*obs.p2[[as.character(n)]]))
    
    pout=sapply(rateout, sum)
    pnull= sum(pout*m.eq.x/m.eq) - pout*m.eq.x/m.eq
    r_nnn=max.person*pout+pnull
    nnn = max(ceiling(r_nnn))
    
    if(nnn<2) nnn=2
    
    pout=pout/nnn
    pnull=pnull/nnn
    ratein=lapply(1:length(ratein), function(n) ratein[[n]]/nnn)
    rateout=lapply(1:length(rateout), function(n) rateout[[n]]/nnn)
    
    nnnk=nnn
    if(i==1440) nnnk=nnn-1
    
    new.t=c(new.t, i+ 1:(nnn-1)/nnn )
    
    for (k in 1:nnnk){
      
      # inner loop, keep gamma and update delta
      t1 = i+(k-1)/nnn; t2 = i+k/nnn; t3=i+(k+1)/nnn
      
      if(k==nnn) tran1=get("transition_forward_int") else tran1=get("transition_forward_fra")
      
      if(k==nnn-1) {tran2=get("transition_backward_int"); obs.p3=obs.p2} else tran2=get("transition_backward_fra")
      
      la1=getSlice(la,t1); lb2=getSlice(lb_old,t2); 
      la2=getSlice(la_old,t2);  lb3=getSlice(lb_old,t3);
      
      lg2=lapply(1:length(locations), function(n) {
        gamma=la2[[n]]*lb2[[n]]
        gamma/sum(gamma)
      } )
      lg2_old=lg2
      
      epsilon=0.5
      
      for ( itera in 0:10 ) {
        
        # update delta over gradient
        for (inner in 0:10){
          #print(inner)
          transition1=tran1(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2)
          lg2_f=transition1$lg2
          la2_tilde=transition1$la2_tilde
          
          transition2=tran2(la2,lb3,ratein, locin, rateout, locout, pout, pnull,obs.p3)
          lg2_b=transition2$lg1
          lb2_tilde=transition2$lb1_tilde
          
          la2=lapply(1:length(locations), function(n) {
            la2[[n]][la2[[n]]<=1e-20]=1e-20
            la2[[n]]
          } )
          lb2=lapply(1:length(locations), function(n) {
            lb2[[n]][lb2[[n]]<=1e-20]=1e-20
            lb2[[n]]
          } )
          delta=lapply(1:length(locations), function(n) {
            la2[[n]]/lb2[[n]]
          } )
          
          la2_tilde=lapply(1:length(locations), function(n) {
            la2_tilde[[n]][la2_tilde[[n]]<=1e-20]=1e-20
            la2_tilde[[n]]
          } )
          lb2_tilde=lapply(1:length(locations), function(n) {
            lb2_tilde[[n]][lb2_tilde[[n]]<=1e-20]=1e-20
            lb2_tilde[[n]]
          } )
          delta_tilde=lapply(1:length(locations), function(n) {
            la2_tilde[[n]] / lb2_tilde[[n]]
          } )
          
          lg2_f=lapply(1:length(locations), function(n) {
            lg2_f[[n]][lg2_f[[n]]<=1e-20]=1e-20
            lg2_f[[n]]
          } )
          lg2_b=lapply(1:length(locations), function(n) {
            lg2_b[[n]][lg2_b[[n]]<=1e-20]=1e-20
            lg2_b[[n]]
          } )
          
          delta_new=lapply(1:length(locations), function(n) {
            #a=delta[[n]]* (exp(lg2_f[[n]]) / exp(lg2_b[[n]]) )^epsilon    ## true gradient
            a=delta[[n]]* ( delta_tilde[[n]] / delta[[n]] )^epsilon    ## approximate gradient
            a[a<=1e-20]=1e-20
            a[a>=1e20]=1e20
            a
          } )
          
          la2=lapply(1:length(locations), function(n) {
            sqrt(lg2[[n]] * delta_new[[n]] )
          } )
          lb2=lapply(1:length(locations), function(n) {
            sqrt(lg2[[n]] / delta_new[[n]] )
          } )
        }
        
        # update gamma over gradient
        for (outer in 0:10){

          transition1=tran1(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2)
          lg_f=transition1$lg2
          
          transition2=tran2(la2,lb3,ratein, locin, rateout, locout, pout, pnull,obs.p3)
          lg_b=transition2$lg1
          
          lg2_old=lg2
          lg2_old=lapply(1:length(locations), function(n) {
            lg2_old[[n]][lg2_old[[n]]<=1e-20]=1e-20
            lg2_old[[n]]
          } )
          
          lg_f=lapply(1:length(locations), function(n) {
            lg_f[[n]][lg_f[[n]]<=1e-20]=1e-20
            lg_f[[n]]
          } )
          lg_b=lapply(1:length(locations), function(n) {
            lg_b[[n]][lg_b[[n]]<=1e-20]=1e-20
            lg_b[[n]]
          } )
          lg2_tilde=lapply(1:length(locations), function(n) 1/2*lg_f[[n]] + 1/2*lg_b[[n]] )
          
          lg2= lapply(1:length(locations), function(n) {
            a=lg2_old[[n]] / ( exp(lg2_old[[n]]) / exp(lg2_tilde[[n]])    )^epsilon  ## true gradient
            #a=lg2_old[[n]] / ( lg2_old[[n]] / lg2_tilde[[n]]    )^epsilon  ## approximate gradient
            a[a<=1e-20]=1e-20
            a[a>=1e20]=1e20
            a
          } )
          la2=lapply(1:length(locations), function(n) sqrt(lg2[[n]] * delta_new[[n]] ) )
          lb2=lapply(1:length(locations), function(n) sqrt(lg2[[n]] / delta_new[[n]] ) )
        }
        
      }

      if(k==nnn){
        la[[i+1]]=la2
        lb[[i+1]]=lb2
        lg[[i+1]]=lg2
      } else{
        
        if(length(attr(lg,'t'))==length.lg){lg = alloc(lg); length.lg = length(lg)}
        if(min(abs(t2-attr(lg,'t')))<1e-6) {
          lg[[which.min(abs(t2-attr(lg,'t')))]] = lg2
        } else {
          attr(lg,'t') = c(attr(lg,'t'),t2);
          lg[[length(attr(lg,'t'))]]=lg2
        }
        
        if(length(attr(la,'t'))==length.la){la = alloc(la); length.la = length(la)}
        if(min(abs(t2-attr(la,'t')))<1e-6) {
          la[[which.min(abs(t2-attr(la,'t')))]] = la2
        } else {
          attr(la,'t') = c(attr(la,'t'),t2);
          la[[length(attr(la,'t'))]]=la2
        }
        
        if(length(attr(lb,'t'))==length.lb){lb = alloc(lb); length.lb = length(lb)}
        if(min(abs(t2-attr(lb,'t')))<1e-6) {
          lb[[which.min(abs(t2-attr(lb,'t')))]] <- lb2
        } else{
          attr(lb,'t') = c(attr(lb,'t'),t2)
          lb[[length(attr(lb,'t'))]]<-lb2
        }
        
      }
      
    } # k
    
  }
  
  new.t=c(1:1441,new.t)
  la = unclass(la)[match(new.t,attr(la,'t'))];  attr(la,'t') = new.t;  attr(la,'c')="a"
  lb = unclass(lb)[match(new.t,attr(lb,'t'))];  attr(lb,'t') = new.t;   attr(lb,'c')="b"
  lg = unclass(lg)[match(new.t,attr(lg,'t'))];   attr(lg,'t') = new.t;   attr(lg,'c')="a"
  
  list(la = la, lb = lb, lg = lg)
}


backward2 = function(la, lb, lg, obs, rate_in_f, rate_out_f, max.person){
  
  la_old=la
  lb_old=lb
  
  new.t = c()
  
  length.la = length(la) 
  length.lb = length(lb)
  length.lg = length(lg)
  
  
  for(i in 1440:1){
    #print(i)
    
    ratein=rate_in_f(i)
    rateout=rate_out_f(i)
    locin=loc_in_f(i)
    locout=loc_out_f(i)
    
    la1=la_old[[i]]
    lb2=lb[[i+1]]
    
    m.eq=numeric(length = length(locations))
    m.eq[unobservable]=sapply(unobservable,function(n) sum(la1[[n]]*lb2[[n]]))
    obs.p2=lapply(observable, function(n) obs.matrix[ 1:(max.person[n]+1) ,observation[[as.character(n)]][i+1]+1 ] )
    names(obs.p2)=observable_nominal
    m.eq[observable]=sapply(observable,function(n) sum(la1[[n]]*lb2[[n]]*obs.p2[[as.character(n)]]))
    m.eq[m.eq==0]=1e-20
    m.eq.x=numeric(length = length(locations))
    m.eq.x[unobservable]=sapply(unobservable,function(n) sum(0:max.person[n]*la1[[n]]*lb2[[n]]))
    m.eq.x[observable]=sapply(observable,function(n) sum(0:max.person[n]*la1[[n]]*lb2[[n]]*obs.p2[[as.character(n)]]))
    
    pout=sapply(rateout, sum)
    pnull= sum(pout*m.eq.x/m.eq) - pout*m.eq.x/m.eq
    r_nnn=max.person*pout+pnull
    nnn = max(ceiling(r_nnn))
    
    if(nnn<2) nnn=2
    
    pout=pout/nnn
    pnull=pnull/nnn
    ratein=lapply(1:length(ratein), function(n) ratein[[n]]/nnn)
    rateout=lapply(1:length(rateout), function(n) rateout[[n]]/nnn)
    
    nnn1=0
    if(i==1) nnn1 = 1
    
    new.t=c(new.t, i+ (nnn-1):1/nnn )
    
    for (k in (nnn-1):nnn1){
      #print(k)
      
      # inner loop, keep gamma and update delta
      t1 = i+(k-1)/nnn; t2 = i+k/nnn; t3=i+(k+1)/nnn
      
      if(k==nnn-1) {tran2=get("transition_backward_int"); obs.p3=obs.p2} else tran2=get("transition_backward_fra")
      
      if(k==0) {
        tran1=get("transition_forward_int")
        obs.p2=lapply(observable, function(n) obs.matrix[ 1:(max.person[n]+1) ,observation[[as.character(n)]][i]+1 ] )
        names(obs.p2)=observable_nominal
      } else tran1=get("transition_forward_fra")
      
      la1=getSlice(la_old,t1); lb2=getSlice(lb_old,t2); 
      la2=getSlice(la_old,t2);  lb3=getSlice(lb,t3);
      
      lg2=lapply(1:length(locations), function(n) {
        gamma=la2[[n]]*lb2[[n]]
        gamma/sum(gamma)
      } )
      lg2_old=lg2
      
      epsilon=0.5
      
      for ( itera in 0:10 ) {
        
        # update delta over gradient
        for (inner in 0:10){
          #print(inner)
          transition1=tran1(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2)
          lg2_f=transition1$lg2
          la2_tilde=transition1$la2_tilde
          
          transition2=tran2(la2,lb3,ratein, locin, rateout, locout, pout, pnull,obs.p3)
          lg2_b=transition2$lg1
          lb2_tilde=transition2$lb1_tilde
          
          la2=lapply(1:length(locations), function(n) {
            la2[[n]][la2[[n]]<=1e-20]=1e-20
            la2[[n]]
          } )
          lb2=lapply(1:length(locations), function(n) {
            lb2[[n]][lb2[[n]]<=1e-20]=1e-20
            lb2[[n]]
          } )
          delta=lapply(1:length(locations), function(n) {
            la2[[n]]/lb2[[n]]
          } )
          
          la2_tilde=lapply(1:length(locations), function(n) {
            la2_tilde[[n]][la2_tilde[[n]]<=1e-20]=1e-20
            la2_tilde[[n]]
          } )
          lb2_tilde=lapply(1:length(locations), function(n) {
            lb2_tilde[[n]][lb2_tilde[[n]]<=1e-20]=1e-20
            lb2_tilde[[n]]
          } )
          delta_tilde=lapply(1:length(locations), function(n) {
            la2_tilde[[n]] / lb2_tilde[[n]]
          } )
          
          lg2_f=lapply(1:length(locations), function(n) {
            lg2_f[[n]][lg2_f[[n]]<=1e-20]=1e-20
            lg2_f[[n]]
          } )
          lg2_b=lapply(1:length(locations), function(n) {
            lg2_b[[n]][lg2_b[[n]]<=1e-20]=1e-20
            lg2_b[[n]]
          } )
          
          delta_new=lapply(1:length(locations), function(n) {
            #a=delta[[n]]* (exp(lg2_f[[n]]) / exp(lg2_b[[n]]) )^epsilon    ## true gradient
            a=delta[[n]]* ( delta_tilde[[n]] / delta[[n]] )^epsilon    ## approximate gradient
            a[a<=1e-20]=1e-20
            a[a>=1e20]=1e20
            a
          } )
          
          la2=lapply(1:length(locations), function(n) {
            sqrt(lg2[[n]] * delta_new[[n]] )
          } )
          lb2=lapply(1:length(locations), function(n) {
            sqrt(lg2[[n]] / delta_new[[n]] )
          } )
        }
        
        # update gamma over gradient
        for (outer in 0:10){
          
          transition1=tran1(la1,lb2,ratein, locin, rateout, locout, pout, pnull,obs.p2)
          lg_f=transition1$lg2
          
          transition2=tran2(la2,lb3,ratein, locin, rateout, locout, pout, pnull,obs.p3)
          lg_b=transition2$lg1
          
          lg2_old=lg2
          lg2_old=lapply(1:length(locations), function(n) {
            lg2_old[[n]][lg2_old[[n]]<=1e-20]=1e-20
            lg2_old[[n]]
          } )
          
          lg_f=lapply(1:length(locations), function(n) {
            lg_f[[n]][lg_f[[n]]<=1e-20]=1e-20
            lg_f[[n]]
          } )
          lg_b=lapply(1:length(locations), function(n) {
            lg_b[[n]][lg_b[[n]]<=1e-20]=1e-20
            lg_b[[n]]
          } )
          lg2_tilde=lapply(1:length(locations), function(n) 1/2*lg_f[[n]] + 1/2*lg_b[[n]] )
          
          lg2= lapply(1:length(locations), function(n) {
            a=lg2_old[[n]] / ( exp(lg2_old[[n]]) / exp(lg2_tilde[[n]])    )^epsilon  ## true gradient
            #a=lg2_old[[n]] / ( lg2_old[[n]] / lg2_tilde[[n]]    )^epsilon  ## approximate gradient
            a[a<=1e-20]=1e-20
            a[a>=1e20]=1e20
            a
          } )
          la2=lapply(1:length(locations), function(n) sqrt(lg2[[n]] * delta_new[[n]] ) )
          lb2=lapply(1:length(locations), function(n) sqrt(lg2[[n]] / delta_new[[n]] ) )
        }
        
      }
      
      if(k==0){
        la[[i]]=la2
        lb[[i]]=lb2
        lg[[i]]=lg2
      } else {
        if(length(attr(lg,'t'))==length.lg){lg = alloc(lg); length.lg = length(lg)}
        if(min(abs(t2-attr(lg,'t')))<1e-6) {
          lg[[which.min(abs(t2-attr(lg,'t')))]] = lg2
        } else {
          attr(lg,'t') = c(attr(lg,'t'),t2);
          lg[[length(attr(lg,'t'))]]=lg2
        }
        
        if(length(attr(la,'t'))==length.la){la = alloc(la); length.la = length(la)}
        if(min(abs(t2-attr(la,'t')))<1e-6) {
          la[[which.min(abs(t2-attr(la,'t')))]] = la2
        } else {
          attr(la,'t') = c(attr(la,'t'),t2);
          la[[length(attr(la,'t'))]]=la2
        }
        
        if(length(attr(lb,'t'))==length.lb){lb = alloc(lb); length.lb = length(lb)}
        if(min(abs(t2-attr(lb,'t')))<1e-6) {
          lb[[which.min(abs(t2-attr(lb,'t')))]] <- lb2
        } else{
          attr(lb,'t') = c(attr(lb,'t'),t2)
          lb[[length(attr(lb,'t'))]]<-lb2
        }
      }
      
    } # i
    
  }
  new.t=c(1:1441, rev(new.t))
  la = unclass(la)[match(new.t,attr(la,'t'))];   attr(la,'t') = new.t;  attr(la,'c')="a"
  lb = unclass(lb)[match(new.t,attr(lb,'t'))];   attr(lb,'t') = new.t;   attr(lb,'c')="b"
  lg = unclass(lg)[match(new.t,attr(lg,'t'))];   attr(lg,'t') = new.t;   attr(lg,'c')="a"
  
  list(la = la,lb=lb,lg=lg)
}


#######################################################################################################
for(iter in 1:100){
  
  la_last=la
  lb_last=lb
  lg_last=lg
  
  aaa = forward2(la, lb, lg, obs, rate_in_f, rate_out_f, max.person)
  
  la=aaa$la
  lb=aaa$lb
  lg=aaa$lg
  
  lg_int = unclass(lg)[match(1:1441,attr(lg,'t'))]; attr(lg_int,'t') = 1:1441; attr(lg_int,'c') ="a"
  
  #  png(paste("The ", i, " forward.png" ))
  #plot the mean of lg
  layout(matrix(1:ncol(obs),ncol=1), heights=pmax(apply(obs,2,max),apply(loc.d,2,max))+1)
  par(mar=c(0,0,0,0),oma=c(5,2,0,1)+.1)
  for(ii in 1:ncol(obs)){
    plot(loc.d[,ii],type='l',col='black',xaxt='n',xlab='',ylab=colnames(obs)[ii],ylim = c(0,max(loc.d[,ii])+1))
    lines(sapply(1:(nrow(obs)), function(n){ 
      gamma=lg_int[[n]][[ii]]
      gamma=gamma/sum(gamma)
      sum(gamma* (0:(length(gamma)-1))) }),col="red",lty=1)
    if(ii==1) text(1300,10,paste(iter, " f lg" ),col = 'red',lwd=3)
  }
  axis(side=1)
  #  dev.off()
  print(sprintf('The %d forward completed',iter))
  
  
  
  
  la_last=la
  lb_last=lb
  lg_last=lg
  
  bbb = backward2(la, lb, lg, obs, rate_in_f, rate_out_f, max.person)
  
  la=bbb$la
  lb=bbb$lb
  lg=bbb$lg
  
  lg_int = unclass(lg)[match(1:1441,attr(lg,'t'))]; attr(lg_int,'t') = 1:1441;  attr(lg_int,'c')="a"
  
  #  png(paste("The ", i, " backward.png" ))
  #plot the mean of lg
  layout(matrix(1:ncol(obs),ncol=1), heights=pmax(apply(obs,2,max),apply(loc.d,2,max))+1)
  par(mar=c(0,0,0,0),oma=c(5,2,0,1)+.1)
  for(ii in 1:ncol(obs)){
    plot(loc.d[,ii],type='l',col='black',xaxt='n',xlab='',ylab=colnames(obs)[ii],ylim = c(0,max(loc.d[,ii])+1))
    lines(sapply(1:(nrow(obs)), function(n){ 
      gamma=lg_int[[n]][[ii]]
      gamma=gamma/sum(gamma)
      sum(gamma* (0:(length(gamma)-1))) }),col="red",lty=1)
    if(ii==1) text(1300,10,paste(iter, "b lg" ),col = 'red',lwd=3)
  }
  axis(side=1)
  #  dev.off()
  print(sprintf('The %d backward completed',iter))
}