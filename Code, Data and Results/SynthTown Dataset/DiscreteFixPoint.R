#Using discrete fix point algorithm
load("inference_1.RData")

upd_forward<-function(v1,inc,out,pn,dec,len){
  v2=v1*pn - 0:(len-1)*v1*out
  v2[2:len]=v2[2:len] + v1[1:(len-1)]*inc
  v2[1:(len-1)]=v2[1:(len-1)] + 1:(len-1)*v1[2:len]*dec
  v2[v2<0]=0
  v2[v2>1e200]=1e200
  v2
}

transition_forward_fra<-function(la1,lb2,ratein, locin, rateout, locout, pout, pn){
  m.inc=sapply(1:length(locations),function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
  m.eq=sapply(1:length(locations),function(n) sum(la1[[n]]*lb2[[n]]))
  m.eq[m.eq==0]=1e-20
  m.dec=sapply(1:length(locations),function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]] ))
  
  fra.inc=m.inc/m.eq
  fra.dec=m.dec/m.eq
  pinc=sapply(1:length(locations),function(n) sum(  ratein[[n]] * fra.dec[locin[[n]] ] )) 
  pdec= sapply(1:length(locations),function(n)  sum(  rateout[[n]] * fra.inc[locout[[n]]  ] ) )

  la2_tilde = lapply(1:length(locations), function(n) {
    v=upd_forward(la1[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1)
    v=v/sum(v*lb2[[n]])
    v
  } )
  list(la2_tilde=la2_tilde)
}

transition_forward_int<-function(la1,lb2,ratein, locin, rateout, locout, pout, pn,obs.p2){
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
  
  la2_tilde = lapply(1:length(locations), function(n) upd_forward(la1[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1) )
  la2_tilde[observable]=lapply(observable, function(n) la2_tilde[[n]] * obs.p2[[as.character(n) ]] )
  la2_tilde = lapply(1:length(locations), function(n) la2_tilde[[n]]/ sum(la2_tilde[[n]]*lb2[[n]]) )
  list(la2_tilde=la2_tilde)
}


upd_backward<-function(v1,inc,out,pn,dec,len){
  v2=v1*pn - 0:(len-1)*v1*out
  v2[1:(len-1)]=v2[1:(len-1)]+v1[2:len]*inc
  v2[2:len]=v2[2:len]+1:(len-1)*v1[1:(len-1)]*dec
  v2[v2<0]=0
  v2[v2>1e200]=1e200
  v2
}

transition_backward_fra<-function(la1,lb2,ratein, locin, rateout, locout, pout, pn){
  m.inc=sapply(1:length(locations),function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
  m.eq=sapply(1:length(locations),function(n) sum(la1[[n]]*lb2[[n]]))
  m.eq[m.eq==0]=1e-20
  m.dec=sapply(1:length(locations),function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]] ))
  
  fra.inc=m.inc/m.eq
  fra.dec=m.dec/m.eq
  pinc=sapply(1:length(locations),function(n) sum(  ratein[[n]] * fra.dec[locin[[n]] ] )) 
  pdec= sapply(1:length(locations),function(n)  sum(  rateout[[n]] * fra.inc[locout[[n]]  ] ) )
  
  lb1_tilde = lapply(1:length(locations), function(n) {
    v=upd_backward(lb2[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1)
    v=v/sum(v*la1[[n]])
    v
  } )

  list(lb1_tilde=lb1_tilde)
}

transition_backward_int<-function(la1,lb2,ratein, locin, rateout, locout, pout, pn,obs.p2){
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
  
  lb2[observable]=lapply(observable, function(n) lb2[[n]] * obs.p2[[as.character(n) ]] )
  lb1_tilde = lapply(1:length(locations), function(n) {
    v=upd_backward(lb2[[n]],pinc[n],pout[n],pn[n],pdec[n],max.person[n]+1)
    v=v/sum(v*la1[[n]])
    v
  } )
  list(lb1_tilde=lb1_tilde)
}



forward2 = function(la, lb, obs, rate_in_f, rate_out_f, max.person){
  new.t = c()
  length.la = length(la)
  
  for(i in 1:1440){
    ratein=rate_in_f(i) # ratein is a list, each element stores the rate constant of the cars moving from its neighbors to the link 
    rateout=rate_out_f(i) # rateout is a list, each element stores the rate constant of the cars moving from the link to its neighbors
    locin=loc_in_f(i)
    locout=loc_out_f(i)
    
    la1=la[[i]]
    lb1=lb[[i]]
    estimates=sapply(1:length(locations), function(n){ 
      gamma=la1[[n]]*lb1[[n]]
      gamma=gamma/sum(gamma)
      sum(gamma* (0:max.person[n]) ) })
    lb2=lb[[i+1]]

    #############  substeps
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

    m.inc=numeric(length = length(locations))
    m.dec=numeric(length = length(locations))
    m.inc[unobservable]=sapply(unobservable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
    m.dec[unobservable]=sapply(unobservable,function(n) sum( 1:max.person[n] *la1[[n]][2:(max.person[n]+1)]* lb2[[n]][1:max.person[n]] ))
    m.inc[observable]=sapply(observable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)]  * obs.p2[[as.character(n) ]][2:(max.person[n]+1)] ))
    m.dec[observable]=sapply(observable,function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]]* obs.p2[[as.character(n) ]][1:max.person[n]] ) )
    
    fra.inc=m.inc/m.eq
    fra.dec=m.dec/m.eq
    #for each link, calculate the prob of transition at all other links 
    tran=lapply(1:length(locations), function(n)  rateout[[n]] * fra.dec[n] * fra.inc[locout[[n]]  ] )
    alltran=sum(unlist(tran))
    trother=numeric(length = length(locations))
    trother[]=alltran
    trother=trother-sapply(tran,sum) # transition at other links = all transition - transition from local link - transition to local link
    for(n in 1:length(locations)) trother[locout[[n]]]=trother[locout[[n]]]-tran[[n]]
    
    r_nnn=estimates*pout+pnull-trother
    nnn = max(ceiling(r_nnn))
    if(nnn<2) nnn=2
    
    pn=1-pnull/nnn+trother/nnn
    pout=pout/nnn
    ratein=lapply(1:length(ratein), function(n) ratein[[n]]/nnn)
    rateout=lapply(1:length(rateout), function(n) rateout[[n]]/nnn)

    new.t=c(new.t,i+1:(nnn-1)/nnn)
    
    for (k in 1:nnn){
      t1 = i+(k-1)/nnn; t2 = i+k/nnn;
      lb2=getSlice(lb,t2);
      
      if(k!=nnn) {
        tran=transition_forward_fra(la1,lb2,ratein, locin, rateout, locout, pout, pn)
        la2=tran$la2_tilde
        
        if(length(attr(la,'t'))==length.la){la = alloc(la); length.la = length(la)}
        if(min(abs(t2-attr(la,'t')))<1e-6) {
          la[[which.min(abs(t2-attr(la,'t')))]] = la2
        } else {
          attr(la,'t') = c(attr(la,'t'),t2);
          la[[length(attr(la,'t'))]]=la2
        }
        
      } else {
        tran=transition_forward_int(la1,lb2,ratein, locin, rateout, locout, pout, pn,obs.p2)
        la2=tran$la2_tilde
        la[[i+1]]=la2
      }
      
      la1=la2
    } # k
  }

  new.t=c(1:1441,new.t)
  la = unclass(la)[match(new.t,attr(la,'t'))]; attr(la,'t') = new.t;  attr(la,'c')="a"
  list(la = la)
}

backward2 = function(la, lb, obs, rate_in_f, rate_out_f, max.person){
  
  new.t = c()
  length.lb = length(lb);
  
  for(i in 1440:1 ){
    ratein=rate_in_f(i) # ratein is a list, each element stores the rate constant of the cars moving from its neighbors to the link 
    rateout=rate_out_f(i) # rateout is a list, each element stores the rate constant of the cars moving from the link to its neighbors
    locin=loc_in_f(i)
    locout=loc_out_f(i)
    
    la1=la[[i]]
    lb1=lb[[i]]
    estimates=sapply(1:length(locations), function(n){ 
      gamma=la1[[n]]*lb1[[n]]
      gamma=gamma/sum(gamma)
      sum(gamma* (0:max.person[n]) ) })
    lb2=lb[[i+1]]
    
    #############  substeps
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
    
    m.inc=numeric(length = length(locations))
    m.dec=numeric(length = length(locations))
    m.inc[unobservable]=sapply(unobservable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)] ) )
    m.dec[unobservable]=sapply(unobservable,function(n) sum( 1:max.person[n] *la1[[n]][2:(max.person[n]+1)]* lb2[[n]][1:max.person[n]] ))
    m.inc[observable]=sapply(observable,function(n) sum( la1[[n]][1:max.person[n]] * lb2[[n]][2:(max.person[n]+1)]  * obs.p2[[as.character(n) ]][2:(max.person[n]+1)] ))
    m.dec[observable]=sapply(observable,function(n) sum( 1:max.person[n] * la1[[n]][2:(max.person[n]+1)] * lb2[[n]][1:max.person[n]]* obs.p2[[as.character(n) ]][1:max.person[n]] ) )
    
    fra.inc=m.inc/m.eq
    fra.dec=m.dec/m.eq
    #for each link, calculate the prob of transition at all other links 
    tran=lapply(1:length(locations), function(n)  rateout[[n]] * fra.dec[n] * fra.inc[locout[[n]]  ] )
    alltran=sum(unlist(tran))
    trother=numeric(length = length(locations))
    trother[]=alltran
    trother=trother-sapply(tran,sum) # transition at other links = all transition - transition from local link - transition to local link
    for(n in 1:length(locations)) trother[locout[[n]]]=trother[locout[[n]]]-tran[[n]]
    
    r_nnn=estimates*pout+pnull-trother
    nnn = max(ceiling(r_nnn))
    if(nnn<2) nnn=2
    
    pn=1-pnull/nnn+trother/nnn
    pout=pout/nnn
    ratein=lapply(1:length(ratein), function(n) ratein[[n]]/nnn)
    rateout=lapply(1:length(rateout), function(n) rateout[[n]]/nnn)
    #############

    if(nnn>1) new.t=c(new.t,i+(nnn-1):1/nnn)
    
    for (k in nnn:1){
      t1 = i+(k-1)/nnn; t2 = i+k/nnn;
      la1=getSlice(la,t1)
      
      if(k!=nnn) {
        tran=transition_backward_fra(la1,lb2,ratein, locin, rateout, locout, pout, pn)        
      } else {
        tran=transition_backward_int(la1,lb2,ratein, locin, rateout, locout, pout, pn,obs.p2)
      }
      
      lb1=tran$lb1_tilde
      lb2 = lb1

      if(k==1){
        lb[[i]]=lb1
      } else{
        if(length(attr(lb,'t'))==length.lb){lb = alloc(lb); length.lb = length(lb)}
        if(min(abs(t1-attr(lb,'t')))<1e-12) {
          lb[[which.min(abs(t1-attr(lb,'t')))]] <- lb1
        } else{
          attr(lb,'t') = c(attr(lb,'t'),t1)
          lb[[length(attr(lb,'t'))]]<-lb1
        }
      }
    }
  }
  
  new.t=c(1:1441,rev(new.t) )
  lb = unclass(lb)[match(new.t,attr(lb,'t'))]; attr(lb,'t') = new.t;  attr(lb,'c')="b"
  list(lb = lb)
}


#######################################################################################################

for(iter in 1:100){
  print(Sys.time())
  aaa = forward2(la, lb, obs, rate_in_f, rate_out_f, max.person)
  print(Sys.time())
  la=aaa$la
  
  la_int = unclass(la)[match(1:1441,attr(la,'t'))]; attr(la_int,'t') = 1:1441;  attr(la_int,'c')="a"
  lb_int = unclass(lb)[match(1:1441,attr(lb,'t'))]; attr(lb_int,'t') = 1:1441; attr(lb_int,'c') ="b"
  
  #plot the mean of la*lb
  layout(matrix(1:ncol(obs),ncol=1), heights=pmax(apply(obs,2,max),apply(loc.d,2,max))+1)
  par(mar=c(0,0,0,0),oma=c(5,2,0,1)+.1)
  for(ii in 1:ncol(obs)){
    plot(loc.d[,ii],type='l',col='black',xaxt='n',xlab='',ylab=colnames(obs)[ii],ylim = c(0,max(loc.d[,ii])+1))
    lines(sapply(1:(nrow(obs)), function(n){ 
      gamma=la_int[[n]][[ii]]*lb_int[[n]][[ii]]
      gamma=gamma/sum(gamma)
      sum(gamma* (0:(length(gamma)-1)) ) }),col="red",lty=1)
    if(ii==1) text(750,1800,paste(iter, " f lg" ),cex=1.5 ,col = 'red',lwd=3)
  }
  axis(side=1)

  print(sprintf('The %d forward completed',iter))

  
  print(Sys.time())
  bbb = backward2(la, lb, obs, rate_in_f, rate_out_f, max.person)
  print(Sys.time())
  lb=bbb$lb
  
  la_int = unclass(la)[match(1:1441,attr(la,'t'))]; attr(la_int,'t') = 1:1441;  attr(la_int,'c')="a"
  lb_int = unclass(lb)[match(1:1441,attr(lb,'t'))]; attr(lb_int,'t') = 1:1441; attr(lb_int,'c') ="b"
  
  #plot the mean of la*lb
  layout(matrix(1:ncol(obs),ncol=1), heights=pmax(apply(obs,2,max),apply(loc.d,2,max))+1)
  par(mar=c(0,0,0,0),oma=c(5,2,0,1)+.1)
  for(ii in 1:ncol(obs)){
    plot(loc.d[,ii],type='l',col='black',xaxt='n',xlab='',ylab=colnames(obs)[ii],ylim = c(0,max(loc.d[,ii])+1))
    lines(sapply(1:(nrow(obs)), function(n){ 
      gamma=la_int[[n]][[ii]]*lb_int[[n]][[ii]]
      gamma=gamma/sum(gamma)
      sum(gamma* (0:(length(gamma)-1))) }),col="red",lty=1)
    if(ii==1) text(750,1800,paste(iter, "b lg"), cex=1.5, col = 'red',lwd=3)
  }
  axis(side=1)

  print(sprintf('The %d backward completed',iter))
}
