load("inference_1.RData")


observable=c(3,22)
unobservable=setdiff(1:length(locations),observable)

obs.prob=sapply(1:ncol(obs.matrix), function(n){
  m= n-1
  v= min( max(m^2/100,1e-6), 10)
  c(m,v)
})


dataempty=lapply(1:nrow(loc.d), function(n){
  rbind(max.person/2,rep(1e10,length(locations)))
})

dataempty=lapply(1:nrow(loc.d), function(n){
  rbind(0,rep(1e10,length(locations)))
})

sliceempty=rbind( rep( 0,length(locations)) ,rep(1e-6,length(locations)))

start=sliceempty
start[1,]=loc.d[1,]

end=sliceempty
end[1,]=loc.d[nrow(loc.d),]

la=dataempty
la[[1]]=start
lb=dataempty
lb[[length(lb)]]=end

alloc = function(x){
  old.t = attr(x,'t')
  old.c = attr(x,'c')
  if(length(attr(x, 't'))==length(x)) length(x) = length(x)*2 #alloc memory
  attr(x,'t') = old.t
  attr(x,'c') = old.c
  x
}

# read a slice from filtration, previous nearest one
getSlice <- function(x, t ){
  tt = attr(x, 't')
  
  if(attr(x,'c')=="a"){
    t0 = which(tt==max(tt[tt<=t]))
    y=x[[t0]]
  }
  if(attr(x,'c')=="b"){
    t0 = which(tt==min(tt[tt>=t]))
    y=x[[t0]]
  }
  y
}


normproduct_mv<-function(m1,v1,m2,v2){
  v=1/(1/v1+1/v2)
  m=v*(m1/v1+m2/v2)
  rbind(m,v)
}
lg=lapply(1:length(la), function(n) {
  a=la[[n]]
  b=lb[[n]]
  normproduct_mv(a[1,],a[2,],b[1,],b[2,])
} )



attr(la,'t') =attr(lb,'t') = attr(lg,'t') = 1:nrow(loc.d)
attr(la,'c')="a"
attr(lb,'c')="b"
attr(lg,'c')="a"

observable_nominal=as.character(observable)
observation=lapply(observable, function(n) obs[,n])
names(observation)=observable_nominal


remove(list = setdiff(ls(),c('observation','obs.matrix','lg','loc.d','rate_in','obs',
                             'rate_out','rate_in_f','rate_out_f','obs.prob','getss',
                             'loc_in','loc_out','loc_in_f','loc_out_f',
                             'la','lb','m.time','max.person','observable_nominal','unobservable','observable','alloc','getSlice','locations')))

save.image(file = "inference_2.RData")
