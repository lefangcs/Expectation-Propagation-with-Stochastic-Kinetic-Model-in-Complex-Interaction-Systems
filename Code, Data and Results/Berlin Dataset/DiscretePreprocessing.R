load("inference_0.RData")


person.state.d=person.state.d[201:1401,]
n.row=nrow(person.state.d)
n.col=ncol(person.state.d)
locations=sort(unique(as.vector(person.state.d)))

loc.d = t(apply(person.state.d, 1, function(x) table(factor(x, levels=locations ))))
n.row=nrow(person.state.d)
n.col=ncol(person.state.d)
person.state.d=matrix(match(as.vector(person.state.d), locations),nrow = n.row,ncol = n.col)


td=seq(from=1,to=nrow(person.state.d),length.out = 21)
m.time=array(0,c(length(locations),length(locations),20))
for(i in 1:20){
  print(i)
  p.s.d=person.state.d[td[i]:td[i+1],]
  nowloc=c(head(p.s.d,-1),1:length(locations))
  nextloc=c(tail(p.s.d,-1),1:length(locations))
  m=table(nowloc, nextloc)
  diag(m)=diag(m)-1
  diag(m)[rowSums(m)==0]=1
  m=sweep(m,1,rowSums(m),"/")
  m.time[,,i]=m
}


#only consider neighbors
rate_in=list()
rate_out=list()
loc_in=list()
loc_out=list()

for(i in 1:dim(m.time)[3]){
  print(i)
  m=m.time[,,i]
  diag(m)=0
  rownames(m)=1:length(locations) 
  colnames(m)=1:length(locations)
  
  rate_in[[i]]=lapply(1:ncol(m), function(n) {
    m[,n][m[,n]!=0]
  })
  loc_in[[i]]=lapply(1:ncol(m), function(n) {
    as.integer(names(rate_in[[i]][[n]]))
  })
  rate_out[[i]]=lapply(1:ncol(m), function(n) {
    m[n,][m[n,]!=0]
  })
  loc_out[[i]]=lapply(1:ncol(m), function(n) {
    as.integer(names(rate_out[[i]][[n]]))
  })
}

rate_in_f=function(i) rate_in[[ceiling((i)/(nrow(loc.d))*length(rate_in))]]
rate_out_f=function(i) rate_out[[ceiling((i)/(nrow(loc.d))*length(rate_out))]]
loc_in_f=function(i) loc_in[[ceiling((i)/(nrow(loc.d))*length(loc_in))]]
loc_out_f=function(i) loc_out[[ceiling((i)/(nrow(loc.d))*length(loc_out))]]


maxloc.d=apply(loc.d,2, max )
max.person=ifelse(maxloc.d<=10,maxloc.d+5,maxloc.d+10)

obs.matrix=obs.matrix[1:1500,1:1500]
obs.matrix[obs.matrix==0]=1e-20
obs.matrix=sweep(obs.matrix,1,rowSums(obs.matrix),FUN = '/')
rownames(obs.matrix)=as.character(1:1500)
colnames(obs.matrix)=as.character(1:1500)

#only main links are observable,colnames=
observable=650:900
unobservable=setdiff(1:length(locations),observable)


dataempty=lapply(1:nrow(loc.d), function(n){
  lapply(1:length(locations), function(m){
    rep(1,max.person[m]+1)
  })
})

sliceempty=lapply(1:length(locations), function(m){
  rep(0,max.person[m]+1)
})

start=sliceempty
for( i in 1:length(locations)) start[[i]][loc.d[1,i]+1]=1

end=sliceempty
for( i in 1:length(locations)) end[[i]][loc.d[nrow(loc.d),i]+1]=1

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


lg=la
for(i in 1:length(lg)){
  lg[[i]]=lapply(1:length(locations), function(n) la[[i]][[n]]*lb[[i]][[n]]/sum(la[[i]][[n]]*lb[[i]][[n]]))
} 


attr(la,'t') =attr(lb,'t') = attr(lg,'t') = 1:nrow(loc.d)
attr(la,'c')="a"
attr(lb,'c')="b"
attr(lg,'c')="a"

obs=loc.d
observable_nominal=as.character(observable)
observation=lapply(observable, function(n) obs[,n])
names(observation)=observable_nominal

remove(list = setdiff(ls(),c('observation','obs.matrix','lg','loc.d','rate_in','obs',
                             'rate_out','rate_in_f','rate_out_f',
                             'loc_in','loc_out','loc_in_f','loc_out_f',
                             'la','lb','m','m.time','max.person','observable_nominal','unobservable','observable','alloc','getSlice','locations')))

save.image(file = "inference_1.RData")