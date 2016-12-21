gsSimulation<-function(N.group=10, N0.m=15, N0.f=45, b=.2, K=100, u.m=0.1, u.f=0.05, 
                       p=0.005, v.m=.5, v.f=.1, d.m=.001, d.f=.1, simtime=100, runs=50, seed=1)
{
  #Parameters
  #N.group: number of groups
  #b: birth rate
  #K: carrying capacity of each group
  #u.m: male natural mortality
  #u.f: female natural mortality
  #p: probably of epidemic outbreak in a group per time period
  #v.m: male mortality due to viral infection
  #v.f: female mortality due to viral infection
  #d.m: male dispersal rate
  #d.f: female dispersal rate
  #simtime: length of simulation
  
  mean.runs<- matrix(nrow=runs, ncol=15)
  
  for (i in 1:runs)
  { #Initial conditions
    N.m.mat<-N.f.mat<-matrix(rep(0, N.group*simtime), nrow = simtime)
    N.m<-rep(N0.m, N.group)
    N.f<-rep(N0.f, N.group)
    N.m.mat[1,]<-N.m
    N.f.mat[1,]<-N.f
    set.seed(seed+i)
    #Updating the group sizes
    for (t in 1:simtime) {
      B.m<-rpois(N.group, b*N.f*pmax(1-(N.f+N.m)/K,0))
      B.f<-rpois(N.group, b*N.f*pmax(1-(N.f+N.m)/K,0))
      N.m<-N.m+B.m
      N.f<-N.f+B.f
      infect<-(runif(N.group)<p*N.m)
      U.m<-rpois(N.group, u.m*N.m+infect*v.m*N.m)
      U.f<-rpois(N.group, u.f*N.f+infect*v.f*N.f)
      N.m<-pmax(N.m-U.m, 1)
      N.f<-pmax(N.f-U.f, 1)
      D.m.out<-pmin(rpois(N.group, d.m*N.m^2/N.f),N.m-1)
      D.f.out<-pmin(rpois(N.group, d.f*N.f^2/N.m),N.f-1)
      N.m<-N.m - D.m.out
      N.f<-N.f - D.f.out
      D.m.in<-rmultinom(1, sum(D.m.out), N.f^2/N.m/(sum(N.f^2/N.m)))
      D.f.in<-rmultinom(1, sum(D.f.out), N.m^2/N.f/(sum(N.m^2/N.f)))
      N.m<-N.m + D.m.in
      N.f<-N.f + D.f.in
      N.m.mat[t,]<-N.m
      N.f.mat[t,]<-N.f
    }
    
    #Compute variances within group
    Var.M<-apply(N.m.mat, 2, var)
    Var.F<-apply(N.f.mat, 2, var)
    Var.G<-apply(N.m.mat+N.f.mat, 2, var)
    mean.N.m<-apply(N.m.mat, 2, mean)
    mean.N.f<-apply(N.f.mat, 2, mean)
    mean.G<-apply(N.f.mat+N.m.mat, 2, mean)
    df.mean<-cbind("mean.m"=mean.N.m, "mean.f"=mean.N.f, "mean.g"=mean.G)
    df.var<-cbind("var.m"=Var.M, "var.f"=Var.F, "var.g"=Var.G)
    df.cov<-cbind("cov.m"=sqrt(Var.M)/mean.N.m, "cov.f"=sqrt(Var.F)/mean.N.f, "cov.g"=sqrt(Var.G)/mean.G)
    df.summary<-cbind(df.mean, df.var, df.cov)
    mean.runs[i,1:9]<-colMeans(df.summary)
    
    #Compute variances between groups
    Var.M.bt<-apply(N.m.mat, 1, var)
    Var.F.bt<-apply(N.f.mat, 1, var)
    Var.G.bt<-apply(N.m.mat+N.f.mat, 1, var)
    mean.m.bt<-apply(N.m.mat, 1, mean)
    mean.f.bt<-apply(N.f.mat, 1, mean)
    mean.g.bt<-apply(N.f.mat+N.m.mat, 2, mean)
    df.mean.bt<-cbind("mean.m"=mean.m.bt, "mean.f"=mean.f.bt, "mean.g"=mean.g.bt)
    df.var.bt<-cbind("var.m"=Var.M.bt, "var.f"=Var.F.bt, "var.g"=Var.G.bt)
    df.cov.bt<-cbind("cov.m"=sqrt(Var.M.bt)/mean.m.bt, "cov.f"=sqrt(Var.F.bt)/mean.f.bt, "cov.g"=sqrt(Var.G.bt)/mean.g.bt)
    df.summary.bt<-cbind(df.var.bt, df.cov.bt)
    mean.runs[i,10:15]<-colMeans(df.summary.bt)
  }
  mean.runs<-as.data.frame(mean.runs)
  names(mean.runs) = c("mean.m", "mean.f", "mean.g", "var.m.w", "var.f.w", "var.g.w", "cov.m.w", "cov.f.w", "cov.g.w",
                       "var.m.b", "var.f.b", "var.g.b", "cov.m.b", "cov.f.b", "cov.g.b")
  return(mean.runs)
}