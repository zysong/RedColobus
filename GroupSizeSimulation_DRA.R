#Parameters
set.seed(100)
N.group<- 10
b<-.02
K<-1000
u.m<-.01
u.f<-.005
p<-.001
v.m<-.5
v.f<-.1
d.m<-.0001
d.f<-.01
simtime<-1000

#Initial conditions
N.m.mat<-N.f.mat<-matrix(rep(0, N.group*simtime), nrow = simtime)
N.m<-rep(20, N.group)
N.f<-rep(50, N.group)
N.m.mat[1,]<-N.m
N.f.mat[1,]<-N.f

#Updating the group sizes
for (t in 1:simtime) {
  B.m<-rpois(N.group, b*N.f*pmax((1-sum(N.m)*(N.f+N.m)/K/N.m),0))
  B.f<-rpois(N.group, b*N.f*pmax((1-sum(N.m)*(N.f+N.m)/K/N.m),0))
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

#Compute variances
Var.M<-apply(N.m.mat, 2, var)
Var.F<-apply(N.f.mat, 2, var)
Var.G<-apply(N.m.mat+N.f.mat, 2, var)
mean.N.m<-apply(N.m.mat, 2, mean)
mean.N.f<-apply(N.f.mat, 2, mean)
mean.G<-apply(N.f.mat+N.m.mat, 2, mean)
df.var<-cbind("Male"=Var.M, "Female"=Var.F, "Group"=Var.G)
df.cov<-cbind("Male"=sqrt(Var.M)/mean.N.m, "Female"=sqrt(Var.F)/mean.N.f, "Group"=sqrt(Var.G)/mean.G)

pdf("../Manuscript/figures/A1B1plot%01d.pdf", onefile=FALSE)
#plot time series
plot.m<-matplot(N.m.mat, ylim = c(0, 100), xlab="Time", ylab="Male group size", type = "l", main="FBD+HVI", cex.lab = 1.5);
plot.f<-matplot(N.f.mat, ylim = c(0, 250), xlab="Time", ylab="Female group size", type = "l", main="FBD+HVI", cex.lab = 1.5);
plot.g<-matplot(N.f.mat+N.m.mat, ylim = c(0, 300), xlab="Time", ylab="Group size", type = "l", main="FBD+HVI", cex.lab = 1.5);
plot.r<-matplot(N.m.mat/N.f.mat, ylim = c(0, 1), xlab="Time", ylab = "Sex ratio", type = "l", main="FBD+HVI", cex.lab = 1.5);

#plot variances
plot.var<-boxplot(df.var, ylim=c(0, 3000), ylab="Variance", main="FBD+HVI", cex.lab = 1.5);
plot.cov<-boxplot(df.cov, ylim=c(0, 1.1), ylab="Coefficient of variation", main="FBD+HVI", cex.lab = 1.5);

dev.off()