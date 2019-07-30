############# libraries ############
library(MASS)
library(plot3D)
library(class)
library(irlba)
############### Functions: Calculating PCA error and SPCA error ########################
SPCA=function(X,d)
{
  dm=dim(as.matrix(X)); n=dm[1]; m=dm[2];

  if(n<=max((d+2),5))
  {  
    cat("Not Enough Samples \n")
    SS=SS_new=0;
    if(min(dim(X))==1)
    {
      c=X;
    }
    else
    {
      c=colMeans(X);
    }
    r=sqrt(sum(X[1,]-c)^2)
    pca=prcomp(X,retx=TRUE,center=TRUE)
    part=pca$x[,1]
    coeff=pca[["rotation"]] # pXp matrix
    if(n>=(d+1))
    {
      V=coeff[,(1:min((d+1),ncol(coeff)))]
    }
    else
    {
      V=cbind(coeff,matrix(data=0,m,(d+1-n)))
    }
    part=c(rep(1,floor(n/2)),rep(-1,ceiling(n/2)))
    mu=c
    c.d=t(V)%*%c
   }

  if(n>max((d+1),5))   # if there are enough samples, fit the data by a sphere or a hyperplane
  {# do d+1 dimensional PCA first
	  mu=colMeans(X);   centeredX=sapply(1:m,function(u){X[,u]-mu[u]}); # substract from each column
	  
	pca=prcomp(X,retx=TRUE,center=TRUE)
	part=pca$x[,1]
	coeff=pca[["rotation"]] # pXp matrix
	explained=100*(pca$sdev)^2/sum((pca$sdev)^2) # explained variance by each principal comp

	V=coeff[,(1:min((d+1),ncol(coeff)))];
	Y_t=pca$x[,(1:min((d+1),m))]  # Y_tride : matrix of V_(d+1)^T (X_i- X_bar)
    	Y=as.matrix(rep(1,n))%*%mu+(Y_t%*%t(V)); # Y: projection of X onto the d+1 dimensional affine space mu+V
    	Y_norm=apply(Y,1,function(u){return(sum(u^2))})  # vector of norm of each Yi
    	meanY_norm=mean(Y_norm);   # mean Y_norm

    	H.list=apply(Y,1,function(u){RR=as.vector(mu-u)%*%t(as.vector(mu-u)); return(list(RR))})  # should be n-length list of mXm matrices of (Yi-X_bar)(Y_i-X_bar)^t for i=1:n
    	H.list=lapply(H.list,data.frame) # transforms each component matrices to data.frames
    	H=Reduce("+",H.list)  # component-wise sum of the data-frame
    	H=as.matrix(H)

    	f=rowSums(sapply(1:n,function(u){rr=(Y_norm[u]-meanY_norm)*(mu-Y[u,]); return(rr) } ))

    	c.d=t(V)%*%(-0.5*ginv(H)%*%f-mu)  # (d+1)-dimensional center
    	c=mu+V%*%c.d;  # center of the sphere in D-dimensional space

    	Riemd=apply(Y,1,function(u){RR=sqrt(sum((c-u)^2)); return(RR)})
    	r=mean(Riemd) # radius
    	SS=sum(sapply(1:n,function(u){rr=sum((X[u,]-Y[u,])^2); return(rr)}))
    	SS=SS+n*var(Riemd); # second part of the error(spherical error)
    	V_new=coeff[,1:d];
    	SS_new=sum(explained[(d+1):length(explained)])
  }
  return(list(c,V,r,SS,SS_new,part,mu,c.d))
}

################################################################################################################################

SS_calc=function(X,mu,c,V,r,d,c.d)
{ 
  n=nrow(X); p=ncol(X);
  d.opt=ncol(V)
  
  Vd=V[,1:d]
  Y.hatd.pc=(X-as.matrix(rep(1,n))%*%mu)%*%Vd
  Y.hatd=cbind(Y.hatd.pc,rep(0,nrow(Y.hatd.pc)))
  Y=as.matrix(rep(1,n))%*%mu+Y.hatd.pc%*%t(Vd); # projection of X onto the d+1 dimensional affine space mu+V
  SS_new=sum(sapply(1:n,function(u){rr=sum((X[u,]-Y[u,])^2); return(rr)}))
  Y.hatD=Y.hatD.pc=Y

  Y.tride=(X-as.matrix(rep(1,n))%*%mu)%*%V    # nX(d+1) matrix of V(d+1)^T (Xi-X.bar) s
  Y.hatd=t(sapply(1:n,function(u){Vy=(Y.tride[u,]-c.d); normVy=as.numeric(sqrt(t(Vy)%*%(Vy)));  rr=c.d+r*(Vy)/normVy; return(rr)}))
  Y.hatD=t(sapply(1:n, function(u) {VVty=V%*%t(V)%*%(X[u,]-c); normVVt=as.numeric(sqrt(t(VVty)%*%(VVty))); rr=c+r*(VVty)/normVVt; return(rr)}))
  SS=sum(sapply(1:n,function(u){rr=sum((X[u,]-Y.hatD[u,])^2); return(rr)}))
  return(list(SS,SS_new,Y.hatD,Y.hatd,Y.hatD.pc,Y.hatd.pc))
} 
 
##################################################################################################################################
SPCA_Manifold_approx=function(x,y,d,epsilon,epsilon2,N.min)
{
  
  n=nrow(x); p=ncol(x); nt=nrow(y);
  k=1 ## number of partitions is 2^k
  center.list=list()
  r.list=list()
  V.list=list()
  mu.list=list()
  indx.list=list()
  num.indx.list=list()
  centerd.list=list()
  
  n.part1=sp.err=sp.err1=NULL
  
  test.error.list=list()
  test.error.list.pc=list()
  y.hatD.list=list()
  y.hatd.list=list()
  y.hatD.pc.list=list()
  y.hatd.pc.list=list()
  
  Yindx.list=list()
  num.Yindx.list=list()
  
  ############# K=1: No partition ###########
  # Train
  mse=NULL; n.part=NULL;
  spcax=SPCA(x,d)
  center.list[[k]]=c=spcax[[1]];  V.list[[k]]=V=as.matrix(spcax[[2]]);  r.list[[k]]=r=spcax[[3]];  
  mu.list[[k]]=mu=spcax[[7]]; centerd.list[[k]]=cd=spcax[[8]];
  mse[1]=pc.error=log10(spcax[[4]]/n); n.part[1]=1
  
  # Test
  n.part1[1]=1
  spcayd=SS_calc(y,mu,c,V,r,d,cd)
  test.error.list[[1]]=spcayd[[1]]; erro=sqrt(test.error.list[[1]]/nt); sp.err[1]=log10(erro)
  test.error.list.pc[[1]]=spcayd[[2]];  erro1=sqrt(test.error.list.pc[[1]]/nt); sp.err1[1]=log10(erro1) # keep this to compare with PCA
  
  y.hatD.list[[1]]=Y.hatD=spcayd[[3]]; y.hatd.list[[1]]=Y.hatd=spcayd[[4]]; 
  y.hatD.pc.list[[1]]=Y.hatD.pc=spcayd[[5]]; y.hatd.pc.list[[1]]=Y.hatd.pc=spcayd[[6]]; 
  
  ############# Partitioning starts ################
  k=k+1; 
  INDX=INDX.new=rep(1,n); max.cls=n;
  YINDX=YINDX.new=rep(1,nt);
  
  while((pc.error>epsilon)&(max.cls>N.min))
  {
    print(k)
    # Train
    Lev_INDX=levels(as.factor(INDX))
    L=length(Lev_INDX)
    INDX.num=array(data=NA,dim=(2*L))
    CENTER=matrix(data=NA,nrow=p,ncol=(2*L))
    CENTERd=matrix(data=NA,nrow=(d+1),ncol=(2*L))
    V.mat=array(data=NA,dim=c(p,(d+1),(2*L))) 
    R.array=array(data=NA,dim=(2*L))
    mu.array=matrix(data=NA,nrow=p,ncol=(2*L))
    error=array(data=NA,dim=(2*L))
    
    # Test
    error.y=array(data=0,dim=(2*L))
    error.y.pc=array(data=0,dim=(2*L))
    ##
    Flag=NULL; count=0;
    num=1
    for(i in 1:L)
    {
      Flag[i]=1
      clust=Lev_INDX[i]
      samps=which(INDX==clust)
      x.use=x[samps,]
      spcaxuse=SPCA(x.use,d)  
      error.use=spcaxuse[[4]]
      center.use=spcaxuse[[1]];  V.use=as.matrix(spcaxuse[[2]]);  
      r.use=spcaxuse[[3]];  mu.use=spcaxuse[[7]]; centerd.use=spcaxuse[[8]];
      part.old=spcaxuse[[6]]
      pc.error.use=log10(error.use/length(samps))
      
      samp.y=which(YINDX==clust)
      if(length(samp.y)>0)
      {
        y.use=as.matrix(y[samp.y,])
        if(min(dim(y.use))==1)
        {
          y.use=t(y.use)
        }
        spcayd=SS_calc(y.use,mu.use,center.use,V.use,r.use,d,centerd.use)         
        
        error.use.y=spcayd[[1]]
        error.use.y.pc=spcayd[[2]]
        Y.hatD[samp.y,]=spcayd[[3]]; Y.hatd[samp.y,]=spcayd[[4]]; 
        Y.hatD.pc[samp.y,]=spcayd[[5]]; Y.hatd.pc[samp.y,]=spcayd[[6]]; 
      }
      else 
      {
        Flag[i]=0
      }
      
      if((pc.error.use>epsilon)&(nrow(x.use)>max(N.min,max(d,p)))) 
      {
        idx=which(part.old>=0)
        z=as.matrix(x.use[idx,])
        idx2=setdiff((1:nrow(x.use)),idx)
        w=as.matrix(x.use[idx2,])
        
        if(length(samp.y)>0)
        {
          y.ct=t(apply(y.use,1,function(u){u=(u-mu.use)}))
          y.part=y.ct%*%V.use[,1]
          idx.y=which(y.part>=0)
          if(length(idx.y)>0)
          { zy=as.matrix(y.use[idx.y,]);         if(min(dim(zy))==1) {zy=t(zy)}  }
          else
          {  count=count+1 }
          
          idx2.y=setdiff((1:nrow(y.use)),idx.y)
          if(length(idx2.y)>0)
          {  wy=as.matrix(y.use[idx2.y,]);         if(min(dim(wy))==1) {wy=t(wy)}  }
          else
          {  count=count+1  }
        }
        
        if(min(nrow(z),nrow(w))<6)
        {
          INDX.new[samps]=as.character(num)
          CENTER[,(2*i-1)]=CENTER[,(2*i)]=center.use; 
          CENTERd[,(2*i-1)]=CENTERd[,(2*i)]=centerd.use; 
          R.array[(2*i-1)]=R.array[(2*i)]=r.use; 
          mu.array[,(2*i-1)]=mu.array[,(2*i)]=mu.use; 
          V.mat[,,(2*i-1)]=V.mat[,,(2*i)]=V.use; 
          error[(2*i-1)]=error[(2*i)]=error.use;
          INDX.num[c(2*i-1,2*i)]=num
          
          if(length(samp.y)>0)
          {
            YINDX.new[samp.y]=as.character(num)
            error.y[(2*i-1)]=error.y[(2*i)]=error.use.y
            error.y.pc[(2*i-1)]=error.y.pc[(2*i)]=error.use.y.pc
          }
          num=num+1
        }
        else
        {
          spcaz=SPCA(z,d); 
          INDX.new[samps[idx]]=as.character(num)
          CENTER[,(2*i-1)]=c=spcaz[[1]]; 
          CENTERd[,(2*i-1)]=cd=spcaz[[8]];
          R.array[(2*i-1)]=r=spcaz[[3]];
          mu.array[,(2*i-1)]=mu=spcaz[[7]]; 
          V.mat[,,(2*i-1)]=V=as.matrix(spcaz[[2]])
          error[(2*i-1)]=spcaz[[4]]; 
          INDX.num[(2*i-1)]=num;
          
          if((length(samp.y)>0)&(length(idx.y)>0))
          {
            YINDX.new[samp.y[idx.y]]=as.character(num)
            spcayd=SS_calc(zy,mu,c,V,r,d,cd)  
            error.y[(2*i-1)]=spcayd[[1]]
            error.y.pc[(2*i-1)]=spcayd[[2]]
            Y.hatD[samp.y[idx.y],]=spcayd[[3]]; Y.hatd[samp.y[idx.y],]=spcayd[[4]]; 
            Y.hatD.pc[samp.y[idx.y],]=spcayd[[5]]; Y.hatd.pc[samp.y[idx.y],]=spcayd[[6]]; 
          }
          if((length(samp.y)>0)&(length(idx.y)==0))
          {  count=count+1 }
          
          num=num+1
          
          spcaw=SPCA(w,d); 
          INDX.new[samps[idx2]]=as.character(num)
          CENTER[,(2*i)]=c=spcaw[[1]];
          CENTERd[,(2*i)]=cd=spcaw[[8]];
          R.array[(2*i)]=r=spcaw[[3]]; 
          mu.array[,(2*i)]=mu=spcaw[[7]]; 
          V.mat[,,(2*i)]=V=as.matrix(spcaw[[2]])
          INDX.num[(2*i)]=num
          error[(2*i)]=spcaw[[4]]; 
          
          if((length(samp.y)>0)&(length(idx2.y)>0))
          {
            YINDX.new[samp.y[idx2.y]]=as.character(num)
            spcayd=SS_calc(wy,mu,c,V,r,d,cd)  
            error.y[2*i]=spcayd[[1]]
            error.y.pc[2*i]=spcayd[[2]]
            Y.hatD[samp.y[idx2.y],]=spcayd[[3]]; Y.hatd[samp.y[idx2.y],]=spcayd[[4]]; 
            Y.hatD.pc[samp.y[idx2.y],]=spcayd[[5]]; Y.hatd.pc[samp.y[idx2.y],]=spcayd[[6]]; 
          }
          if((length(samp.y)>0)&(length(idx2.y)==0))
          {  count=count+1 }
          
          num=num+1
        }
      }
      else
      {
        INDX.new[samps]=as.character(num)
        CENTER[,(2*i-1)]=CENTER[,(2*i)]=center.use; 
        CENTERd[,(2*i-1)]=CENTERd[,(2*i)]=centerd.use; 
        R.array[(2*i-1)]=R.array[(2*i)]=r.use; 
        mu.array[,(2*i-1)]=mu.array[,(2*i)]=mu.use; 
        V.mat[,,(2*i-1)]=V.mat[,,(2*i)]=V.use; 
        error[(2*i-1)]=error[(2*i)]=error.use;
        INDX.num[c(2*i-1,2*i)]=num
        
        if(length(samp.y)>0)
        {
          YINDX.new[samp.y]=as.character(num)
          error.y[(2*i-1)]=error.y[(2*i)]=error.use.y
          error.y.pc[(2*i-1)]=error.y.pc[(2*i)]=error.use.y.pc
        }
        num=num+1
      }
    }
    drop.indx=which(duplicated(INDX.num)==TRUE)
    if(length(drop.indx)>0)
    {
      INDX.num=INDX.num[-drop.indx]
      R.array=R.array[-drop.indx]
      CENTER=CENTER[,-drop.indx]
      CENTERd=CENTERd[,-drop.indx]
      mu.array=mu.array[,-drop.indx]
      V.mat=V.mat[,,-drop.indx]
      error=error[-drop.indx]
      
      error.y=error.y[-drop.indx]
      error.y.pc=error.y.pc[-drop.indx]
    }
    num=(num-1)
    
    pc.error=mse[k]=log10(sqrt((sum(error))/n))
    n.part[k]=length(INDX.num)
    center.list[[k]]=CENTER
    centerd.list[[k]]=CENTERd
    mu.list[[k]]=mu.array
    r.list[[k]]=R.array
    V.list[[k]]=V.mat
    indx.list[[k]]=INDX.new
    num.indx.list[[k]]=INDX.num
    
    test.error.list[[k]]=error.y
    test.error.list.pc[[k]]=error.y.pc
    y.hatd.list[[k]]=Y.hatd
    y.hatD.list[[k]]=Y.hatD
    y.hatd.pc.list[[k]]=Y.hatd.pc
    y.hatD.pc.list[[k]]=Y.hatD.pc
    Yindx.list[[k]]=YINDX.new
    
    n.part1[k]=length(unique(YINDX.new))
    sp.err[k]=log10(sqrt(sum(error.y)/nt))
    sp.err1[k]=log10(sqrt(sum(error.y.pc)/nt))
    k=k+1
    INDX=INDX.new
    YINDX=YINDX.new
    max.cls=max(unname(table(INDX)))
    if(k>4)
    {
      if(n.part[k-1]==n.part[k-2])
      {
        print("Further partitions are not possible!")
        break
      }
      if((abs(mse[k-2]-mse[k-1])<epsilon2))
      {
        print("There is little decrease in errors over further partitions!")
        break
      }
    }
  }
  K=k-1 # keep the k
  
  SP1=sp.err; NP1=log10(n.part1)
  Result=cbind(NP1,SP1)
  colnames(Result)=c("log10(#partitions)","log10(MSE)")
  print(Result)
  
  return_res=list(y.hatD.list,Result)
}



