##################################################################################################################
################################## Efficient Manifold Learning Using Spherelets ##################################
##################################################################################################################

This project provides codes for manifold learning using sphrelets (see https://arxiv.org/abs/1706.08263 )

The R-file SPCA_functions.R contains all the functions required for manifold estimation using Spherelets.
The three main functions are:

(1) SPCA: Given the data and the intrinsic dimension d, this function estimates the center, radius of the best 
sphere of dimension d using SPCA and calculates the error in estimation. The inputs to this function are:

(a) X: the data matrix in no. of data-points times dimension format
(b) d: The intrinsic dimension

The output to this function is a list which including:

(a) c: The center of the fitted sphere (dimension D) 
(b) V: (d+1)-dimensional PCA rotation matrix
(c) r: Radius of the fitted sphere
(d) SS: Sum of squared error of Spherical-PCA projection
(e) SS_new: Sum of squared error of PCA projection
(f) c.d: The center of the fitted sphere in the d+1-dimensional affine space 


(2) SS_calc: Given the center c, radius r, V, intrinsic dimension d and a test data matrix X, this function finds 
the SPCA-projection Proj(X) of X to the sphere determined by (V,c,r), and calculates the sum of squared difference
between X and Proj(X). Inputs to this function are:

(a) X: the test data matrix in nXD format
(b) mu: component-wise mean of the training data
(c) c: center of the best fitting d-dimensional sphere obtained from the training data
(d) r: radius of the best fitting d-dimensional sphere obtained from the training data
(e) V: (d+1)-dimensional PCA rotation matrix obtained from the training data
(f) d: intrinsic dimension
(g) c.d: The center of the best fitting sphere in the d+1-dimensional affine space obtained from the training data

The output of this function includes:
(a) SS: Sum of squared difference between X and Proj(X) when Spherical-PCA projection is applied
(b) SS_new: Sum of squared difference between X and Proj(X) when d-dimensional PCA projection is applied
(c) Y.hatD: Estimated values of the test points


(3) SPCA_Manifold_approx: Given a data-matrix, intrinsic dimension d, and an error threshold epsilon this function
recursively partitions the data using PCA and fits SPCA in each partition until the average error in estimation is 
less than the threshold value. Given a test data-matrix y, this projection operator is applied on y without 
retraining to obtain predictive values of y, and mean squared error is calculated. The inputs to this function 
are:
(a) x: training data matrix in nXD format
(b) y: test data matrix in nXD format
(c) d: intrinsic dimension
(e) epsilon: log10(error) threshold, e.g., epsilon=-5
(f) N.min: minimum on of samples required to partition a subspace further, e.g., N.min=10
(g) epsilon2: To avoid over-parametrization, the algorithm stops if there is negligible decrease in error for 
	      further partitioning. The epsilon2 threshold is set to control this. If the decrease in error for 
	      further partitioning is less than epsilon2, then the algorithm stops.

The output of the function includes:
(a) y.hatD: A list of length equal to the level of partitions. The i-th component of the list contains estimated 
	    values of y at the i-th level of partitioning.
(b) Result: A table providing the number of partitions and mean squared error in log10 scale.

##################################################################################################################
We demonstrate the use of the functions using the following example:  

## Generating data:

EulerSpiral=function(n)
{
  t=seq(0,3,by=(2/n))
  len=length(t)
  X=matrix(data=0,len,2)
  sint2=function(u){ sin(u^2) }
  cost2=function(u){ cos(u^2) }
  for(i in 2:len)
  {
    X[i,1]=X[(i-1),1]+(integrate(sint2,t[i-1],t[i])$value)
    X[i,2]=X[(i-1),2]+(integrate(cost2,t[i-1],t[i])$value)
  }
  mu=c(0,0);
  Sigma=0.0001*diag(2);
  noise=mvrnorm(nrow(X),mu,Sigma);
  Y=X+noise;
  return(list(X,Y))
}

X=EulerSpiral(1000)
x=X[[1]]; y=X[[2]]

## Manifold approximation using spherelets 
source('SPCA_functions.R')
d=1
N.min=10
epsilon=-5
epsilon2=0.0001

SPCA_part=SPCA_Manifold_approx(x,y,d,epsilon,epsilon2,N.min)
Y_hat=SPCA_part[[1]]

par(mfrow=c(2,3))
for(i in 1:(length(Y_hat)-1))
{
plot(Y_hat[[i+1]],xlab=" ", ylab=" ")
}
