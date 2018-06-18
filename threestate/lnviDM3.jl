
dat = readdlm("mds_G7-1950.csv",';')
cnames = dat[1,:]
dat    = dat[2:end,:]

nrow,ncol = size(dat)
pdat      = hcat([log.(dat[:,i]) - log.(dat[:,j]) for j=1:7, i=1:7 if i<j]...)
pnames    = [cnames[i] * " - " * cnames[j] for j=1:7, i=1:7 if i<j]

w = pdat[:,1]
b = [.4 ,.8,   1.1,    .3,    .3,.3,.3,.3,.3,.4,.1,.2,.3]
function lnviD3(b,w)
  
  # b      = [.4 ,.8,   1.1,    .3,    .3,.3,.3,.3,.3,.4,.1,.2,.3]
  # b = inits
  # w = randn(1000)
  datsize = length(w)

  d1hat, d2hat, d3hat  = b[1:3];
  p11hat, p22hat, p33hat = b[4:6];
  p12hat, p21hat, p31hat = b[7:9]; 
  sigma1hat, sigma2hat, sigma3hat =  b[10], b[10], b[10];
  mu1hat, mu2hat, mu3hat = b[11:13];
  
  rho1hat=0; theta1hat=0;
  # tranmax = Array(Float64,3,3);
  # [tranmax[i,i] = b[i+3]  for i = 1:3]
  tranmax = diagm(b[4:6])
  
  tranmax[[4 2 3]] = [p12hat,p21hat,p31hat]
  tranmax[[7 8 6]] = 1-sum(tranmax, 2)
  
  mu     = [mu1hat,mu2hat,mu3hat]
  sigma  = [sigma1hat,sigma2hat,sigma3hat]
  coefd1 = cumprod((d1hat+(0:(datsize-2)))./(1:(datsize-1)))
  coefd2 = cumprod((d2hat+(0:(datsize-2)))./(1:(datsize-1)))
  coefd3 = cumprod((d3hat+(0:(datsize-2)))./(1:(datsize-1)))
  coeff  = [coefd1 coefd2 coefd3]
  
  # muhat is to store the selected path of state: from date 2 */
  eta1  = zeros(9,datsize); # row 1:3 for path1, 4:6 for path2 */
  eta   = zeros(9,datsize); # same as above */
  xihat = zeros(9,datsize); # same as above */
  ehat  = zeros(3,datsize); # row 1, for path1, row 2 for path2 */
  
  # to keep the long memory part, not iid */
  zsigmahat = zeros(3,datsize) # row 1, for path1, row 2 for path2 */
  zhat      = zeros(3,datsize) # row 1, for path1, row 2 for path2 */
  lnlk      = zeros(3,datsize)
  
  path1 = zeros(Int8,3,datsize)
  path2 = zeros(Int8,3,datsize)
  path3 = zeros(Int8,3,datsize)
  path1[:,1] = [1,0,0]
  path2[:,1] = [0,1,0]
  path3[:,1] = [0,0,1]
  
  zsigmahat[:,1] = w[1] - mu 
  # zsigmahat is de-drifted series for each given path; or it is z*sigmahat
  
  zhat[:,1] = zsigmahat[:,1] ./ sigma;
  
  ehat[:,1] = zhat[:,1];
  
  limprobs  = (tranmax^7)[1,:]
  lnlk[:,1] = log.(exp.(-ehat[:,1].^2./2 + 1e-50 ) ./ (sqrt(2*pi)).*limprobs);
  
  j = 2;
  
  
  ####################!!!!!!!!!!!!!!!!!!!!!
  while j <= datsize
    # tempsMUZ and temps are needed for likelihood possibility
    
    tempsMUZ = vec(transpose([(zsigmahat[1,1:(j-1)])[end:-1:1] (zsigmahat[2,1:(j-1)])[end:-1:1] (zsigmahat[3,1:(j-1)])[end:-1:1]]) * coeff[1:(j-1),:])
    
    temps    = repeat(mu; inner=[3]) + rho1hat * repeat(sigma;inner=[3]) .* repeat(zhat[:,j-1];outer=[3]) +
      theta1hat*repeat(sigma,inner=[3]) .* repeat(ehat[:,j-1], outer=[3])
    
    # likelihood possibility
    eta1[:,j]  = exp.(-(w[j] - temps - tempsMUZ).^2 ./ (2*repeat(sigma;inner=[3]).^2+1e-50)) ./
      (sqrt.(2*pi*repeat(sigma,inner=[3]).^2))
    
    eta[:,j] = eta1[:,j] .* vec(tranmax)
    
    # Case is 1 for 1->1; 2 for 2->1; 3 for 1->2; 4 for 2->2
    case  = reshape([log(eta[x,j]+1e-50) + lnlk[(x-1) % 3 + 1,j-1] for x in 1:9],3,3)
    case  = findmax(case,1)[2]
    case2 = (case-1) % 3 + 1
    
    # holding previous paths in temp memory
    oldlnlk = lnlk
    immutable oldPath2
      path
    end
    
    oldzsigmahat = zsigmahat
    oldzhat = zhat
    oldehat = ehat
    
    # recalculating values
    # path ends up with state 1:
    lnlk[1,j]    = log(eta[case[1],j]+1e-50) + oldlnlk[case2[1],j-1]
    path1[:,1:j] = [oldPath[((case2[1]-1)*3+1):(case2[1]*3),1:(j-1)] [1,0,0]]
    zsigmahat[1,1:j] = append!(oldzsigmahat[case2[1],1:(j-1)], w[j] - mu[1] - tempsMUZ[case[1]])
    zhat[1,1:j]  = append!(oldzhat[case2[1],1:(j-1)],zsigmahat[1,j]/(sigma1hat));
    ehat[1,1:j]  = append!(oldehat[case2[1],1:(j-1)],(w[j] - tempsMUZ[case[1]] -temps[case[1]])/sigma1hat);
    
    # path ends up with state 2:
    lnlk[2,j]    = log(eta[case[2],j]+1e-50) + oldlnlk[case2[2],j-1]
    path2[:,1:j] = [oldPath[((case2[2]-1)*3+1):(case2[2]*3),1:(j-1)] [0,1,0]];
    zsigmahat[2,1:j] = append!(oldzsigmahat[case2[2],1:(j-1)], w[j] - mu[2] - tempsMUZ[case[2]])
    zhat[2,1:j]  = append!(oldzhat[case2[2],1:(j-1)],zsigmahat[2,j]/(sigma2hat));
    ehat[2,1:j]  = append!(oldehat[case2[2],1:(j-1)],(w[j] - tempsMUZ[case[2]] -temps[case[2]])/sigma2hat)
    
    # path ends up with state 3:
    lnlk[3,j]    = log(eta[case[3],j]+1e-50) + oldlnlk[case2[3],j-1]
    path3[:,1:j] = [oldPath[((case2[3]-1)*3+1):(case2[3]*3),1:(j-1)] [0,0,1]]
    zsigmahat[3,1:j] = append!(oldzsigmahat[case2[3],1:(j-1)], w[j] - mu[3] - tempsMUZ[case[3]])
    zhat[3,1:j]  = append!(oldzhat[case2[3],1:(j-1)],zsigmahat[3,j]/(sigma3hat));
    ehat[3,1:j]  = append!(oldehat[case2[3],1:(j-1)],(w[j] - tempsMUZ[case[3]] -temps[case[3]])/sigma3hat)
    
    # for(x in 0:2){
    #   lnlk[x+1,j]   = log(eta[case[x+1],j]+1e-50) + oldlnlk[case2[x+1],j-1]
    #   paths[(3*x+1):(3*x+3),1:j] = cbind(oldPath[(3*x+1):(3*x+3),1:(j-1)],(1:3 == (x+1)) * 1)
    #   zsigmahat[x+1,1:j] = c(oldzsigmahat[case2[x+1],1:(j-1)],  w[j] - mu[case2[x+1]] - tempsMUZ[case[x+1]])
    #   zhat[x+1,1:j] = c(oldzhat[case2[x+1],1:(j-1)], zsigmahat[x+1,j]/(sigma[x+1]))
    #   ehat[x+1,1:j] = c(oldehat[case2[x+1],1:(j-1)], (w[j] - tempsMUZ[case[x+1]] - temps[case[x+1]])/sigma[x+1]) # ONEMLI!!!
    # }
    
    j=j+1

  end
  #   lnlipath1 = t(lnlk1);
  #   lnlipath2 = t(lnlk2);
  lnlkRes = findmax(lnlk[:,datsize])[1][1]
  # if(lnlk[1,datsize] > lnlk[2,datsize]){
  #   lnlkRes = lnlk[1,datsize]} else {lnlkRes = lnlk[2,datsize]}
  # 
  return lnlkRes
end

dat = readdlm("mds_G7-1950.csv",';')
cnames = dat[1,:]
dat    = dat[2:end,:]

nrow,ncol = size(dat)
pdat      = hcat([log.(dat[:,i]) - log.(dat[:,j]) for j=1:7, i=1:7 if i<j]...)
pnames    = [cnames[i] * " - " * cnames[j] for j=1:7, i=1:7 if i<j]

source("~/Documents/threestate/arfimaSimH.R")

genseries2 = function(n){
  n1 = ceil(n/3)
  n2 = ceil(n/3)
  n3 = n - n1 - n2
  ser1 = arfima.simH(n1, model = list(dfrac = .9))
  ser2 = arfima.simH(n2, model = list(dfrac = .6)) + tail(ser1,1)
  ser3 = arfima.simH(n3, model = list(dfrac = .2)) + tail(ser2,1)
  c(ser1,ser2,ser3)
}


# inits      = c(.4 ,.8,1.1,.8,.8,.8,.1,.1,.1,.01,.1)
# inits      = c(runif(3,-2,2),runif(3,.3,.6))
# inits      = c(inits,.8-inits[4:6],.01,.1)                
# lnviD3(inits, ww)
# 
# ww = log(read.csv("../Downloads/g7_1950.csv",sep=";")[:,-1])
# ww = ww[:,1] - ww[:,2]
# lowerV = c(-2,-2,-2,.7,.7,.7,0,0,0,0.001,-5,-5,-5)
# upperV = c(2,2,2,.999,.999,.999,.3,.3,.3,5,5,5,5)
# inits  = c(.4 ,   .8,   1.1,    .3,    .3,.3,.3,.3,.3,.01,.1)
# optim(inits, function(p) -lnviD3(p,ww),method="L-BFGS-B",control = list(trace=3))
# 
# aa =  c(2,-2,2,0.201,0.201,0.201,0.999,0.999,0.201,0.001,-5 )
# lnviD3(aa, ww)

ser = genseries2(200)
plot(ser)

lowerV = c(-2,-2,-2, .6 , .6 , .6 , 0, 0, 0,0.001,-5,-5,-5)
inits  = c( 1,.5, 0, .61, .61, .61,.1,.1,.1,   .4,.1,.1,.1) 
upperV = c( 2, 2, 2,.999,.999,.999,.3,.3,.3,    5, 5, 5, 5)

const_mat = matrix(0,13,13)
diag(const_mat) = 1
const_mat = rbind(const_mat,-const_mat)
dmmvec = -c(0,0,0,1,0,0,1,0,0,0,0,0,0)
dmmvec = c(dmmvec,0,dmmvec,0,dmmvec)[1:39]
const_mat = rbind(const_mat,matrix(dmmvec,3,,T))
const_mat = cbind(const_mat,c(lowerV,-upperV,rep(-1,3)))

res = constrOptim(inits, function(p) -lnviD3(p,ser), NULL, ui = const_mat[:,-14], const_mat[:,14])
thepath = return_path(res$par,ser)
plot(apply(res$par[1:3]*thepath,2,sum),type="l")
