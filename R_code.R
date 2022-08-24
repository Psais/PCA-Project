library(jpeg)#We should have jpeg already installed
#------
loadImages = function() {
  name = c("male.pacole","male.pspliu","male.sjbeck","male.skumar",
           "male.rsanti","female.anpage","female.asamma","female.klclar",
           "female.ekavaz","female.drbost") 
  #------  
  x = matrix(0,nrow=10*5,ncol=200*180)
  #------
  k = 0
  
  for(i in 1:10) {
    for(j in 1:5) {
      k = k + 1
      #The first argument below should be the location of the folder.
      filename = paste("C:/Users/Manas/Downloads/grayfaces/",name[i],".",j,".jpg",sep="")
      x[k,] = as.vector(readJPEG(filename))
    }
  }
  #------
  return(x)
}
#---------------
#Next, we perform PCA. This is exactly what I mentioned above!

process = function(x) {
  meanx = apply(x,2,mean)
  
  y = scale(x,scale=F) #Rows of y are the cases
  A = y %*% t(y) 
  
  #To calculate eigenvalues of Y'Y, we calculate eigenvalues of YY'!
  #Because by a very standard result,
  # non zero eigenvalues of AB and BA are equivalent --- (Fact 1)
  
  
  eig = eigen(A)    
  
  P = t(y) %*% eig$vec[,-50] 
  
  #Columns of P are e vectors of Y'Y
  
  Q = apply(P,2,function(x) x/sqrt(sum(x*x))) 
  
  #Columns of Q form onb for rowspace of Y
  
  
  scores = y %*% Q
  #scores[i,] is for i-th image
  return(list(centre=meanx,onb=Q,scores=scores,values=eig$values))
}
#------------------
#Each principal component is again a 200*180 vector, and hence an image 
#It might be instructive to take a look at these.
#The following function helps you to just that.

showFace = function(newCoord, i) {
  plot(1:2,ty='n',main="0")
  y = abs(newCoord$onb[,i])
  extreme = range(y)
  y = (y-extreme[1])/(extreme[2]-extreme[1])
  dim(y) = c(200,180)
  rasterImage(as.raster(y),1,1,2,2)
}
#----------------

showSteps = function(newCoord,i) {
  meanx = newCoord$centre
  Q = newCoord$onb
  scores = newCoord$scores
  values = newCoord$values
  
  coeff = as.vector(scores[i,])
  
  plot(1:2,ty='n',main="0")
  y = meanx
  dim(y) = c(200,180)
  plot(as.raster(y))
  readline()
  for(k in 1:49) {
    if(k==1)
      temp = Q[,1]*coeff[1]
    else
      temp=Q[,1:k] %*% as.vector(coeff[1:k])
    recons = meanx + temp
    recons[recons<0]=0
    recons[recons>1]=1
    dim(recons) = c(200,180)
    plot(as.raster(recons))
    readline()
  }
}

filename = "C:/Users/Manas/Downloads/grayfaces/male.pacole.1.jpg"
tmp=readJPEG(filename)

rasterImage(tmp[90:130,70:115],-3,0,3,25)

x = loadImages()
newcoord = process(x)
showSteps(newcoord,26)

pc1 = newcoord$onb[,1]
showFace(newcoord,1)
expl = 100*cumsum(newcoord$values)/sum(newCoord$values) 
# This checks the explanation
plot(expl,type = 'l')
expl[10] #96.77246
print(expl[9]) #95.67942
print(expl[8]) #93.60762

#Hence in the spirit of alpha = 0.05, we choose 9 principal components!

princ_scores = newcoord$scores[,1:9] #9 principal components matrix

#--------
clust_size = 5
clust_num = 10
Within_Matrices = list() 
#List of all W_i matrices 
#(Which calculate variance within cluster of photos)
for (i in 1:clust_num) {
  j = clust_size*(i-1)
  Within_Matrices[[i]] = cov(princ_scores[(j+1):(j+clust_size),])  
}
W = (clust_size-1)*Reduce('+',Within_Matrices)/(clust_num*clust_size-1)

#W is our wanted within matrix!

Between_Vectors = list()
#List of all means of vectors of each cluster
for (i in 1:clust_num) {
  j = clust_size*(i-1)
  Between_Vectors[[i]] = apply(princ_scores[(j+1):(j+clust_size),],2,mean)
}
mean_mat = Between_Vectors[[1]] #Mean matrix
for (i in 1:clust_num) {
  mean_mat = rbind(mean_mat, Between_Vectors[[i]])
}
B = cov(mean_mat) 

#B is our between matrix!

W_inv = solve(W) 

#W inverse. Non singular since ONB with non zero eigenvalues
#We write next function to do PCA in essence with matrice W-1B!

cluster_ana = function(x) {
  A = W_inv %*% B
  eig = eigen(A)
  Q = apply(eig$vec,2,function(x) x/sqrt(sum(x*x))) 
  #Columns of Q form onb for rowspace of Y
  clust_scores = x %*% Q
  return(clust_scores)
  #clust_scores is the score matrix for cluster positions!
}

#This is the function with stores ranges!
rl = list() #List with range vectors as entries
for (i in 1:clust_num) {
  rl[[i]] = list(r1 = range(scores_clust[i,1]), r2 = range(scores_clust[i,2]))
}

scores_clust = cluster_ana(princ_scores)



#Final function to check if scores lie in a certain area!
cluster_check = function(x){
  scores = cluster_ana(x)
  cv = scores[1:2]
  for (i in 1:clust_num){
    if 
    ((rl[[i]]$r1[1]<=cv[1])&&
    (cv[1]<=rl[[i]]$r1[2])&&
    (rl[[i]]$r2[1]<=cv[2])&&
    (cv[2]<=rl[[i]]$r2[2])){
      print(i)}
    else {
      next}
  }
}

#Checks for first image! We are done!
cluster_check(princ_scores[1,])
