DrNMF  <- function(networks,ppi,k,numIter,alpha,beta)
{
networks[is.na(networks)] = 0;
ppi[is.na(ppi)] = 0

numT = dim(networks)[3]
numobj = dim(networks)[1]

  
  degNet = diag(rowSums(ppi))
  Lppi = degNet-ppi


for(i in 1:numT)
{

  #initial solutions 
  networks[,,i] = (networks[,,i]+t(networks[,,i]))/2
  diag(networks[,,i])=0

  

  bmatrix = matrix(abs(rnorm(numobj*k[1],mean=0,sd=1)),byrow=T,nrow=numobj)
  fmatrix = matrix(abs(rnorm(numobj*k[1],mean=0,sd=1)),byrow=T,nrow=k[1])
 
  if(i==1)
   {
     #update rules:
    for(index in 1:numIter)
     {
    #print(i)
      #1 update basis matrix
     #print(dim(networks))
     #print(dim(fmatrix_1))
     #print(dim(fmatrix_2)) 
   Bnum = networks[,,i]%*%t(fmatrix)
   Bdenom = bmatrix%*%fmatrix%*%t(fmatrix)
   Bdenom = Bdenom+1e-32
   bmatrix = Bnum/Bdenom*bmatrix;
   for(t in 1:dim(bmatrix)[2])
       {
          tempvalue = sqrt(t(bmatrix[,t])%*%bmatrix[,t]);
           if(tempvalue>1e-16)
           {
           bmatrix[,t] = bmatrix[,t]/tempvalue;
            }
          }
   #2 update feature matrix
   fnum = t(bmatrix)%*%networks[,,i]
   fdenom = t(bmatrix)%*%bmatrix%*%fmatrix+beta*fmatrix%*%Lppi
   fdenom = fdenom+1e-32
   fmatrix = fnum/fdenom*fmatrix;
   for(t in 1:dim(fmatrix)[2])
       {
          tempvalue = sqrt(t(fmatrix[,t])%*%fmatrix[,t]);
           if(tempvalue>1e-16)
           {
           fmatrix[,t] = fmatrix[,t]/tempvalue;
            }
          }

        }
     


    }else
    {
    #other time points

     
  degTime = diag(rowSums(networks[,,i-1]))
  Lnet = degTime-networks[,,i-1]

       #update rules:
    for(index in 1:numIter)
     {
    #print(i)
      #1 update basis matrix
     #print(dim(networks))
     #print(dim(fmatrix_1))
     #print(dim(fmatrix_2)) 
   Bnum = networks[,,i]%*%t(fmatrix)
   Bdenom = bmatrix%*%fmatrix%*%t(fmatrix)+alpha*Lnet%*%bmatrix
   Bdenom = Bdenom+1e-32
   bmatrix = Bnum/Bdenom*bmatrix;
   for(t in 1:dim(bmatrix)[2])
       {
          tempvalue = sqrt(t(bmatrix[,t])%*%bmatrix[,t]);
           if(tempvalue>1e-16)
           {
           bmatrix[,t] = bmatrix[,t]/tempvalue;
            }
          }
   #2 update feature matrix
   fnum = t(bmatrix)%*%networks[,,i]
   fdenom = t(bmatrix)%*%bmatrix%*%fmatrix+alpha*fmatrix%*%Lppi
   fdenom = fdenom+1e-32
   fmatrix = fnum/fdenom*fmatrix;
   for(t in 1:dim(fmatrix)[2])
       {
          tempvalue = sqrt(t(fmatrix[,t])%*%fmatrix[,t]);
           if(tempvalue>1e-16)
           {
           fmatrix[,t] = fmatrix[,t]/tempvalue;
            }
          }

        }





    }


   filename = paste(i,"RData",sep=".")
   save(bmatrix,fmatrix,file=filename)

    } 



 

}
