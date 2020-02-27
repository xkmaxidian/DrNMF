#****************************************************************
#*****                       DrNMF                  *************
#*****        for dynamic community detection       *************
#*****          during cancer progression           *************
#*****                   MA Xiaoke                  *************
#*****                   July-2018                  *************
#*****               Xidian University              *************
#*****             Email: xkma@xidian.edu.cn        *************
#****************************************************************
#library(NMI)
#library(NMF)



#******************* step1: parameter setting and data loading ***
 

lname = load("Data.RData"); print(lname);  print(dim(networks))
numcluster = c(); #number of evolving communities at each time
numT = dim(networks)[3]-1; # number of time steps
alpha = 1; #parameter for regularization
beta = 10;
numIter = 50
numnode = dim(networks)[1]
numcluster = c(214,198,215,228)
print(numcluster)

#*****************************************************************

#************ step 2 ********************************************



results = array(0,dim=c(numnode,numT))

  #the first time we use the NMF to obtain 
  #index = which(community[,1]>0)
  #tempNet = networks[index,index,1]
  #res = nmf(tempNet,numcluster[1],"lee",maxIter=75);
  
  #temresult = c()
  #Prevsolution = array(0,dim=c(numnode,numcluster[1]))
  # for(i in 1:dim(basis(res))[1])
   #{
   #temp1 = which.max(basis(res)[i,])
   # temresult= c(temresult,temp1[1])
   # Prevsolution[i,temp1[1]] = 1
   # }
  #results[index,1] = temresult; 
  #end of the first time
  #From the second time point
  source("./DrNMF.R");
  DrNMF(networks[,,1:numT],networks[,,numT+1],numcluster,numIter,alpha,beta)
#  for(i in 2:numT)
#  {
#   index = which(community[,i]>0)
#   regNet = Prevsolution%*%t(Prevsolution)
#   regNet[regNet>0] = 1;
#   diag(regNet) = 0;
#   regNet = regNet[index,index]
#   TempNet = networks[index,index,c((i-1):i)]
#   tempResult=RTLECDN(TempNet,regNet,numcluster[(i-1):i],numIter,gamma);
#   print("we are here")
#    Prevsolution = array(0,dim=c(numnode,numcluster[i]))
#    temresult = c()
#   for(t in 1:dim(tempResult)[1])
#   {
#   temp1 = which.max(tempResult[t,])
#    temresult= c(temresult,temp1[1])
#    Prevsolution[t,temp1[1]] = 1
#    }
#   results[index,i] = temresult; 
#   }





#**************** step3: accuracy of algorithm ******************
#using the normalized mutual information

#accuracy = c()

#for(i in 1:numT)
#{
# print(i)
# index = which(community[,i]>0)
# trueY = cbind(c(1:length(index)),community[index,i])
# PredY = cbind(c(1:length(index)),results[index,i])
#save(trueY,PredY,file="test.RData")
# tempa =NMI(trueY,PredY);
# accuracy = c(accuracy,tempa$value)
#print(tempa)

#}

#print(accuracy)
#save(results,file="expand.RData")

