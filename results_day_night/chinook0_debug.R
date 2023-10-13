# Goal
# run marss model with just temp diff and flow and check aicc


num_years = 15
num_rows = num_years*2
list_combinations <- get_covariate_combinations(1:5)
out.tab.hatchery<- NULL
fits.hatchery <- list()
for(i in 1:length(list_combinations)){
  
  covariate_number <- length(list_combinations[[i]])
  covariates <- list_combinations[[i]]
  print(covariates)
  c = NULL
  name = NULL
  if(covariates[1] == 1){
    
    for(kk in c(1,3,4,6,9)){
      c = NULL
      name = NULL
      for(j in covariates){
        if(j == 1){
          k = kk
        }
        else if(j==2){
          k = 7
        }
        else if(j==3){
          k = 2
        }
        else if(j==4){
          k = 5
        }
        else if(j==5){
          k = 10
        }
        
        c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
        name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
        name_individual = substr(name_long,1,nchar(name_long)-9)
        name = paste(name, substr(name_long,1,nchar(name_long)-9))
        
      }
      # print(c)
      print(name)
      c_num <- length(covariates)
      if(k==10){
        print("has hatchery")
        has_hatchery = TRUE
        c_num <- length(covariates)
        if(c_num == 1){
          
          fit.model = c(list(c= c), mod.list_1_0_h)
        }
        
        else if(c_num == 2){
          
          fit.model = c(list(c= c), mod.list_2_0_h)
        }
        else if(c_num == 3){
          
          fit.model = c(list(c= c), mod.list_3_0_h)
        }
        else if(c_num == 4){
          
          fit.model = c(list(c= c), mod.list_4_0_h)
        }
        else if(c_num == 5){
          
          fit.model = c(list(c= c), mod.list_5_0_h)
        }
      }
      else{
        has_hatchery = FALSE
        print("no hatchery")
        c_num <- length(covariates)
        if(c_num == 1){
          
          fit.model = c(list(c= c), mod.list_1_0)
        }
        
        else if(c_num == 2){
          
          fit.model = c(list(c= c), mod.list_2_0)
        }
        else if(c_num == 3){
          
          fit.model = c(list(c= c), mod.list_3_0)
        }
        else if(c_num == 4){
          
          fit.model = c(list(c= c), mod.list_4_0)
        }
        else if(c_num == 5){
          
          fit.model = c(list(c= c), mod.list_5_0)
        }
      }
      
      # fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
      #              control=list(maxit=2000))
      # 
      # 
      # out=data.frame(c=name, d = "None",
      #                logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
      #                num.iter=fit$numIter, converged=!fit$convergence,
      #                stringsAsFactors = FALSE)
      # out.tab.hatchery=rbind(out.tab.hatchery,out)
      # fits.hatchery=c(fits.hatchery,list(fit))
    }
    
  }
  else{
    for(j in covariates){
      if(j == 1){
        k = 1
      }
      else if(j==2){
        k = 7
      }
      else if(j==3){
        k = 2
      }
      else if(j==4){
        k = 5
      }
      else if(j==5){
        k = 10
      }
      
      c = rbind(c,covariates_chinook0[((1+(k-1)*num_rows):(k*num_rows)),])
      name_long = rownames(covariates_chinook0)[1+(k-1)*num_rows]
      num_individual = substr(name_long,1,nchar(name_long)-9)
      name = paste(name, substr(name_long,1,nchar(name_long)-9))
      
    }
    # print(c)
    print(name)
    c_num <- length(covariates)
    if(k==10){
      print("has hatchery")
      has_hatchery = TRUE
      c_num <- length(covariates)
      if(c_num == 1){
        
        fit.model = c(list(c= c), mod.list_1_0_h)
      }
      
      else if(c_num == 2){
        
        fit.model = c(list(c= c), mod.list_2_0_h)
      }
      else if(c_num == 3){
        
        fit.model = c(list(c= c), mod.list_3_0_h)
      }
      else if(c_num == 4){
        
        fit.model = c(list(c= c), mod.list_4_0_h)
      }
      else if(c_num == 5){
        
        fit.model = c(list(c= c), mod.list_5_0_h)
      }
    }
    else{
      has_hatchery = FALSE
      print("no hatchery")
      c_num <- length(covariates)
      if(c_num == 1){
        
        fit.model = c(list(c= c), mod.list_1_0)
      }
      
      else if(c_num == 2){
        
        fit.model = c(list(c= c), mod.list_2_0)
      }
      else if(c_num == 3){
        
        fit.model = c(list(c= c), mod.list_3_0)
      }
      else if(c_num == 4){
        
        fit.model = c(list(c= c), mod.list_4_0)
      }
      else if(c_num == 5){
        
        fit.model = c(list(c= c), mod.list_5_0)
      }
    }
    
    
    # fit <- MARSS(subset_chinook_summer_perhour, model=fit.model, silent = TRUE, method = "BFGS",
    #              control=list(maxit=2000))
    # 
    # 
    # out=data.frame(c=name, d = "None",
    #                logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
    #                num.iter=fit$numIter, converged=!fit$convergence,
    #                stringsAsFactors = FALSE)
    # out.tab.hatchery=rbind(out.tab.hatchery,out)
    # fits.hatchery=c(fits.hatchery,list(fit))
    # 
  }
  
  
}

# the problem is that there is a space in from of name_individual that was not detecting
#the presence of hatchery covariate correctly
# if I change the condition to be based on k insteado f name_individual, it should work