arfima.simH<-function(n, model, seed=NULL,...){
  require("arfima")
  seed<-ifelse(is.null(seed),as.double( Sys.time())/runif(1,1,56984),seed)  
  if( model$dfrac <.5 && model$dfrac> -1){set.seed(seed);
                                          return(arfima::arfima.sim(n=n, model=model,...))
  } else {
    if(model$dfrac<0) {model$dfrac<-model$dfrac+1;
                       return(diff((arfima.simH(n=n+1,  model=model,seed=seed,... ) ) ) )}
    else{ model$dfrac<-model$dfrac-1;
          return(cumsum((arfima.simH(n=n,  model=model,seed=seed,... ) ) ) )
    }
  }
}

