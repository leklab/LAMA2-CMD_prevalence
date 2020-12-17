#!/usr/bin/Rscript

library(dplyr)


# get the posterior params 
get_posterior_w_v <- function(x,n,w,v){
  a <- x+v
  b <- 2*n-x+w
  return(list(a,b))
}

# estimate the mean of AF posterior distribution
pre_w_v_plugged <- function(v,w,n,x){
 # tmp <- ((x+v-1)^2)/(2*n+v+w-1)^2
  tmp.1 <- (x+v)/(2*n+v+w)
  
  return(tmp.1)
}

#######################################################################################
############# use normal distribution to approximate the sum of beta r.v.s ############
#######################################################################################

#### get the mean of the approximated normal distribution ####
get_mean_normal <- function(vs,ws){
  tmp <- vs/(vs+ws)
  return(sum(tmp))
}

#### get the variance of the approximated normal distribution ####
get_var_normal <- function(vs,ws){
  vars <- vs*ws/((vs+ws+1)*(vs+ws)^2)
  return(sum(vars))
}

get_prevelance_estimate <- function(pop){

  pop_ac_col = paste(pop,"AC",sep="_")
  pop_an_col = paste(pop,"AN",sep="_")


  #refiltered_data = filtered_data
  refiltered_data = filtered_data %>% filter(.[[pop_an_col]] > 10000)

  AC_interested = pull(refiltered_data,pop_ac_col)
  AN_interested = pull(refiltered_data,pop_an_col)


  types = c("frameshift_variant","splice_acceptor_variant","splice_donor_variant","missense_variant","stop_gained","exon_variant","UTR_variant","other_variant")

  af_changed <- rep(0,nrow(refiltered_data))

  posterior_param <- data.frame(v=rep(0,nrow(refiltered_data)),w=rep(0,nrow(refiltered_data)))
  inds_updated <- c()

  for(i in 1:length(types)){

    ind.tmp <- grep(types[i],refiltered_data$Annotation)
    inds_updated <- c(inds_updated,ind.tmp)

    v = params$v[which(params$type==types[i])]
    w = params$w[which(params$type==types[i])]
    ac = AC_interested[ind.tmp]
    an = AN_interested[ind.tmp]

    post_params <- get_posterior_w_v(ac,an/2,w,v)

    posterior_param[ind.tmp,1]<- post_params[[1]]
    posterior_param[ind.tmp,2]<- post_params[[2]]

    af_changed_tmp <- pre_w_v_plugged(v,w,an/2,ac)
    af_changed[ind.tmp] <- af_changed_tmp

  }


  mu <- get_mean_normal(posterior_param$v[which(posterior_param$w!=0)],posterior_param$w[which(posterior_param$w!=0)])
  sigma.2 <- get_var_normal(posterior_param$v[which(posterior_param$w!=0)],posterior_param$w[which(posterior_param$w!=0)])

  # get the estimates
  E.s.2 = mu^2 + sigma.2
  prev.estimated =  E.s.2


  # get the confidence interval
  confidence = 0.95
  alpha = 1 - as.numeric(confidence)
  # lower bound
  lb <- qchisq(alpha/2,1,ncp = mu^2/sigma.2)*sigma.2
  # upper bound
  up <- qchisq(1-alpha/2,1,ncp = mu^2/sigma.2)*sigma.2


  cat(pop,"estimated prevalence (per million): ",prev.estimated*1e6,"\n")
  cat(pop,"confidence interval with ",as.numeric(confidence)*100,"% confidence: ",lb*1e6,"-",up*1e6,"\n")

}


data = read.csv("LAMA2_known_pathogenic_novel_lof_gnomad_af.tsv",header=T,sep="\t")
params = read.csv("data/beta_parameter_prior_ExAC.txt",header=T,sep="\t")

#filtered_data = subset(data, ALL_AC < 30 & Annotation != 'start_lost' & gnomAD_Source != 'Genome_r3')
filtered_data = subset(data, ALL_AC < 30 & Annotation != 'start_lost')
#filtered_data = subset(data, ALL_AC < 30 & INFO != 'Novel_gnomAD_LoF' & Annotation != 'start_lost')

get_prevelance_estimate("ALL")
get_prevelance_estimate("NFE")
get_prevelance_estimate("AFR")
get_prevelance_estimate("AMR")
get_prevelance_estimate("EAS")











