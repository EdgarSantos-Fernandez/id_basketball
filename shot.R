#---------------------------------------------------------------------
# Edgar Santos Fernandez
# 2019-07-24
# Description: Takes the game data, computes the ID and plots the clusters 
#---------------------------------------------------------------------

.libPaths('C:\\1\\R')
#install.packages('Rcpp', dependencies = T)
library(Rcpp)
library(ggplot2)
library(mvtnorm)
library(plyr)
library(dplyr)
library(rgl)
library(reshape2)
library(plot3D)
library(Matrix)
library(knitr)
library(superheat)
library(mcclust)
library(RCurl)
library(jsonlite)
library(sp)
library(zoo)
library(lubridate)
library(stringr)
library(gridExtra)
library(grid)
library(xtable)

Rcpp::sourceCpp("BNPHidalgoSB___plain_and_rep.cpp")
Rcpp::sourceCpp("posterior_analysis_of_likelihood_and_Nq.cpp")


#-----------------------------------------------------
#functions:

subset_data_both <- function(m4 = m4, team = team){ #subset the data
  possession <- ifelse(team == '_h', 'home', 'away')
  vars <- c(names(m4)[grepl('_h', names(m4))] , names(m4)[grepl('_a', names(m4))] , 'x_loc_mean','z_loc_ball', 'desc', 'event.id' , 'ball_possession')  
  DATA <- m4[, vars] 
  DATA$repeated <- ifelse(DATA$x_loc_mean - lag(DATA$x_loc_mean) == 0,1,0) #repeated row with different actions
  DATA <- DATA %>% group_by( x_loc_mean ) %>%  slice(n()) %>% ungroup() %>% arrange(event.id)
  DATA <- dplyr::filter(DATA, ball_possession == possession)
  
  ball_possession <- DATA$ball_possession
  categ <- DATA$desc # event_team outcome
  event_id <- DATA$event.id
  
  ball_names <- c(names(DATA)[grepl('_ball', names(DATA))]) # removing the coordinates of the ball
  home <- DATA %>% dplyr::select(-repeated, -x_loc_mean, -ball_names)
  
  DATA <- DATA %>% dplyr::select(-desc, -repeated, -x_loc_mean, -event.id, -ball_possession, -ball_names)
  list(DATA, ball_possession, categ, event_id)
}


# Just a wrapper around to F.Denti ID functions  
# returns a list with the output of HIDALGO, q, mu and the processing time
int_dim <- function( DATA = DATA, 
                     categ = categ,
                     q = 3,
                     nsim = 40000, # Iterations
                     BI = 20000,   # Burn In
                     L = 50,      # Truncation Threshold for the Stick-Breaking process
                     CSI = .85){
  
  # Computing mu = 2nd NN / 1st NN
  distDD <- dist(DATA)
  distDD <- as.matrix(distDD)

  mu <- numeric(nrow(distDD))
  for(i in 1:nrow(distDD)){
    sD <- sort(distDD[i,])
    mu[i] <-sD[3]/sD[2]
  }

  ###########  compute Nq - functions ############################
  f_which.min <- function(vec, idx) sort(vec, index.return = TRUE)$ix[idx]
  Nq_x2 <- function(q,dist_mat){
    n <- nrow(dist_mat)
    # negli q non รจ contato i==j!
    # sort(c,index.return=T) !!!!!
    Nq <- t(apply(
      dist_mat, 1, function(x) { 
        index <- f_which.min(x,idx =2:(q+1) ) 
        l <- numeric(n)
        l[index] <- 1
        l
      }
    ) )
    return(Nq)
  }
  
  # actually compute Nq, q=3
  #q=q #3#   
  Nq=Nq_x2(q, as.matrix(distDD))
  N=length(mu)
  
  ########### RUNNING THE MCMC #########################################################################################
  ########
  # the prior on d: high variance
  #curve(dgamma(x,1.25,.25),0,10) # NB: selection of the prior? Impact of posterior
  #alpha=.2
  #sum(alpha/(alpha+(1:N)-1))
  #########
  start <- Sys.time()
  out <- BNP_HIDALGO_SB_Rcpp_repulsive_faster(mu_obser = mu,
                                              L = L, 
                                              NSIM = nsim, 
                                              burn_in = BI,
                                              thinning = 1,
                                              verbose = T,
                                              verbose_step = 500,
                                              # hyperparameters
                                              a_alpha = 1, b_alpha = 1,
                                              # hyperparameters for Gamma on d
                                              a0_d = 1.25, b0_d = .25, 
                                              # fixed alpha for the SB process?
                                              fixedAB = F, fixedALPHA = 1,
                                              Nq = Nq,q = q,CSI = CSI,
                                              # use a repulsive prior? (I would run the algorithm without it first, and then re-run it with the rep prior id the d overlaps to great extent)
                                              REPULSIVE = F,  
                                              t1 = .5,  t2 = .001)
  
  end <- Sys.time()
  time <- end - start
  list(out = out, q = q, mu = mu, time = time)
}



# reading the game data
m4 <- readRDS('game.12.25.2015.CLE.at.GSW.RDS')
#m4 <- dplyr::filter(df_games, game == games[i], desc %in% c("ShotMade", "ShotMissed"))

DATA_home <- subset_data_both(m4 = m4, team = '_h')[[1]] #home team attacking 
categ_home <- subset_data_both(m4 = m4, team = '_h')[[3]]

DATA_away <- subset_data_both(m4 = m4, team = '_a')[[1]] #away team attacking 
categ_away <- subset_data_both(m4 = m4, team = '_a')[[3]]

nsim <- 40000

id_home0 <- int_dim(DATA = DATA_home, categ = categ_home, q = 3, nsim = nsim)
glimpse(id_home0)

id_away0 <- int_dim(DATA = DATA_away, categ = categ_away, q = 3, nsim = nsim)
glimpse(id_home0)

out_home <- id_home0[[1]] 
out_away <- id_away0[[1]]


# Home team. following the latent labels Z, we rebuild the chains (to avoid label switching)
N = length(id_home0[[3]]) #length of mu
nsim=40000
tracked_d_home=Tracking(out_home$AC, out_home$AD,N,nsim)

# Away team.following the latent labels Z, we rebuild the chains (to avoid label switching)
N=length(id_away0[[3]]) #length of mu
tracked_d_away=Tracking(out_away$AC, out_away$AD,N,nsim)


# posterior mean of the ID
df_home <- data.frame(team = 'home',
                      mean=apply(tracked_d_home,1,mean),
                      sd=apply(tracked_d_home,1,sd),
                      ball_possession = subset_data_both(m4 = m4, team = '_h')[[2]],
                      event_id = subset_data_both(m4 = m4, team = '_h')[[4]],
                      categ = subset_data_both(m4 = m4, team = '_h')[[3]])

df_away <- data.frame(team = 'away',
                      mean=apply(tracked_d_away,1,mean),
                      sd=apply(tracked_d_away,1,sd),
                      ball_possession =subset_data_both(m4 = m4, team = '_a')[[2]],
                      event_id = subset_data_both(m4 = m4, team = '_a')[[4]],
                      categ = subset_data_both(m4 = m4, team = '_a')[[3]])

#---------------------------------------------------------------------

# plotting
M_home <- mcclust::comp.psm(t(out_home$AC)) ## computes the pairwise coclustering matrix 
CL_home <- mcclust::minbinder(cls = t(out_home$AC), psm = M_home) # Compute Posterior Expectation of Binders Loss Function

labels <- subset_data_both(m4 = m4, team = '_h')[[4]] # event.id
#ball_possession <- subset_data_both(m4 = m4, team = '_h')[[2]] # ball_possesion

clus <- data.frame(cluster = apply(t(out_home$AC),2,median),  # most likely cluster group
                   event.id = labels)

m4$desc <- droplevels(m4$desc) 

M_home <- data.frame(M_home)
names(M_home) <- labels
row.names(M_home) <- labels

ball_possession <- ifelse(subset_data_both(m4 = m4, team = '_h')[[2]]=='home',1,0)
yrcol <- ifelse(ball_possession==1,"#0072B2", "#F0E442" )

desc <- ifelse(subset_data_both(m4 = m4, team = '_h')[[3]]=='ShotMade',1,0)
ytcol <- ifelse(desc== 1, "#009E73", "#D55E00")

if(nrow(M_home) > 100){M_home <- M_home[1:100, 1:100]
ball_possession <- ball_possession[1:100]
yrcol <- yrcol[1:100]
desc <- desc[1:100]
ytcol <- ytcol[1:100]
}

set.seed(2019)
x11(); gsw_attack <- superheat(M_home,
                        yr= desc, yr.obs.col = ytcol, yr.axis.name = 'missed or scored', yr.point.size = 1,
                        pretty.order.rows = T,
                        pretty.order.cols = T,
                        # generate column clusters
                        n.clusters.rows = 4,
                        n.clusters.cols = 4,
                        left.label = 'variable',
                        print.plot = T,
                        bottom.label.text.size = 2.5,
                        left.label.text.size = 1.5,
                        clustering.method = 'hierarchical') 


# away team

#out_away$AC[out_away$AC > nrow(out_away)] <- nrow(out_away) ## NB: All elements of cls must be integers in 1:nobs. had to set cluster to max ncol
M_away <- mcclust::comp.psm(t(out_away$AC)) ## computes the pairwise coclustering matrix 
CL_away <- mcclust::minbinder(cls = t(out_away$AC),psm = M_away) 

labels <- subset_data_both(m4 = m4, team = '_a')[[4]] # event.id
clus2 <- data.frame( cluster2 = apply(t(out_away$AC),2,median),  # most likely cluster group
                     event.id = labels)

M_away <- data.frame(M_away)
names(M_away) <- labels
row.names(M_away) <- labels

ball_possession2 <- ifelse(subset_data_both(m4 = m4, team = '_a')[[2]]=='home',1,0)
yrcol2 <- ifelse(ball_possession2==1,"#0072B2", "#F0E442" )

desc2 <- ifelse(subset_data_both(m4 = m4, team = '_a')[[3]]=='ShotMade',1,0)
ytcol2 <- ifelse(desc2== 1, "#009E73", "#D55E00")

if(nrow(M_away) > 100){M_away <- M_away[1:100, 1:100]
ball_possession2 <- ball_possession2[1:100]
yrcol2 <- yrcol2[1:100]
desc2 <- desc2[1:100]
ytcol2 <- ytcol2[1:100]
}

set.seed(2019)
x11(); cle_attack <-superheat(M_away,
                        yr= desc2, yr.obs.col = ytcol2, yr.axis.name = 'missed or scored', yr.point.size = 1,
                        pretty.order.rows = T,
                        pretty.order.cols = T,
                        # generate column clusters
                        n.clusters.rows = 4,
                        n.clusters.cols = 4,
                        left.label = 'variable',
                        print.plot = T,
                        bottom.label.text.size = 2.5,
                        left.label.text.size = 1.5,
                        clustering.method = 'hierarchical')



###############################################

# cluster GSW (home team)

cluster <- data.frame(clus = gsw_attack$membership.cols) #cluster membership
cluster$id <- as.numeric(as.character(rownames(cluster)))
cluster <- cluster %>% arrange(id)
cluster$id2 = colnames(M_home)
cluster$descr = desc
cluster$possess = ball_possession
cluster <- cluster %>% left_join(m4[,c('event.id','desc')], by = c( 'id' = 'event.id'))
cluster <- cluster %>% left_join(df_home[,c('event_id','mean')], by = c( 'id' = 'event_id'))

cluster_gsw <- cluster %>% group_by(clus) %>% dplyr::summarise(ns = n(), 
                                                               prop = mean(descr),
                                                               meanID = mean(mean))

xtable(data.frame(cluster_gsw[,c(2:4)]), type = "latex", digits = 3)


# cluster CLE (away team)
cluster2 <- data.frame(clus = cle_attack$membership.cols) #cluster2 membership
cluster2$id <- as.numeric(as.character(rownames(cluster2)))
cluster2 <- cluster2 %>% arrange(id)
cluster2$id2 = colnames(M_away)
cluster2$descr = desc2
cluster2$possess = ball_possession2
cluster2 <- cluster2 %>% left_join(m4[,c('event.id','desc')], by = c( 'id' = 'event.id'))
cluster2 <- cluster2 %>% left_join(df_away[,c('event_id','mean')], by = c( 'id' = 'event_id'))

cluster_cle <- cluster2 %>% group_by(clus) %>% dplyr::summarise(ns = n(),
                                                                prop = mean(descr),
                                                                meanID = mean(mean))
cluster_cle

xtable(data.frame(cluster_cle[,c(2:4)]), type = "latex", digits = 3)
