library(INLA)
library(plyr)
library(ggplot2)
library(reshape)
library(latex2exp)
library(ggforce)
library(xtable)
library(scales)
library(gridExtra)

setwd("/u/s35he/algal_sim_results/scenario0")
true_param <- c(2, 0.25, -1, -1, 50, 20, 0, 0, 0.7, 7, 0.25)
true_param2 <- c(2, 0.25, -1, -1, 0.02, 0.05, 0, 0, 0.7, 7, 0.25)

zz <- file("./sim_result_1_12_100.txt", open = "wt")
sink(zz ,type = "output")
sink(zz, type = "message")

summary_mean <- summary_rmse <- rmse_mean <- rmse_list <- time_mean <- error_term <- list()
insitu_true = c(true_param[1], true_param[3:4], true_param[6:length(true_param)])
satellite_true = c(true_param[1], true_param[3:5], true_param[7:length(true_param)])
fusion_true = true_param

for (k in 1:12){

    summary_insitu <- summary_satellite <- summary_fusion <- list()
    sample_insitu <- sample_satellite <- sample_fusion <- list()
    rmse_insitu <- rmse_satellite <- rmse_fusion <- list()
    time_insitu <- time_satellite <- time_fusion <- list()
    error_term <- c()
    
    for (seed in 1:100){
    tryCatch({
    print(c(k, seed))
    load(file=paste0("./files/spde_simulation_", seed, "_", k, ".rda"))
    
    summary_insitu[[seed]] = insitu_list$summary$summary
    summary_satellite[[seed]] = satellite_list$summary$summary
    summary_fusion[[seed]] = fusion_list$summary$summary
    
    sample_insitu[[seed]] = sapply(1:ncol(insitu_list$sample.par), function(i) sqrt(mean((insitu_list$sample.par[,i]-insitu_true[i])^2)))
    sample_satellite[[seed]] = sapply(1:ncol(satellite_list$sample.par), function(i) sqrt(mean((satellite_list$sample.par[,i]-satellite_true[i])^2)))
    sample_fusion[[seed]] = sapply(1:ncol(fusion_list$sample.par), function(i) sqrt(mean((fusion_list$sample.par[,i]-fusion_true[i])^2)))
    
    rmse_insitu[[seed]] = insitu_list$summary$rmse$pred_rd
    rmse_satellite[[seed]] = satellite_list$summary$rmse$pred_rd
    rmse_fusion[[seed]] = fusion_list$summary$rmse$pred_rd
    
    time_insitu[[seed]] = insitu_list$summary$time
    time_satellite[[seed]] = satellite_list$summary$time
    time_fusion[[seed]] = fusion_list$summary$time
    
    rm(rspde_list, insitu_list, satellite_list, fusion_list)
    
    }, error= function(e){
      message(paste("An error occurred for item", seed,":\n"), e)
    })

    }
    
    error_term[[k]] <- 0
    for (seed in 1:length(summary_insitu)){
      if(is.null(summary_insitu[[seed]])){
        error_term[[k]] = c(error_term[[k]], seed)
      }
    }
    
    if (length(summary_insitu) < 100) error_term[[k]] = c(error_term[[k]], seq(length(summary_insitu)+1, 100))
    
    error_term[[k]] <- error_term[[k]][-1]
    
    if(length(which(sapply(summary_insitu, is.null))) > 0){

    summary_insitu = summary_insitu[-which(sapply(summary_insitu, is.null))]
    summary_satellite = summary_satellite[-which(sapply(summary_satellite, is.null))]
    summary_fusion = summary_fusion[-which(sapply(summary_fusion, is.null))]

    rmse_insitu = rmse_insitu[-which(sapply(rmse_insitu, is.null))]
    rmse_satellite = rmse_satellite[-which(sapply(rmse_satellite, is.null))]
    rmse_fusion = rmse_fusion[-which(sapply(rmse_fusion, is.null))]
    }
    
    summary_insitu_mean = do.call(cbind, lapply(summary_insitu, function(x) unlist(x[,1])))
    summary_satellite_mean = do.call(cbind, lapply(summary_satellite, function(x) unlist(x[,1])))
    summary_fusion_mean = do.call(cbind, lapply(summary_fusion, function(x) unlist(x[,1])))
    
    summary_mean[[k]] = list(rbind(summary_insitu_mean[1,], NA, summary_insitu_mean[2:3, ], NA, summary_insitu_mean[4:nrow(summary_insitu_mean), ]),
                              rbind(summary_satellite_mean[1,],NA, summary_satellite_mean[2:4, ], NA, summary_satellite_mean[5:nrow(summary_satellite_mean), ]),
                              summary_fusion_mean)
    
    summary_insitu_rmse = do.call(cbind, sample_insitu)
    summary_satellite_rmse = do.call(cbind, sample_satellite)
    summary_fusion_rmse = do.call(cbind, sample_fusion)

    summary_rmse[[k]] = list(rbind(summary_insitu_rmse[1,], NA, summary_insitu_rmse[2:3, ], NA, summary_insitu_rmse[4:nrow(summary_insitu_rmse), ]),
                              rbind(summary_satellite_rmse[1,],NA, summary_satellite_rmse[2:4, ], NA, summary_satellite_rmse[5:nrow(summary_satellite_rmse), ]),
                              summary_fusion_rmse)
    
    rmse_mean[[k]] = cbind(apply(laply(rmse_insitu, as.matrix), 2, mean), 
                           apply(laply(rmse_satellite, as.matrix), 2, mean), 
                           apply(laply(rmse_fusion, as.matrix), 2, mean))
    
    rmse_list[[k]] = list(laply(rmse_insitu, as.matrix), 
                          laply(rmse_satellite, as.matrix), 
                          laply(rmse_fusion, as.matrix))
    
    time_mean[[k]] = cbind(unlist(time_insitu), unlist(time_satellite), unlist(time_fusion))
    
    rm(summary_insitu, summary_satellite, summary_fusion, rmse_insitu, rmse_satellite, rmse_fusion)
      
} 

save(summary_mean, summary_rmse, rmse_mean, rmse_list, time_mean, error_term, file="./sim_comb_1_12_100.rda")

sink()

load("./sim_comb_1_12_100.rda")
k <- 12
num_grd_pts = c(5, 30)
missingness = c(0.5, 0.8)
mesh_size = c(0.15, 0.1, 0.05)
true_param = c(2, 0.25, -1, -1, 50, 0, 0, 0.7, 7, 0.25)
true_param2 = c(2, 0.25, -1, -1, 0.02, 0, 0, 0.7, 7, 0.25)
param=c("intercept", "bias", "beta1", "beta2", "precision", "theta1", "theta2", "rho", "kappa", "variance")
param_levels=c("intercept", "bias", "beta1", "beta2", "rho", "theta1", "theta2", "precision", "kappa", "variance")
length.pts = length(num_grd_pts)
length.missingness = length(missingness)
length.mesh = length(mesh_size)
length.param = length(true_param)

########## param bias plot

summary_df <- lapply(1:k, function(x) data.frame(rbind(summary_mean[[x]][[1]], summary_mean[[x]][[2]], summary_mean[[x]][[3]])))
summary_df <- lapply(1:k, function(x) cbind(summary_df[[x]], type=rep(c("Insitu", "Satellite", "Fusion"), each=length.param), 
                                            k = x, param=rep(param,3), true=rep(true_param, 3)))
summary_df <- lapply(1:k, function(x) melt(summary_df[[x]], id.vars=c("type", "k", "param", "true")))
summary_df <- data.frame(do.call(rbind, summary_df))
summary_df$variable <- as.numeric(summary_df$variable)
summary_df$mesh_size <- as.factor(mesh_size[ceiling(summary_df$k/(k/length.mesh))])
summary_df$num_grd_pts <- as.factor(ifelse(summary_df$k %% 2 ==1, 5, 30))
summary_df$missingness <- as.factor(ifelse(summary_df$k %% 4 %in% c(1,2), 0.5, 0.8))
summary_df$k <- as.factor(summary_df$k)
summary_df$type <- factor(summary_df$type, levels=c("Insitu", "Satellite", "Fusion"),
                          labels = c('Insitu' = parse(text='In~situ'), 
                                     'Satellite'='Satellite', 
                                     'Fusion'='Fusion'))
summary_df$bias <- summary_df$value - summary_df$true
summary_df$param_name <- summary_df$param
summary_df$param <- factor(summary_df$param, levels=param_levels,
                           labels=c('intercept'=parse(text=TeX('$\\beta_{0}$')),
                                               'bias'=parse(text=TeX('$a$')),
                                               'beta1'=parse(text=TeX('$\\beta_{1}$')),
                                               'beta2'=parse(text=TeX('$\\beta_{2}$')),
                                               'rho'=parse(text=TeX('$\\rho$')),
                                               'theta1'=parse(text=TeX('$\\theta_1$')),
                                               'theta2'=parse(text=TeX('$\\theta_2$')),
                                               'precision'=parse(text=TeX('$\\tau_z$')),
                                               'kappa'=parse(text=TeX('$\\kappa$')),
                                               'variance'=parse(text=TeX('$\\sigma_\\omega^2$'))))
summary_df$yintercept <- ifelse(is.na(summary_df$bias), NA, 0)

tex <- data.frame("param"=NA, "k"=NA, "Insitu"=NA,
                  "Satellite"=NA,
                  "Fusion"=NA)[numeric(0), ]
for (par in param[c(1:5,8,9,10)]){
  for (i in 1:k){
    tmp1 <- summary_df$bias[summary_df$param_name == par & summary_df$k == i & summary_df$type=="In ~ situ"]
    tmp2 <- summary_df$bias[summary_df$param_name == par & summary_df$k == i & summary_df$type=="Satellite"]
    tmp3 <- summary_df$bias[summary_df$param_name == par & summary_df$k == i & summary_df$type=="Fusion"]
    tex <- rbind(tex, as.vector(c(par, i, round(mean(tmp1, na.rm=T), 4),
                                  round(mean(tmp2, na.rm=T),4),
                                  round(mean(tmp3, na.rm=T),4))))
  }
}
colnames(tex) = c("param", "k", "Insitu","Satellite","Fusion")
print(xtable(tex), include.rownames=FALSE)

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  #l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("1e\\+", "10^", l)
  l <- gsub("1e", "10^", l)
  # return this as an expression
  parse(text=l)
}

ann_text <- data.frame(k=6, bias=0, lab = "NA",
                       param = factor("bias", levels = param_levels,
                                      labels=c('intercept'=parse(text=TeX('$\\beta_{0}$')),
                                                        'bias'=parse(text=TeX('$a$')),
                                                        'beta1'=parse(text=TeX('$\\beta_{1}$')),
                                                        'beta2'=parse(text=TeX('$\\beta_{2}$')),
                                                        'rho'=parse(text=TeX('$\\rho$')),
                                                        'theta1'=parse(text=TeX('$\\theta_1$')),
                                                        'theta2'=parse(text=TeX('$\\theta_2$')),
                                                        'precision'=parse(text=TeX('$\\tau_z$')),
                                                        'kappa'=parse(text=TeX('$\\kappa$')),
                                                        'variance'=parse(text=TeX('$\\sigma_\\omega^2$')))),
                       type=factor("In situ", levels = c("In situ", "Satellite", "Fusion"),
                                   labels = c('Insitu' = parse(text='In~situ'), 
                                              'Satellite'='Satellite', 
                                              'Fusion'='Fusion')))
ann_text2 <- data.frame(k=6, bias=0, lab = "NA",
                        param = factor("bias", levels = param_levels,
                                       labels=c('intercept'=parse(text=TeX('$\\beta_{0}$')),
                                                'bias'=parse(text=TeX('$a$')),
                                                'beta1'=parse(text=TeX('$\\beta_{1}$')),
                                                'beta2'=parse(text=TeX('$\\beta_{2}$')),
                                                'rho'=parse(text=TeX('$\\rho$')),
                                                'theta1'=parse(text=TeX('$\\theta_1$')),
                                                'theta2'=parse(text=TeX('$\\theta_2$')),
                                                'precision'=parse(text=TeX('$\\tau_z$')),
                                                'kappa'=parse(text=TeX('$\\kappa$')),
                                                'variance'=parse(text=TeX('$\\sigma_\\omega^2$')))),
                        type=factor("Satellite", levels = c("In situ", "Satellite", "Fusion"),
                                    labels = c('Insitu' = parse(text='In~situ'), 
                                               'Satellite'='Satellite', 
                                               'Fusion'='Fusion')))
trans_ <- scales::trans_new("power_trans_neg",
                           transform=function(x) sign(x)*(abs(x)^(1/3)),
                           inverse=function(x) sign(x)*(abs(x)^3))
trans_2 <- scales::trans_new("power_trans_neg",
                            transform=function(x) sign(x)*(abs(x)^(1/9)),
                            inverse=function(x) sign(x)*(abs(x)^9))

pdf(paste0("./plot_param_bias.pdf"), width=7, height=10)
p1 <- ggplot(summary_df[summary_df$param_name %in% param[c(1:4,8)], ]) +
  geom_boxplot(aes(k, bias, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","type"), scales = "free_y", labeller=label_parsed) + 
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("Bias") + xlab("") + 
  geom_text(data = ann_text, aes(k, bias, label = lab)) +
  geom_text(data = ann_text2, aes(k, bias, label = lab)) +
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust=0.15),
        strip.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.background=element_rect(fill="white"),
        plot.margin = unit(c(0,1,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))

p2 <- ggplot(summary_df[summary_df$param_name %in% param[c(5)], ]) +
  geom_boxplot(aes(k, bias, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","type"), scales = "free_y", labeller=label_parsed) + 
  scale_y_continuous(trans=trans_2, breaks=c(-100,-10,-1,0,1,10,100,1000,10^4), labels=fancy_scientific) +
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("") + ggtitle("") + scale_fill_discrete(name = "Mesh max. \ninner length") +
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust=0.15),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.text.y = element_text(size = 14),
        strip.background=element_rect(fill="white"),
        plot.margin = unit(c(-1,1.18,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))

p3 <- ggplot(summary_df[summary_df$param_name %in% param[c(9,10)], ]) +
  geom_boxplot(aes(k, bias, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","type"), scales = "free_y", labeller=label_parsed) + 
  scale_y_continuous(trans=trans_, breaks=c(-100,-10,-1, 0,1,10,100,1000,10^4), labels=fancy_scientific) +
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("") + ggtitle("") + scale_fill_discrete(name = "Mesh max. \ninner length") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(-1,0.9,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))

grid.arrange(p1,p2,p3, nrow=3, layout_matrix=rbind(1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3))

dev.off()

pdf(paste0("./plot_param_bias_by_scenario.pdf"), width=7, height=10)
p1 <- ggplot(summary_df[summary_df$param_name %in% param[c(1:4,8)], ]) +
  geom_boxplot(aes(type, bias, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","k"), scales = "free_y", labeller=label_parsed) + 
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("Bias") + xlab("") + 
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust=0.15),
        strip.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.background=element_rect(fill="white"),
        plot.margin = unit(c(0,1,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))

p2 <- ggplot(summary_df[summary_df$param_name %in% param[c(5)], ]) +
  geom_boxplot(aes(type, bias, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","k"), scales = "free_y", labeller=label_parsed) + 
  scale_y_continuous(trans=trans_2, breaks=c(-100,-10,-1,0,1,10,100,1000,10^4), labels=fancy_scientific) +
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("") + ggtitle("") + scale_fill_discrete(name = "Mesh max. \ninner length") +
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust=0.15),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.text.y = element_text(size = 14),
        strip.background=element_rect(fill="white"),
        plot.margin = unit(c(-1,1.18,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))

p3 <- ggplot(summary_df[summary_df$param_name %in% param[c(9,10)], ]) +
  geom_boxplot(aes(type, bias, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","k"), scales = "free_y", labeller=label_parsed) + 
  scale_y_continuous(trans=trans_, breaks=c(-100,-10,-1, 0,1,10,100,1000,10^4), labels=fancy_scientific) +
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("") + ggtitle("") + scale_fill_discrete(name = "Mesh max. \ninner length") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust=1),
        plot.margin = unit(c(-1,0.9,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))

grid.arrange(p1,p2,p3, nrow=3, layout_matrix=rbind(1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3))

dev.off()

########## param rmse plot

summary_df_rmse <- lapply(1:k, function(x) data.frame(rbind(summary_rmse[[x]][[1]], summary_rmse[[x]][[2]], summary_rmse[[x]][[3]])))
summary_df_rmse <- lapply(1:k, function(x) cbind(summary_df_rmse[[x]], type=rep(c("Insitu", "Satellite", "Fusion"), each=length.param), 
                                            k = x, param=rep(param,3), true=rep(true_param2, 3)))
summary_df_rmse <- lapply(1:k, function(x) melt(summary_df_rmse[[x]], id.vars=c("type", "k", "param", "true")))
summary_df_rmse <- data.frame(do.call(rbind, summary_df_rmse))
summary_df_rmse$variable <- as.numeric(summary_df_rmse$variable)
summary_df_rmse$mesh_size <- as.factor(mesh_size[ceiling(summary_df_rmse$k/(k/length.mesh))])
summary_df_rmse$k <- as.factor(summary_df_rmse$k)
summary_df_rmse$type <- factor(summary_df_rmse$type, levels=c("Insitu", "Satellite", "Fusion"),
                               labels = c('Insitu' = parse(text='In~situ'), 
                                          'Satellite'='Satellite', 
                                          'Fusion'='Fusion'))
summary_df_rmse$param_name <- summary_df_rmse$param
summary_df_rmse$param <- factor(summary_df_rmse$param, levels=param, 
                           labels=c('intercept'=parse(text=TeX('$\\beta_{0}$')),
                                    'bias'=parse(text=TeX('$a$')),
                                    'beta1'=parse(text=TeX('$\\beta_{1}$')),
                                    'beta2'=parse(text=TeX('$\\beta_{2}$')),
                                    'precision'=parse(text=TeX('$\\tau_z$')),
                                    'theta1'=parse(text=TeX('$\\theta_1$')),
                                    'theta2'=parse(text=TeX('$\\theta_2$')),
                                    'rho'=parse(text=TeX('$\\rho$')),
                                    'kappa'=parse(text=TeX('$\\kappa$')),
                                    'variance'=parse(text=TeX('$\\sigma_\\omega^2$'))))

summary_df_rmse$yintercept <- ifelse(is.na(summary_df_rmse$value), NA, ifelse(summary_df_rmse$param_name %in% param[c(5,9,10)], -Inf, 0))

tex <- data.frame("param"=NA, "k"=NA, "In situ"=NA,
                  "Satellite"=NA,
                  "Fusion"=NA)[numeric(0), ]
for (par in param_levels[c(1:5,8,9,10)]){
  for (i in 1:k){
    tmp1 <- summary_df_rmse$value[summary_df_rmse$param_name == par & summary_df_rmse$k == i & summary_df_rmse$type=="In ~ situ"]
    tmp2 <- summary_df_rmse$value[summary_df_rmse$param_name == par & summary_df_rmse$k == i & summary_df_rmse$type=="Satellite"]
    tmp3 <- summary_df_rmse$value[summary_df_rmse$param_name == par & summary_df_rmse$k == i & summary_df_rmse$type=="Fusion"]
    tex <- rbind(tex, as.vector(c(par, i, round(mean(tmp1, na.rm=T), 4),
                                  round(mean(tmp2, na.rm=T),4),
                                  round(mean(tmp3, na.rm=T),4))))
  }
}

colnames(tex) = c("param", "k", "In situ","Satellite","Fusion")
print(xtable(tex), include.rownames=FALSE)


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e\\+", "10^", l)
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

ann_text <- data.frame(k=6, value=0.075, lab = "NA",
                       param = factor("bias", levels = param_levels,
                                      labels=c('intercept'=parse(text=TeX('$\\beta_{0}$')),
                                               'bias'=parse(text=TeX('$a$')),
                                               'beta1'=parse(text=TeX('$\\beta_{1}$')),
                                               'beta2'=parse(text=TeX('$\\beta_{2}$')),
                                               'rho'=parse(text=TeX('$\\rho$')),
                                               'theta1'=parse(text=TeX('$\\theta_1$')),
                                               'theta2'=parse(text=TeX('$\\theta_2$')),
                                               'precision'=parse(text=TeX('$\\tau_z$')),
                                               'kappa'=parse(text=TeX('$\\kappa$')),
                                               'variance'=parse(text=TeX('$\\sigma_\\omega^2$')))),
                       type=factor("In situ", levels = c("In situ", "Satellite", "Fusion"),
                                   labels = c('Insitu' = parse(text='In~situ'), 
                                              'Satellite'='Satellite', 
                                              'Fusion'='Fusion')))
ann_text2 <- data.frame(k=6, value=0.075, lab = "NA",
                        param = factor("bias", levels = param_levels,
                                       labels=c('intercept'=parse(text=TeX('$\\beta_{0}$')),
                                                'bias'=parse(text=TeX('$a$')),
                                                'beta1'=parse(text=TeX('$\\beta_{1}$')),
                                                'beta2'=parse(text=TeX('$\\beta_{2}$')),
                                                'rho'=parse(text=TeX('$\\rho$')),
                                                'theta1'=parse(text=TeX('$\\theta_1$')),
                                                'theta2'=parse(text=TeX('$\\theta_2$')),
                                                'precision'=parse(text=TeX('$\\tau_z$')),
                                                'kappa'=parse(text=TeX('$\\kappa$')),
                                                'variance'=parse(text=TeX('$\\sigma_\\omega^2$')))),
                        type=factor("Satellite", levels = c("In situ", "Satellite", "Fusion"),
                                    labels = c('Insitu' = parse(text='In~situ'), 
                                               'Satellite'='Satellite', 
                                               'Fusion'='Fusion')))

pdf(paste0("./plot_param_rmse.pdf"), width=7, height=10)
p1 <- ggplot(summary_df_rmse[summary_df_rmse$param_name %in% param[c(1:4,8)], ]) +
          geom_boxplot(aes(k, value, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
          facet_grid(c("param","type"), scales = "free_y", labeller=label_parsed) + 
          geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
          ylab("RMSE") + xlab("") + 
          geom_text(data = ann_text, aes(k, value, label = lab)) +
          geom_text(data = ann_text2, aes(k, value, label = lab)) +
          theme(legend.position = "none", 
                axis.title.y = element_text(hjust=0.15),
                strip.text = element_text(size = 14),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 8),
                strip.background=element_rect(fill="white"),
                plot.margin = unit(c(0,1,0,1), "cm"),
                panel.background = element_rect(fill="white", colour = 'lightgrey'))
p2 <- ggplot(summary_df_rmse[summary_df_rmse$param_name %in% param[c(5,9,10)], ]) +
        geom_boxplot(aes(k, value, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
        facet_grid(c("param","type"), scales = "free_y", labeller=label_parsed) + 
        scale_y_continuous(trans=log_trans(), breaks=c(0.01,0.1,1,10,100,1000,10^4,10^5,10^6, 10^7,10^8), labels=fancy_scientific) +
        geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
        ylab("") + ggtitle("") + scale_fill_discrete(name = "Mesh max. \ninner length") +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              strip.text.y = element_text(size = 14),
              axis.text = element_text(size = 8),
              plot.margin = unit(c(-1,0.8,0,1), "cm"),
              panel.background = element_rect(fill="white", colour = 'lightgrey'))
grid.arrange(p1,p2, nrow=2, layout_matrix=rbind(1,1,1,2,2))
dev.off()

pdf(paste0("./plot_param_rmse_by_scenario.pdf"), width=7, height=10)
p1 <- ggplot(summary_df_rmse[summary_df_rmse$param_name %in% param[c(1:4,8)], ]) +
  geom_boxplot(aes(type, value, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","k"), scales = "free_y", labeller=label_parsed) + 
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("RMSE") + xlab("") + 
  theme(legend.position = "none", 
        axis.title.y = element_text(hjust=0.15),
        strip.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.background=element_rect(fill="white"),
        plot.margin = unit(c(0,1,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))
p2 <- ggplot(summary_df_rmse[summary_df_rmse$param_name %in% param[c(5,9,10)], ]) +
  geom_boxplot(aes(type, value, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
  facet_grid(c("param","k"), scales = "free_y", labeller=label_parsed) + 
  scale_y_continuous(trans=log_trans(), breaks=c(0.01,0.1,1,10,100,1000,10^4,10^5,10^6, 10^7,10^8), labels=fancy_scientific) +
  geom_hline(aes(yintercept=yintercept), linetype=2, linewidth=0.3) + 
  ylab("") + ggtitle("") + scale_fill_discrete(name = "Mesh max. \ninner length") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust=1),
        plot.margin = unit(c(-1,0.8,0,1), "cm"),
        panel.background = element_rect(fill="white", colour = 'lightgrey'))
grid.arrange(p1,p2, nrow=2, layout_matrix=rbind(1,1,1,2,2))
dev.off()

########## rmse plot

rmse_df <- lapply(1:k, function(x) data.frame(rbind(rmse_list[[x]][[1]], rmse_list[[x]][[2]], rmse_list[[x]][[3]])))
rmse_df <- lapply(1:k, function(x) cbind(rmse_df[[x]], type=rep(c("Insitu", "Satellite", "Fusion"), each=nrow(rmse_list[[x]][[1]])), 
                                                 k = x))
rmse_df <- lapply(1:k, function(x) melt(rmse_df[[x]], id.vars=c("type", "k")))
rmse_df <- data.frame(do.call(rbind, rmse_df))
rmse_df$day <- as.numeric(rmse_df$variable)
rmse_df$mesh_size <- as.factor(mesh_size[ceiling(rmse_df$k/(k/length.mesh))])
rmse_df$k <- as.factor(rmse_df$k)
rmse_df$type <- factor(rmse_df$type, levels=c("Insitu", "Satellite", "Fusion"),
                          labels = c('Insitu' = parse(text='In~situ'), 
                                     'Satellite'='Satellite', 
                                     'Fusion'='Fusion'))

pdf(paste0("./plot_pred_rmse.pdf"), height=6, width=7)
print(ggplot(rmse_df[rmse_df$day %in% c(1,5,10,15,19),]) +
        geom_boxplot(aes(k, value, fill=mesh_size), outlier.size=0.1, lwd=0.01) + 
        facet_grid(c("day","type"), scales = "free_y", labeller=label_parsed) + 
        ylab("RMSE") + scale_fill_discrete(name = "Mesh max. \ninner length") +
        theme(legend.position = "bottom", 
              strip.text = element_text(size = 14),
              axis.text = element_text(size = 8),
              strip.background=element_rect(fill="white"),
              panel.background = element_rect(fill="white", colour = 'lightgrey')))

dev.off()

rmse_df <- lapply(1:k, function(x) rmse_mean[[x]] <- cbind(rmse_mean[[x]], 
                                                           rep(mesh_size[ceiling(x/(k/length.mesh))], nrow(rmse_mean[[x]]))))
rmse_df <- data.frame(do.call(rbind, lapply(rmse_df, function(x) unlist(x))))
colnames(rmse_df) <- c("Insitu", "Satellite", "Fusion", "mesh_size")
rmse_df$k <- as.factor(rep(1:k, each=nrow(rmse_mean[[1]])))
rmse_df$day <- as.factor(rep(1:nrow(rmse_mean[[1]]), k))
rmse_df$num_grd_pts <- as.factor(rep(rep(num_grd_pts, each=nrow(rmse_mean[[1]])), k/length.pts))
rmse_df$missingness <- as.factor(rep(rep(missingness, each=length.pts*nrow(rmse_mean[[1]])), k/length.pts/length.missingness))
rmse_df$mesh_size <- as.factor(rmse_df$mesh_size)
rmse_df_agg <- melt(rmse_df, id.vars=c("mesh_size","k", "day", "num_grd_pts", "missingness"), variable_name="type")


rects <- data.frame(xstart = c(1,14), xend = c(14,19), Type = c("Train","Test"))
rects$Type <- factor(rects$Type, levels=c("Train","Test"))

pdf(paste0("./plot_pred_rmse_ts.pdf"), height=7, width=10)
print(ggplot(rmse_df_agg) +
        geom_line(aes(day, value, group=type, col=type)) +
        geom_point(aes(day, value, group=type, col=type), size=1) +
        geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = Type), alpha = 0.2) +
        facet_wrap(.~k, ncol=4) +
        scale_fill_manual(values=c('Train' = 'white', 'Test'='grey'))+
        # scale_color_brewer(palette = "Set2") +
        guides(color=guide_legend(title="Model")) +
        xlab("Day") + ylab("RMSE") + 
        theme(legend.position = "bottom", 
              strip.text = element_text(size = 12),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 8),
              strip.background=element_rect(fill="white"),
              panel.background = element_rect(fill="white", colour = 'lightgrey')))
dev.off()


