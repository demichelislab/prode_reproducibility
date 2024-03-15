# Analysis on robustness to noise in input data --------------------------------

# Imports ======================================================================

library(data.table)
library(prodeTool)
library(ggplot2)
library(dplyr)

set.seed(100)

# Load data ====================================================================

sg_data <- readRDS('./data/data_22Q1_sgRNA.rds')

sh_data <- readRDS('./data/data_22Q1_shRNA.rds')

ppi_data <- fread('./data/edge_list_combined_400_genes.txt', data.table=F)

# Run analysis on sgRNA dataset ================================================

## For each standard dev. thr compute the profile of essentiality scores .......

ge_dt <- t(sg_data$crispr) # input Gene Knock-Out Effcts
co_dt <- sg_data$meta_data # input column data

# Fill out NA with mean scores .................................................

mm <- matrix(
    rep(
        rowMeans(ge_dt, na.rm=T), 
        ncol(ge_dt)
    ), 
    ncol = ncol(ge_dt)
)

ge_dt[which(is.na(ge_dt))] <- mm[which(is.na(ge_dt))]

# Run analysis .................................................................

sg_out <- lapply(c(0, .5, 1, 2, 3, 4, 5), function(sd_thr){
    
    print(sd_thr)
    
    noise_m <- matrix(
        rnorm(
            n = nrow(ge_dt)*ncol(ge_dt),
            mean = 0, sd = sd_thr
        ), 
        ncol = ncol(ge_dt), 
        nrow = nrow(ge_dt)
    )
    
    prode_input <- getProdeInput(
        score_matrix = ge_dt + noise_m, # add noise to Gene Eff. 
        col_data     = co_dt,
        edge_table   = ppi_data
    )
    
    prode_results <- runProde(
        prodeInput = prode_input, 
        scaledEst  = F, 
        filterCtrl = F
    )
    
    results(prode_results)
    
})

# Run analysis on shRNA dataset ================================================

## For each standard dev. thr compute the profile of essentiality scores .......

ge_dt <- t(sh_data$crispr) # input Gene Knock-Out Effects
co_dt <- sh_data$meta_data # input column data

# Fill out NA with mean scores .................................................

mm <- matrix(rep(rowMeans(ge_dt, na.rm=T), ncol(ge_dt)), ncol = ncol(ge_dt))
ge_dt[which(is.na(ge_dt))] <- mm[which(is.na(ge_dt))]

# ..............................................................................

sh_out <- lapply(c(0, .5, 1, 2, 3, 4, 5), function(sd_thr){
    
    print(sd_thr)
    
    noise_m <- matrix(
        rnorm(
            n = nrow(ge_dt)*ncol(ge_dt),
            mean = 0, sd = sd_thr
        ), 
        ncol = ncol(ge_dt), 
        nrow = nrow(ge_dt)
    )
    
    prode_input <- getProdeInput(
        score_matrix = ge_dt + noise_m, # add noise to Gene Eff. 
        col_data     = co_dt,
        edge_table   = ppi_data
    )
    
    prode_results <- runProde(
        prodeInput = prode_input, 
        scaledEst  = F, 
        filterCtrl = F
    )
    
    results(prode_results)
    
})

# Compute correlations between scores ==========================================

dt_prode_sg <- do.call(cbind, lapply(sg_out, function(x) x$NIE_score))
dt_prode_sh <- do.call(cbind, lapply(sh_out, function(x) x$NIE_score))
dt_avgge_sg <- do.call(cbind, lapply(sg_out, function(x) x$Estimate))
dt_avgge_sh <- do.call(cbind, lapply(sh_out, function(x) x$Estimate))

cor_pro_sg <- apply(dt_prode_sg, 2, function(x){
    cor(x, dt_prode_sg[,1], method = 'spearman')
})

cor_avg_sg <- apply(dt_avgge_sg, 2, function(x){
    cor(x, dt_avgge_sg[,1], method = 'spearman')
})

cor_pro_sh <- apply(dt_prode_sh, 2, function(x){
    cor(x, dt_prode_sh[,1], method = 'spearman')
})

cor_avg_sh <- apply(dt_avgge_sh, 2, function(x){
    cor(x, dt_avgge_sh[,1], method = 'spearman')
})

# Plot =========================================================================

data.frame(
    'cor' = c(cor_pro_sg, cor_avg_sg, cor_pro_sh, cor_avg_sh), 
    'met' = rep(rep(c('PRODE', 'Average Gene Eff.'), each = c(7)), 2), 
    'nms' = rep(paste0(rep(c(0, 0.5, 1, 2, 3, 4, 5), 2), 'sd'), 2), 
    'dat' = rep(c('sgRNA', 'shRNA'), each = 14)
) %>%
    filter(nms != '0.5sd') %>% 
    ggplot() + 
    geom_line(
        aes(
            x = nms, 
            y = cor,
            group = paste0(met, dat)
        ), 
        col = 'black', 
        linetype = 'dashed'
    ) + 
    geom_point(
        aes(
            x = nms, 
            y = cor, 
            fill = met, 
            shape = dat
        ), 
        size = 2.5
    )  + 
    theme_light() + 
    labs(
        y = 'Spearman Corr. \n (noise vs. no-noise)', 
        x = 'Gaussian-noise', 
        shape = '', 
        fill = ''
    ) + 
    theme(
        legend.position = 'top'
    ) + 
    ylim(0, 1) + 
    scale_fill_manual(
        values = c(
            "#eb9a8a", "#ffd261"
        )
    ) + 
    scale_shape_manual(
        values = c(21, 24)
    ) + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE), 
           shape=guide_legend(nrow=2,byrow=TRUE))


















