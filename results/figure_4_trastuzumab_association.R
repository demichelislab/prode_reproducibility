#_______________________________________________________________________________
#                       Figure 4 - association OXPHOS / Trastuzumab
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)

set.seed(100)

# Load data ====================================================================

## PRODE top oxphos 
oxph_pro <- fread('./tables/supp_table_7.txt', data.table=F)

dt1 <- readRDS('./data/cohort_transneo.rds')
dt2 <- readRDS('./data/cohort_GSE66399.rds')
dt3 <- readRDS('./data/cohort_GSE50948.rds')
dt4 <- readRDS('./data/cohort_GSE42822.rds')
dt5 <- readRDS('./data/cohort_GSE37946.rds')

val_dts <- list(dt1, dt2, dt3, dt4, dt5)

# Pre-process expression data though scaling ----------------------------------

val_dts[[1]]$mRNA2 <- t(apply(val_dts[[1]]$mRNA, 1, function(x) scale(x)))
val_dts[[2]]$mRNA2 <- t(apply(val_dts[[2]]$mRNA, 1, function(x) scale(x)))
val_dts[[3]]$mRNA2 <- t(apply(val_dts[[3]]$mRNA, 1, function(x) scale(x)))
val_dts[[4]]$mRNA2 <- t(apply(val_dts[[4]]$mRNA, 1, function(x) scale(x)))
val_dts[[5]]$mRNA2 <- t(apply(val_dts[[5]]$mRNA, 1, function(x) scale(x)))

# Dataset names 
dt1 <- 'transneo'
dt2 <- 'GSE66399'
dt3 <- 'GSE50948'
dt4 <- 'GSE42822'
dt5 <- 'GSE37946'


# Compute association coefficients to response ---------------------------------

oxph_pro <- oxph_pro %>% filter(top100_PRODE == 'yes') %>% pull(oxphos_gene)

all_pvals <- list()
all_coefs <- list()

for (id in 1:5){
    
    oxph.tmp <- oxph_pro[which(oxph_pro %in% rownames(val_dts[[id]]$mRNA2))]
    
    coefs <- vapply(oxph.tmp, function(x){
        
        # prepare vectors
        v1 <- rank(val_dts[[id]]$mRNA[x,], na.last = 'keep')/sum(!is.na(val_dts[[id]]$mRNA[x,]))
        v2 <- rank(val_dts[[id]]$mRNA['ERBB2',], na.last = 'keep')/sum(!is.na(val_dts[[id]]$mRNA['ERBB2',]))
        
        fit <- glm(val_dts[[id]]$outcome~v1+v2)
        coefficients(summary(fit))[2,3]
        
    }, .1)
    
    all_coefs[[id]] <- coefs
    
    pvals <- vapply(oxph.tmp, function(x){
        
        # prepare vectors
        v1 <- rank(val_dts[[id]]$mRNA[x,], na.last = 'keep')/sum(!is.na(val_dts[[id]]$mRNA[x,]))
        v2 <- rank(val_dts[[id]]$mRNA['ERBB2',], na.last = 'keep')/sum(!is.na(val_dts[[id]]$mRNA['ERBB2',]))
        
        fit <- glm(val_dts[[id]]$outcome~v1+v2)
        coefficients(summary(fit))[2,4]
        
    }, .1)
    
    all_pvals[[id]] <- pvals
    
}

all_coefs <- as.data.frame(bind_rows(all_coefs))
rownames(all_coefs) <- c(dt1, dt2, dt3, dt4, dt5)

all_pvals <- as.data.frame(bind_rows(all_pvals))
rownames(all_pvals) <- c(dt1, dt2, dt3, dt4, dt5)

#all_coefs <- t(apply(all_coefs, 1, function(x) (x - m1) / s1))

all_coefs %>%
    as.matrix() %>%
    reshape2::melt() %>%
    ggplot() +
    geom_point(
        aes(
            x = Var2,
            y = Var1,
            fill = value
        ), shape = 21,
        size = 3.5
    ) +
    geom_text(
        aes(
            x = Var2,
            y = Var1,
            label = ifelse(value <= 0.05, '*', '')
        ), data = all_pvals %>% as.matrix() %>% reshape2::melt(),
        vjust=0.75,
        size = 5
    ) +
    scale_fill_gradient2(
        low = 'steelblue',
        high = '#f55142',
        mid = 'white',
        limits=c(-3,3)
    ) + theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position='top') +
    labs(x = 'OXPHOS genes', y='', fill = 'Association coefficent to\n Trastuzumab response')

# Comparing to background distribution =========================================

mms <- rowMeans(all_coefs, na.rm=T)
sds <- apply(all_coefs,1,sd,na.rm=T)

all_random <- list(
    'a' = c(),
    'b' = c(),
    'c' = c(),
    'd' = c(),
    'e' = c()
)

# NOTE: this takes few min. - it computes association coefficients 
# by taking randomly samples genes instead of the 29 OXPHOS selected by PRODe. 
for (id in 1:5){ 
    
    print(id)
    
    oo <- lapply(1:1000, function(i){
        
        coefs <- unlist(lapply(sample(1:nrow(val_dts[[id]]$mRNA), length(oxph_pro)), function(x){
            
            # prepare vectors
            v1 <- rank(val_dts[[id]]$mRNA[x,], na.last = 'keep')/sum(!is.na(val_dts[[id]]$mRNA[x,]))
            v2 <- rank(val_dts[[id]]$mRNA['ERBB2',], na.last = 'keep')/sum(!is.na(val_dts[[id]]$mRNA['ERBB2',]))
            
            fit <- glm(val_dts[[id]]$outcome~v1+v2)
            coefficients(summary(fit))[2,3]
            
        }))
        
        mean(coefs)
        
    })
    
    all_random[[id]] <- unlist(oo)
    
}

random <- as.data.frame(bind_cols(all_random))
colnames(random) <- c(dt1, dt2, dt3, dt4, dt5)

m1 <- colMeans(random)
s1 <- apply(random, 2, sd)

random <- t(apply(random, 1, function(x) (x - m1)/s1))

scores1 <- (mms - m1) / s1

mean(random[,1] <= scores1[1]) # "GSE37946"
mean(random[,2] <= scores1[2]) # "GSE42822"
mean(random[,3] <= scores1[3]) # "GSE50948"
mean(random[,4] <= scores1[4]) # "GSE66399_trastuzumab"
mean(random[,5] <= scores1[5]) # "transneo"

dtt <- data.frame(
    variable = c(dt1, dt2, dt3, dt4, dt5),
    value = scores1
)

random %>%
    as.data.frame() %>%
    reshape2::melt() %>%
    ggplot() +
    geom_point(
        aes(
            x = value,
            y = variable
        ), data = dtt,
        shape=21,
        size = 3,
        fill = 'gray'
    ) +
    geom_hline(yintercept=c(1:5), linetype='dashed', col='gray')+
    ggridges::geom_density_ridges(
        aes(
            x = value,
            y = variable
        ), alpha = 0, height=0.2,
        col = 'gray',
    )  +
    theme_light()+
    labs(
        x = 'Average association coefficients to Trastuzumab response',
        y = ''
    )
 ,






