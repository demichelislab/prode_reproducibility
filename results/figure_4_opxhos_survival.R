#-------------------------------------------------------------------------------
#                            Figure 5 - OXPHOS 
#_______________________________________________________________________________

library(patchwork)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(dplyr)

# Load data ====================================================================

prode_res <- fread('./tables/supp_table_6.txt', 
                   data.table=F)

oxphos_gns <- fread('./tables/supp_table_7.txt', 
                   data.table=F)

## TCGA data 
tcga_dat <- readRDS('./data/brca_tcga_data.rds')

# OXPHOS activity expession ====================================================

## Compute activity scores by fitting log-norm reg. for each BRCA pat. 

{
    bin_v   <- (rownames(tcga_dat$expr) %in% oxphos_gns$oxphos_gene)*1 # binary vector
    
    mod_fit <- lm(log(tcga_dat$expr+1)~bin_v) # fit all models 
    
    scores  <- vapply(summary(mod_fit), function(x){ # extract t.values 
        coefficients(x)[2,3]
    }, .1)
}

# Plot -------------------------------------------------------------------------

p <- data.frame(
    'score'   = c(scores), 
    'subtype' = factor(tcga_dat$surv$subtype)
) %>% 
    ggplot(   
        aes(
            x = reorder(subtype, score), 
            y = score
        )
    ) +
    geom_violin(aes(fill = subtype == 'Her2')) + 
    geom_boxplot(width = 0.2, outlier.size = .5) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c('Her2', 'Other')
        ), hide.ns = F,label = 'p.signif'
    ) +
    theme_light() +
    labs(x = 'Molecular Subtype') +
    ggpubr::stat_compare_means(
        comparisons = list(
            c('Her2', 'Basal'),
            c('Her2', 'LumB'),
            c('Her2', 'LumA')
        ), label = 'p.signif', hide.ns = T
    ) +
    labs(y='OXPHOS genes activity\n coefficient ') +
    theme(legend.position='none') +
    scale_fill_manual(
        values = c('lightgray', '#996662')
    ) 

ggsave(plot=p, './figures/fig4_oxphos_ass_coefs.pdf', height = 3.2, width = 3)

# TCGA survival analysis =======================================================

## Retrieve PRODE private OXPHOS genes .........................................

oxph_pro <- oxphos_gns %>% filter(top100_PRODE == 'yes') %>% pull(oxphos_gene)

## Compute score (averge oxph_pro expression) ..................................

surv <- tcga_dat$surv
expr <- tcga_dat$expr

score <- rowMeans(apply(expr[oxph_pro,], 1, scale))

ttype <- ifelse(grepl('Her2', surv$subtype), 1, 0) # subtype labeling 

## Compose dataset for survival analysis .......................................

dt <- data.frame(
    sub  = ttype,  
    sco  = score, 
    age  = surv$age,
    sex  = surv$sex,
    gii  = surv$gii,
    sur  = surv$survival,
    cen  = surv$survival_censored
) %>% 
    group_by(sub) %>% 
    mutate( # Divide into two groups based on median score
        'group'  = ifelse(
            sco <= quantile(sco, .5), 
            paste0('OXPHOS', ' low'), 
            paste0('OXPHOS', ' high')
        ),
    ) 

## Fit and plot ................................................................

{
    ## Her2+
    fit <- survfit(Surv(sur, cen)~group, data=dt %>% filter(sub == 1))
    
    pa <- ggsurvplot(fit, pval = T, cumevents = T)
    
    p1 <- pa[[1]] + theme_light() + 
        scale_color_manual(values = 
                               c('gray', 'steelblue', 'red')) + 
        labs(col = '') + 
        ggtitle('HER2+') +
        theme(legend.position='none') 
    
    p3 <- pa[[2]]$data %>% 
        ggplot() +
        geom_text(
            aes(
                x = time,
                y = group,
                label = n.risk,
                col = group
            ), size = 2.5
        )+
        theme_light() +
        scale_color_manual(
            values = c('gray', 'steelblue', 'red')
        ) +
        theme(legend.position='none') +
        ylab('') +
        xlab('Time') 
    
    ## not Her2+
    fit <- survfit(Surv(sur, cen)~group, data=dt %>% filter(sub == 0))
    
    pb <- ggsurvplot(fit, pval = T, risk.table = T)
    
    p2 <- pb[[1]] + theme_light() + 
        scale_color_manual(
            values =  c('gray', 'steelblue', 'red')) + 
        labs(col = '') + 
        ggtitle('Others') +
        theme(legend.position='none')+
        ylab('')
    
    p4 <- pb[[2]]$data %>% 
        ggplot() +
        geom_text(
            aes(
                x = time,
                y = group,
                label = n.risk,
                col = group
            ),  size = 2.5
        )+
        theme_light() +
        scale_color_manual(
            values = c('gray', 'steelblue', 'red')
        ) +
        theme(legend.position='none') +
        ylab('') +
        xlab('Time')+
        theme(axis.text.y = element_blank())
    
    }


p1 + p2 + p3  + p4 +patchwork::plot_layout(ncol = 2, nrow = 2, heights = c(3,1))






















