# ==============================================================================
#                Fig 3 - Performance of SL datasets 
# ==============================================================================

library(data.table)
library(patchwork)
library(dplyr)
library(ggplot2)

# Helper functions =============================================================

auc.conf.int.expit <- function(estimate,num.pos,num.neg,conf.level=0.95) {
    ## Calculates confidence interval for an AUC estimate using expit.
    
    ## convert to logit scale
    est.logit = log(estimate/(1-estimate))
    ## standard error (from Kevin Eng)
    se.logit = sqrt(estimate*(1-estimate)/num.pos)*(1/estimate + 1/(1-estimate))
    ## confidence interval in logit
    ci.logit = est.logit+qnorm(c((1-conf.level)/2,(1+conf.level)/2))*se.logit
    
    ## back to original scale
    ci = exp(ci.logit)/(1+exp(ci.logit))
    attr(ci,"conf.level") = conf.level
    attr(ci,"method") = "expit"
    return(ci)
}

# Loading data =================================================================

scores <- fread('./tables/supp_table_4.txt', data.table=F)

# Converting the scores so that 0 is neg. ctrl, 1 is pos. ctrl and scores are 
# directly proportional 

scores$label <- ifelse(scores$label == 'SL', 1, 0)

# ROC AUC when taking all the datasets combined 

roc1 <- pROC::roc(scores$label, (-1)*scores$PRODE_NICE, direction = '<')
roc2 <- pROC::roc(scores$label, (-1)*scores$diff_Gene_Effect, direction = '<')
roc3 <- pROC::roc(scores$label, (-1)*scores$MUFFINN_sum, direction = '<')
roc4 <- pROC::roc(scores$label, (-1)*scores$MUFFINN_max, direction = '<')
roc5 <- pROC::roc(scores$label, (-1)*scores$RWR, direction = '<')
roc6 <- pROC::roc(scores$label, (-1)*scores$NetSig, direction = '<')

pROC::roc.test(roc1, roc2)
pROC::roc.test(roc1, roc3)
pROC::roc.test(roc1, roc4)
pROC::roc.test(roc1, roc5)
pROC::roc.test(roc1, roc6)

# Computing perfs ==============================================================

methods  <- colnames(scores[-c(1:3, ncol(scores))])
datasets <- unique(scores$source) 

perfs <- bind_rows(
    
    o <- lapply(datasets, function(dd){
        
        print(dd)
        
        all_dd <- list()
        
        for (meth in methods){
            
            pred <- ROCR::prediction(
                labels = scores %>% filter(source == dd) %>% pull(label),
                predictions = scores %>% filter(source == dd) %>% pull(meth)*(-1)            
            )
            
            perf1 <- ROCR::performance(pred, 'prec', 'rec')
            perf2 <- ROCR::performance(pred, 'sens', 'spec')
            
            all_dd[[meth]] <- data.frame(
                'rec' =  perf1@x.values[[1]], 
                'pre' =  perf1@y.values[[1]], 
                'sen' =  perf2@x.values[[1]], 
                'spe' =  perf2@y.values[[1]],
                'met' = meth, 
                'dat' = dd
            )
            
        } 
        
        return(bind_rows(all_dd))
    
    })
    
)

aucs <- bind_rows(
    
    o <- lapply(datasets, function(dd){
        
        print(dd)
        
        all_dd <- list()
        
        for (meth in methods){
            
            dt <- ROCR::prediction(
                labels = scores %>% filter(source == dd) %>% pull(label),
                predictions = scores %>% filter(source == dd) %>% pull(meth)*(-1)
            )
            
            perf1 <- ROCR::performance(dt, 'aucpr')
            perf2 <- ROCR::performance(dt, 'auc')
            perf3 <- ROCR::performance(dt, 'prec', 'rec')
            
            conf1 <- auc.conf.int.expit(
                estimate = perf1@y.values[[1]], 
                num.pos  = sum(scores %>% filter(source == dd) %>% pull(label) == 1), 
                num.neg = sum(scores %>% filter(source == dd) %>% pull(label) == 0)
            )
            
            conf2 <- auc.conf.int.expit(
                estimate = perf2@y.values[[1]], 
                num.pos  = sum(scores %>% filter(source == dd) %>% pull(label) == 1), 
                num.neg = sum(scores %>% filter(source == dd) %>% pull(label) == 0)
            )
            
            all_dd[[meth]] <- data.frame(
                'aucpr'     =  perf1@y.values[[1]], 
                'aucpr_low' =  conf1[[1]], 
                'aucpr_high'=  conf1[[2]], 
                'avg_prec'  =  mean(perf3@y.values[[1]], na.rm=T), 
                'auroc'     =  perf2@y.values[[1]], 
                'auroc_low' =  conf2[[1]], 
                'auroc_high'=  conf2[[2]], 
                'met'       = meth, 
                'dat'       = dd
            )
            
        }
        
        return(bind_rows(all_dd))
        
    })
    
)

p1 <- perfs %>%
    mutate(met = factor(met, levels = methods)) %>% 
    ggplot() +
    geom_path(
        aes(
            y=sen,
            x=spe,
            col=met
        )
    ) +
    theme_light() +
    scale_color_manual(
        values = c("#eb9a8a", "#ffd261", "steelblue","#936da3","#a83281","#62c3e3")
    ) +
    scale_x_reverse() +
    #coord_flip() +
    geom_abline(
        linetype = 'dashed',
        color = 'gray',
        slope = 1,
        intercept = 1
    ) +
    theme(legend.position = 'none') +
    labs(y = 'Sensitivity', x = 'Specificity') +
    facet_grid(.~dat)+
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black'))

p2 <- perfs %>%
    mutate(met = factor(met, levels = methods)) %>% 
    ggplot() +
    geom_path(
        aes(
            x=rec,
            y=pre,
            col=met
        )
    ) +
    ylim(0.5, 1) +
    theme_light() +
    scale_color_manual(
        values = c("#eb9a8a", "#ffd261", "steelblue","#936da3","#a83281","#62c3e3")
    )+
    theme(legend.position = 'none') +
    labs(x = 'Recall', y = 'Precision') +
    facet_grid(.~dat) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black'))


p3 <- aucs %>% 
    mutate(met = factor(met, levels = methods)) %>% 
    mutate('PPIs' = rep(c('Gene Eff.', rep('PPIs +\n Gene Eff.', 5)), 3)) %>% 
    ggplot() +
    facet_grid(.~dat) +
    geom_bar(aes(     x = PPIs, 
                      fill = met, 
                      y = auroc,#, 
                      width =rep(c(0.2, 1, 1, 1, 1, 1), 3)), 
             stat="identity", position = position_dodge(preserve = 'total', width = 1), size=0.8, col ='white') +
    geom_errorbar(aes(x=PPIs, fill = met, y=auroc, ymin=auroc_low, ymax=auroc_high, 
                      width =rep(c(0.03, 0.1, 0.1, 0.1, 0.1, 0.1), 3)), 
                  position = position_dodge(preserve = "total", width=1), col="black",
                  stat="identity", size=0.3) +
    geom_text(
        aes(
            x = `PPIs`, 
            fill = met, 
            label = round(auroc, 2),
            y = 0.1,
        ), position = position_dodge(preserve = "total", width=1), 
        angle = 90, 
        col = 'black', 
        size = 2.5
    ) +
    theme_light() +
    theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, 'cm'), legend.position='none') +
    scale_fill_manual(
        values = c("#eb9a8a", "#ffd261", "steelblue","#936da3","#a83281","#62c3e3")
    ) +
    labs(y = 'ROC AUC', x = '')+
    coord_cartesian(ylim = c(0, 1)) +
    facet_grid(.~dat) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black'))

p4 <- aucs %>% 
    mutate(met = factor(met, levels = methods)) %>% 
    mutate('PPIs' = rep(c('Gene Eff.', rep('PPIs +\n Gene Eff.', 5)), 3)) %>% 
    ggplot() +
    facet_grid(.~dat) +
    geom_bar(aes(     x = PPIs, 
                      fill = met, 
                      y = aucpr,#, 
                      width =rep(c(0.2, 1, 1, 1, 1, 1), 3)), 
             stat="identity", position = position_dodge(preserve = 'total', width = 1), size=0.8, col ='white') +
    geom_errorbar(aes(x=PPIs, fill = met, y=aucpr, ymin=aucpr_low, ymax=aucpr_high, 
                      width =rep(c(0.03, 0.1, 0.1, 0.1, 0.1, 0.1), 3)), 
                  position = position_dodge(preserve = "total", width=1), col="black",
                  stat="identity", size=0.3) +
    geom_text(
        aes(
            x = `PPIs`, 
            fill = met, 
            label = round(aucpr, 2),
            y = 0.1,
        ), position = position_dodge(preserve = "total", width=1), 
        angle = 90, 
        col = 'black', 
        size = 2.5
    ) +
    theme_light() +
    theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, 'cm'), legend.position='none') +
    scale_fill_manual(
        values = c("#eb9a8a", "#ffd261", "steelblue","#936da3","#a83281","#62c3e3")
    ) +
    labs(y = 'PR AUC', x = '')+
    coord_cartesian(ylim = c(0, 1)) +
    facet_grid(.~dat) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black'))

(p1+p2)/(p3+p4)




