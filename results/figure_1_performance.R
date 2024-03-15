# Computing ROC AUC and PR AUC discriminating Essential vs. Not Essential ......

# Import libraries -------------------------------------------------------------

library(ggplot2)
library(data.table)
library(dplyr)
library(ROCR)

# Helper functions -------------------------------------------------------------

auc.conf.int.expit <- function(estimate, num.pos, num.neg, conf.level=0.95) {
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

# Load datasets ----------------------------------------------------------------

val_dt <- fread('./paper_tables/supp_table_1.txt', data.table=F)

val_dt$label <- ifelse(val_dt$label == 'essential', 1, 0)

val_dt_sg <- val_dt %>% filter(!is.na(avg_Gene_Effect_sgRNA))
val_dt_sh <- val_dt %>% filter(!is.na(avg_Gene_Effect_shRNA))

# Compute ROC AUC --------------------------------------------------------------

## `direction = '>'` --> Controls (not-essential genes) have higher values

## sgRNA (knock-out) dataset ...................................................
sg_roc0 <- pROC::roc(val_dt_sg$label, val_dt_sg$avg_Gene_Effect_sgRNA, direction = '>')
sg_roc1 <- pROC::roc(val_dt_sg$label, val_dt_sg$PRODE_NIE_sgRNA, direction = '>')
sg_roc2 <- pROC::roc(val_dt_sg$label, val_dt_sg$MUFFINN_sum_sgRNA, direction = '>')
sg_roc3 <- pROC::roc(val_dt_sg$label, val_dt_sg$MUFFINN_max_sgRNA, direction = '>')
sg_roc4 <- pROC::roc(val_dt_sg$label, val_dt_sg$RWR_sgRNA, direction = '>')

sg_bar_roc0 <- pROC::ci.auc(val_dt_sg$label, val_dt_sg$avg_Gene_Effect_sgRNA, direction = '>')
sg_bar_roc1 <- pROC::ci.auc(val_dt_sg$label, val_dt_sg$PRODE_NIE_sgRNA, direction = '>')
sg_bar_roc2 <- pROC::ci.auc(val_dt_sg$label, val_dt_sg$MUFFINN_sum_sgRNA, direction = '>')
sg_bar_roc3 <- pROC::ci.auc(val_dt_sg$label, val_dt_sg$MUFFINN_max_sgRNA, direction = '>')
sg_bar_roc4 <- pROC::ci.auc(val_dt_sg$label, val_dt_sg$RWR_sgRNA, direction = '>')

## shRNA (knock-down) dataset ..................................................
sh_roc0 <- pROC::roc(val_dt_sh$label, val_dt_sh$avg_Gene_Effect_shRNA, direction = '>')
sh_roc1 <- pROC::roc(val_dt_sh$label, val_dt_sh$PRODE_NIE_shRNA, direction = '>')
sh_roc2 <- pROC::roc(val_dt_sh$label, val_dt_sh$MUFFINN_sum_shRNA, direction = '>')
sh_roc3 <- pROC::roc(val_dt_sh$label, val_dt_sh$MUFFINN_max_shRNA, direction = '>')
sh_roc4 <- pROC::roc(val_dt_sh$label, val_dt_sh$RWR_shRNA, direction = '>')

sh_bar_roc0 <- pROC::ci.auc(val_dt_sh$label, val_dt_sh$avg_Gene_Effect_shRNA, direction = '>')
sh_bar_roc1 <- pROC::ci.auc(val_dt_sh$label, val_dt_sh$PRODE_NIE_shRNA, direction = '>')
sh_bar_roc2 <- pROC::ci.auc(val_dt_sh$label, val_dt_sh$MUFFINN_sum_shRNA, direction = '>')
sh_bar_roc3 <- pROC::ci.auc(val_dt_sh$label, val_dt_sh$MUFFINN_max_shRNA, direction = '>')
sh_bar_roc4 <- pROC::ci.auc(val_dt_sh$label, val_dt_sh$RWR_shRNA, direction = '>')

# Compute PR AUC ---------------------------------------------------------------

# ROCR::prediction() expects the lower level corresponding to the negative class (0), 
# the upper level to the positive class (1), so changing direction to scores 

val_dt_sg <- val_dt %>% filter(!is.na(avg_Gene_Effect_sgRNA))
val_dt_sh <- val_dt %>% filter(!is.na(avg_Gene_Effect_shRNA))

## sgRNA (knock-out) dataset ...................................................
sg_pred0 <- ROCR::prediction(val_dt_sg$avg_Gene_Effect_sgRNA*(-1), val_dt_sg$label)
sg_pred1 <- ROCR::prediction(val_dt_sg$PRODE_NIE_sgRNA*(-1), val_dt_sg$label)
sg_pred2 <- ROCR::prediction(val_dt_sg$MUFFINN_sum_sgRNA*(-1), val_dt_sg$label)
sg_pred3 <- ROCR::prediction(val_dt_sg$MUFFINN_max_sgRNA*(-1),val_dt_sg$label)
sg_pred4 <- ROCR::prediction(val_dt_sg$RWR_sgRNA*(-1),val_dt_sg$label)

sg_rocPr0 <- ROCR::performance(sg_pred0, 'aucpr')@y.values[[1]]
sg_rocPr1 <- ROCR::performance(sg_pred1, 'aucpr')@y.values[[1]]
sg_rocPr2 <- ROCR::performance(sg_pred2, 'aucpr')@y.values[[1]]
sg_rocPr3 <- ROCR::performance(sg_pred3, 'aucpr')@y.values[[1]]
sg_rocPr4 <- ROCR::performance(sg_pred4, 'aucpr')@y.values[[1]]

sg_bar_rocPr0_low <- auc.conf.int.expit(sg_rocPr0, sum(val_dt_sg$label), sum(!val_dt_sg$label))[1]
sg_bar_rocPr1_low <- auc.conf.int.expit(sg_rocPr1, sum(val_dt_sg$label), sum(!val_dt_sg$label))[1]
sg_bar_rocPr2_low <- auc.conf.int.expit(sg_rocPr2, sum(val_dt_sg$label), sum(!val_dt_sg$label))[1]
sg_bar_rocPr3_low <- auc.conf.int.expit(sg_rocPr3, sum(val_dt_sg$label), sum(!val_dt_sg$label))[1]
sg_bar_rocPr4_low <- auc.conf.int.expit(sg_rocPr4, sum(val_dt_sg$label), sum(!val_dt_sg$label))[1]

sg_bar_rocPr0_hig <- auc.conf.int.expit(sg_rocPr0, sum(val_dt_sg$label), sum(!val_dt_sg$label))[2]
sg_bar_rocPr1_hig <- auc.conf.int.expit(sg_rocPr1, sum(val_dt_sg$label), sum(!val_dt_sg$label))[2]
sg_bar_rocPr2_hig <- auc.conf.int.expit(sg_rocPr2, sum(val_dt_sg$label), sum(!val_dt_sg$label))[2]
sg_bar_rocPr3_hig <- auc.conf.int.expit(sg_rocPr3, sum(val_dt_sg$label), sum(!val_dt_sg$label))[2]
sg_bar_rocPr4_hig <- auc.conf.int.expit(sg_rocPr4, sum(val_dt_sg$label), sum(!val_dt_sg$label))[2]

## shRNA (knock-down) dataset ..................................................
sh_pred0 <- ROCR::prediction(val_dt_sh$avg_Gene_Effect_shRNA*(-1), val_dt_sh$label)
sh_pred1 <- ROCR::prediction(val_dt_sh$PRODE_NIE_shRNA*(-1), val_dt_sh$label)
sh_pred2 <- ROCR::prediction(val_dt_sh$MUFFINN_sum_shRNA*(-1), val_dt_sh$label)
sh_pred3 <- ROCR::prediction(val_dt_sh$MUFFINN_max_shRNA*(-1),val_dt_sh$label)
sh_pred4 <- ROCR::prediction(val_dt_sh$RWR_shRNA*(-1),val_dt_sh$label)

sh_rocPr0 <- ROCR::performance(sh_pred0, 'aucpr')@y.values[[1]]
sh_rocPr1 <- ROCR::performance(sh_pred1, 'aucpr')@y.values[[1]]
sh_rocPr2 <- ROCR::performance(sh_pred2, 'aucpr')@y.values[[1]]
sh_rocPr3 <- ROCR::performance(sh_pred3, 'aucpr')@y.values[[1]]
sh_rocPr4 <- ROCR::performance(sh_pred4, 'aucpr')@y.values[[1]]

sh_bar_rocPr0_low <- auc.conf.int.expit(sh_rocPr0, sum(val_dt_sh$label), sum(!val_dt_sh$label))[1]
sh_bar_rocPr1_low <- auc.conf.int.expit(sh_rocPr1, sum(val_dt_sh$label), sum(!val_dt_sh$label))[1]
sh_bar_rocPr2_low <- auc.conf.int.expit(sh_rocPr2, sum(val_dt_sh$label), sum(!val_dt_sh$label))[1]
sh_bar_rocPr3_low <- auc.conf.int.expit(sh_rocPr3, sum(val_dt_sh$label), sum(!val_dt_sh$label))[1]
sh_bar_rocPr4_low <- auc.conf.int.expit(sh_rocPr4, sum(val_dt_sh$label), sum(!val_dt_sh$label))[1]

sh_bar_rocPr0_hig <- auc.conf.int.expit(sh_rocPr0, sum(val_dt_sh$label), sum(!val_dt_sh$label))[2]
sh_bar_rocPr1_hig <- auc.conf.int.expit(sh_rocPr1, sum(val_dt_sh$label), sum(!val_dt_sh$label))[2]
sh_bar_rocPr2_hig <- auc.conf.int.expit(sh_rocPr2, sum(val_dt_sh$label), sum(!val_dt_sh$label))[2]
sh_bar_rocPr3_hig <- auc.conf.int.expit(sh_rocPr3, sum(val_dt_sh$label), sum(!val_dt_sh$label))[2]
sh_bar_rocPr4_hig <- auc.conf.int.expit(sh_rocPr4, sum(val_dt_sh$label), sum(!val_dt_sh$label))[2]

# Statistical tests on ROC curves ----------------------------------------------

## Of note: PRODE roc curve is `sg_roc0`
pROC::roc.test(sg_roc0, sg_roc1, alternative = 'less') 
pROC::roc.test(sg_roc2, sg_roc1, alternative = 'less')
pROC::roc.test(sg_roc3, sg_roc1, alternative = 'less')
pROC::roc.test(sg_roc4, sg_roc1, alternative = 'less')

pROC::roc.test(sh_roc0, sh_roc1, alternative = 'less') 
pROC::roc.test(sh_roc2, sh_roc1, alternative = 'less')
pROC::roc.test(sh_roc3, sh_roc1, alternative = 'less')
pROC::roc.test(sh_roc4, sh_roc1, alternative = 'less')
# Plot ==================h======================================================

# Plot ROC AUC -----------------------------------------------------------------

data.frame(
    'sgRNA input dataset' = c(
        sg_roc1$auc,
        sg_roc0$auc,
        sg_roc2$auc,
        sg_roc3$auc,
        sg_roc4$auc
    ),
    'shRNA input dataset' = c(
        sh_roc1$auc,
        sh_roc0$auc,
        sh_roc2$auc,
        sh_roc3$auc,
        sh_roc4$auc
    ),
    'method' = c(
        'PRODE',
        'Average\nGene Effect',
        'MUFFINN (Sum)',
        'MUFFINN (Max)',
        'RWR'
    )
) %>% 
    reshape2::melt() %>% 
    mutate(variable = gsub('[.]', ' ', variable)) %>% 
    mutate(
        method = factor(
            method,
            levels=c(
                'PRODE',
                'MUFFINN (Sum)',
                'MUFFINN (Max)',
                'RWR',
                'Average\nGene Effect'
            )
        )
    ) %>% 
    mutate(
        'Use of PPIs' = rep(c('PPIs + \nGene Eff.', 'Gene Eff.','PPIs + \nGene Eff.', 'PPIs + \nGene Eff.', 'PPIs + \nGene Eff.'), 2)
    ) %>%  
    mutate(
        CI_low = c(
            c(sg_bar_roc1)[1],
            c(sg_bar_roc0)[1],
            c(sg_bar_roc2)[1],
            c(sg_bar_roc3)[1],
            c(sg_bar_roc4)[1],
            c(sh_bar_roc1)[1],
            c(sh_bar_roc0)[1],
            c(sh_bar_roc2)[1],
            c(sh_bar_roc3)[1],
            c(sh_bar_roc4)[1]
        ), 
        
        CI_high = c(
            c(sg_bar_roc1)[3],
            c(sg_bar_roc0)[3],
            c(sg_bar_roc2)[3],
            c(sg_bar_roc3)[3],
            c(sg_bar_roc4)[3],
            c(sh_bar_roc1)[3],
            c(sh_bar_roc0)[3],
            c(sh_bar_roc2)[3],
            c(sh_bar_roc3)[3],
            c(sh_bar_roc4)[3]
        )
    ) %>% 
    ggplot() +
    geom_bar(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            y = value
        ), 
        stat='identity', 
        col = 'white', 
        size = 0.8, 
        position = position_dodge2(preserve = "single", width=1), 
        #width = 1
    ) +
    geom_errorbar(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            y = value, 
            ymin = CI_low, 
            ymax = CI_high, 
            width = rep(c(0.1, 0.03, 0.1, 0.1, 0.1), 2)
        ), position = position_dodge(preserve = "total", width=0.9), 
        #width=0.3,
        size=0.3 
    ) +
    theme_light() +
    facet_grid(.~variable, space = 'fixed', scales = 'free', 
               labeller=as_labeller(list(
                   'sgRNA input dataset' = 'sgRNA input dataset', 
                   'shRNA input dataset' = 'shRNA input dataset' 
               ))) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black')) +
    ylab('ROC AUC \n (Essential vs. Non essential genes)') +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
    scale_fill_manual(values = c("#ffd261", "steelblue","#936da3","#a83281","#eb9a8a", "lightgray")) +
    xlab('') +
    coord_cartesian(ylim=c(0.8, 1)) +
    scale_x_discrete(expand=c(0.1,0.3,0.3,0.1)) +
    geom_text(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            label = round(value, 3),
            y = 0.81,
        ), position = position_dodge(preserve = "total", width=0.9), 
        angle = 90, 
        col = 'black', 
        size = 3
    ) +
    theme(legend.position = 'bottom') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    labs(fill='')

# Plot PR AUC ------------------------------------------------------------------

data.frame(
    'sgRNA input dataset' = c(
        sg_rocPr1,
        sg_rocPr0,
        sg_rocPr2,
        sg_rocPr3,
        sg_rocPr4
        
    ),
    'shRNA input dataset' = c(
        sh_rocPr1,
        sh_rocPr0,
        sh_rocPr2,
        sh_rocPr3,
        sh_rocPr4
    ),
    'method' = c(
        'PRODE',
        'Average\nGene Effect',
        'MUFFINN (Sum)',
        'MUFFINN (Max)',
        'RWR'
        # 'Node Degree'
    )
) %>% 
    reshape2::melt() %>% 
    mutate(variable = gsub('[.]', ' ', variable)) %>% 
    mutate(
        method = factor(
            method,
            levels=c(
                'PRODE',
                'MUFFINN (Sum)',
                'MUFFINN (Max)',
                'RWR',
                'Average\nGene Effect'
            )
        )
    ) %>% 
    mutate(
        'Use of PPIs' = rep(c('PPIs + \nGene Eff.', 'Gene Eff.','PPIs + \nGene Eff.', 'PPIs + \nGene Eff.', 'PPIs + \nGene Eff.'), 2)
    ) %>%  
    mutate(
        CI_low = c(
            c(sg_bar_rocPr1_low),
            c(sg_bar_rocPr0_low),
            c(sg_bar_rocPr2_low),
            c(sg_bar_rocPr3_low),
            c(sg_bar_rocPr4_low),
            c(sh_bar_rocPr1_low),
            c(sh_bar_rocPr0_low),
            c(sh_bar_rocPr2_low),
            c(sh_bar_rocPr3_low),
            c(sh_bar_rocPr4_low)
        ), 
        
        CI_high = c(
            c(sg_bar_rocPr1_hig),
            c(sg_bar_rocPr0_hig),
            c(sg_bar_rocPr2_hig),
            c(sg_bar_rocPr3_hig),
            c(sg_bar_rocPr4_hig),
            c(sh_bar_rocPr1_hig),
            c(sh_bar_rocPr0_hig),
            c(sh_bar_rocPr2_hig),
            c(sh_bar_rocPr3_hig),
            c(sh_bar_rocPr4_hig)
        )
    ) %>% 
    ggplot() +
    geom_bar(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            y = value
        ), 
        stat='identity', 
        col = 'white', 
        size = 0.8, 
        position = position_dodge2(preserve = "single", width=1), 
        #width = 1
    ) +
    geom_errorbar(
        aes(
            x = `Use of PPIs`,
            fill = method,
            y = value,
            ymin = CI_low,
            ymax = CI_high,
            width = rep(c(0.1, 0.03, 0.1, 0.1, 0.1), 2)
        ), position = position_dodge(preserve = "total", width=0.9),
        # width=0.3,
        size=0.3
    ) +
    theme_light() +
    facet_grid(.~variable, space = 'fixed', scales = 'free', 
               labeller=as_labeller(list(
                   'sgRNA input dataset' = 'sgRNA input dataset', 
                   'shRNA input dataset' = 'shRNA input dataset' 
               ))) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black')) +
    ylab('PR AUC \n (Essential vs. Non essential genes)') +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
    scale_fill_manual(values = c("#ffd261", "steelblue","#936da3","#a83281","#eb9a8a", "lightgray")) +
    xlab('') +
    coord_cartesian(ylim=c(0.8, 1)) +
    scale_x_discrete(expand=c(0.1,0.3,0.3,0.1)) +
    geom_text(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            label = round(value, 3),
            y = 0.81,
        ), position = position_dodge(preserve = "total", width=0.9), 
        angle = 90, 
        col = 'black', 
        size = 3
    ) +
    theme(legend.position = 'bottom') +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    labs(fill='') 

