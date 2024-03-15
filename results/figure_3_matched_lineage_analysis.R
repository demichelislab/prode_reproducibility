#-------------------------------------------------------------------------------
#                     Figure 3 - lineage-matched analyses
#_______________________________________________________________________________

library(data.table)
library(patchwork)
library(dplyr)
library(ggplot2)

# Load data ====================================================================

scores <- fread('./data/lineage_matched_iso_data.txt', 
                data.table=F)


scores %>% 
    group_by(GeneA, lineage) %>% 
    summarise(
        ccNIE = cor(
            PRODE_NIE, Experiment_zscore, use='complete.obs', 
            method = 'spearman'
        ), 
        ccDGE = cor(
            Diff_Gene_Effect, Experiment_zscore, use='complete.obs', 
            method = 'spearman'
        )
    )


# Check if top scores in top PRODE / NIE ---------------------------------------

out <- list()

for (th in c(0.01, 0.05, 0.1, 0.25)){
    
    scores %>% 
        group_by(GeneA, lineage) %>% 
        mutate(
            q1 = rank(PRODE_NIE, na.last = 'keep')/sum(!is.na(PRODE_NIE)), 
            q2 = rank(Diff_Gene_Effect, na.last = 'keep')/sum(!is.na(Diff_Gene_Effect)), 
            q3 = rank(Experiment_zscore, na.last = 'keep') 
        ) %>% 
        summarise(
            top_100_NIE = sum(q3 <= 100 & q1 <= th, na.rm=T)/100, 
            top_100_DGE = sum(q3 <= 100 & q2 <= th, na.rm=T)/100
        ) %>% 
        mutate(
            thr = th
        ) -> out[[as.character(th)]]
    
}

out <- bind_rows(out)

# Plot proportion of nominated top hits ========================================

## statistical tests ...........................................................

toAst <- function(pv){
    
    ast <- rep('', length(pv))
    ast[which(pv <= 0.1)] <- '*'
    ast[which(pv <= 0.05)] <- '**'
    ast[which(pv <= 0.01)] <- '***'
    return(ast)
    
}

pvals <- vapply(c(.01, .05, .1, .25), function(th){
    out %>% filter(thr == th) -> tmp 
    wilcox.test(tmp$top_100_DGE, tmp$top_100_NIE, paired = T)$p.value 
}, .1) %>% toAst()


## Plot data ...................................................................

out %>% 
    group_by(thr) %>% 
    summarise(
        avg_NIE = mean(top_100_NIE), 
        avg_DGE = mean(top_100_DGE), 
        sd_NIE  = sd(top_100_NIE)/sqrt(length(top_100_NIE)), 
        sd_DGE  = sd(top_100_DGE)/sqrt(length(top_100_NIE))
    ) -> tmp 

data.frame(
    thr = rep(tmp$thr, 2), 
    avg = c(tmp$avg_NIE, tmp$avg_DGE), 
    sds = c(tmp$sd_NIE, tmp$sd_DGE), 
    met = rep(
        c('PRODE', 'Diff. Gene Eff'), 
        each = 4
    )
) %>% 
    mutate(thr = factor(thr)) %>% 
    ggplot() +
        geom_line(
            aes(
                x = thr, 
                y = avg, 
                group = met
            )
        ) +
    geom_errorbar(
        aes(
            x = thr, 
            ymin = avg - sds, 
            ymax = avg + sds, 
            y = avg, 
            group = met
        ), width = 0.1
    ) +
    geom_point(
        aes(
            x = thr, 
            y = avg, 
            fill = met
        ), size = 3, 
        shape = 21
    ) + 
    annotate(
        'text', 
        x = 1:4, 
        y = tmp$avg_NIE + 0.1, 
        label = pvals
    ) + 
    ylim(0, 0.6) +
    theme_light() +
    scale_fill_manual(
        values = c("#eb9a8a", "#ffd261")
    )  +
    labs(fill = 'Method', 
         x = 'Decision threshold on \n public data analysis (top %)', 
         y = 'Mean fraction of hits confirmed \n by public data analysis'
    )

# Plot AUCs top / bottom hits ==================================================

scores %>% 
    group_by(GeneA, lineage) %>% 
    mutate(keep = 
        rank(+Experiment_zscore, na.last = 'keep') <= 100 | 
        rank(-Experiment_zscore, na.last = 'keep') <= 100
    ) %>% 
    mutate(
        labl = ifelse(Experiment_zscore < 0, 1, 0)
    ) %>% 
    filter(keep) %>% 
    summarise(
        NIE_AUC = pROC::auc(labl, PRODE_NIE, direction = '>')[[1]], 
        DGE_AUC = pROC::auc(labl, Diff_Gene_Effect, direction = '>')[[1]]
    ) -> auc_perfs

auc_perfs %>% 
    reshape2::melt() %>% 
    mutate(
        variable = factor(variable, levels = c(
            'DGE_AUC', 'NIE_AUC'
        ))
    ) %>% 
    ggplot() +
    geom_boxplot(
        aes(
            x = variable, 
            y = value, 
        ), width = 0.3
    ) +
    geom_point(
        aes(
            x = variable, 
            y = value, 
            fill = variable
        ), shape = 21
    ) +
    geom_line(
        aes(
            x = variable, 
            y = value, 
            group = paste0(GeneA, '-', lineage)
        ), 
        col = 'lightgray', 
        linetype = 'dashed'
    ) +
    scale_fill_manual(
        values = c("#eb9a8a", "#ffd261")
    )  +
    theme_light() +
    geom_hline(
        yintercept = 0.5, 
        linetype = 'dashed'
    ) +
    labs(y = 'AUC\n(Lineage-Matched analysis)', 
         x = '') +
    theme(legend.position = 'none') +
    ylim(0.2, 1) +
    ggpubr::stat_compare_means(
        aes(
            x = variable, 
            y = value
        ), label = 'p.signif', 
        paired = T
    )


































