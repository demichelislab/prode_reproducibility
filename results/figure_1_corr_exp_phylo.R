# Compute correlation between scores, expression CV and gene conservation scores

# Libraries ====================================================================

library(data.table)
library(ggplot2)
library(dplyr)

# Load data ====================================================================

exp_data <- readRDS('./data/data_22Q1_sgRNA.rds')$exp

phylo_dat <- readRDS('./data/phylo_corr.rds')

gn2NIE_sg    <- readRDS('./data/gene2NIE_sgRNA.rds')
gn2AvGe_sg   <- readRDS('./data/gene2AvGE_sgRNA.rds')
gn2MUFS_sg   <- readRDS('./data/gn2MufSum_sgRNA.rds')
gn2MUFM_sg   <- readRDS('./data/gn2MufMax_sgRNA.rds')
gn2RWR_sg    <- readRDS('./data/gn2RWR_sgRNA.rds')

gn2NIE_sh    <- readRDS('./data/gene2NIE_shRNA.rds')
gn2AvGe_sh   <- readRDS('./data/gene2AvGE_shRNA.rds')
gn2MUFS_sh   <- readRDS('./data/gn2MufSum_shRNA.rds')
gn2MUFM_sh   <- readRDS('./data/gn2MufMax_shRNA.rds')
gn2RWR_sh    <- readRDS('./data/gn2RWR_shRNA.rds')

# Compute Expression CV and phylo scores =======================================

gn2cv <- apply(exp_data, 2, function(x){
    var(x, na.rm=T) / mean(x, na.rm=T)
})

gn2ex <- colMeans(exp_data)

gn2ph <- rowMeans(phylo_dat$scores)

# Plot cor with Gene Expression ------------------------------------------------

# all_gns <- Reduce(intersect, list(names(gn2NIE_sg), names(gn2NIE_sh), names(gn2ex)))
# 
# {
#     c1.all_sg <- cor(gn2NIE_sg [all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman') *(-1)
#     c2.all_sg <- cor(gn2AvGe_sg[all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c3.all_sg <- cor(gn2MUFS_sg[all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c4.all_sg <- cor(gn2MUFM_sg[all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c5.all_sg <- cor(gn2RWR_sg [all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     
#     c1.all_sh <- cor(gn2NIE_sh [all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c2.all_sh <- cor(gn2AvGe_sh[all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c3.all_sh <- cor(gn2MUFS_sh[all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c4.all_sh <- cor(gn2MUFM_sh[all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
#     c5.all_sh <- cor(gn2RWR_sh [all_gns], gn2ex[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
# }
# 
# data.frame(
#     'sgRNA input dataset' = c(c1.all_sg, c2.all_sg, c3.all_sg, c4.all_sg, c5.all_sg),
#     'shRNA input dataset' = c(c1.all_sh, c2.all_sh, c3.all_sh, c4.all_sh, c5.all_sh),
#     'method' = c(
#         'PRODE',
#         'Average\nGene Effect',
#         'MUFFINN (Sum)',
#         'MUFFINN (Max)',
#         'RWR'
#         # 'Node Degree'
#     )
# ) %>% 
#     reshape2::melt() %>% 
#     mutate(
#         method = factor(
#             method,
#             levels=c(
#                 'PRODE',
#                 'MUFFINN (Sum)',
#                 'MUFFINN (Max)',
#                 'RWR',
#                 # 'Node Degree',
#                 'Average\nGene Effect'
#             )
#         )
#     ) %>% 
#     mutate(
#         'Use of PPIs' = rep(c('PPIs + \nGene Eff.', 'Gene Eff.','PPIs + \nGene Eff.', 'PPIs + \nGene Eff.', 'PPIs + \nGene Eff.'), 2)
#     ) %>% 
#     ggplot()  +
#     geom_linerange(
#         aes(
#             #xmin = `Use of PPIs`, 
#             x = `Use of PPIs`,
#             fill = method, 
#             ymin = 0,
#             ymax = value
#         ), 
#         #stat='identity', 
#         size = 0.5, 
#         col = 'black',
#         linetype='dashed',
#         position = position_dodge2(preserve = "single", width = 1), 
#         width = 1
#     ) +
#     geom_point(
#         aes(
#             x = `Use of PPIs`, 
#             fill = method, 
#             y = value
#         ), 
#         #stat='identity', 
#         col = 'black', 
#         shape = 21, 
#         size = 3.5, 
#         position = position_dodge2(preserve = "single", width = 1), 
#         linewidth=0.1
#     ) +
#     theme_light() +
#     facet_grid(.~variable, space = 'fixed', scales = 'free', 
#                labeller=as_labeller(list(
#                    'sgRNA input dataset' = 'sgRNA.input.dataset', 
#                    'shRNA input dataset' = 'shRNA.input.dataset' 
#                ))) +
#     ylim(0, .8) +
#     theme(strip.background =element_rect(fill="white"))+
#     theme(strip.text = element_text(colour = 'black')) +
#     ylab('Cor. between Essentiality Scores \n and average Gene Expression') +
#     theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
#     scale_fill_manual(values = c("#ffd261", "steelblue","#936da3","#a83281","#eb9a8a", "lightgray")) +
#     xlab('') +
#     theme(legend.position= 'bottom')

# Plot cor with Gene CV --------------------------------------------------------

all_gns <- Reduce(intersect, list(names(gn2NIE_sg), names(gn2NIE_sh), names(gn2ex)))

{
    c1.all_sg <- cor(gn2NIE_sg [all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c2.all_sg <- cor(gn2AvGe_sg[all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c3.all_sg <- cor(gn2MUFS_sg[all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c4.all_sg <- cor(gn2MUFM_sg[all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c5.all_sg <- cor(gn2RWR_sg [all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    
    c1.all_sh <- cor(gn2NIE_sh [all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c2.all_sh <- cor(gn2AvGe_sh[all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c3.all_sh <- cor(gn2MUFS_sh[all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c4.all_sh <- cor(gn2MUFM_sh[all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
    c5.all_sh <- cor(gn2RWR_sh [all_gns], gn2cv[all_gns], use = 'complete.obs', method = 'spearman')
}

# Reported in table 2 - sgRNA 
print(
    round(c(c2.all_sg, c1.all_sg, c3.all_sg, c4.all_sg, c5.all_sg), 2)
)

# Reported in table 2 - shRNA 
print(
    round(c(c2.all_sh, c1.all_sh, c3.all_sh, c4.all_sh, c5.all_sh), 2)
)


data.frame(
    'sgRNA input dataset' = c(c1.all_sg, c2.all_sg, c3.all_sg, c4.all_sg, c5.all_sg),
    'shRNA input dataset' = c(c1.all_sh, c2.all_sh, c3.all_sh, c4.all_sh, c5.all_sh),
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
    mutate(
        method = factor(
            method,
            levels=c(
                'PRODE',
                'MUFFINN (Sum)',
                'MUFFINN (Max)',
                'RWR',
                # 'Node Degree',
                'Average\nGene Effect'
            )
        )
    ) %>% 
    mutate(
        'Use of PPIs' = rep(c('PPIs + \nGene Eff.', 'Gene Eff.','PPIs + \nGene Eff.', 'PPIs + \nGene Eff.', 'PPIs + \nGene Eff.'), 2)
    ) %>% 
    ggplot()  +
    geom_linerange(
        aes(
            #xmin = `Use of PPIs`, 
            x = `Use of PPIs`,
            fill = method, 
            ymin = 0,
            ymax = value
        ), 
        #stat='identity', 
        size = 0.5, 
        col = 'black',
        linetype='dashed',
        position = position_dodge2(preserve = "single", width = 1), 
        width = 1
    ) +
    geom_point(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            y = value
        ), 
        #stat='identity', 
        col = 'black', 
        shape = 21, 
        size = 3.5, 
        position = position_dodge2(preserve = "single", width = 1), 
        linewidth=0.1
    ) +
    theme_light() +
    facet_grid(.~variable, space = 'fixed', scales = 'free', 
               labeller=as_labeller(list(
                   'sgRNA input dataset' = 'sgRNA.input.dataset', 
                   'shRNA input dataset' = 'shRNA.input.dataset' 
               ))) +
    ylim(0, .8) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black')) +
    ylab('Cor. between Essentiality Scores \n and average Gene Expression CV') +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
    scale_fill_manual(values = c("#ffd261", "steelblue","#936da3","#a83281","#eb9a8a", "lightgray")) +
    xlab('') +
    theme(legend.position= 'bottom')

ggsave('./figures/fig1_corr_expcv_nie.pdf', height = 4, width = 6)

# Plot cor with cons. scores ---------------------------------------------------

all_gns <- Reduce(intersect, list(names(gn2NIE_sg), names(gn2NIE_sh), names(gn2ex)))

{
    c1.all_sg <- cor(gn2NIE_sg [all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c2.all_sg <- cor(gn2AvGe_sg[all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c3.all_sg <- cor(gn2MUFS_sg[all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c4.all_sg <- cor(gn2MUFM_sg[all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c5.all_sg <- cor(gn2RWR_sg [all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    
    c1.all_sh <- cor(gn2NIE_sh [all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c2.all_sh <- cor(gn2AvGe_sh[all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c3.all_sh <- cor(gn2MUFS_sh[all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c4.all_sh <- cor(gn2MUFM_sh[all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
    c5.all_sh <- cor(gn2RWR_sh [all_gns], gn2ph[all_gns], use = 'complete.obs', method = 'spearman')*(-1)
}

# Reported in table 2 - sgRNA 
print(
    round(c(c2.all_sg, c1.all_sg, c3.all_sg, c4.all_sg, c5.all_sg), 2)
)

# Reported in table 2 - shRNA 
print(
    round(c(c2.all_sh, c1.all_sh, c3.all_sh, c4.all_sh, c5.all_sh), 2)
)

data.frame(
    'sgRNA input dataset' = c(c1.all_sg, c2.all_sg, c3.all_sg, c4.all_sg, c5.all_sg),
    'shRNA input dataset' = c(c1.all_sh, c2.all_sh, c3.all_sh, c4.all_sh, c5.all_sh),
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
    mutate(
        method = factor(
            method,
            levels=c(
                'PRODE',
                'MUFFINN (Sum)',
                'MUFFINN (Max)',
                'RWR',
                # 'Node Degree',
                'Average\nGene Effect'
            )
        )
    ) %>% 
    mutate(
        'Use of PPIs' = rep(c('PPIs + \nGene Eff.', 'Gene Eff.','PPIs + \nGene Eff.', 'PPIs + \nGene Eff.', 'PPIs + \nGene Eff.'), 2)
    ) %>% 
    ggplot()  +
    geom_linerange(
        aes(
            #xmin = `Use of PPIs`, 
            x = `Use of PPIs`,
            fill = method, 
            ymin = 0,
            ymax = value
        ), 
        #stat='identity', 
        size = 0.5, 
        col = 'black',
        linetype='dashed',
        position = position_dodge2(preserve = "single", width = 1), 
        width = 1
    ) +
    geom_point(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            y = value
        ), 
        #stat='identity', 
        col = 'black', 
        shape = 21, 
        size = 3.5, 
        position = position_dodge2(preserve = "single", width = 1), 
        linewidth=0.1
    ) +
    theme_light() +
    facet_grid(.~variable, space = 'fixed', scales = 'free', 
               labeller=as_labeller(list(
                   'sgRNA input dataset' = 'sgRNA.input.dataset', 
                   'shRNA input dataset' = 'shRNA.input.dataset' 
               ))) +
    ylim(0, .8) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black')) +
    ylab('Cor. between Essentiality Scores \n and Phylogenetic Conservation Scores') +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
    scale_fill_manual(values = c("#ffd261", "steelblue","#936da3","#a83281","#eb9a8a", "lightgray")) +
    xlab('') +
    theme(legend.position= 'bottom')

ggsave('./figures/fig1_corr_phylo_nie.pdf', height = 4, width = 6)

# Plot correlation shRNA and sgRNA ---------------------------------------------

all_gns <- Reduce(intersect, list(names(gn2NIE_sg), names(gn2NIE_sh), names(gn2ex)))

{
    c1 <- cor(gn2NIE_sg [all_gns], gn2NIE_sh [all_gns], use = 'complete.obs', method = 'spearman')
    c2 <- cor(gn2AvGe_sg[all_gns], gn2AvGe_sh[all_gns], use = 'complete.obs', method = 'spearman')
    c3 <- cor(gn2MUFS_sg[all_gns], gn2MUFS_sh[all_gns], use = 'complete.obs', method = 'spearman')
    c4 <- cor(gn2MUFM_sg[all_gns], gn2MUFM_sh[all_gns], use = 'complete.obs', method = 'spearman')
    c5 <- cor(gn2RWR_sg [all_gns], gn2RWR_sh [all_gns], use = 'complete.obs', method = 'spearman')
}

data.frame(
    'sgRNA input dataset' = c(c1, c2, c3, c4, c5),
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
    mutate(
        method = factor(
            method,
            levels=c(
                'PRODE',
                'MUFFINN (Sum)',
                'MUFFINN (Max)',
                'RWR',
                # 'Node Degree',
                'Average\nGene Effect'
            )
        )
    ) %>% 
    mutate(variable = gsub('[.]', ' ', variable)) %>% 
    mutate(
        'Use of PPIs' = rep(c('PPIs + \nGene Eff.', 'Gene Eff.','PPIs + \nGene Eff.', 'PPIs + \nGene Eff.', 'PPIs + \nGene Eff.'))
    ) %>% 
    ggplot()  +
    geom_linerange(
        aes(
            #xmin = `Use of PPIs`, 
            x = `Use of PPIs`,
            fill = method, 
            ymin = 0,
            ymax = value
        ), 
        #stat='identity', 
        size = 0.5, 
        col = 'black',
        linetype='dashed',
        position = position_dodge2(preserve = "single", width = 1), 
        width = 1
    ) +
    geom_point(
        aes(
            x = `Use of PPIs`, 
            fill = method, 
            y = value
        ), 
        #stat='identity', 
        col = 'black', 
        shape = 21, 
        size = 3.5, 
        position = position_dodge2(preserve = "single", width = 1), 
        linewidth=0.1
    ) +
    theme_light() +
    facet_grid(.~variable, space = 'fixed', scales = 'free', 
               labeller=as_labeller(list(
                   'sgRNA input dataset' = 'sgRNA.input.dataset', 
                   'shRNA input dataset' = 'shRNA.input.dataset' 
               ))) +
    ylim(0, .8) +
    theme(strip.background =element_rect(fill="white"))+
    theme(strip.text = element_text(colour = 'black')) +
    ylab('Cor. between Essentiality Scores \n and average Gene Expression') +
    theme(axis.text.x = element_text(angle = 60,  hjust=1)) +
    scale_fill_manual(values = c("#ffd261", "steelblue","#936da3","#a83281","#eb9a8a", "lightgray")) +
    xlab('') +
    theme(legend.position= 'bottom')


