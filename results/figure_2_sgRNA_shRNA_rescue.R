# ==============================================================================
#                               Figure 2 code 
# ==============================================================================
#
# ProDe rescues sgRNA essential genes undetected by shRNA dependency effects. 
#
# Libraries ====================================================================

library(data.table)
library(ggplot2)
library(prodeTool)
library(dplyr)

# Load data ====================================================================

pr_output_sg <- readRDS('./data/HT1197_prode_output_sgRNA.rds')
pr_output_sh <- readRDS('./data/HT1197_prode_output_shRNA.rds')

val_dt <- fread('./tables/supp_table_1.txt', data.table=F)

# Check the intersection of ess. genes =========================================

## Thresholds obtained from max. sens / spec (supp. analyses)
ess_pro_sg  <- pr_output_sg %>% as.data.frame() %>% filter(NIE_score <= -4.14) %>% pull(gene)
ess_pro_sh  <- pr_output_sh %>% as.data.frame() %>% filter(NIE_score <= -3.48) %>% pull(gene)
ess_avge_sg <- pr_output_sg %>% as.data.frame() %>% filter(Estimate  <= -0.43) %>% pull(gene)
ess_avge_sh <- pr_output_sh %>% as.data.frame() %>% filter(Estimate  <= -0.22) %>% pull(gene)

## Consider only genes that are tested in both shRNA and sgRNA datasets
cmn <- intersect(pr_output_sg$gene, pr_output_sh$gene) 

pos_gns <- val_dt %>% filter(label == 'essential') %>% pull(gene)
neg_gns <- val_dt %>% filter(label != 'essential') %>% pull(gene)

pos_ctrls <- pos_gns[which(pos_gns %in% cmn)]
neg_ctrls <- neg_gns[which(neg_gns %in% cmn)]

## Retrieve private essential (sgRNA or shRNA)
priv_sg <- ess_avge_sg[which(ess_avge_sg %in% pos_ctrls & !ess_avge_sg %in% ess_avge_sh)]
priv_sh <- ess_avge_sh[which(ess_avge_sh %in% pos_ctrls & !ess_avge_sh %in% ess_avge_sg)]
cmn_sgsh <-  intersect(ess_avge_sg[which(ess_avge_sg %in% pos_ctrls)], ess_avge_sh[which(ess_avge_sh %in% pos_ctrls)])

# Found by PRODE - shRNA 
mean(priv_sg %in% ess_pro_sh)
sum(priv_sg %in% ess_pro_sh)

library(VennDiagram)

## Create Venn diagram
venn_data <- list(
    " " = pos_ctrls[which(pos_ctrls %in% ess_avge_sh)],
    " " = pos_ctrls[which(pos_ctrls %in% ess_avge_sg)]
)

# Run only if you want to produce the venn diagram in fig. 2 \\\\\\\\\\\\\\\\\\\
# venn.diagram(
#     x = venn_data,
# 
#     filename = './figures/fig2_shRNA_sgRNA_venn_diagram.png',
#     output=T, 
# 
#     # Output features
#     imagetype="png" ,
#     height = 380 ,
#     width = 380 ,
#     resolution = 300,
#     compression = "lzw",
# 
#     # Circles
#     lwd = 0.5,
#     lty = c(1,2),
#     fill = c('#3285a8', '#32a891'),
#     alpha = .8,
#     # 
#     # # Numbers
#     cex = .6,
#     #fontface = "bold",
#     fontfamily = "Arial", # may be not installed on your device 
# 
#     # Set names
#     cat.cex = 0.5,
#     cat.fontface = "bold",
#     cat.default.pos = "outer",
#     cat.pos = c(-27, 27),
#     cat.dist = c(0.055, 0.055),
#     cat.fontfamily = "Arial",
#     margin = .05,
#     disable.logging = T
# )
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# Violin plots -----------------------------------------------------------------

# Plot the PRODE NIE shRNA scores for genes found only by sgRNA 
pr_output_sh %>% 
    as.data.frame() %>% 
    # mutate(
    #   prode_score = rank(prode_score)  
    # ) %>% 
    filter(
        gene %in% c(priv_sg, neg_ctrls)
    ) %>% 
    mutate(
        'ctrl' = ifelse(gene %in% priv_sg, 'Private sgRNA\n essentials', 'Non-expressed\n genes')
    ) %>% 
    ggplot(
        aes(
            x = ctrl, 
            y = NIE_score
        )
    ) +
    geom_violin(
        aes(
            x = ctrl, 
            y = NIE_score, 
            fill = ctrl
        )
    ) +
    geom_boxplot(
        aes(
            x = ctrl, 
            y = NIE_score
        ), width=.1, 
        outlier.size = 0
    ) +
    geom_hline(
        yintercept = -3.48, 
        linetype = 'dotted'
    ) +
    theme_light() +
    ggpubr::stat_compare_means(
        comparisons = list(
            c('Non-expressed\n genes', 'Private sgRNA\n essentials')
        ), label = 'p.signif'
    ) +
    theme(legend.position='none', 
          axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(
        values = c(
            'lightgray', '#32a891'
        )
    ) +
    labs(
        x = '', 
        y = 'PRODE NIE score\n HT1197 - (shRNA)'
    ) +
    annotate('text', x = 2.5, y = -2, label='Not Essential', angle=90, size = 3) + 
    annotate('text', x = 2.5, y = -6, label='Essential', angle=90, size = 3)+
    ylim(-7.1, 0.5)

ggsave('./figures/fig2_violin_plot_private_sgRNA.pdf', 
       height = 4, width = 2)

# Scatterplot ------------------------------------------------------------------

high_gns <- c(
    'DIS3', 'RIOK1'
)

low_gns <- c(
    'CINP'
)

dt <- pr_output_sh %>% 
    as.data.frame() %>% 
    mutate(
        r_y = rank(NIE_score), 
        r_x = rank(Estimate)
    ) %>% 
    filter(gene %in% c(priv_sg, neg_ctrls))

ys <- sum(pr_output_sh$NIE_score <= -3.48)
xs <- sum(pr_output_sh$Estimate <= -.22)

dt %>% 
    filter(gene%in%priv_sg) %>% 
    ggplot() +
    geom_hline(
        yintercept = ys,
        linetype='dotted'
    ) +
    geom_vline(
        xintercept = xs,
        linetype='dotted'
    ) +
    geom_abline(
        linetype='dashed'
    ) +
    geom_point(
        aes(
            x = r_x,
            y = r_y,
            col = factor(
                ifelse(gene %in% priv_sg, 'Private sgRNA essential', 'Non-expressed genes'), 
                levels=c('Private sgRNA essential', 'Non-expressed genes')
            )
        ), data = dt  %>% filter(gene%in%priv_sg) 
    ) +
    geom_point(
        aes(
            x = r_x,
            y = r_y,
        ), data = dt %>% filter(gene %in% c(high_gns, low_gns)),
        shape = 23,
        fill = '#32a891',#'#ff5c33',
        size = 4
    ) +
    ggrepel::geom_label_repel(
        aes(
            x = r_x,
            y = r_y,
            label = as.character(gene)
        ),
        data = dt %>% filter(gene %in% c(high_gns, low_gns)),
        box.padding = 1,
        color = 'black'#'#ff5c33'
    ) +
    labs(color = '') + 
    theme_light() + 
    theme(legend.position ='bottom') +
    scale_color_manual(
        values = c(
            '#32a891', 'lightgray'
        )
    ) +
    labs(
        x = 'Ranked Avg. Gene Eff. \nHT1197 - (shRNA)', 
        y = 'Ranked PRODE NIE score \nHT1197 - (shRNA)'
    ) +
    ylim(0, 13500) +
    xlim(0, 13500) 

ggsave('./figures/fig2_ranked_dependency_effects_scatterplot.pdf', 
       height = 3.7, width = 3.5)


# Number of connections --------------------------------------------------------

adj_m <- pr_input_sh@adjMatrix # Prode object stores the adjacency matrix 
adj_m <- as.matrix(adj_m) # from Matrix to matrix (R base)

ess_m <- adj_m[pos_ctrls, c('CINP', 'DIS3', 'RIOK1')]
ess_s <- colSums(ess_m) # number of connections per gene 

bk_d <- colSums(adj_m[pos_ctrls, -which(colnames(adj_m)%in%pos_ctrls)]) # every gene n. connection with positive controls (essential genes)

ggplot() + 
    geom_density(
        aes(
            x = bk_d + 1
        )
    ) + 
    geom_segment(
        aes(
            x = ess_s,
            y = 2,
            yend = 0,
            xend = ess_s
        ),
        col = '#32a891',
        linetype = 'dashed'
    ) +
    geom_point(
        aes(
            x = ess_s,
            y = 2
        ), 
        col = '#32a891', 
        #linetype = 'dashed', 
        size = 2
    ) + 
    annotate(
        'text', 
        x = ess_s + c(0, -10, +10), 
        y = 2, 
        label = c('CINP', 'DIS3', 'RIOK1'), 
        col = '#32a891', 
        angle = 60, 
        vjust = 0, 
        hjust = -0.2, 
        size = 3
    ) + 
    scale_x_log10() + 
    theme_light() + 
    labs(
        x = 'Number of connections (+1) with \n reference essential genes', 
        y = 'density'
    ) +
    ylim(0, 3)


ggsave('./paper_figures/fig2_ess_connection_density.pdf', height = 2, width = 4)


