# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: R [conda env:anaconda-ukbb-trait-analysis-R]
#     language: R
#     name: conda-env-anaconda-ukbb-trait-analysis-R-r
# ---

# %%
snakefile = normalizePath("../Snakefile")


# %%
# rm(snakemake)

# %%
if (! exists("snakemake")) {
    rule = "compose_paper_figure"
    python = "/opt/modules/i12g/anaconda/envs/florian4/bin/python"
    wildcards = paste(
        # 'comparison=all',
        sep=','
    )

    cmd=c(
        "-m snakemk_util.main",
        "--rule", rule,
        "--snakefile", normalizePath(snakefile),
        "--root", dirname(normalizePath(snakefile)),
        "--wildcards", paste0('"', wildcards, '"'),
        "--gen-preamble RScript",
        "--create_dirs"
    )
    eval(parse(text=system2(python, cmd, stdout=TRUE, stderr="")))
}


# %%
snakemake@input

# %%
snakemake@output

# %%
snakemake@params

# %% [markdown]
# # load libraries

# %%
library(IRdisplay)
library(Cairo)

library(data.table)
library(arrow)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(ggtext)
library(patchwork)
library(grid)
library(gridExtra)

# %%
# CairoSVG function for ggplot2 -------------------------------------------
# Solves some issues with standard svg device - mainly usage of symbols
cairosvg <- function(filename, width, height){
  library(Cairo)
  CairoSVG(file = filename, width = width, height = height)
}

# %%
# Set display size of plots
options(repr.plot.width=12, repr.plot.height=8)

# %%
THEME = theme_bw(base_size=12, base_family = 'Helvetica')

# %% [markdown]
# # subplots

# %% [markdown]
# ## number of significant genes

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$compare_associations_dir, "/num_significants.scatter_plot.parquet")))
plot_df

# %%
plot_df$ratio = plot_df$AbExp_all_tissues / plot_df$LOFTEE_pLoF
plot_df$significant = (plot_df$ratio >= 2) | (plot_df$ratio <= 0.5)
plot_df$phenotype = ifelse(plot_df$significant, plot_df$phenotype_col, "")
plot_df

# %%
plot_df[order(`ratio`, decreasing = TRUE)]

# %%
plot_df[`LOFTEE_pLoF` > `AbExp_all_tissues`]

# %%
plot_df$phenotype_col

# %%
traits_to_show = c(
    # "LDL direct",
    #"Albumin",
    #"Alkaline\nphosphatase",
    # "Aspartate\naminotransferase",
    # "Basophill\ncount",
    # "Phosphate",
    # "IGF1",
    # "Testosterone",
    # "Aspartate\naminotransferase",
    # "Direct\nbilirubin",
    # "SHBG",
    # "Urate",
    "HDL\ncholesterol",
    'Mean sphered\ncell volume',
    "Triglycerides",
    #"Alanine\naminotransferase",
    #"Apolipoprotein\nA",
    # "c reactive\nprotein"
    ""
)
print(traits_to_show)

# %%
plot_1 = (
    ggplot(plot_df, aes(x=`LOFTEE_pLoF`, y=`AbExp_all_tissues`))
    + geom_point(size=3)
    + geom_abline(slope=1, color="black", linetype="dashed")
    + labs(
        title="Number of significantly associating genes\n(p-values, alpha=0.05)",
        x="Genes discovered by LOFTEE pLoF",
        y=" \nGenes discovered by AbExp"
    )
    + scale_x_continuous(limits = c(0, 30))
    + scale_y_continuous(limits = c(0, 30))
    + theme(
        # figure.size=c(6, 4),
        axis.text.x=element_text(
            angle=30,
            # hjust=1
        ),
        strip.text.y=element_text(
            angle=0,
        ),
        # title=element_text(linespacing=1.4),
    )
    + geom_text_repel(
        # data = plot_df[(significant == TRUE)], 
        data = plot_df[(`phenotype_col` %in% traits_to_show)],
        min.segment.length = 0,
        point.size = 5,
        aes(label=`phenotype_col`)
    )
    # + facet_wrap("covariates")
    + coord_equal()
    + THEME
)

plot_1

# %%
plot_1 + geom_text_repel(aes(label=`phenotype_col`))

# %%
path = paste0(snakemake@params$output_basedir, "/num_significants")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot_1, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), plot_1, width = 8, height = 6, dpi=600, device=cairo_pdf)

display_png(file=paste0(path, ".png"))

# %% [markdown]
# ## common-variant vs rare-variant

# %%
snakemake@params$alanine_aminotransferase_dir

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$alanine_aminotransferase_dir, "/predictions_vs_common_PRS.parquet")))
plot_df

# %%
phenotype_col = "Alanine Aminotransferase"

lm_coef = coef(lm(plot_df$`Age+Sex+PC+PRS+AbExp_all_tissues`~plot_df$`Age+Sex+PC+PRS`))
lm_coef

x_intercept = 1.25
sd_shift = 1

# %%
plot_2 <- (
    ggplot(plot_df, aes(
        x=`Age+Sex+PC+PRS+AbExp_all_tissues`, 
        y=`Age+Sex+PC+PRS`,
        color=`full_model_new_risk`
    ))
    + geom_bin2d(bins=100, linewidth=0.3)
    + geom_smooth(method="lm", color="red")
    + geom_abline(slope = lm_coef[2], intercept=lm_coef[1] + sd_shift, linetype="dashed")
    + geom_abline(slope = lm_coef[2], intercept=lm_coef[1] - sd_shift, linetype="dashed")
    + annotate(
        "segment",
        x = x_intercept,
        xend = x_intercept + lm_coef[1] + sd_shift * lm_coef[2],
        y = lm_coef[1] + x_intercept * lm_coef[2] + lm_coef[1] + sd_shift * lm_coef[2],
        yend = lm_coef[1] + x_intercept * lm_coef[2],
        arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm")),
        # label = "my line"
    )
    + annotate(
        geom="richtext",
        x = (2*x_intercept + lm_coef[1] + sd_shift * lm_coef[2]) / 2,
        y = lm_coef[1] + x_intercept * lm_coef[2] + (lm_coef[1] + sd_shift * lm_coef[2]) / 2,
        angle=atan(-lm_coef[2]) * 180 / pi,
        vjust=-0.5,
        # hjust=1,
        color="black",
        fill="white",
        label = paste0("+/- ", sd_shift, " SD")
    )
    # + geom_text(
    # geom_label(aes(label=str_wrap(Period,12), x=(StartDate + EndDate)/2), size=3)
    + scale_fill_gradient(name = "Individuals", trans = "log", breaks=c(1, 10, 100, 1000), low="lightgrey", high="black")
    + scale_color_manual(
        name="",
        values=c(`TRUE`="orange"),
        labels=c(`TRUE`="Prediction differs by 1 SD"),
        na.value = "#00000000"
    )
    + coord_fixed()
    + labs(
        x="Prediction based on common variants and AbExp",
        y="Prediction based on common variants",
        title=paste0("Predictions for '", phenotype_col, "' of\ncommon-variant model vs. AbExp-DNA model")
    )
    + THEME
    + theme(plot.title = element_text(hjust = 0.5))
)
plot_2

# %%
path = paste0(snakemake@params$output_basedir, "/common_vs_rare_variant_scatter")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot_2, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), plot_2, width = 8, height = 6, dpi=600, device=cairo_pdf)

display_png(file=paste0(path, ".png"))

# %% [markdown]
# ## rare-variant vs. y_true

# %%
snakemake@params$alanine_aminotransferase_dir

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$alanine_aminotransferase_dir, "/predictions_cloud.parquet")))
plot_df = plot_df[`variable` == 'Age+Sex+PC+PRS+AbExp_all_tissues \n r²=0.300']
plot_df

# %%
unique(plot_df$variable)

# %%
phenotype_col = "Alanine Aminotransferase"

plot_3 <- (
    ggplot(plot_df, aes(
        x=`value`,
        y=`measurement`,
        color=`full_model_new_risk`,
    ))
    + geom_bin2d(bins=100, linewidth=0.3)
    + geom_smooth(method="lm", color="red")
    + scale_fill_gradient(name = "Individuals", trans = "log", breaks=c(1, 10, 100, 1000), low="lightgrey", high="black")
    + scale_color_manual(
        name="",
        values=c(`TRUE`="orange"),
        labels=c(`TRUE`="Prediction differs by 1 SD"),
        na.value = "#00000000"
    )
    # + coord_fixed()
    + labs(
        x="Prediction based on common variants and AbExp",
        y=paste0("'", phenotype_col, "'\nmeasurement"),
        title=paste0("Predictions for '", phenotype_col, "' of\ncommon-variant model vs. AbExp-DNA model")
    )
    + THEME
    + theme(plot.title = element_text(hjust = 0.5))
)
plot_3

# %%
path = paste0(snakemake@params$output_basedir, "/rare_variant_vs_y_true")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot_3, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), plot_3, width = 8, height = 6, dpi=600, device=cairo_pdf)

display_png(file=paste0(path, ".png"))

# %% [markdown]
# ## r² bar plot proportional difference

# %%
phenotype_label_order = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_dir, "/r2_bar_plot_proportional_difference.LOFTEE_pLoF__vs__AbExp_all_tissues.parquet")))
phenotype_label_order = phenotype_label_order[,.(`phenotype_label_order`=mean(`difference_to_LOFTEE_pLoF`)), by=`phenotype_col`]
phenotype_label_order = phenotype_label_order[order(`phenotype_label_order`, decreasing = TRUE)]$phenotype_col
phenotype_label_order

# %%
snakemake@params$compare_risk_scores_dir

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_dir, "/r2_bar_plot_proportional_difference.LOFTEE_pLoF__vs__AbExp_all_tissues.parquet")))
plot_df

# %%
plot_4 = (
    ggplot(plot_df, aes(
        x=reorder(`phenotype_col`, `proportional_difference_to_LOFTEE_pLoF`),
        y=`proportional_difference_to_LOFTEE_pLoF`,
        color=`significant`,
    ))
    # + geom_boxplot()
    # + geom_bar(
    #     stat=stat_summary(fun_y=np.mean),
    # )
    + scale_x_discrete(breaks=phenotype_label_order)
    + stat_summary(
        fun.ymin=function(x){ mean(x) - sd(x) },
        fun.ymax=function(x){ mean(x) + sd(x) },
        geom = "errorbar",
        color="black"
    )
    + stat_mean(
        size=3,
        geom = "point"
    )
    + scale_color_manual(
        name="",
        values=c(`TRUE`="red"),
        labels=c(`TRUE`="significant"),
        na.value = "black"
    )
     + scale_y_continuous(
         labels=scales::percent
     )
    + labs(
        x="phenotype",
        y="relative increase in r² between 'AbExp' and 'LOFTEE pLoF'",
        title="Comparison of phenotype prediction models using different feature sets",
    )
    + THEME
    + theme(
        # legend_text=element_text(linespacing=1.4),
        # figure_size=(8, 12),
        # axis_text_x=element_text(
        # #     rotation=45,
        # #     hjust=1
        #     # vjust=10,
        # ),
        # strip_text_y=element_text(
        #     rotation=0,
        # ),
        # title=element_text(linespacing=1.4, vjust=-10),
        # axis_title_x=element_text(linespacing=1.4, vjust=-10),
    )
    # + coord_equal()
    + coord_flip()
)

plot_4

# %%
path = paste0(snakemake@params$output_basedir, "/r2_bar_plot_proportional_difference")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot_4, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), plot_4, width = 8, height = 6, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))
display_png(file=paste0(path, ".png"))

# %% [markdown]
# ## number of individuals where the absolute error increased/decreased

# %%
snakemake@params$compare_risk_scores_dir

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_dir, "/num_individuals_with_changed_abserr.diverging_barplot.parquet")))
plot_df = plot_df[`sd_cutoff_label` %in% c("0.5 SD", "0.75 SD", "1.0 SD")]
plot_df

# %%
plot_5 = (
    ggplot(subset(plot_df, sd_cutoff == '1'), aes(x=reorder(`phenotype_col`, `reduce`), fill = `feature_set`, width=.8))
    + geom_col(aes(y=-`increase`), position=position_dodge(width=0.8, preserve='single'), alpha=1)
    + geom_col(aes(y=`reduce`), position=position_dodge(width=0.8, preserve='single'))
    # + geom_col(aes(y="-increase"), position=positions.position_dodge(preserve='single'), alpha=0.5)
    # + geom_col(aes(y="reduce"), position=positions.position_dodge(preserve='single'))
    + geom_hline(aes(yintercept = 0))
    + scale_x_discrete(breaks=phenotype_label_order)
    + scale_fill_manual(values=c("orange", "#619CFF"),labels=c(
                `AbExp_all_tissues` = "AbExp",
                `LOFTEE_pLoF` = "LOFTEE pLoF"
    ))
    + theme(
        # figure_size=(8, 8),
    )
    #+ facet_wrap("sd_cutoff_label", scales="free_x")
    # + scale_x_discrete(limits=plot_df.query("feature_set=='AbExp_all_tissues' and sd_cutoff==1").sort_values("reduce")["phenotype_col"].to_list())
    # + scale_fill_manual(["red", "blue"],breaks=reversed(["LOFTEE_pLoF", "AbExp_all_tissues"]))
    + coord_flip()
    + labs(
        y='Nr. of individuals with:\n◀---- increased prediction error ---- ┃ ---- reduced prediction error ----▶',
        x="phenotype",
        title="Number of individuals where\nthe absolute error compared to the common-variant model\nchanges by more than 1.0 standard deviation"
    )
    + THEME
    + theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust=1, vjust=1, angle = 45)
    )
)

plot_5

# %%
path = paste0(snakemake@params$output_basedir, "/num_individuals_with_changed_abserr")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot_5, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), plot_5, width = 8, height = 6, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))
display_png(file=paste0(path, ".png"))

# %% [markdown]
# # common plot

# %%
lower_lhs = (
    (
        plot_4 
        + ggtitle(NULL)
        + xlab("")
        + ylab("Relative increase in r² between\n AbExp and LOFTEE pLoF")
        + theme(
            # axis.title.y = element_blank(),
            legend.position = "bottom"
            # legend.title = element_text("Nr. of individuals")
            # legend.title = element_blank(),
        )
    ) 
    + (
        plot_5
        # + facet_wrap("sd_cutoff_label")
        + scale_fill_discrete(
            name="",
            labels=c(
                `AbExp_all_tissues` = "AbExp-DNA",
                `LOFTEE_pLoF` = "LOFTEE pLoF"
            )
        )
        + ggtitle(NULL)
        + ylab("Individuals with:\n◀-- increased error -- ┃ -- reduced error --▶")
        + theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            # legend.title = element_text("Nr. of individuals")
            # legend.title = element_blank(),
            legend.position = "bottom"
        )
        + scale_fill_manual(values=c("orange", "#619CFF"),labels=c(
                `AbExp_all_tissues` = "AbExp",
                `LOFTEE_pLoF` = "LOFTEE pLoF"
    ))
    )
    #+ ylab("proportional difference in r² between\n AbExp-DNA and LOFTEE pLoF")
    + plot_layout(widths = c(1, 1))
    & theme(plot.margin = margin(0,0,0,0))
)
lower_lhs

# %%
as_ggplot(cowplot::get_legend(
    plot_2
    + guides(fill="none")
))

# %%
rhs = (
    (
        plot_2
        + ggtitle(NULL)
        + guides(
            color="none",
            fill=guide_colourbar(
                direction = "horizontal",
                title.position = "top"
            )
        )
        + theme(
            legend.position = c(0.22, 0.87),
            legend.background = element_rect(fill="#00000000"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            # axis.ticks.x = element_blank()
        )
    )
    # / plot_spacer()
    / (
        plot_3
        + ggtitle(NULL)
        + guides(
            color="none",
            fill=guide_colourbar(
                direction = "horizontal",
                title.position = "top"
            )
        )
        + theme(
            legend.position = c(0.22, 0.87),
            legend.background = element_rect(fill="#00000000")
            # legend.position = c(0.3, 0.65)
            # legend.title = element_text("Nr. of individuals")
            # legend.title = element_blank(),
        )
    )
    + plot_layout(nrow=2)
    & theme(plot.margin = margin(0,0,0,0))
)
rhs = (
    wrap_plots(rhs) / as_ggplot(cowplot::get_legend(
        plot_2
        + guides(fill="none")
        + guides(color = guide_legend(override.aes = list(colour = "orange", linewidth = 10)))
        + theme(
            legend.background = element_rect(fill="#00000000")
        )
    ))
    + plot_layout(nrow=2, heights = c(15, 1))
)
rhs

# %%
commonplot <- ggarrange(
  ncol = 2, nrow = 1, widths = c(1,2),
  ggarrange(nrow = 2, ncol = 1, heights= c(1,2), labels = c('a', 'b'), legend="bottom", align="v",
            plot_1 + ggtitle(NULL),
            rhs
  ),
  ggarrange(labels=c("c"),
      lower_lhs
  )
  # ggarrange(nrow = 2, labels = c('b', 'c'), heights = c(1,1), common.legend = TRUE, legend="bottom", align = "h",
  #           plot_2 + ggtitle("") + theme(
  #               legend.margin = margin(6, 6, 6, 6)
  #           ),
  #           plot_3 + ggtitle("") + theme(
  #               legend.margin = margin(6, 6, 6, 6)
  #           )
  # )
)
commonplot

# %%
path = paste0(snakemake@params$output_basedir, "/paper_figure")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), commonplot, width = 16, height = 12, dpi=600, type = "cairo")
ggsave(paste0(path, ".svg"), commonplot, width = 16, height = 12, dpi=600, device=svg)
ggsave(paste0(path, ".pdf"), commonplot, width = 16, height = 12, dpi=600, device=cairo_pdf)

display_png(file=paste0(path, ".png"))

# %%

# %%

# %%
