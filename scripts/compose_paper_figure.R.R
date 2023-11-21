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
library(dplyr)
library(stringr)

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
THEME = theme_bw(base_size=12, base_family = 'Helvetica') + theme(plot.tag = element_text(face = "bold"))

# %% [markdown]
# # subplots

# %% [markdown]
# ## number of significant genes

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$compare_associations_dir, "/num_significants.scatter_plot.parquet")))
plot_df

# %%
plot_df$ratio = plot_df$AbExp_all_tissues / plot_df$LOFTEE
plot_df$significant = (plot_df$ratio >= 2) | (plot_df$ratio <= 0.5)
# plot_df$phenotype = ifelse(plot_df$significant, plot_df$phenotype_col, "")
plot_df

# %%
plot_df[order(`ratio`, decreasing = TRUE)]

# %%
plot_df[`LOFTEE` > `AbExp_all_tissues`]

# %%
plot_df[`LOFTEE` < `AbExp_all_tissues`]

# %%
sum(plot_df[,`LOFTEE`])

# %%
sum(plot_df[,`AbExp_all_tissues`])

# %%
sum(plot_df[,`AbExp_all_tissues`]) / sum(plot_df[,`LOFTEE`])

# %%
plot_df$phenotype_col

# %%
traits_to_show = c(
    # "LDL direct",
    #"Albumin",
    "Alkaline\nphosphatase",
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
max_n = max(plot_df$`LOFTEE`, plot_df$`AbExp_all_tissues`)
max_n

# %%
plot_1 = (
    ggplot(plot_df, aes(x=`LOFTEE`, y=`AbExp_all_tissues`))
    + geom_point(size=2)
    + geom_abline(slope=1, color="black", linetype="dashed")
    + labs(
        title="Number of significantly associating genes\n(p-values, alpha=0.05)",
        x="Genes discovered using LOFTEE",
        y=" \nGenes discovered using AbExp"
    )
    + scale_x_continuous(limits = c(0, max_n))
    + scale_y_continuous(limits = c(0, max_n))
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
    # + coord_equal()
    + THEME
    + labs(tag="a")
)

plot_1

# %%
options(repr.plot.width=20, repr.plot.height=15)

plot_1 + geom_text_repel(aes(label=`phenotype_col`), max.overlaps=50)

# %%
w=4
h=4
path = paste0(snakemake@params$output_basedir, "/num_significants")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot_1, width = w, height = h, dpi=600, type = "cairo")
ggsave(paste0(path, ".svg"), plot_1, width = w, height = h, dpi=600, device=svg)
ggsave(paste0(path, ".pdf"), plot_1, width = w, height = h, dpi=600, device=cairo_pdf)

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
        x="Prediction based on\ncommon variants and AbExp",
        y="Prediction based on\ncommon variants",
        title=paste0(phenotype_col, " level")
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
plot_df = plot_df[`variable` == 'Age+Sex+PC+PRS+AbExp_all_tissues \n r²=0.295']
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
        x="Prediction based on\ncommon variants and AbExp",
        y=paste0("Measurement"),
        title=paste0(phenotype_col, " level")
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
phenotype_label_order = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_dir, "/r2_bar_plot_proportional_difference.LOFTEE__vs__AbExp_all_tissues.parquet")))
phenotype_label_order = phenotype_label_order[,.(`phenotype_label_order`=mean(`difference_to_LOFTEE`)), by=`phenotype_col`]
phenotype_label_order = phenotype_label_order[order(`phenotype_label_order`, decreasing = TRUE)]$phenotype_col
phenotype_label_order

# %%
snakemake@params$compare_risk_scores_dir

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_dir, "/r2_bar_plot_proportional_difference.LOFTEE__vs__AbExp_all_tissues.parquet")))
plot_df

# %%
plot_4 = (
    ggplot(plot_df, aes(
        x=reorder(`phenotype_col`, `proportional_difference_to_LOFTEE`),
        y=`proportional_difference_to_LOFTEE`,
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
        y="relative difference in R² between 'AbExp' and 'LOFTEE pLoF'",
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
plot_df[,.(`error_increase` = sum(`increase`), `error_reduce` = sum(`reduce`)), by=c("feature_set", "sd_cutoff", "sd_cutoff_label")]

# %%
dodge_width=0.9
plot_5 = (
    ggplot(subset(plot_df, sd_cutoff == '1'), aes(x=reorder(`phenotype_col`, `reduce`), fill = `feature_set`, width=.8))
    + geom_col(aes(y=-`increase`), position=position_dodge(width=dodge_width, preserve='single'), alpha=1)
    + geom_col(aes(y=`reduce`), position=position_dodge(width=dodge_width, preserve='single'))
    # + geom_col(aes(y="-increase"), position=positions.position_dodge(preserve='single'), alpha=0.5)
    # + geom_col(aes(y="reduce"), position=positions.position_dodge(preserve='single'))
    + geom_hline(aes(yintercept = 0))
    + scale_x_discrete(breaks=phenotype_label_order)
    + scale_fill_manual(
        labels=c(
            `AbExp_all_tissues` = "AbExp",
            `LOFTEE` = "LOFTEE"
        ),
        # values=c("orange", "#619CFF"),
        values=c(
            `AbExp_all_tissues` = "#48a462",
            `LOFTEE` = "#999999"
            # `LOFTEE` = "#619CFF"
        ),
    )
    + theme(
        # figure_size=(8, 8),
    )
    #+ facet_wrap("sd_cutoff_label", scales="free_x")
    # + scale_x_discrete(limits=plot_df.query("feature_set=='AbExp_all_tissues' and sd_cutoff==1").sort_values("reduce")["phenotype_col"].to_list())
    # + scale_fill_manual(["red", "blue"],breaks=reversed(["LOFTEE", "AbExp_all_tissues"]))
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
    )
    #+ ylab("proportional difference in r² between\n AbExp-DNA and LOFTEE pLoF")
    + plot_layout(widths = c(1, 1))
    & theme(plot.margin = margin(0,0,0,0))
)
lower_lhs

# %%
lower_lhs = ggarrange(
    (
        plot_4 
        + ggtitle(NULL)
        + xlab("")
        + ylab("Relative increase in R² between\n AbExp and LOFTEE")
        + theme(
            # axis.title.y = element_blank(),
            legend.position = "bottom"
            # legend.title = element_text("Nr. of individuals")
            # legend.title = element_blank(),
        )
        + labs(tag = "d")
    ), (
        plot_5
        # + facet_wrap("sd_cutoff_label")
        # + scale_fill_discrete(
        #     name="",
        #     labels=c(
        #         `AbExp_all_tissues` = "AbExp-DNA",
        #         `LOFTEE` = "LOFTEE pLoF"
        #     )
        # )
        + ggtitle(NULL)
        + ylab("Individuals with changed error:\n◀ increased ┃ reduced ▶")
        # + ylab("Individuals with:\n◀-- increased error -- ┃ -- reduced error --▶")
        + theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            # legend.title = element_text("Nr. of individuals")
            legend.title = element_blank(),
            legend.position = "bottom"
        )
        # + scale_fill_manual(values=c("orange", "#619CFF"),labels=c(
        #         `AbExp_all_tissues` = "AbExp",
        #         `LOFTEE` = "LOFTEE pLoF"
        # ))
        + labs(tag = "e")
    ),
    align="h",
    # labels=c("c", "d"),
    widths=c(1.7, 1)
    # #+ ylab("proportional difference in r² between\n AbExp-DNA and LOFTEE pLoF")
    # + plot_layout(widths = c(1, 1))
    # & theme(plot.margin = margin(0,0,0,0))
)
lower_lhs

# %%
path = paste0(snakemake@params$output_basedir, "/combined_r2_num_indiv")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), lower_lhs, width = 10, height = 8, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), lower_lhs, width = 10, height = 8, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))
display_png(file=paste0(path, ".png"))

# %%
w=8
h=8

path = paste0(snakemake@params$output_basedir, "/combined_r2_num_indiv.shallow")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), lower_lhs, width = w, height = h, dpi=600, type = "cairo")
ggsave(paste0(path, ".svg"), lower_lhs, width = w, height = h, dpi=600, device=svg)
ggsave(paste0(path, ".pdf"), lower_lhs, width = w, height = h, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))
display_png(file=paste0(path, ".png"))

# %%
# as_ggplot(cowplot::get_legend(
#     plot_2
#     + guides(fill="none")
# ))

# %%
apply_consistent_x_lims <- function(this_plot){
    num_plots <- length(this_plot$layers)
    x_lims <- lapply(1:num_plots, function(x) ggplot_build(this_plot[[x]])$layout$panel_scales_x[[1]]$range$range)
    min_x <- min(unlist(x_lims))
    max_x <- max(unlist(x_lims))
    this_plot & xlim(min_x, max_x)
}
apply_consistent_y_lims <- function(this_plot){
    num_plots <- length(this_plot$layers)
    y_lims <- lapply(1:num_plots, function(x) ggplot_build(this_plot[[x]])$layout$panel_scales_y[[1]]$range$range)
    min_y <- min(unlist(y_lims))
    max_y <- max(unlist(y_lims))
    this_plot & ylim(min_y, max_y)
}

# %%
rhs = (
    (
        plot_2
        # + ggtitle(NULL)
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
            # axis.title.x = element_blank(),
            # axis.text.x = element_blank(),
            # axis.ticks.x = element_blank()
        )
        + labs(tag = "b")
    )
    # / plot_spacer()
    / (
        plot_3
        # + ggtitle(NULL)
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
        + labs(tag = "c")
    )
    + plot_layout(nrow=2)
    # & theme(plot.margin = margin(0,0,0,0))
)
rhs = apply_consistent_x_lims(rhs)
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
w=4
h=8
path = paste0(snakemake@params$output_basedir, "/alanine_aminotransferase")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), rhs, width = w, height = h, dpi=600, type = "cairo")
ggsave(paste0(path, ".svg"), rhs, width = w, height = h, dpi=600, device=svg)
ggsave(paste0(path, ".pdf"), rhs, width = w, height = h, dpi=600, device=cairo_pdf)

display_png(file=paste0(path, ".png"))

# %%
commonplot <- ggarrange(
    ncol = 2, nrow = 1, widths = c(1.1,2),
    ggarrange(nrow = 2, ncol = 1, heights= c(1, 2.4), legend="bottom", align="v",
            (
                plot_1
                + ggtitle(NULL)
                + coord_fixed(ratio=1)
                + labs(tag = "a")
            ),
            rhs
    ),
    ggarrange(
        lower_lhs
        + ggtitle(NULL)
        + xlab("")
        # + ylab("Relative increase in r² between\n AbExp and LOFTEE pLoF")
        + theme(
            # axis.title.y = element_blank(),
            legend.position = "bottom"
            # legend.title = element_text("Nr. of individuals")
            # legend.title = element_blank(),
        )
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
w=11.5
h=12
ggsave(paste0(path, ".png"), commonplot, width = w, height = h, dpi=600, type = "cairo")
ggsave(paste0(path, ".svg"), commonplot, width = w, height = h, dpi=600, device=svg)
ggsave(paste0(path, ".pdf"), commonplot, width = w, height = h, dpi=600, device=cairo_pdf)

display_png(file=paste0(path, ".png"))

# %% [markdown]
# # number of significant genes for AbExp aggregations

# %%
rename_models = c(
    `LOFTEE` = "LOFTEE",
    `AbExp_all_tissues` = "AbExp all tissues",
    `minimum_AbExp` = "Minimum AbExp",
    `median_AbExp` = "Median AbExp"
)

# %%
plot_df = as.data.table(read_parquet(paste0(snakemake@params$compare_associations_dir, "/num_significants.scatter_plot.parquet")))
plot_df

# %%
plot_df[, .(sum(`AbExp_all_tissues`), sum(`LOFTEE`), sum(`median_AbExp`), sum(`minimum_AbExp`)), by="covariates"]

# %%
plot_df$ratio = plot_df$AbExp_all_tissues / plot_df$LOFTEE
plot_df$significant = (plot_df$ratio >= 2) | (plot_df$ratio <= 0.5)
# plot_df$phenotype = ifelse(plot_df$significant, plot_df$phenotype_col, "")
plot_df

# %%
traits_to_show = c(
    # "LDL direct",
    #"Albumin",
    "Alkaline\nphosphatase",
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
max_n = max(plot_df$`minimum_AbExp`, plot_df$`median_AbExp`, plot_df$`AbExp_all_tissues`)
max_n

# %%
melted_plot_df = melt(
    plot_df,
    id.vars=c("phenotype_col", "covariates"),
    measure.vars=c("minimum_AbExp", "median_AbExp")
)
melted_plot_df = merge(
    melted_plot_df,
    plot_df[,.(`phenotype_col`, `covariates`, `AbExp_all_tissues`)],
    on=c("phenotype_col", "covariates"),
    how="outer"
)
melted_plot_df[, variable:=recode(melted_plot_df$variable, `minimum_AbExp` = "Minimum", `median_AbExp` = "Median")]
melted_plot_df

# %%
plot = (
    ggplot(melted_plot_df, aes(x=`AbExp_all_tissues`, y=`value`, color=`phenotype_col`))
    + geom_point(size=2)
    + geom_abline(slope=1, color="black", linetype="dashed")
    + labs(
        # title="Number of significantly associating genes\n(p-values, alpha=0.05)",
        y="Genes discovered by AbExp aggregated across tissues",
        x="Genes discovered by AbExp using all tissues"
    )
    + scale_x_continuous(limits = c(0, max_n))
    + scale_y_continuous(limits = c(0, max_n))
    # + geom_text_repel(
    #     # data = plot_df[(significant == TRUE)], 
    #     data = plot_df[(`phenotype_col` %in% traits_to_show)],
    #     min.segment.length = 0,
    #     point.size = 5,
    #     aes(label=`phenotype_col`)
    # )
    # + facet_wrap("covariates")
    # + coord_equal()
    + coord_flip()
    + facet_wrap("variable", nrow=1)
    + THEME
    + theme(
        legend.position="bottom",
        legend.title = element_blank()
    )
    # + labs(tag="a")
)

plot

# %%
w=8
h=8
path = paste0(snakemake@params$output_basedir, "/num_significants_on_aggregates")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), plot, width = w, height = h, dpi=600, device=cairo_pdf)
ggsave(paste0(path, ".svg"), plot, width = w, height = h, dpi=600, device=svg)

display_png(file=paste0(path, ".png"))

# %% [markdown]
# # QQ-Plot

# %%
combined_qqplot_df = as.data.table(read_parquet(paste0(snakemake@input$`combined_qqplot_pq`)))
combined_qqplot_df

# %%
plot_df = combined_qqplot_df[`subsample` & `covariates` == 'randomized_sex_age_genPC_CLMP_PRS']
plot_df$model = recode(factor(plot_df$feature_set, levels=names(rename_models)), !!!rename_models)
plot_df

# %%
qq_plot = (
    ggplot(plot_df, aes(x=`theoretical`, y=`sample`))
    + geom_point()
    + geom_abline(slope=1, linetype="dashed", color="red")
    # + pn.scale_x_log10(limits=(10**-20, 1))
    # + pn.scale_y_log10(limits=(10**-20, 1))
    + labs(
        # title="\n".join([
        #     f"""Q-Q plot of randomized p-values vs. random uniform distribution""",
        #     # f"""phenotype: '{snakemake.wildcards["phenotype_col"]}'""",
        #     # f"""feature set: '{snakemake.wildcards["feature_set"]}'""",
        #     # f"""covariates: '{snakemake.wildcards["covariates"]}'""",
        # ]),
        x="-log10(p) theoretical",
        y="-log10(p) sample",
    )
    + facet_wrap("model", scales="free")
    + THEME
    + theme(
        legend.position="bottom",
        legend.title = element_blank()
    )
)
qq_plot

# %%
w=8
h=8
path = paste0(snakemake@params$output_basedir, "/qq_plot")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), qq_plot, width = w, height = h, dpi=600, type = "cairo")
display_png(file=paste0(path, ".png"))
ggsave(paste0(path, ".pdf"), qq_plot, width = w, height = h, dpi=600, device=cairo_pdf)
ggsave(paste0(path, ".svg"), qq_plot, width = w, height = h, dpi=600, device=svg)

# %% [markdown]
# # Linear risk scores

# %% [markdown]
# ## r² bar plot proportional difference

# %%
linear_phenotype_label_order = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_linear_dir, "/r2_bar_plot_proportional_difference.LOFTEE__vs__AbExp_all_tissues.parquet")))
linear_phenotype_label_order = linear_phenotype_label_order[,.(`phenotype_label_order`=mean(`difference_to_LOFTEE`)), by=`phenotype_col`]
linear_phenotype_label_order = linear_phenotype_label_order[order(`phenotype_label_order`, decreasing = TRUE)]$phenotype_col
linear_phenotype_label_order

# %%
linear_pheno_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_linear_dir, "/r2_bar_plot_proportional_difference.LOFTEE__vs__AbExp_all_tissues.parquet")))
linear_pheno_df

# %%
linear_plot_4 = (
    ggplot(linear_pheno_df, aes(
        x=reorder(`phenotype_col`, `proportional_difference_to_LOFTEE`),
        y=`proportional_difference_to_LOFTEE`,
        color=`significant`,
    ))
    # + geom_boxplot()
    # + geom_bar(
    #     stat=stat_summary(fun_y=np.mean),
    # )
    + scale_x_discrete(breaks=linear_phenotype_label_order)
    + stat_summary(
        fun.min=function(x){ mean(x) - sd(x) },
        fun.max=function(x){ mean(x) + sd(x) },
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
        y="relative difference in R² between 'AbExp' and 'LOFTEE pLoF'",
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

linear_plot_4

# %%
path = paste0(snakemake@params$output_basedir, "/r2_bar_plot_proportional_difference__linear")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), linear_plot_4, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), linear_plot_4, width = 8, height = 6, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))
display_png(file=paste0(path, ".png"))

# %% [markdown]
# ## number of individuals where the absolute error increased/decreased

# %%
snakemake@params$compare_risk_scores_dir

# %%
linear_pheno_indivs_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_linear_dir, "/num_individuals_with_changed_abserr.diverging_barplot.parquet")))
linear_pheno_indivs_df = linear_pheno_indivs_df[`sd_cutoff_label` %in% c("0.5 SD", "0.75 SD", "1.0 SD")]
linear_pheno_indivs_df

# %%
linear_pheno_indivs_df[,.(`error_increase` = sum(`increase`), `error_reduce` = sum(`reduce`)), by=c("feature_set", "sd_cutoff", "sd_cutoff_label")]

# %%
dodge_width=0.9
linear_plot_5 = (
    ggplot(subset(linear_pheno_indivs_df, sd_cutoff == '1'), aes(x=reorder(`phenotype_col`, `reduce`), fill = `feature_set`, width=.8))
    + geom_col(aes(y=-`increase`), position=position_dodge(width=dodge_width, preserve='single'), alpha=1)
    + geom_col(aes(y=`reduce`), position=position_dodge(width=dodge_width, preserve='single'))
    # + geom_col(aes(y="-increase"), position=positions.position_dodge(preserve='single'), alpha=0.5)
    # + geom_col(aes(y="reduce"), position=positions.position_dodge(preserve='single'))
    + geom_hline(aes(yintercept = 0))
    + scale_x_discrete(breaks=linear_phenotype_label_order)
    + scale_fill_manual(
        labels=c(
            `AbExp_all_tissues` = "AbExp",
            `LOFTEE` = "LOFTEE"
        ),
        # values=c("orange", "#619CFF"),
        values=c(
            `AbExp_all_tissues` = "#48a462",
            `LOFTEE` = "#999999"
            # `LOFTEE` = "#619CFF"
        ),
    )
    + theme(
        # figure_size=(8, 8),
    )
    #+ facet_wrap("sd_cutoff_label", scales="free_x")
    # + scale_x_discrete(limits=plot_df.query("feature_set=='AbExp_all_tissues' and sd_cutoff==1").sort_values("reduce")["phenotype_col"].to_list())
    # + scale_fill_manual(["red", "blue"],breaks=reversed(["LOFTEE", "AbExp_all_tissues"]))
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

linear_plot_5

# %%
path = paste0(snakemake@params$output_basedir, "/num_individuals_with_changed_abserr__linear")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), linear_plot_5, width = 8, height = 6, dpi=600, type = "cairo")
ggsave(paste0(path, ".pdf"), linear_plot_5, width = 8, height = 6, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))
display_png(file=paste0(path, ".png"))

# %% [markdown]
# ## rsquared

# %%
linear_rsquared_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_linear_dir, "/rsquared.parquet")))
linear_rsquared_df = melt(linear_rsquared_df, id.vars=c("phenotype_col", "feature_set", "covariates"), variable="predictor_variables", value.name = "linear model R²")
linear_rsquared_df[, `predictor_variables` := recode(as.factor(`predictor_variables`), !!!c(`full_model_r2`= "full", `restricted_model_r2` = "restricted"))]
linear_rsquared_df = linear_rsquared_df[`predictor_variables` != "basic_model_r2"]
linear_rsquared_df

# %%
nonlinear_rsquared_df = as.data.table(read_parquet(paste0(snakemake@params$compare_risk_scores_dir, "/rsquared.parquet")))
nonlinear_rsquared_df = melt(nonlinear_rsquared_df, id.vars=c("phenotype_col", "feature_set", "covariates"), variable="predictor_variables", value.name = "gradient boosted trees model R²")
nonlinear_rsquared_df[, `predictor_variables` := recode(as.factor(`predictor_variables`), !!!c(`full_model_r2`= "full", `restricted_model_r2` = "restricted"))]
nonlinear_rsquared_df = nonlinear_rsquared_df[`predictor_variables` != "basic_model_r2"]
nonlinear_rsquared_df

# %%
rsquared_df = merge(linear_rsquared_df, nonlinear_rsquared_df, by=c("phenotype_col", "feature_set", "covariates", "predictor_variables"), how="inner")
rsquared_df[, `feature_set` := recode(`feature_set`, `AbExp_all_tissues` = "AbExp all tissues", `LOFTEE` = "LOFTEE")]
rsquared_df$phenotype_col = str_replace_all(rsquared_df$phenotype_col, "_", " ")
rsquared_df

# %%
plot = (
    ggplot(rsquared_df, aes(x=`linear model R²`, y=`gradient boosted trees model R²`, color=`phenotype_col`))
    + geom_point(size=2)
    + geom_abline(slope=1, color="black", linetype="dashed")
    + labs(
        # title="Number of significantly associating genes\n(p-values, alpha=0.05)",
        # y="Genes discovered by AbExp aggregated across tissues",
        # x="Genes discovered by AbExp using all tissues"
    )
    # + geom_text_repel(
    #     # data = plot_df[(significant == TRUE)], 
    #     data = plot_df[(`phenotype_col` %in% traits_to_show)],
    #     min.segment.length = 0,
    #     point.size = 5,
    #     aes(label=`phenotype_col`)
    # )
    + facet_grid("feature_set ~ predictor_variables")
    # + coord_equal()
    # + coord_flip()
    # + facet_wrap("variable", nrow=1)
    + THEME
    + theme(
        legend.position="bottom",
        legend.title = element_blank()
    )
    # + labs(tag="a")
)

plot

# %%
plot_df = subset(rsquared_df, (`predictor_variables` == 'full') & (`feature_set` %in% c("AbExp all tissues", "LOFTEE")))

# %%
plot_df$phenotype_col

# %%
traits_to_show = c(
    # "LDL direct",
    #"Albumin",
    "Apolipoprotein B",
    "Lipoprotein A",
    # "Alkaline phosphatase",
    # "Aspartate\naminotransferase",
    # "Basophill\ncount",
    # "Phosphate",
    # "IGF1",
    # "Testosterone",
    # "Aspartate\naminotransferase",
    # "Direct\nbilirubin",
    # "SHBG",
    # "Urate",
    "HDL cholesterol",
    # 'Mean sphered cell volume',
    "Triglycerides",
    #"Alanine\naminotransferase",
    #"Apolipoprotein\nA",
    # "c reactive\nprotein"
    ""
)
print(traits_to_show)

# %%
rsquared_plot = (
    ggplot(
        plot_df,
        aes(
            x=`linear model R²`,
            y=`gradient boosted trees model R²`
            # color=`phenotype_col`
        )
    )
    + geom_point(size=2)
    + geom_abline(slope=1, color="black", linetype="dashed")
    + geom_text_repel(
        # data = plot_df[(significant == TRUE)], 
        data = plot_df[(`phenotype_col` %in% traits_to_show)],
        min.segment.length = 0,
        point.size = 5,
        nudge_x = -0.02,
        nudge_y = 0.1,
        aes(label=`phenotype_col`)
    )
    + xlim(0, 0.7)
    + ylim(0, 0.7)
    + labs(
        x="Elastic net model R²",
        y="Gradient boosted trees model R²"
        # title="Number of significantly associating genes\n(p-values, alpha=0.05)",
        # y="Genes discovered by AbExp aggregated across tissues",
        # x="Genes discovered by AbExp using all tissues"
    )
    # + geom_text_repel(
    #     # data = plot_df[(significant == TRUE)], 
    #     data = plot_df[(`phenotype_col` %in% traits_to_show)],
    #     min.segment.length = 0,
    #     point.size = 5,
    #     aes(label=`phenotype_col`)
    # )
    + facet_wrap("feature_set")
    # + coord_equal()
    # + coord_flip()
    # + facet_wrap("variable", nrow=1)
    + THEME
    + theme(
        legend.position="bottom",
        legend.title = element_blank()
    )
    # + labs(tag="a")
)

rsquared_plot

# %% [markdown]
# ## common plot

# %%
supplementary_s4_plot = ggarrange(
    ggarrange(
        (
            linear_plot_4 
            + ggtitle(NULL)
            + xlab("")
            + ylab("Relative increase in R² between\n AbExp and LOFTEE")
            + theme(
                # axis.title.y = element_blank(),
                legend.position = "bottom"
                # legend.title = element_text("Nr. of individuals")
                # legend.title = element_blank(),
            )
            + labs(tag = "a")
        ), (
            linear_plot_5
            # + facet_wrap("sd_cutoff_label")
            # + scale_fill_discrete(
            #     name="",
            #     labels=c(
            #         `AbExp_all_tissues` = "AbExp-DNA",
            #         `LOFTEE` = "LOFTEE pLoF"
            #     )
            # )
            + ggtitle(NULL)
            + ylab("Individuals with changed error:\n◀ increased ┃ reduced ▶")
            # + ylab("Individuals with:\n◀-- increased error -- ┃ -- reduced error --▶")
            + theme(
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                # legend.title = element_text("Nr. of individuals")
                legend.title = element_blank(),
                legend.position = "bottom"
            )
            # + scale_fill_manual(values=c("orange", "#619CFF"),labels=c(
            #         `AbExp_all_tissues` = "AbExp",
            #         `LOFTEE` = "LOFTEE pLoF"
            # ))
            + labs(tag = "b")
        ),
        align="h",
        # labels=c("c", "d"),
        widths=c(1.7, 1)
        # #+ ylab("proportional difference in r² between\n AbExp-DNA and LOFTEE pLoF")
        # + plot_layout(widths = c(1, 1))
        # & theme(plot.margin = margin(0,0,0,0))
    ),
    (
        rsquared_plot
        + labs(tag = "c")
    ),
    nrow=2,
    heights=c(2, 1)
)
supplementary_s4_plot

# %%
supplementary_s4_plot = ggarrange(
    (
        (
            linear_plot_4 
            + ggtitle(NULL)
            + xlab("")
            + ylab("Relative increase in R² between\n AbExp and LOFTEE")
            + theme(
                # axis.title.y = element_blank(),
                legend.position = "bottom"
                # legend.title = element_text("Nr. of individuals")
                # legend.title = element_blank(),
            )
            + labs(tag = "a")
        ) + (
            linear_plot_5
            # + facet_wrap("sd_cutoff_label")
            # + scale_fill_discrete(
            #     name="",
            #     labels=c(
            #         `AbExp_all_tissues` = "AbExp-DNA",
            #         `LOFTEE` = "LOFTEE pLoF"
            #     )
            # )
            + ggtitle(NULL)
            + ylab("Individuals with changed error:\n◀ increased ┃ reduced ▶")
            # + ylab("Individuals with:\n◀-- increased error -- ┃ -- reduced error --▶")
            + theme(
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                # legend.title = element_text("Nr. of individuals")
                legend.title = element_blank(),
                legend.position = "bottom"
            )
            # + scale_fill_manual(values=c("orange", "#619CFF"),labels=c(
            #         `AbExp_all_tissues` = "AbExp",
            #         `LOFTEE` = "LOFTEE pLoF"
            # ))
            + labs(tag = "b")
        )
    ) /
    (
        rsquared_plot
        + labs(tag = "c")
    )
    + plot_layout(heights=c(2, 1))
)
supplementary_s4_plot

# %%
path = paste0(snakemake@params$output_basedir, "/supplementary_s4")
print(paste0("Saving to ", path, "..."))
ggsave(paste0(path, ".png"), supplementary_s4_plot, width = 10, height = 12, dpi=600, type = "cairo")
display_png(file=paste0(path, ".png"))
ggsave(paste0(path, ".svg"), supplementary_s4_plot, width = 10, height = 12, dpi=600, device=svg)
ggsave(paste0(path, ".pdf"), supplementary_s4_plot, width = 10, height = 12, dpi=600, device=cairo_pdf)

# display_pdf(file=paste0(path, ".pdf"))

# %%

# %%
