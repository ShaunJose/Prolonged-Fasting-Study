
# Color set from Lancet journal standard (ggisci package)+ #4 orange
colors_custom <- c("#00468BFF", "#ED0000FF", "#42B540FF",
                   "#925E9FFF", "#0099B4FF", "#FDAF91FF",
                   "#AD002AFF", "#ADB6B6FF", "#1B1919FF"
)

# Color set as above but 1st color (blue) switched to the last one
colors_custom_2 <- c("#ED0000FF", "#42B540FF",
                   "#925E9FFF", "#0099B4FF", "#FDAF91FF",
                   "#AD002AFF", "#ADB6B6FF", "#1B1919FF",
                   "#00468BFF"
)

# Color set for correlation scatter plots, few colors
colors_custom_3 <- c("blue3", "green4")


# 04/2021 AM
# Histogram plots - by group
# Expects df, a list of variable list and grouping variable as input
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotHistogramsByGroup(df, var_list, Timepoint)
func_PlotHistogramsByGroup <- function(df, var_list, by_group) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Histogram # ", var, ":", sep = ""))
    print(target[var])
    
    plot_hist <-
      ggplot(df, aes_string(target[var])) +
      geom_histogram(aes_string(fill = by_group), alpha = 0.8, bins = 10) +
      facet_grid(rows=by_group, drop = TRUE, scales = "free_y") +
      scale_fill_manual(
        values=colors_custom,
        aesthetics = c("fill"), name = "Timepoint") +
      xlab(target[var])
    
    # Print each plot to HTML
    show(plot_hist)
  
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_histogram_", target[var], ".png", sep="")
    figure_svg <- paste("svg_fig_histogram_", target[var], ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
  }
    
}

# 04/2021 AM
# Histogram plots - by group
# Expects df, a list of variable list and grouping variable as input
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotHistogramsByGroup(df, var_list, Timepoint)
func_PlotHistogramsByGroup.out <- function(df, var_list, by_group) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  # Remove outliers & NAs (NAs rm need to select the var column)
  df_out <- df
  df_out <- func_removeOutliersDF(df_out)
  # df_out <- df_out[complete.cases(df_out),]  #remove NAs
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Histogram # ", var, ":", sep = ""))
    print(target[var])

    plot_hist <-
      ggplot(df_out, aes_string(target[var])) +
      geom_histogram(aes_string(fill = by_group), alpha = 0.8, bins = 10) +
      facet_grid(rows=by_group, drop = TRUE, scales = "free_y") +
      scale_fill_manual(
        values=colors_custom,
        aesthetics = c("fill"), name = "Timepoint") +
      xlab(target[var])
    
    # Print each plot to HTML
    show(plot_hist)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_histogram_", target[var], "_out", ".png", sep="")
    figure_svg <- paste("svg_fig_histogram_", target[var], "_out", ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
  }
  
}



# 03/2021 AM
# Box plots with dots and lines - by group
# Expects a list of variable list as input
# Group and color vars are hard coded! 
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotBoxByGroup(df, var_list)
func_PlotBoxByGroup <- function(df, var_list) {
  
  # Libraries
  require(cowplot)
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot # ", var, ":", sep = ""))
    print(target[var])
    
    plot_bar <-
      
      ggplot(df, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = ID, fill = ID), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = ID), size = 1.2, alpha = 0.5) +
      theme_cowplot(14) + # sets basic white theme, font size 14
      # labs(y = labels) + # replacing var names with labels - does not work
      labs(x = NULL) + # removes "Timepoint" x axis title
      theme(axis.title.y = element_text(size = rel(1.5))) + # increase y size
      theme(legend.position = "none") # removes legend
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot__", target[var], ".png", sep="")
    figure_svg <- paste("svg_fig_boxplot__", target[var], ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
  }
}
# 04/2021 AM
# Box plots with points and lines gated to cut-offs (case_group)
# Expects a list of variable list as input
# Group and color vars are hard coded! 

func_PlotBoxByCutoff <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot # ", var, ":", sep = ""))
    print(target[var])
    
    plot_bar <-
      
      ggplot(df, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = case_group, fill = case_group), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = case_group), size = 1.2, alpha = 0.5) +
      scale_fill_manual(
        values = colors_custom,
        aesthetics = c("colour", "fill"), name = "Grouped By:")
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot__", target[var], ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot__", target[var], ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
  }
}



# 04/2021 AM
# Box plots with dots and lines - by group - outliers removed
# Expects a list of variable list as input
# Group and color vars are hard coded! 
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotBoxByGroup(df, var_list)
func_PlotBoxByGroup.out <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  # Remove outliers & NAs (NAs rm need to select the var column)
  df_out <- df
  df_out <- func_removeOutliersDF(df_out)
  # df_out <- df_out[complete.cases(df_out),]  #remove NAs  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot # ", var, ":", sep = ""))
    print(target[var])

    plot_bar <-
      
      ggplot(df_out, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = ID, fill = ID), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = ID), size = 1.2, alpha = 0.5) +
      theme(legend.position = "none") # removes legend
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot_", target[var], "_out", ".png", sep="")
    figure_svg <- paste("svg_fig_boxplot_", target[var], "_out", ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
  }
}


# 04/2021 AM
# Box plots with points and lines gated to cut-offs (case_group)
# Expects a list of variable list as input
# Group and color vars are hard coded! 

func_PlotBoxByCutoff.out <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  # Remove outliers & NAs (NAs rm need to select the var column)
  df_out <- df
  df_out <- func_removeOutliersDF(df_out)
  # df_out <- df_out[complete.cases(df_out),]  #remove NAs  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot # ", var, ":", sep = ""))
    print(target[var])
    
    plot_bar <-
      
      ggplot(df_out, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = case_group, fill = case_group), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = case_group), size = 1.2, alpha = 0.5) +
      scale_fill_manual(
        values = colors_custom,
        aesthetics = c("colour", "fill"), name = "Grouped By:")
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot_", target[var], "_out", ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot_", target[var], "_out", ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
  }
}




# 03/2021 AM - same as above, saving names are different
# Box plots with dots and lines - by group
# Expects a list of variable list as input
# Group and color vars are hard coded! 
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotBoxByGroup.f(df_fr, var_list)
func_PlotBoxByGroup.f <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot (fraction of Baseline) # ", var, ":", sep = ""))
    print(target[var])
    
    plot_bar <-
      
      ggplot(df, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = ID, fill = ID), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = ID), size = 1.2, alpha = 0.5) +
      theme(legend.position = "none") # removes legend
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot_F_", target[var], ".png", sep="")
    figure_svg <- paste("svg_fig_boxplot_F_", target[var], ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
  }
}


# 04/2021 AM
# Box plots with points and lines gated to cut-offs (case_group)
# Expects a list of variable list as input
# Group and color vars are hard coded! 

func_PlotBoxByCutoff.f <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot (fraction of Baseline) # ", var, ":", sep = ""))
    print(target[var])
    
    plot_bar <-
      
      ggplot(df, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = case_group, fill = case_group), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = case_group), size = 1.2, alpha = 0.5) +
      scale_fill_manual(
        values = colors_custom,
        aesthetics = c("colour", "fill"), name = "Grouped By:")
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot_F_", target[var], ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot_F_", target[var], ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
  }
}



# 04/2021 AM - version for removing outliers, NAs before plot
# Box plots with dots and lines - by group
# Expects a list of variable list as input
# Group and color vars are hard coded! 
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotBoxByGroup.f(df_fr, var_list)
func_PlotBoxByGroup.f.out <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
    
  # Remove outliers & NAs (NAs rm need to select the var column)
  df_out <- df
  df_out <- func_removeOutliersDF(df_out)
  # df_out <- df_out[complete.cases(df_out),]  #remove NAs
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot (fraction of Baseline) # ", var, ":", sep = ""))
    print(target[var])

    plot_bar <-
      
      ggplot(df_out, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = ID, fill = ID), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = ID), size = 1.2, alpha = 0.5) +
      theme(legend.position = "none") # removes legend
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot_F_", target[var], "_out", ".png", sep="")
    figure_svg <- paste("svg_fig_boxplot_F_", target[var], "_out", ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
  }
}




# 04/2021 AM - version for removing outliers, NAs before plot
# Box plots with points and lines gated to cut-offs (case_group)
# Expects a list of variable list as input
# Group and color vars are hard coded! 

func_PlotBoxByCutoff.f.out <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  # Remove outliers & NAs (NAs rm need to select the var column)
  df_out <- df
  df_out <- func_removeOutliersDF(df_out)
  # df_out <- df_out[complete.cases(df_out),]  #remove NAs
  
  for(var in seq(1,len_vars)){
    
    # Print plot tile and number for HTML output
    print(paste("Boxplot (fraction of Baseline) # ", var, ":", sep = ""))
    print(target[var])
    
    plot_bar <-
      
      ggplot(df_out, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = case_group, fill = case_group), size = 2, alpha = 0.5) +
      geom_line(aes(group = ID, color = case_group), size = 1.2, alpha = 0.5) +
      scale_fill_manual(
        values = colors_custom,
        aesthetics = c("colour", "fill"), name = "Grouped By:")
    
    # Print each plot to HTML
    show(plot_bar)
    
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot_F_", target[var], "_out", ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot_F_", target[var], "_out", ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
  }
}


# 08/2021 AM
# Box plots with dots and lines - by group - interactive via Plotly
# Expects a list of variable list as input
# Group and color vars are hard coded! 
# Usage: 
# var_list <- names(df)
# var_list <- var_list[-c(1,2)]
# func_PlotBoxByGroup(df, var_list)
func_PlotBoxByGroupPlotly <- function(df, var_list) {
  
  target<- var_list
  # get the number of variables
  len_vars <- length(target)  
  
  plotlist = list()
  
  for(var in seq(1,len_vars)){
    
    # Print plot title and number for HTML output
    # print(paste("Boxplot # ", var, ":", sep = ""))
    # print(target[var])
    
    plot_bar <-
      
      ggplot(df, aes_string("Timepoint", target[var])) +
      geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour = "grey30") +
      geom_point(aes(group = ID, color = ID), size = 2, alpha = 0.6) +
      geom_line(aes(group = ID, color = ID), size = 1.2, alpha = 0.4)
      
    # make paths and file names for figures
    figure_path <- paste("../DATA/Figures")
    figure_png <- paste("fig_boxplot__", target[var], ".png", sep="")
    figure_svg <- paste("svg_fig_boxplot__", target[var], ".svg", sep="")
    
    # save .png and .svg figures
    ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    ggsave(filename = figure_svg, path = figure_path)
    
    # Print each plot to HTML via plotly (interactive)
    # show(ggplotly(plot_bar, tooltip = c("group")))
    plotlist[[var]] = ggplotly(plot_bar,
                               width = 672,
                               height= 480,
                               tooltip = c("group", "y"))
  }
  
  htmltools::tagList(setNames(plotlist, NULL))
}



# 09/2021 AM 
# Plot scatter and linear model
# Arguments: df with values, 1st column is ID
# Arguments: cor table generated with cor table function
# Usage: 
# func_PlotCorScatter(df, df_table)
func_PlotCorScatter <- function(df, df_table) {
  
  # Libraries
  require(cowplot)
  require(ggpmisc) 
  require(ggrepel)
  require(ggpp)
  
  # Set the first column (PF ID) to factor, if not set previously
  df[1] <- lapply(df[1], factor)
  
  # notation for linear model formula to display on plot
  my_formula <- y ~ x

  for(i in 1:nrow(df_table)){
    # Get information out of row into vars
    nr <- unlist(df_table[i,1])
    x <- unlist(df_table[i,2])
    y <- unlist(df_table[i,3])
    r <- unlist(df_table[i,4])
    p <- unlist(df_table[i,5])
    cor_meth <- unlist(df_table[i,6])
    N <- unlist(df_table[i,7])
    
    # Get ymax xmin for placing text label on graph
    # xmin <- df %>% select(x) %>% min()
    # ymax <- df %>% select(y) %>% max()

    # p value condition
    p_value = ifelse(p < 0.001,
                     ", p < 0.001, ",
                     paste(", p = ", round(p, digits = 3), ", ", sep = ""))
    
    # Paste plot title, number and info for HTML/plot output
    text1 <- paste("Plot # ", nr, ":", "  \n", sep = "")
    text2 <- paste(x, " vs ", y,  ", r =", round(r, digits = 2),
                   p_value,
                   "  \n", cor_meth, ", N = ", N, "  \n", sep = "")
    text3 <- paste("r = ", round(r, digits = 2),
                   p_value,
                   cor_meth, ",  N = ", N, sep = "")

    # Print to console
    cat(text1) 
    cat(text2)
    
    # ggplot
    p <-
      ggplot(df, aes_string(x, y)) +
      geom_smooth(method = "lm", color="red3", formula = my_formula) +
      stat_poly_eq(formula = my_formula, 
                   aes(label = paste(
                     ..eq.label.., ..rr.label.., ..p.value.label..,
                     sep = "*plain(\",\")~~~")),
                   parse = TRUE) + 
      geom_text_repel(label = df[,1], max.overlaps = 2,
                      size = 2.4, color = "grey50") +
      geom_text_npc(aes(npcx = 0.055, npcy = 1, label = text3, 
                        check_overlap = T)) +
      geom_point(size = 2.5, alpha = 0.6, color = "blue3") +
      theme_cowplot(12)
    
    
      # theme(legend.position = "none") # removes legend
    
    # # make paths and file names for figures - save turned off
    # figure_path <- paste("../DATA/Figures")
    # figure_png <- paste("fig_boxplot__", target[var], ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot__", target[var], ".svg", sep="")
    # 
    # # save .png and .svg figures
    # ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
    # Print each plot to HTML
    show(p)
  }

}



# 09/2021 AM 
# Plot scatter and linear model
# Arguments: df with values, 1st column is ID
# Arguments: cor table generated with cor table function
# Usage: 
# func_PlotCorScatter(df, df_table)
func_PlotCorScatterGroup <- function(df, df_table) {
  
  # Libraries
  require(cowplot)
  require(ggpmisc) 
  require(ggrepel)
  require(ggpp)
  
  # Set the first two columns (PF ID, Group) to factor, if not set previously
  df[1:2] <- lapply(df[1:2], factor)
  
  # notation for linear model formula to display on plot
  my_formula <- y ~ x
  
  for(i in 1:nrow(df_table)){
    # Get information out of row into vars
    nr <- unlist(df_table[i,1])
    x <- unlist(df_table[i,2])
    y <- unlist(df_table[i,3])
    r <- unlist(df_table[i,4])
    p <- unlist(df_table[i,5])
    cor_meth <- unlist(df_table[i,6])
    N <- unlist(df_table[i,7])
    
    # Get ymax xmin for placing text label on graph
    # xmin <- df %>% select(x) %>% min()
    # ymax <- df %>% select(y) %>% max()
    
    # p value condition
    p_value = ifelse(p < 0.001,
                     ", p < 0.001, ",
                     paste(", p = ", round(p, digits = 3), ", ", sep = ""))
    
    # Paste plot title, number and info for HTML/plot output
    text1 <- paste("Plot # ", nr, ":", "  \n", sep = "")
    text2 <- paste(x, " vs ", y,  ", r =", round(r, digits = 2),
                   p_value,
                   "  \n", cor_meth, ", N = ", N, "  \n", sep = "")
    text3 <- paste("r = ", round(r, digits = 2),
                   p_value,
                   cor_meth, ",  N = ", N, sep = "")
    
    # Print to console
    cat(text1) 
    cat(text2)
    
    # ggplot
    p1 <-
      ggplot(df, aes_string(x, y)) +
      facet_wrap(. ~ Group, scales = "free") +
      geom_smooth(method = "lm", color="red3", formula = my_formula) +
      stat_poly_eq(formula = my_formula, 
                   aes(label = paste(
                     ..rr.label.., ..p.value.label..,
                     sep = "*plain(\",\")~~~")),
                   parse = TRUE) + 
      geom_text_repel(label = df[,1], max.overlaps = 2,
                      size = 2.4, color = "grey50") +
      # geom_text_npc(aes(npcx = 0.055, npcy = 1, label = text3, 
      # check_overlap = T)) +
      geom_point(aes(group = Group, color = Group, fill = Group),
                 size = 2.5, alpha = 0.6) +
      scale_fill_manual(
        values = colors_custom_3, aesthetics = c("colour", "fill")) +
      theme_cowplot(12) +
      theme(legend.position = "none") # removes legend
    
    
    
    p2 <-
      ggplot(df, aes_string(x, y)) +
      geom_smooth(method = "lm", color="red3", formula = my_formula) +
      stat_poly_eq(formula = my_formula, 
                   aes(label = paste(
                     ..eq.label.., ..rr.label.., ..p.value.label..,
                     sep = "*plain(\",\")~~~")),
                   parse = TRUE) + 
      geom_text_repel(label = df[,1], max.overlaps = 2,
                      size = 2.4, color = "grey50") +
      # geom_text_npc(aes(npcx = 0.055, npcy = 1, label = text3,
      # check_overlap = T)) +
      geom_point(aes(group = Group, color = Group, fill = Group),
                 size = 2.5, alpha = 0.6) +
      scale_fill_manual(
        values = colors_custom_3, aesthetics = c("colour", "fill")) +
      theme_cowplot(12)
    
    p <- cowplot::plot_grid(p1, NULL, p2,  nrow = 2,
                            align = "h",
                            rel_widths = c(4, 1, 1),
                            greedy = TRUE,
                            axis = "b")
    
    
    # theme(legend.position = "none") # removes legend
    
    # # make paths and file names for figures - save turned off
    # figure_path <- paste("../DATA/Figures")
    # figure_png <- paste("fig_boxplot__", target[var], ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot__", target[var], ".svg", sep="")
    # 
    # # save .png and .svg figures
    # ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
    # Print each plot to HTML
    show(p)
  }
  
}



# 09/2021 AM 
# Plot scatter and linear model - PLOTLY - plot labels do not work!
# Arguments: df with values, 1st column is ID
# Arguments: cor table generated with cor table function
# Usage: 
# func_PlotCorScatterPlotly(df, df_table)
func_PlotCorScatterPlotly <- function(df, df_table) {
  
  # Set firt column (PF ID) to factor, if not set previously
  df[1] <- lapply(df[1], factor)
  
  # notation for linear model formula to display on plot
  my_formula <- y ~ x
  
  # Plotly plotlist
  plotlist = list()
  
  for(i in 1:nrow(df_table)){
    # Get information out of row into vars
    x <- unlist(df_table[i,2])
    y <- unlist(df_table[i,3])
    r <- unlist(df_table[i,4])
    p <- unlist(df_table[i,5])
    cor_meth <- unlist(df_table[i,6])
    N <- unlist(df_table[i,7])
    
    # Paste Info for plot - does not work with plotly!
    text1 <- paste("Plot # ", i, ":", sep = "")
    text2 <- paste(x, " vs ", y,  ", r =", round(r, digits = 2),
                   ", p = ", round(p, digits = 4),
                   ", ", cor_meth, ", N = ", N, sep = "")
    text3 <- paste("r = ", round(r, digits = 2),
                   ",  p = ", round(p, digits = 4),
                   ", ", cor_meth, ",  N = ", N, sep = "")
    # NOTE: Printing does not work with plotly
    
    # ggplot
    p <-
      ggplot(df, aes_string(x, y)) +
      geom_smooth(method = "lm", formula = my_formula) +
      # stat_poly_eq(formula = my_formula, 
      #              aes(label = paste(
      #                ..eq.label.., ..rr.label.., ..p.value.label..,
      #                sep = "*plain(\",\")~~~")),
      #              parse = TRUE) + 
      # geom_text_npc(aes(npcx = 0.055, npcy = 1, label = text3, 
      #                   check_overlap = T)) +
      theme_minimal() +
      geom_point(aes(group = ID, color = ID), size = 2, alpha = 0.8) +
      theme(legend.position = "none") # removes legend
    
    # # make paths and file names for figures - save turned off
    # figure_path <- paste("../DATA/Figures")
    # figure_png <- paste("fig_boxplot__", target[var], ".png", sep="")
    # figure_svg <- paste("svg_fig_boxplot__", target[var], ".svg", sep="")
    # 
    # # save .png and .svg figures
    # ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
    # ggsave(filename = figure_svg, path = figure_path)
    
    # Print each plot to HTML via plotly (interactive)
    plotlist[[i]] = ggplotly(p,
                               width = 480,
                               height= 350,
                               tooltip = c("group"))
  }
  
  htmltools::tagList(setNames(plotlist, NULL))
}


# # 03/2021 AM - recent IaPL version - sets of var lists
# # Box plots with dots and lines - by group
# # Expects a list of variable list as input
# # Usage: pass a list of vars with two grouping vars as elements 1-2
# # func_PlotBoxLineByGroup(df, boxplot_str_name)
# func_PlotBoxLineByGroup <- function(df, var_list) {
#   
#   # get the number of total variables
#   len_vars <- length(var_list)
#   # get timepoint group variable name, it's the 1st element of var_list
#   timepoint <- var_list[2]   
#   # get subject_id group variable name, it's the 2nd element of var_list
#   subject_id <- var_list[1]
#   
#   # select continuos vars for plotting 
#   target<- var_list[3:len_vars]
#   # get the number of variables in target
#   len_vars <- length(target)  
#   
#   for(var in seq(1,len_vars)){
#     
#     # create new subset data frame for plot, to remove NAs
#     df1 <- df[,c(target[var], timepoint, subject_id)]
#     df1 <- df1[complete.cases(df1),]  #remove NAs
#     
#     # Print table tile and number for HTML output
#     print(paste("Boxplot # ", var, ":", sep = ""))
#     print(target[var])
#     
#     # Plot and asign each plot to plot_bar
#     plot_bar <-
#       
#       ggplot(df1, aes_string(timepoint, target[var])) +
#       geom_boxplot(width = 0.3, size = 0.5, fatten = 1.5, colour="grey30") +
#       geom_point(aes_string(group = subject_id, color = subject_id, fill = subject_id), size=2, alpha=0.5) +
#       geom_line(aes_string(group = subject_id, color = subject_id), size = 1.2, alpha=0.5) +
#       theme(legend.position = "none") # removes legend 
#     
#     # Place each plot image into HTML
#     show(plot_bar)
#     
#     # make paths and file names for figures
#     figure_path <- paste("../DATA/Figures")
#     figure_png <- paste("fig_boxplot_", var, "_", target[var], ".png", sep="")
#     figure_svg <- paste("svg_fig_boxplot_", var, "_", target[var], ".svg", sep="")
#     
#     # save .png and .svg figures
#     ggsave(filename = figure_png, path = figure_path, type = "cairo-png")
#     ggsave(filename = figure_svg, path = figure_path)
#     
#   }
# }




# Custom adjustment of ggpairs scatter plots (normally used in the lower triangle) 
# Note: not generalised - fill variable is hard-coded!
func_ggpairsScatter <- function(data, mapping, ...){
  ggplot(data = data, mapping=mapping) +
    geom_point(
      size = 1, 
      alpha=0.6,
      mapping = aes_string(fill="Timepoint"))
}


# Custom adjustment of ggpairs density plots (normally used in the diagonal) 
# Note: not generalised - fill variable is hard-coded!
func_ggpairsDensity <- function(data, mapping, ...){
  ggplot(data = data, mapping=mapping) +
    geom_density(mapping = aes_string(color="Timepoint"), fill=NA)
}



# 04/2021
# Outlier removal functions
# Remove outliers from a column
func_removeOutliersCol <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.45 * IQR(x, na.rm = na.rm) # change it i back to 1.5 later!!
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Removes all outliers from a data set
func_removeOutliersDF <- function(df){
  df[] <- lapply(df, function(x) if (is.numeric(x))
    func_removeOutliersCol(x) else x)
  df
}