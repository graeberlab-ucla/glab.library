require("plotly")
require("htmlwidgets")

# pca = an example dataframe with plotting components (e.g. PCA) and additional info
# df created as a temp copy that will be manipulated to facilitate plotting
df = pca   # df = readRDS("~/Dropbox/df.rds")
# samples.plus = samples df with some additional custom columns
# this next step adds info from samples df to the plotting df
df$plot.type = samples.plus$plot.type[match(df$sample_key, samples.plus$sample_key)]

colors1 = c()
colors1 = c('#BF382A', '#0C4B8E', "grey", "grey", "green", "orange", "yellow", "red")
plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~plot.type
        , colors =colors1) %>%
  add_markers() %>% layout(scene = list(xaxis = list(title = 'PC1'), yaxis = list(title = 'PC2'), zaxis = list(title = 'PC3')))

#save as an html widget
fig3d <- plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~plot.type,
                 colors = colors1)  %>%
  add_markers() %>% layout(scene = list(xaxis = list(title = 'PC1'), yaxis = list(title = 'PC2'), zaxis = list(title = 'PC3')))
fig3d
saveWidget(widget = fig3d, file = "./plotly3d.html", selfcontained = T)


