

VizCor <- function(x, title, method){
  # Function which creates a ggplot object of the correlation matrix based on a data.frame on tidy format!
  #
  # Args:
  #     x: A data frame where data should be on tidy format. One observation per row and each column is a variable.
  #     title: The title of the ggplot figure
  #     method: What type of method should be used when calculating the correlation matrix? 
  #             - Method options: "spearman","pearson","kendall" 
  #
  # Returns:
  #     A ggplot2 object where the variables are ordered. 
  
  # Order columns
  x <- x[,order(colnames(x))]
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  upper <- get_upper_tri(cor(x, method = method))
  # reorder to get the label and lower part right.
  p1 <- ggplot(
    within(
      melt(upper, na.rm=TRUE),
      Var2 <- ordered(Var2, levels=rev(sort(unique(Var2))) 
      )),
    aes(Var1,Var2,fill=value)) +
    geom_raster() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Correlation") +
    coord_fixed() +
    xlab("") +
    ylab("")
  return(p1)
}
