measure=function(x, y, measure_type) {
  #' @param x vector
  #' @param y vector
  #' @param measure_type choice between "Correlation" or "Euclidean" distance in clustering
  #' @return actual distance measure between x and y

  if (measure_type=="Euclidean"){
    output=sqrt(sum((x-y)^2))}
  if (measure_type=="Correlation"){
    output=1-cor(x,y)}
  return(output)
}