setGeneric("compute", signature="object",
           function(object,
                    ...)
           standardGeneric("compute")
           )


setGeneric("exclude", signature="x",
           function(x,
                    ...)
           standardGeneric("exclude")
           )

setGeneric("report", signature="x",
           function(x,
                    ...)
           standardGeneric("report")
           )

setGeneric("plot", 
           function(x, y,
                    ...)
           standardGeneric("plot")
           )