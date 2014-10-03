setGeneric("compute", signature="object",
           function(object, ...)
           standardGeneric("compute")
           )

setGeneric("exclude", signature="x",
           function(x, ...)
           standardGeneric("exclude")
           )

setGeneric("report", 
           function(x, y, ...)
           standardGeneric("report")
           )

setGeneric("plot", 
           function(x, y, ...)
           standardGeneric("plot")
           )

setGeneric("project", signature="x",
           function(x, ...)
           standardGeneric("project")
           )

setGeneric(".validity",
           function(object)
           standardGeneric(".validity")
           )
