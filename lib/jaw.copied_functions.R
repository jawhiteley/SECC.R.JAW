##################################################
# functions written by people other than me, which I copied
# (copied / acquired / adopted / etc.)
# Jonathan Whiteley		R v2.12		2011-01-26
##################################################

cleanVarName <- function(variable.name)
{
  # clean.var.name by John Myles White, package "ProjectTemplate"
  variable.name <- gsub('_', '.', variable.name, perl = TRUE)
  variable.name <- gsub('-', '.', variable.name, perl = TRUE)
  variable.name <- gsub('\\s+', '.', variable.name, perl = TRUE)
  return(variable.name)
}
