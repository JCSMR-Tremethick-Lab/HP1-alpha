loadConfigFile <- function(file_name = NULL, file_type = c("YAML", "JSON")){
  stopifnot(!is.null(file_name) & file_type %in% c("YAML", "JSON"))
  if (file_type == "YAML"){
    library(yaml)
    config <- yaml::yaml.load_file(file_name)
  } else if (file_type == "JSON"){
    library(jsonlite)
    config <- jsonlite::fromJSON(file_name)
  }
  return(config)
}
