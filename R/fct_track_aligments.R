#' get_name
#' @description
#' A short description...
#' @param name description
#' @import JBrowseR
#' @return return a CNV plot for selected gene and sample.
#'
#' @noRd
#' @export
is_url <- function(string) {
  grepl("www.|http:|https:", string)
}

#' get_name
#' @description
#' A short description...
#' @param name description
#' @import JBrowseR
#' @return return a CNV plot for selected gene and sample.
#'
#' @noRd
#' @export
get_name <- function(string) {
  if (is_url(string)) {
    # get assembly name from URL
    name <- parse_url_name(string)
  } else {
    # get assembly name from file
    name <- parse_file_name(string)
  }
}

parse_url_name <- function(url) {
  stringr::str_split(url, "/", simplify = TRUE) %>%
    last_value() %>%
    stringr::str_split("[.]", simplify = TRUE) %>%
    first_value()
}

parse_file_name <- function(path) {
  stringr::str_split(basename(path), "[.]", simplify = TRUE) %>%
    first_value()
}

last_value <- function(x) {
  x[length(x)]
}

first_value <- function(x) {
  x[1]
}


#' get_aligment
#' @description
#' A short description...
#' @param name description
#' @import JBrowseR
#' @return return a CNV plot for selected gene and sample.
#'
#' @noRd
#' @export
get_alignment_adapter <- function(track_data, assembly) {
  # strip the .gz extension before checking for BAM or CRAM
  track_non_gz <- strip_gz(track_data)
  
  if (stringr::str_ends(track_non_gz, ".bam")) {
    index <- stringr::str_c(track_data, ".bai")
    as.character(
      stringr::str_glue(
        '"adapter": {{ ',
        '"type": "BamAdapter", ',
        '"maxHeight" : 4000, ', # Not sure it changes something
        '"style": {{ "height" : 1 }},', # Not sure it changes something
        '"bamLocation": {{ ',
        '"uri": "{track_data}" ',
        "}}, ",
        '"index": {{ "location": {{ "uri": "{index}" }} }} ',
        #"}}"
        ', "fetchSizeLimit" : 1000}}' # Not sure it changes something
      )
    )
  } else if (stringr::str_ends(track_non_gz, ".cram")) {
    index <- stringr::str_c(track_data, ".crai")
    sequence_adapter <- get_assembly_adapter(assembly)
    as.character(
      stringr::str_glue(
        '"adapter": {{ ',
        '"type": "CramAdapter", ',
        '"cramLocation": {{ ',
        '"uri": "{track_data}" ',
        "}}, ",
        '"craiLocation": {{ "uri": "{index}" }}, ',
        '"sequenceAdapter": {sequence_adapter} ',
        "}}"
      )
    )
  } else {
    stop("alignment data must be either BAM or CRAM")
  }
}

#' strip
#' @description
#' A short description...
#' @param name description
#' @import JBrowseR
#' @return return a CNV plot for selected gene and sample.
#'
#' @noRd
#' @export
strip_gz <- function(track_data) {
  if (stringr::str_ends(track_data, ".gz")) {
    stringr::str_trunc(track_data, nchar(track_data) - 3, "right", "")
  } else {
    track_data
  }
}

#' get_assembly_name
#' @description
#' A short description...
#' @param name description
#' @import JBrowseR
#' @return return a CNV plot for selected gene and sample.
#'
#' @noRd
#' @export
get_assembly_name <- function(assembly) {
  assembly_list <- jsonlite::fromJSON(assembly)
  
  assembly_list$name
}

#' track_alignments
#' @description
#' A short description...
#' @param name description
#' @import JBrowseR
#' @return return a CNV plot for selected gene and sample.
#'
#' @noRd
#' @export
track_alignments <- function (track_data, assembly) 
{
  type <- "AlignmentsTrack"
  name <- get_name(track_data)
  assembly_name <- get_assembly_name(assembly)
  track_id <- stringr::str_c(assembly_name, "_", name)
  adapter <- GermlineVarDB::get_alignment_adapter(track_data, assembly)
  as.character(stringr::str_glue("{{", "\"type\": \"{type}\", ", 
                                 "\"name\": \"{name}\", ", "\"assemblyNames\": [\"{assembly_name}\"], ", 
                                 "\"trackId\": \"{track_id}\", ", "{adapter} ", "}}"))
}