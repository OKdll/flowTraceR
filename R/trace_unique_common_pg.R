#' Trace unique_common categorization for proteinGroup level
#'
#' Unique_common categorizations are analyzed on proteinGroup level
#'
#' For each submitted dataframe the unique_common proteinGroup_precursor connection is analyzed to highlight potential differences in proteinGroup denotations for common precursors.
#'
#' @param input_df1 A tibble with flowTraceR's unique_common categorization for the proteinGroup_precursor connection.
#' @param input_df2 A tibble which is the counter part for input_df1 - which was used to generate the unique_common categorization for the proteinGroup_precursor connection.
#' @param analysis_name1 output tibble name for input_df1 - default is \code{"input_df1"}.
#' @param analysis_name2 output tibble name for input_df2 - default is \code{"input_df2"}.
#' @param string_analysis Logical value, default is \code{FALSE}. If TRUE, only keeps proteinGroup identifications of input_df1 in which protein denotations are not present in the counterpart - the proteinGroups of input_df2 - and vice versa.
#'
#' @author Oliver Kardell
#'
#' @import dplyr
#' @import stringr
#'
#' @return This function returns a list with two \code{tibbles} - input_df1 and input_df2 - and with the following columns :
#' \itemize{
#'  \item traceR_connected_pg_prec - categorization of the connection between proteinGroup and precursor - unique_common
#'  \item traceR_proteinGroups_input_df1 - proteinGroup denotations of input_df1 for common precursor between input_df1 and input_df2
#'  \item traceR_precursor - common precursor between input_df1 and input_df2
#'  \item traceR_proteinGroups_input_df2 - proteinGroup denotations of input_df2 for common precursor between input_df1 and input_df2
#' }
#'
#' @export
#'
#' @examples
#' # Load libraries
#' library(dplyr)
#' library(stringr)
#' library(tibble)
#'
#' # DIA-NN example data
#' diann <- tibble::tibble(
#'   "traceR_connected_pg_prec" = c("common_common", "common_unique",
#'   "unique_common", "unique_common"),
#'   "traceR_proteinGroups" = c("P02768", "P02671", "Q92496", "P04433"),
#'   "traceR_precursor" = c("AAC(UniMod:4)LLPK1", "RLEVDIDIK2",
#'   "EGIVEYPR2", "ASQSVSSYLAWYQQK2"),
#' )
#'
#' # Spectronaut example data
#' spectronaut <- tibble::tibble(
#'   "traceR_connected_pg_prec" = c("common_common", "common_unique",
#'   "unique_common", "unique_common"),
#'   "traceR_proteinGroups" = c("P02768", "P02671", "Q02985", "A0A0A0MRZ8;P04433"),
#'   "traceR_precursor" = c("AAC(UniMod:4)LLPK1", "M(UniMod:35)KPVPDLVPGNFK2",
#'   "EGIVEYPR2", "ASQSVSSYLAWYQQK2"),
#' )
#'
#' # Find difference in pg denotation
#' # string_analysis = TRUE
#' resultA <- trace_unique_common_pg(input_df1 = diann,
#'  input_df2 = spectronaut,
#'  analysis_name1 = "DIA-NN",
#'  analysis_name2 = "Spectronaut",
#'  string_analysis = TRUE)
#'
#' # Find difference in pg denotation
#' # string_analysis = FALSE
#' # compare with resultA
#' resultB <- trace_unique_common_pg(input_df1 = diann,
#'  input_df2 = spectronaut,
#'  analysis_name1 = "DIA-NN",
#'  analysis_name2 = "Spectronaut",
#'  string_analysis = FALSE)

trace_unique_common_pg <- function(input_df1,
                                input_df2,
                                analysis_name1 = "input_df1",
                                analysis_name2 = "input_df2",
                                string_analysis = FALSE) {

  #dependency##
  dependency <- list(input_df1, input_df2)
  for (i in seq_len(length(dependency))) {

      if ("traceR_connected_pg_prec" %in% colnames(dependency[[i]]) == FALSE) {
        stop("For connected_levels proteinGroup_precursor: traceR_connected_pg_prec column must be present in submitted data.")
      }

    if ("traceR_proteinGroups" %in% colnames(dependency[[i]]) == FALSE | "traceR_precursor" %in% colnames(dependency[[i]]) == FALSE) {
      stop("For tracing unique_common: traceR_proteinGroups and traceR_precursor column must be present in submitted data.")
    }

      if (sum(str_detect(string = unique(dependency[[i]]$traceR_connected_pg_prec), pattern = "unique_common")) == 0) {
        stop("No unique_common entries detected.")
      }
    }
  ###

    connecting_col <- "traceR_connected_pg_prec"
    level_col <- "traceR_proteinGroups"

  #reduce data frames
  reduced_df1 <- input_df1 %>%
    dplyr::filter(!!ensym(connecting_col) == "unique_common") %>%
    dplyr::select(all_of(connecting_col), all_of(level_col), .data$traceR_precursor) %>%
    dplyr::distinct(.data$traceR_precursor, .keep_all = TRUE)

  reduced_df2 <- input_df2 %>%
    dplyr::filter(!!ensym(connecting_col) == "unique_common") %>%
    dplyr::select(all_of(connecting_col), all_of(level_col), .data$traceR_precursor) %>%
    dplyr::distinct(.data$traceR_precursor, .keep_all = TRUE)

  #Get these entries of other input_df
  input_df1_prepared <- input_df1 %>%
    dplyr::select(all_of(level_col), .data$traceR_precursor) %>%
    dplyr::distinct(.data$traceR_precursor, .keep_all = TRUE)

  input_df2_prepared <- input_df2 %>%
    dplyr::select(all_of(level_col), .data$traceR_precursor) %>%
    dplyr::distinct(.data$traceR_precursor, .keep_all = TRUE)

  input_df1_joined <- dplyr::left_join(
    reduced_df1,
    input_df2_prepared,
    suffix = c(paste0("_", analysis_name1), paste0("_", analysis_name2)),
    by = "traceR_precursor"
  )

  input_df2_joined <- dplyr::left_join(
    reduced_df2,
    input_df1_prepared,
    suffix = c(paste0("_", analysis_name2), paste0("_", analysis_name1)),
    by = "traceR_precursor"
  )

  #start string analysis
  if (string_analysis == TRUE) {

    input_df1_string <- analyze_string_pg(input_df = input_df1_joined)
    input_df2_string <- analyze_string_pg(input_df = input_df2_joined)

   output_list <- list(
      A = input_df1_string,
      B = input_df2_string
    )

    names(output_list) <- c(analysis_name1, analysis_name2)

  if(nrow(output_list[[1]]) == 0 & nrow(output_list[[1]]) == 0){
    message("Message: No complete difference in proteinGroup denotation detected. At least one ProteinID is always part of the other proteinGoupID.")
  }

    return(output_list)

  } else if (string_analysis == FALSE) {

    output_list <- list(
      A = input_df1_joined,
      B = input_df2_joined
    )

    names(output_list) <- c(analysis_name1, analysis_name2)

    return(output_list)
  }

}
