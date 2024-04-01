# read .narrowPeak files into a GRanges object
narrow_peaks_to_granges <- function(peak_dir, srr_id) {
  peak_file <- paste0(peak_dir, srr_id, "_peaks.narrowPeak")
  narrow_granges <- readNarrowPeak(peak_file)
  mcols(narrow_granges) <- cbind(mcols(narrow_granges), 
                                 DataFrame(sample = srr_id))
  narrow_granges
}