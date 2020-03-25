#
#********************************************************************
#
#  Functions to modify grids that failed `screenmill_calibrate`
#
#  Almost all of these borrowed from Eric Edward Bryant
#
# ----------------------------------------------------------
#
#                   read_plate function
#
# reads screenmill-annotation and -calibration data for a specific `plate_id`
#
# Input
#   plate : a `plate_id` from screenmill annotations
#   dir : a directory containing screenmill data files
#   view : whether to display an image of the plate with grid overlay - if available
#
# Returns a lst of path, anno, crop, grid
#

#' read a single screenmill plate and annotation data by plate_id and display
#'
#' @param plate plate id from a screenmill annotation
#' @param dir directory containing screenmill data files
#' @param view logical display plate with grid defaults to TRUE
#' @export

read_plate <- function(plate = "2016-03-02-001-001-003", # or a previously read plate
                       dir = "data",
                       view = TRUE) {
  anno <- screenmill:::read_annotations(dir) %>% filter(.data$plate_id == !!plate)

  crop <-
    screenmill:::read_calibration_crop(dir) %>%
    semi_join(anno, by = c("template", "position"))

  grid <-
    read_csv(file.path(dir, 'screenmill-calibration-grid.csv'), col_types = cols()) %>%
    semi_join(anno, by = c("template", "position"))

  path <- file.path(dir, anno$file)

  img <- apply_calibration(path, anno, crop, grid)

  result <- lst(path, img, anno, crop, grid)
  if (view) view_plate(result)
  return(invisible(result))
}

apply_calibration <- function(path, anno, crop, grid) {
  img     <- screenmill:::read_greyscale(path, crop$invert)
  rough   <- with(crop, img[ rough_l:rough_r, rough_t:rough_b ])
  rotated <- EBImage::rotate(rough, crop$rotate)
  cropped <- with(crop, rotated[ fine_l:fine_r, fine_t:fine_b ])
  return(cropped)
}

update_plate <- function(plate, view = TRUE) {
  plate$img <- with(plate, apply_calibration(path, anno, crop, grid))
  if (view) view_plate(plate)
  return(invisible(plate))
}

#--------------------------------------------------------------------
#
#      view_plate
#
#' display single palte image with grid overlay
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param grid_color eponymous
#' @export
view_plate <- function(plate,grid_color='blue') {
  EBImage::display(plate$img, method = 'raster')
  with(plate$grid, segments(l, t, r, t, col = grid_color))
  with(plate$grid, segments(l, b, r, b, col = grid_color))
  with(plate$grid, segments(l, t, l, b, col = grid_color))
  with(plate$grid, segments(r, t, r, b, col = grid_color))
  return(invisible(plate))
}

#-------------------------------------------------------
#
#'   shifts plate grid as a unit
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param left number of pixels to move grid to the left
#' @param right number of pixels to move grid to the right
#' @param up number of pixels to move grid up
#' @param down number of pixels to move grid down
#' @param view logical display plate with grid defaults to TRUE
#' @export

shift_grid <- function(plate,
                       left = 0, right = 0, up = 0, down = 0,
                       view = TRUE) {
  plate$grid <-
    plate$grid %>%
    mutate(
      x = x + (right - left),
      y = y + (down - up),
      l = l + (right - left),
      r = r + (right - left),
      t = t + (down - up),
      b = b + (down - up)
    )
  if (view) view_plate(plate)
  return(invisible(plate))
}

#-------------------------------------------------------
#
#' shift_grid_cell
#'
#'   Adds a conditional to shift grid so that grid can be moved by colony_row,
#'    colony_col or both to specify a specific cell. The conditional should
#'    be captured with `quo` as such:
#'      \code{cond = quo(colony_row == 16 & colony_col%in%c(14:22))}
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param cond a quoted conditional statement to select range of rows or cols
#' @param left number of pixels to move grid to the left
#' @param right number of pixels to move grid to the right
#' @param up number of pixels to move grid up
#' @param down number of pixels to move grid down
#' @param view logical display plate with grid defaults to TRUE
#' @export
#
#
shift_grid_cell <- function(plate,cond,
                            left = 0, right = 0, up = 0, down = 0,
                            view = TRUE) {
  plate$grid <-
    plate$grid %>%
    mutate(
      x = if_else (!!cond,x + (right - left),x),
      y = if_else (!!cond,y + (down - up),y),
      l = if_else (!!cond,l + (right - left),l),
      r = if_else (!!cond,r + (right - left),r),
      t = if_else (!!cond,t + (down - up),t),
      b = if_else (!!cond,b + (down - up),b)
    )
  if (view) view_plate(plate)
  return(invisible(plate))
}


#
# ------------------------------------------------------------------------
#
#'              justify_grid_edge
#'
#'
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param cond a quoted conditional statement to select range of rows or cols
#' @param side one of "t", "b", "l", "r"
#' @param fcn function by which to align selected edge, median, max or min suggested
#' @param view logical display plate with grid defaults to TRUE
#' @export

justify_grid_edge <- function(plate,cond,side="t",fcn,view=TRUE) {

  target_value <- as.integer(fcn(pull(filter(plate$grid,!!cond)[,side])))
  # top or bottom
  if (side == "t" | side == "b") {
    if (side == "t") {
      plate$grid <-
        plate$grid %>%
        mutate(
          t = if_else (!!cond,target_value,t),
          y = if_else (!!cond,t + (b - t)%/%2,y),
        )
    } else { # else `b`
      plate$grid <-
        plate$grid %>%
        mutate(
          b = if_else (!!cond,target_value,b),
          y = if_else (!!cond,t + (b - t)%/%2,y),
        )
    }
 # left or right
  } else if (side == "l" | side == "r") {
    if (side == "l") {
      plate$grid <-
        plate$grid %>%
        mutate(
          l = if_else (!!cond,target_value,t),
          x = if_else (!!cond,t + (r - l)%/%2,x),
        )
    } else { # else `r`
        plate$grid <-
          plate$grid %>%
          mutate(
            r = if_else (!!cond,target_value,t),
            x = if_else (!!cond,t + (r - l)%/%2,x),
          )
      }
  }
  if (view) view_plate(plate)
  return(invisible(plate))
}

#
# ------------------------------------------------------------------------
#'
#'             standardize_grid
#'
#' Throwing some kind of bug
#'
#' @param plate plate id from a screenmill annotation
#' @param view logical display plate with grid defaults to TRUE
#' @export

standardize_grid <- function(plate,view=TRUE) {
  plate$grid <-
    plate$grid %>%
    group_by(colony_row) %>%
    mutate(
      t = median(t),
      b = median(b),
      y = median(y)
    ) %>%
    ungroup() %>%
    group_by(colony_col) %>%
    mutate(
      l = median(l),
      r = median(r),
      x = median(x)
    ) %>%
    ungroup()
  if (view) view_plate(plate)
  return(invisible(plate))
}

#
# ------------------------------------------------------------------------
#
#'              align_grid_row
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param cond a quoted conditional statement to select range of rows or cols
#' @param side one of "t", "b"
#' @param view logical display plate with grid defaults to TRUE
#' @export

align_grid_row <- function(plate,cond,side="t", view=TRUE) {

  if (side == "t") {
    # compute highest cells in selected row
    a_top <- plate$grid %>%
      filter(!!cond) %>%
      pull(t) %>%
      min()
    #message (paste0("a_top = ",a_top))
    plate$grid <-
      plate$grid %>%
      mutate(
        y = if_else (!!cond,y + (a_top - t),y),
        t = if_else (!!cond,t + (a_top - t),t),
        b = if_else (!!cond,b + (a_top - t),b)
      )
  } else if (side == "b") {
    a_bot <- plate$grid %>%
      filter(!!cond) %>%
      pull(b) %>%
      max()
    message (paste0("a_bot = ",a_bot))
    plate$grid <-
      plate$grid %>%
      mutate(
        y = if_else (!!cond,y + (b - a_bot),y),
        t = if_else (!!cond,t + (b - a_bot),t),
        b = if_else (!!cond,b + (b - a_bot),b)
      )
  } else warning ('parameter side must be one of "t" or "b"')
  if (view) view_plate(plate)
  return(invisible(plate))
}

#------------------------------------------------------------------
#
#'    morph_grid_cell
#'
#'   Modifies grid by specific edges in order to shrink or expand as necessary
#'     to capture data or avoid plate anomalies (shadows/reflections/contaminants)
#'     Uses a conditional as in shift_grid_cell so that specific cells/rows/cols can be
#'     modified. x and y are recalculated when cells are updated. Positive values for
#'     left,right,up,down move named grid edges in that direction, negative values reverse this
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param cond a quoted conditional statement to select range of rows or cols
#' @param left number of pixels to move grid to the left
#' @param right number of pixels to move grid to the right
#' @param up number of pixels to move grid up
#' @param down number of pixels to move grid down
#' @param view logical display plate with grid defaults to TRUE
#' @export

morph_grid_cell <- function(plate,cond,
                            left = 0, right = 0, up = 0, down = 0,
                            view = TRUE) {
  plate$grid <-
    plate$grid %>%
    mutate(
      x = if_else (!!cond,x + (right - left)%/%2,x),
      y = if_else (!!cond,y + (down - up)%/%2,y),
      l = if_else (!!cond,l - left,l),
      r = if_else (!!cond,r + right,r),
      t = if_else (!!cond,t - up,t),
      b = if_else (!!cond,b + down,b)
    )
  if (view) view_plate(plate)
  return(invisible(plate))
}

#---------------------------------------------------------------------------
#
#' shift_crop
#'
#' shifts the crop of a plate by the specified pixel amount in the specified
#' direction
#'
#' @param plate a screenmill plate object loaded by the read_plate function
#' @param cond a quoted conditional statement to select range of rows or cols
#' @param left number of pixels to move grid to the left
#' @param right number of pixels to move grid to the right
#' @param up number of pixels to move grid up
#' @param down number of pixels to move grid down
#' @param view logical display plate with grid defaults to TRUE
#' @export

shift_crop <- function(plate,
                       left = 0, right = 0, up = 0, down = 0,
                       view = TRUE) {
  plate$crop <-
    plate$crop %>%
    mutate(
      fine_l = fine_l + (right - left),
      fine_r = fine_r + (right - left),
      fine_t = fine_t + (down - up),
      fine_b = fine_b + (down - up)
    )
  result <- update_plate(plate, view = view)
  return(invisible(result))
}

#-----------------------------------------------------------
#
#' save_plate_calibration
#'
#' saves modified annotation data to the appropriate screenmill
#' files after modifications
#'
#' @param plate
#' @export

save_plate_calibration <- function(plate) {
  new_anno <- plate$anno
  new_crop <- plate$crop
  new_grid <- plate$grid
  dir      <- dirname(plate$path)

  anno <-
    screenmill::read_annotations(dir) %>%
    anti_join(new_anno, by = "plate_id") %>%
    bind_rows(new_anno) %>%
    write_csv(file.path(dir, "screenmill-annotations.csv"))

  crop <-
    screenmill::read_calibration_crop(dir) %>%
    anti_join(new_crop, by = c("template", "position")) %>%
    bind_rows(new_crop) %>%
    write_csv(file.path(dir, "screenmill-calibration-crop.csv"))

  grid <-
    read_csv(file.path(dir, "screenmill-calibration-grid.csv"), col_types = cols()) %>%
    anti_join(new_grid, by = c("template", "position")) %>%
    bind_rows(new_grid) %>%
    write_csv(file.path(dir, "screenmill-calibration-grid.csv"))

  return(invisible(plate))
}

#
# ------------------------------------------------------------
#
#'                 use_calibration_from
#'
#' INPUT
#'    acceptor - a single screenmill plate image with associated metadata produced by read_plate
#'                'accepts' the donated grid
#'    donor    - a single screenmill plate image with associated metadata produced by read_plate
#'                'donates' a grid
#'
#'      progamatically, the acceptor gives the identity elements of the plate (image-file-path, group, etc).
#'      Identity data is copied onto the donor object with the desired crop and grid data
#'      The the 'donor' is updated, which reloads the image using the new image path etc.
#' Copies grid annotation data from the donor onto the acceptor image
#'
#' @param acceptor A screenmill plate image object produced by `read_plate`
#' @param donor A screenmill plate image object produced by `read_plate`
#' @export

use_calibration_from <- function(acceptor, donor) {
  donor$path          <- acceptor$path
  donor$anno          <- acceptor$anno
  donor$grid$template <- acceptor$anno$template
  donor$grid$position <- acceptor$anno$position
  donor$grid$group    <- acceptor$anno$group
  donor$crop$template <- acceptor$anno$template
  donor$crop$position <- acceptor$anno$position
  result <- update_plate(donor)
  return(invisible(result))
}
