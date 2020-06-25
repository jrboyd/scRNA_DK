#this takes a subsetted monocle object and initializes it as a new processed monocle object.
#this procedure increases local resolution of trajectories
source("Step0_setup.R")

select_dir = shiny_choose_branch(output_path = output_root_dir, allow_new = FALSE)
select_name = sub(paste0(prefix_subset, "."), "", basename(select_dir))
select_file = file.path(select_dir, file_subset)
stopifnot(file.exists(select_file))

new_name = shiny_new_name(select_name, output_path = output_root_dir)



