function usc_rigid_reg(moving_filename, static_filename, output_filename, similarity, moving_mask)
%A wrapper for USC's rigid registration.
opts.similarity = similarity;
opts.moving_mask = moving_mask;
register_files_affine(moving_filename, static_filename, output_filename, opts);

