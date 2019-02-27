function usc_rigid_reg(moving_filename, static_filename, output_filename, similarity)
%A wrapper for USC's rigid registration.
opts.similarity = similarity;
register_files_affine(moving_filename, static_filename, output_filename, opts);

