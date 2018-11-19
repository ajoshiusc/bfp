function CFn_opts = regOptsMI(opts)
% returns options for NMI-registration as structure 

CFn_opts.nbins = opts.nbins;
CFn_opts.win_width = opts.parzen_width;
CFn_opts.nthreads = opts.nthreads;
CFn_opts.log_lookup = opts.log_lookup;
CFn_opts.log_thresh = opts.log_thresh;

if opts.log_lookup
   load('log_lookup.mat');
   CFn_opts.log_lookup1 = log_lookup1;
   CFn_opts.log_lookup2 = log_lookup2;
   CFn_opts.range1 = range1;
   CFn_opts.range2 = range2;
   CFn_opts.scale1 = ones(opts.nbins*opts.nbins,1)/range1(1);
   CFn_opts.scale2 = ones(2*opts.nbins, 1)/range2(1);
   clear log_lookup1 log_lookup2 range1 range2
end
end
