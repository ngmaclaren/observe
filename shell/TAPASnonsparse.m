function TAPASnonsparse(infile)

addpath(genpath("~/tapas-6.1.0/rDCM"));
datadir = "/projects/academic/naokimas/neil/brains-ns50/";

pID = strrep(infile, ".txt", "");
pID = strrep(pID, datadir, "");
disp(pID)
structoutfile = pID + "output.mat";
matrixoutfile = pID + "matrix.csv";
% outdir = "~/Documents/reduction/data/TAPASnonsparse/";
outdir = "/projects/academic/naokimas/neil/TAPASnonsparse/";

mat = load(infile); % infile is time series data from HCP

Y = struct("y", mat, "dt", 0.72); % HCP docs says TR = 0.72
dcm = tapas_rdcm_model_specification(Y,[],[]); % empty args specify resting state fMRI

output = tapas_rdcm_estimate(dcm, "r", [], 1); % "r" for empirical, 2 for sparse, 1 for original, empty arg is for options

save(outdir + structoutfile, "-struct", "output");
writematrix(output.Ep.A, outdir + matrixoutfile);

end