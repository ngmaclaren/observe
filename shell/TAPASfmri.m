function TAPASfmri(infile)

addpath(genpath("~/tapas-6.1.0/rDCM"));
datadir = "/projects/academic/naokimas/neil/brains-ns50/";

pID = strrep(infile, ".txt", "");
pID = strrep(pID, datadir, "");
disp(pID)
structoutfile = pID + "output.mat";
matrixoutfile = pID + "matrix.csv";
outdir = "~/Documents/reduction/data/TAPAS/";

mat = load(infile);

Y = struct("y", mat, "dt", 0.72);
dcm = tapas_rdcm_model_specification(Y,[],[]);
output = tapas_rdcm_estimate(dcm, "r", [], 2);

save(outdir + structoutfile, "-struct", "output");
writematrix(output.Ep.A, outdir + matrixoutfile);

end
