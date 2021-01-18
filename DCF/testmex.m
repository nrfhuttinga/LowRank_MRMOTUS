% Author: Nick Zwart
% Date: 2011 aug 18
% Rev: 2011 aug 21
% A simple test of the sdc mex compilation.

'Test SDC:'
'    read crds'
% nVec: 3
% dim1: 1812
% dim2: 3
% dim3: 79
fid = fopen('spi_crds.raw');
    tmp = squeeze(fread(fid,inf, 'float32'));
    fclose(fid);
    size(tmp)
    tmp = reshape(tmp,[3,1812,3,79]);
    size(tmp) 
    crds = tmp;

'   read comparison set'
% nVec: 1
% dim1: 1812
% dim2: 3
% dim3: 79
fid = fopen('soln_DCF.raw');
    tmp = squeeze(fread(fid,inf, 'float32'));
    fclose(fid);
    size(tmp)
    tmp = reshape(tmp,[1812,3,79]);
    size(tmp) 
    comp = tmp;


'   SDC params:'
numIter = 25
effMtx  = 50
osf     = 2.1
verbose = 1

'   start SDC calc'
DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf);
% DCF=iterative_dcf_estimation(crds,25);
DCF = single(DCF); % float32



'   find difference between calc and soln'
diff = (DCF-comp).^2;
sum(sum(sum(diff)))

'   write output DCF'
tmp = DCF;
    size(tmp)
    fid = fopen('DCF.raw','w');
    fwrite(fid,tmp,'float32');
    fclose(fid);



