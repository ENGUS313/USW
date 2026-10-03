function S = recon(cfg, rf)
% cfg: acquisition-matched Resource, Trans, PData, TW, TX, Receive, TGC;
% retain RcvFilter/RcvProfile when used. PData coordinates are in wavelengths.
% rf: int16 [samples_per_angle * angles, channels, compound_frames].
% S = recon(cfg, rf); save('recon.mat','-struct','S','-v7.3');
% In the Verasonics workspace: filename = 'recon.mat'; VSX
% Image output: ImgData{1}(:,:,1,:); requires the licensed Verasonics SDK.
required = {'Resource','Trans','PData','TW','TX','Receive','TGC'};
assert(all(isfield(cfg, required)), 'Missing acquisition metadata.');
assert(isa(rf,'int16') && ~isempty(rf) && ndims(rf) <= 3, 'Invalid RF array.');
[nr, nc, nf] = size(rf);
na = numel(cfg.TX);
assert(na >= 1 && mod(nr,na) == 0 && numel(cfg.Receive) >= na, ...
    'RF layout must contain all angles in acquisition order.');
assert(nc == cfg.Resource.Parameters.numRcvChannels && ...
    nr == cfg.Resource.RcvBuffer(1).rowsPerFrame, 'RF dimensions disagree with metadata.');
assert(isscalar(cfg.PData) && cfg.PData.Size(3) == 1, 'Only 2-D reconstruction is supported.');
for k = 1:numel(required), S.(required{k}) = cfg.(required{k}); end
for name = {'RcvFilter','RcvProfile','TPC'}
    if isfield(cfg,name{1}), S.(name{1}) = cfg.(name{1}); end
end
S.Resource = struct('Parameters',cfg.Resource.Parameters);
S.Resource.Parameters.simulateMode = 2;
S.Resource.Parameters.initializeOnly = 0;
S.Resource.Parameters.startEvent = 1;
S.Resource.RcvBuffer = struct('datatype','int16', 'rowsPerFrame',nr, ...
    'colsPerFrame',nc, 'numFrames',nf);
S.Resource.InterBuffer = struct('numFrames',1);
S.Resource.ImageBuffer = struct('numFrames',nf);
S.RcvData = {rf};
S.Receive = repmat(reshape(cfg.Receive(1:na),1,[]), 1, nf);
S.Recon = repmat(struct('senscutoff',0.6, 'pdatanum',1, ...
    'rcvBufFrame',1, 'IntBufDest',[1 1], 'ImgBufDest',[1 1], 'RINums',1:na), 1, nf);
S.ReconInfo = repmat(struct('mode','accumIQ', 'txnum',1, ...
    'rcvnum',1, 'regionnum',1), 1, nf*na);
S.Event = repmat(struct('info','Reconstruct', 'tx',0, 'rcv',0, ...
    'recon',0, 'process',0, 'seqControl',0), 1, nf+1);
for f = 1:nf
    indices = (f-1)*na + (1:na);
    S.Recon(f).rcvBufFrame = f;
    S.Recon(f).ImgBufDest = [1 f];
    S.Recon(f).RINums = indices;
    for angle = 1:na
        k = indices(angle);
        S.Receive(k).bufnum = 1;
        S.Receive(k).framenum = f;
        S.Receive(k).acqNum = angle;
        S.ReconInfo(k).txnum = angle;
        S.ReconInfo(k).rcvnum = k;
    end
    S.ReconInfo(indices(1)).mode = 'replaceIQ';
    S.ReconInfo(indices(end)).mode = 'accumIQ_replaceIntensity';
    if na == 1, S.ReconInfo(indices(1)).mode = 'replaceIntensity'; end
    S.Event(f).recon = f;
end
S.SeqControl = struct('command','stop');
S.Event(end).info = 'Stop';
S.Event(end).seqControl = 1;
lambda_m = cfg.Resource.Parameters.speedOfSound / (cfg.Trans.frequency*1e6);
S.x_m = (cfg.PData.Origin(1) + (0:cfg.PData.Size(2)-1)*cfg.PData.PDelta(1))*lambda_m;
S.z_m = (cfg.PData.Origin(3) + (0:cfg.PData.Size(1)-1)*cfg.PData.PDelta(3))*lambda_m;
end
