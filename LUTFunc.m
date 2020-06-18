function output = LUTFunc(param)

    % local parameter
    lambda = param.Core.lambda; % wave-length
    Nsamp = param.Core.NsampPerBlock;  % data sample number
    Nch = param.Trans.numelements; % data channel (transducer element) number
    zPNum = param.PData.Size(1); % image axial size
    xPNum = param.PData.Size(2); % image lateral size
    xOrg = param.PData.Origin(1); % image origin (x axis) [lambda]
    zOrg = param.PData.Origin(3); % image origin (z axis) [lambda]
    pdeltaX = param.PData.PDelta(1); % image resolution (x axis)
    pdeltaZ = param.PData.PDelta(3); % image resolution (z axis)
    sampPerWavL = param.Core.sampPerWavL; % samples per wave-length          
    elemPoX = param.Trans.ElementPos(:,1); % transducer position (x axis) [lambda]
    elemPoZ = param.Trans.ElementPos(:,3); % transducer position (z axis) [lambda]
    alpha = param.Core.soundSpeed.system/param.Core.soundSpeed.assum; %Speed Correction Factor
    addDelay = param.Core.PAaddDelay*param.Core.soundSpeed.system/lambda; % system delay [lambda]

    % indices for a sparse matrix
    imgPo = zeros(1,zPNum*xPNum*Nch*2); % approximately assign memory
    datPo = zeros(1,zPNum*xPNum*Nch*2);
    poVal = zeros(1,zPNum*xPNum*Nch*2);
    flag = 0;
    
    for zInx = 1:zPNum
        for xInx = 1:xPNum
            xPos = xOrg+pdeltaX*(xInx-1);
            zPos = zOrg+pdeltaZ*(zInx-1);

            % PathLength (one-way)
            ci=1:Nch;
            txPath = 0;
            rxPath = sqrt((elemPoX(ci)-xPos).^2 + (elemPoZ(ci)-zPos).^2);
            totalPath = alpha*(txPath+rxPath)+addDelay; %[lambda]
            inxR = totalPath*sampPerWavL;
            effCh = find(inxR>1 & inxR+2<Nsamp);

            inxN = floor(inxR); % for interpolation
            wei = inxR-inxN;

            zzInx = zInx*ones(size(effCh));
            xxInx = xInx*ones(size(effCh));

            if(~isempty(effCh))
                imgPo(flag+(1:length(effCh))) = sub2ind([Nsamp,Nch],inxN(effCh),effCh);
                datPo(flag+(1:length(effCh))) = sub2ind([zPNum,xPNum,Nch],zzInx,xxInx, effCh);
                poVal(flag+(1:length(effCh))) = 1-wei(effCh);
                flag = flag+length(effCh);

                imgPo(flag+(1:length(effCh))) = sub2ind([Nsamp,Nch],inxN(effCh)+1,effCh);
                datPo(flag+(1:length(effCh))) = sub2ind([zPNum,xPNum,Nch],zzInx,xxInx, effCh);
                poVal(flag+(1:length(effCh))) = wei(effCh);
                flag = flag+length(effCh);
            end
        end
    end
    
    imgPo = imgPo(1:flag);
    datPo = datPo(1:flag);
    poVal = poVal(1:flag);

    dataDim_ip = Nsamp*Nch;
    dataDim_op = zPNum*xPNum*Nch;
    alignmentMat = sparse(datPo,imgPo,poVal,dataDim_op,dataDim_ip);
    
    output = alignmentMat;
end
    