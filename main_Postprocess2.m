%%%
%%% matlab -nodisplay -nojvm -r main_Postprocess2
%%%

function main_Postprocess2

    DataPrefix = '170626_2136';
    DAT = {'Data1','Data2','Data3','Data4','Data5', 'Data6'};
    addpath('/home/urakubo-h/npy_matlab/');

    tmp = sprintf('./Data/%s/NUM.npy',DAT{1});
    NUM = readNPY(tmp)
    All_Tot_Len  = [];
    All_Head_Len = [];
    All_Neck_Len = [];
    for ID = NUM'; %'
          fprintf('ID: %g',ID);
          [Tot_Len, Head_Len, Neck_Len] = HeadNeckLen(ID, DataPrefix)
	  All_Tot_Len  = [ All_Tot_Len , Tot_Len ];
          All_Head_Len = [ All_Head_Len, Head_Len ];
          All_Neck_Len = [ All_Neck_Len, Neck_Len ];
    end;
    for j = 1:numel(DAT)
       OutputDir = sprintf('./Data/%s/',DAT{j});
       FILENAME = sprintf('%sHead_Len.dat',OutputDir);
       dlmwrite(FILENAME, All_Head_Len,'precision','%g','delimiter','\n');
       FILENAME = sprintf('%sNeck_Len.dat',OutputDir);
       dlmwrite(FILENAME, All_Neck_Len,'precision','%g','delimiter','\n');
       FILENAME = sprintf('%sTot_Len.dat',OutputDir);
       dlmwrite(FILENAME, All_Tot_Len,'precision','%g','delimiter','\n');
    end


    exit;
%%%
%%%
function  [Tot_Len, Head_Len, Neck_Len] =  HeadNeckLen(ID, DataPrefix)
%%%
%%%

%%%
%%% Load data
%%%

	fname = sprintf('./Morph/%g/%sFinal.mat',ID, DataPrefix);
	fprintf(fname);
	load(fname);

%%%
%%% Dilution of areas of SpineHead and SpineNeck
%%%

	SpineHeadI  = logical(SpineHeadI > 0.5);
	SpineNeckI  = logical(SpineNeckI > 0.5);

%%%
%%% Detect lines crossing a target spine
%%%


	skel = SelectBranchCrossingDomain(skel, SpineI);
	fprintf('Number of branches for a centerline: %g\n\n', numel(skel));
	CenterLine = [];
	for i=1:length(skel);
		CenterLine = [CenterLine; flipud(skel{i})];
	end;
        disp(skel)
	if (length(skel) == 1);
		CenterLine = skel{1};
	elseif (length(skel) == 2)
		CenterLine =  [skel{1}(1:end-1,:); flipud(skel{2})] % just for No.16
        else
	        disp(skel{1})
	        disp(skel{2})
	        disp(skel{3})
		fprintf('Please manually set the centerline. \n');
                % 18
		return;
	end;
	% elseif (ID == 18)
	% 	CenterLine =  [skel{1}(1:end-1,:); skel{2}]

%%%
%%% "PSD-near end" will be set as the first index of center line.
%%%

	BC_PSD = SurfTriBaryCenter(node, PSDface);
	SS     = norm(CenterLine(1,:)   - BC_PSD);
	SG     = norm(CenterLine(end,:) - BC_PSD);
	if (SG < SS)
		CenterLine = flipud(CenterLine);
	end;
%%%
%%%
	CenterLine = SelectSegmentCrossingDomain(CenterLine, SpineI);
	HeadLine   = SelectSegmentCrossingDomain(CenterLine, SpineHeadI);
	NeckLine   = SelectSegmentCrossingDomain(CenterLine, SpineNeckI);

	CenterLine = CenterLine * xypitch;
	HeadLine   = HeadLine   * xypitch;
	NeckLine   = NeckLine   * xypitch;
	
	Tot_Len    = TotalDistance(CenterLine);
	Head_Len   = TotalDistance(HeadLine);
	Neck_Len   = TotalDistance(NeckLine);

%%%
%%%
%%% Functions
%%%
%%%

function TotalDist = TotalDistance(skeleton)

	if numel(skeleton) == 0
	  TotalDist = 0;
	  return;
        end;

	tmpd1 = skeleton;
	tmpd2 = skeleton;

	tmpd1(1,:)    = [];
	tmpd2(end,:)  = [];
	TotD          = tmpd1 - tmpd2;
	TotalDist     = 0;
	for i = 1:numel(TotD(:,1));
		tmp = sqrt(TotD(i,:)*TotD(i,:)'); %'
		TotalDist = TotalDist + tmp;
	end;


function PSDcenter = SurfTriBaryCenter(node,PSDface)
	PSDcenter = [];
	for i = 1:numel(PSDface(:,1));
		tmp = node(PSDface(i,1),:) + node(PSDface(i,2),:) + node(PSDface(i,3), :);
		PSDcenter = [PSDcenter; tmp/3];
	end;
	PSDcenter = sum(PSDcenter)/length(PSDcenter(:,1));
	PSDcenter = PSDcenter(:,1:3);


function CenterL = SelectSegmentCrossingDomain(CenterLine, SpineI)
	L = CenterLine;
	len = numel(L(:,1));
	flag = zeros(len,1);
	for j=2:len;
		SS = SpineI(floor(L(j-1,1)), floor(L(j-1,2)), floor(L(j-1,3)));
		SG = SpineI(floor(L(j,1)), floor(L(j,2)), floor(L(j,3)));
		if (SS > 0.5)|(SG > 0.5);
			flag(j,1) = 1;
		end;
	end;
	CenterL = CenterLine(find(flag),:); %% This may cause trouble


function skel = SelectBranchCrossingDomain(skel, SpineI)
	len = numel(skel);
	flag = zeros(len,1);
	for i=1:len;
		L=skel{i};
		for j = 1:numel(L(:,1))
			SS = SpineI(floor(L(j,1)), floor(L(j,2)), floor(L(j,3)));
			if (SS > 0);
				flag(i,1) = 1;
				break;
			end;
		end;
	end;
	skel = skel(find(flag));



