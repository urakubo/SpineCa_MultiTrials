%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mod by H Urakubo                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
function  main_HeadNeckLen

	Addpaths();
	%%
	All_Tot_Len  = [];
        All_Head_Len = [];
        All_Neck_Len = [];
        for ID = [1, 3, 5, 6, 8, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 25, 26];
% 27
          [Tot_Len, Head_Len, Neck_Len] = HeadNeckLen(ID)
	  All_Tot_Len  = [ All_Tot_Len , Tot_Len ];
          All_Head_Len = [ All_Head_Len, Head_Len ];
          All_Neck_Len = [ All_Neck_Len, Neck_Len ];	
	end;
	%%
        OutputDir = 'Data';
	FILENAME = sprintf('./%s/Head_Len.dat',OutputDir);
        dlmwrite(FILENAME, All_Head_Len,'precision','%g','delimiter','\n');
	FILENAME = sprintf('./%s/Neck_Len.dat',OutputDir);
        dlmwrite(FILENAME, All_Neck_Len,'precision','%g','delimiter','\n');
	FILENAME = sprintf('./%s/Tot_Len.dat',OutputDir);
        dlmwrite(FILENAME, All_Tot_Len,'precision','%g','delimiter','\n');
%%%
%%%
function  [Tot_Len, Head_Len, Neck_Len] =  HeadNeckLen(ID)
%%%
%%%
	grey = [0.6,0.6,0.6];

%%%
%%% Load data
%%%

	fname = sprintf('./Morph/%02g/170706_2220Final.mat',ID);
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
	fprintf('Number of branches for a centerline: %g\n', numel(skel));
	CenterLine = [];
	for i=1:length(skel);
		CenterLine = [CenterLine; flipud(skel{i})];
	end;
	if (length(skel) == 1);
		CenterLine = skel{1};
	elseif (length(skel) == 2)
		CenterLine =  [skel{1}(1:end-1,:); flipud(skel{2})] % just for No.16
	else
		fprintf('Please manually set the centerline. \n');
		return;
	end;

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
%%% Graph plot for confirmation
%%%

% GraphPlot1(YY,XX,ZZ,SpineHeadI,SpineNeckI,DendriteI,grey,az,el, HeadLine, NeckLine);

%%%
%%%
%%% Functions
%%%
%%%


function Addpaths()

	set(groot,'defaultAxesFontName',    'Arial');
	set(groot,'defaultTextFontName',    'Arial');
	set(groot,'defaultLegendFontName',  'Arial');
	set(groot,'defaultColorbarFontName','Arial');

function TotalDist = TotalDistance(skeleton)

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



function PSDcenter = SurfTriBaryCenter(node,PSDface);
	PSDcenter = [];
	for i = 1:numel(PSDface(:,1));
		tmp = node(PSDface(i,1),:) + node(PSDface(i,2),:) + node(PSDface(i,3), :);
		PSDcenter = [PSDcenter; tmp/3];
	end;
	PSDcenter = sum(PSDcenter)/length(PSDcenter(:,1));
	PSDcenter = PSDcenter(:,1:3);




function CenterL = SelectSegmentCrossingDomain(CenterLine, SpineI);
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


function skel = SelectBranchCrossingDomain(skel, SpineI);
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


function GraphPlot1(YY,XX,ZZ,SpineHeadI,SpineNeckI,DendriteI,grey,az,el, HeadLine, NeckLine);

	fig = figure;
	view(az,el);
	fv1   = isosurface(YY,XX,ZZ,DendriteI,0.5);
	p1    = patch(fv1,'FaceColor',grey,'EdgeColor','none','FaceAlpha',.5);
	hold on;
	fv2   = isosurface(YY,XX,ZZ,SpineNeckI,0.5);
	p2    = patch(fv2,'FaceColor',[0.3,0.3,1],'EdgeColor','none','FaceAlpha',.5);
	fv3   = isosurface(YY,XX,ZZ,SpineHeadI,0.5);
	p3    = patch(fv3,'FaceColor','b','EdgeColor','none','FaceAlpha',.5);
	plot3(	HeadLine(:,2), ...
			HeadLine(:,1), ...
			HeadLine(:,3), ...
		'o-','Color','r','MarkerSize',2,'LineWidth',0.5);
	plot3(	NeckLine(:,2), ...
			NeckLine(:,1), ...
			NeckLine(:,3), ...
		'o-','Color','g','MarkerSize',2,'LineWidth',0.5);
	view(az,el);
	axis([4.5, 8, 4, 6, 0, 1.5]);
	axis equal;



