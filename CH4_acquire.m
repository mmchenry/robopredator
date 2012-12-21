function CH4_acquire
%a = raw untransformed coords of prey and pred nose in pixels
%a.preyx    [mouthtip COM  tailtip]
%a.preyy
%a.preyz

%b = transformed coords of prey in cm (0,0,0) = pred nose
%[mouthtip COM  tailtip]


%s =  conversions, raw frames numbers, and conditions

%f = flow conditions [mouthtip COM  tailtip]

%f.V mag of velocity
%f.S shear deformation
%f.Rdist    response distance of prey



%location of data files on dropbox



warning off all;
clear;


path = '/Users/WilliamStewart/Dropbox/Robopredator/behavior';

%location of raw videos
vidlocation = '/Volumes/RED/Ch4 Data';





CollectData         = 1;

%gets raw pixel locations of animals in x,y,z coords
GetPredNose         = 1;
GetPreyPosition     = 1;
GetStage2Position   = 1;


TransformCoords     = 1;

CheckResults        = 1;


batch               = 1; %run Transform coords in batch mode
clf('reset');
figure;





%% Collect Data

if CollectData
    

    %allows user to pick sequence from folder
    [seq,lit,speed,namestub] = GetSequence;
    
    
    %Gets name of file (sequence number) containing raw data
    datafile = sprintf('%03d',seq);
    datafile = ['S' datafile '.mat'];
    
    load([path filesep 'single coords/calibration' filesep...
        datafile]);
    

    prompt = {'Sequence number',...
        'Dorsal Stage 2 frame'...
        'Lateral Stage 2 frame'...
        'Lights? (0/1)',...
        'Speed',...
        'Dorsal FS frame',...
        'Lateral FS frame',...
        'Dorsal conversion? (0 = determine new)',...
        'Lateral conversion? (same for Dorsal)'};

    defaults = {num2str(seq),...
        num2str(s.DorS2frame),...
        num2str(s.LatS2frame),...
        num2str(lit),...
        num2str(speed),...
        num2str(s.DorFSframe),...
        num2str(s.LatFSframe),...
        num2str(s.Dconv),...
        num2str(s.Lconv)};
    answer = inputdlg(prompt,'',1,defaults);


    %if calibration unkown
    if str2num(cell2mat(answer(8))) == 0

        %open up calibration photos
        [fname, pname] = uigetfile('*.*',...
            'select DORSAL calibration image');

        cd(pname);
        txt = 'Pick two points that are 20 mm apart';
        imshow(fname);
        [x,y] = ChoosePoints(0,2,txt);

        distpx = ((y(2) - y(1))^2 + (x(2) - x(1))^2)^0.5;

        %pix/mm
        Dconv = distpx / 20;


        [fname, pname] = uigetfile('*.*',...
            'select LATERAL calibration image');

        txt = 'Pick two points that are 20 mm apart';
        imshow(fname);
        [x,y] = ChoosePoints(0,2,txt);

        distpx = ((y(2) - y(1))^2 + (x(2) - x(1))^2)^0.5;

        %pix/mm
        Lconv = distpx / 20;

        answer{8} = num2str(Dconv);
        answer{9} = num2str(Lconv);
        
        
    end
    
    s.DorS2frame   = str2num(answer{2});
    s.LatS2frame   = str2num(answer{3});
    s.lit          = str2num(answer{4});
    s.speed        = str2num(answer{5});
    s.DorFSframe   = str2num(answer{6});
    s.LatFSframe   = str2num(answer{7});
    s.Dconv        = str2num(answer{8});
    s.Lconv        = str2num(answer{9});

 
    save([path filesep 'single coords/calibration' filesep...
        datafile],'s');
   

    commandwindow;


%% %GetPredNose
%determines pixel coordinates of predator nose in Lat and Dor images

if GetPredNose
    
    load([path filesep 'single coords/rawdata_FS' filesep...
        datafile]);
    
    %open up dorsal FS frame
    img = [namestub num2str(s.DorFSframe) '.tif'];
    
    imshow(img);
    xlabel(img,'Interpreter','none');
    
    %pick off coordinates of predator from dorsal image 
    txt = 'select a point on each side of predator head';
    [x,y] = ChoosePoints(0,2,txt);

    [x,y] = GetPredPosition(x,y);
    
    predy = y;
    
    %Now do lateral image
    img = ['Lat' namestub(4:end) num2str(s.LatFSframe) '.tif'];
    
    imshow(img);
    xlabel(img,'Interpreter','none');
    
    txt = 'Select tip of pred nose';
    [x,y] = ChoosePoints(0,1,txt);
    hold on;
    plot(x,y,'g.');
    
    predx = x;
    predz = y;
    
    load([path filesep 'single coords/rawdata_FS' filesep...
        datafile]);
    
    a.predx = predx;
    a.predy = predy;
    a.predz = predz;
    
    save([path filesep 'single coords/rawdata_FS' filesep...
        datafile],'a');
    
    clear a;
end

%%GetPreyPosition
%determines pixel coordinates of prey [head COM tail] in Lat and Dor
%images

if GetPreyPosition
    
   
    close;
    
    %Call up dorsal FS image to specify y coords of prey
    img = [namestub num2str(s.DorFSframe) '.tif'];
    imshow(img);
    xlabel(img,'Interpreter','none');
    
    txt = 'select 3 points along prey: head tip, post. of sb, &tail tip';
    [x,y] = ChoosePoints(0,3,txt);
    
    preyy = zeros(1,3);
    preyy(1) = y(1);
    preyy(3) = y(3);
    
    %determine COM as midpoint between head and swim bladder
    preyy(2) = (y(1) + y(2))/2;
    
    
    
    %Now call up lateral FSimage to specify x z coords
    img = ['Lat' namestub(4:end) num2str(s.LatFSframe) '.tif'];
    
    close;
     
    imshow(img);
    xlabel(img,'Interpreter','none');
    
    txt = 'select 3 points along prey: head tip, post. of sb, &tail tip';
    [x,y] = ChoosePoints(0,3,txt);
    
    
    preyx = zeros(1,3);
    preyz = zeros(1,3);
    
    preyx(1) = x(1);
    preyx(3) = x(3);
    
    preyx(2) = (x(1) + x(2))/2;
    
    preyz(1) = y(1);
    preyz(3) = y(3);
    
    preyz(2) = (y(1) + y(2))/2;
    
    
    load([path filesep 'single coords/rawdata_FS' filesep...
        datafile]);
    
    a.preyz(:) = preyz;
    a.preyx(:) = preyx;
    a.preyy(:) = preyy;
    
    
    save([path filesep 'single coords/rawdata_FS' filesep...
        datafile],'a');
       
    clear a;
end

%% GetStage2Position
%Measures pixel coordinates of prey at end of stage 2 response

if GetStage2Position
    
    close;
    
    %Call up dorsal Stage 2 image to specify y coords of prey
    img = [namestub num2str(s.DorS2frame) '.tif'];
    imshow(img);
    xlabel(img,'Interpreter','none');
    
    txt = 'select 3 points along prey: head tip, post. of sb, &tail tip';
    [x,y] = ChoosePoints(0,3,txt);
    
    preyy = zeros(1,3);
    preyy(1) = y(1);
    preyy(3) = y(3);
    
    %determine COM as midpoint between head and swim bladder
    preyy(2) = (y(1) + y(2))/2;
    
    
    
    %Now call up lateral Stage 2 image to specify x z coords
    img = ['Lat' namestub(4:end) num2str(s.LatS2frame) '.tif'];
    
    close;
     
    imshow(img);
    xlabel(img,'Interpreter','none');
    
    txt = 'select 3 points along prey: head tip, post. of sb, &tail tip';
    [x,y] = ChoosePoints(0,3,txt);
    
    
    preyx = zeros(1,3);
    preyz = zeros(1,3);
    
    preyx(1) = x(1);
    preyx(3) = x(3);
    
    preyx(2) = (x(1) + x(2))/2;
    
    preyz(1) = y(1);
    preyz(3) = y(3);
    
    preyz(2) = (y(1) + y(2))/2;
    
    load([path filesep 'single coords/rawdata_stage2' filesep...
        datafile]);
    
    %save raw coords as a.preyz2 etc. 
    a2.prey2z(:) = preyz;
    a2.prey2x(:) = preyx;
    a2.prey2y(:) = preyy;
    
    
    save([path filesep 'single coords/rawdata_stage2' filesep...
        datafile],'a2');
    
    clear a2;
    
    
    
end




end











%% Transform Coords

if TransformCoords
%Transfroms raw prey coords to coords in predframe of reference
%pred nose at 0,0,0 (x,y,z), z = vertical, y = lateral
%position is in  cm


    %seq = 1;
    
    %[seq,lit,speed,namestub] = GetSequence;

    while 1
    
    
        %Gets name of file (sequence number) containing raw data
        datafile = sprintf('%03d',seq);
        datafile = ['S' datafile '.mat'];

        %load s
        load([path filesep 'single coords/calibration' filesep...
            datafile]);

        %load a
        load([path filesep 'single coords/rawdata_FS' filesep...
            datafile]);

        %load a2
        load([path filesep 'single coords/rawdata_stage2' filesep...
            datafile]);


        

        %Get raw coords of predator nose
        predxraw    = a.predx;
        predyraw    = a.predy;
        predzraw    = a.predz;

        %Transform coords of FS first

        preyxraw    = a.preyx;
        preyyraw    = a.preyy;
        preyzraw    = a.preyz;

        preyx       = preyxraw - predxraw;
        preyy       = (preyyraw - predyraw) * -1;
        preyz       = (preyzraw - predzraw) * -1; %-1 b/c pixel axis. 

        %convert to CM
        preyx       = (preyx ./ s.Lconv) / 10; 
        preyy       = (preyy ./ s.Dconv) / 10;
        preyz       = (preyz ./ s.Lconv) / 10;

        
        %Now transform coords of Stage 2 response
        
        preyxraw2    = a2.prey2x;
        preyyraw2    = a2.prey2y;
        preyzraw2    = a2.prey2z;

        preyx2       = preyxraw2 - predxraw;
        preyy2       = (preyyraw2 - predyraw) * -1;
        preyz2       = (preyzraw2 - predzraw) * -1; %-1 b/c pixel axis. 

        %convert to CM
        preyx2       = (preyx2 ./ s.Lconv) / 10; 
        preyy2       = (preyy2 ./ s.Dconv) / 10;
        preyz2       = (preyz2 ./ s.Lconv) / 10;
        
        
        %adjust x position because predator was moving!
        tdiff = s.LatS2frame - s.LatFSframe; %in frames
        tdiff = tdiff / 500;  % convert to seconds given 500 fps
        x_corr = s.speed * tdiff;
        preyx2 = preyx2 + x_corr;
        
        
        %now save b
        load([path filesep 'Transformed_Prey_Coords.mat']);

        
        b.preyx(seq,:)  = preyx;
        b.preyy(seq,:)  = preyy;
        b.preyz(seq,:)  = preyz;
        
        b.preyx2(seq,:)  = preyx2;
        b.preyy2(seq,:)  = preyy2;
        b.preyz2(seq,:)  = preyz2;
        
        b.speed(seq)    = s.speed;
        b.lit(seq)      = s.lit;
        b.LL(seq)       = s.LL;
        

        save([path filesep 'Transformed_Prey_Coords.mat'],'b');
        
 
        
        
        
        
        
        %visualize raw and transformed coords
        if CheckResults
            
            clf;
            
          
            
            
            subplot(4,2,1),
            
            [namestub] = NameStubber(s,seq,1,1);
            imshow([vidlocation filesep sprintf('%03d',seq) ...
                filesep namestub]);
            hold on;
            for i = 1:3
                plot([0 640],[a.preyy(i) a.preyy(i)],'b');
            end
            
            subplot(4,2,2),
            
            [namestub] = NameStubber(s,seq,1,2);
            imshow([vidlocation filesep sprintf('%03d',seq) ...
                filesep namestub]);
            hold on;
            plot(a.preyx,a.preyz,'bo');
            
            
            subplot(4,2,3),
            
            [namestub] = NameStubber(s,seq,2,1);
            imshow([vidlocation filesep sprintf('%03d',seq) ...
                filesep namestub]);
            hold on;
            for i = 1:3
                plot([0 640],[a2.prey2y(i) a2.prey2y(i)],'b');
            end
            
            subplot(4,2,4),
            
            [namestub] = NameStubber(s,seq,2,2);

            imshow([vidlocation filesep sprintf('%03d',seq) ...
                filesep namestub]);
            hold on;
            plot(a2.prey2x,a2.prey2z,'bo');

            
            subplot(4,2,5:8);
            
            plot3(b.preyx(seq,:),b.preyy(seq,:),b.preyz(seq,:),'b');
            hold on;
            plot3(b.preyx(seq,1),b.preyy(seq,1),b.preyz(seq,1),'b*');
            
            plot3(b.preyx2(seq,:),b.preyy2(seq,:),...
                b.preyz2(seq,:),'m');
            plot3(b.preyx2(seq,1),b.preyy2(seq,1),...
                b.preyz2(seq,1),'m*');
            
            arrow([b.preyx(seq,2) b.preyy(seq,2) b.preyz(seq,2)],...
                [b.preyx2(seq,2) b.preyy2(seq,2) ...
                b.preyz2(seq,2)],20,45,60,10);
            
            
            
            [xp, yp, zp] = ellipsoid(-2, 0, 0, 2, 0.5, 0.5);
            c = ones(21,21);
            who = xp < -2;
            xp(who) = nan;
            surf(xp,yp,zp,c);
            axis equal;
            hold on;
            clear xp yp zp c;
            
            set(gca,'view',[90 0]);
            rotate3d;
            
            

            if batch == 1
                pause;
            end
            
            
            
            
            
            
            
            
            
            
        end
        
        if batch == 0
            break
        else
            disp(num2str(seq));
            seq = seq+1;
        end

        
        clear b s a a2;

        
    end
    
end
 



commandwindow;

end






%% NameStubber

function [namestub] = NameStubber(s,seq,stage,orientation)
%creates namestub for calling up individual image files
% stage...1 = response, 2 = stage 2
% orientation...1 = Dorsal, 2 = Lateral

if orientation == 1
    view = 'Dor';
    if stage == 1
        frame = num2str(s.DorFSframe);
    else
        frame = num2str(s.DorS2frame);
    end
    
else
    view = 'Lat';
    if stage == 1
        frame = num2str(s.LatFSframe);
    else
        frame = num2str(s.LatS2frame);
    end 
end

if s.lit == 0
    lights = 'Da';
else 
    lights = 'Li';
end

speed = sprintf('%02d',s.speed);

seqs = sprintf('%03d',seq);

namestub = [view '_' lights '_' speed '_seq' seqs '-'...
    num2str(frame) '.tif'];





end










%% GET SEQUENCE

function [seq,lit,speed,namestub] = GetSequence


DataLocation  = '/Volumes/RED/Ch4 Data';


d = uigetdir;
seq = str2num(d(end-2:end));
cd(d);

D = dir(d);

if strcmp(D(4).name(1:2),'._')
    namestub = D(4).name(3:19);
else
    namestub = D(4).name(1:17);  %change if change pad!
end

lit = strcmp(namestub(5:6),'Li');
speed = str2num(namestub(8:9));


end



%% CHOOSE POINTS

function [x,y] = ChoosePoints(img,numpoints,txt)
%Used for finding coordinate points on a static image 'img'.
warning off all
%if size(img,1) > 1
%    imshow(img);

%end


title(txt)
hold on;
set(gcf,'DoubleBuffer','on');
%disp(' '); disp(' ');
%disp('Left mouse button picks points.');disp(' ');
%disp('Right mouse button removes last point.');disp(' ');
%disp('Press return to stop.')
n = 0;
but = 1;
while 1
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==1 || but==92
        n = n+1;
        x(n) = xi;
        y(n) = yi;
        
        h = plot(x,y,'ro');
        if numpoints > 0 && n >= numpoints
            break
        end
    
    elseif but==8
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        xrange = xlim;
        yrange = ylim;
       
        
        hold off
        %imshow(img);
        title(txt);
        set(gca,'xlim',xrange);
        set(gca,'ylim',yrange);
        hold on
        
        h = plot(x,y,'ro');
        
    elseif but == 3 || but == 61
       zoomcenter(xi,yi,2);
    elseif but == 31
       zoom out;
    end
end

%delete(h)

%make column vectors!!!!
x = x'; y = y';  
warning on all
end


%%
function [xpred,ypred] = GetPredPosition(x,y)
    xrange = xlim;
    m = (y(2)-y(1)) / (x(2)-x(1));
    m2 = -1 * (1/m);
    
    
    midy = mean(y);
    midx = mean(x);
    
    b = midy - (m2*midx);
    
    hold on;
    plot(midx,midy,'ro');
    plot(xrange,polyval([m2,b],xrange));
    
    title('select tip of rostrum');
    [xpred(1),ypred(1)] = ginput(1);
    hold on;
    
    plot(xpred(1),ypred(1),'ro');
    
   
    
end