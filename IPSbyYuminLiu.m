function varargout = IPSbyYuminLiu(varargin)
% IPSBYYUMINLIU MATLAB code for IPSbyYuminLiu.fig
%      IPSBYYUMINLIU, by itself, creates a new IPSBYYUMINLIU or raises the existing
%      singleton*.
%
%      H = IPSBYYUMINLIU returns the handle to a new IPSBYYUMINLIU or the handle to
%      the existing singleton*.
%
%      IPSBYYUMINLIU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPSBYYUMINLIU.M with the given input arguments.
%
%      IPSBYYUMINLIU('Property','Value',...) creates a new IPSBYYUMINLIU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IPSbyYuminLiu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IPSbyYuminLiu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IPSbyYuminLiu

% Last Modified by GUIDE v2.5 15-Dec-2015 02:00:09

%% Initialization Procedure
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IPSbyYuminLiu_OpeningFcn, ...
                   'gui_OutputFcn',  @IPSbyYuminLiu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%%%% --- Executes just before IPSbyYuminLiu is made visible.
function IPSbyYuminLiu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IPSbyYuminLiu (see VARARGIN)

% Choose default command line output for IPSbyYuminLiu
handles.output = hObject;

%%%% show welcome image
axes(handles.axes1);
welcomeFig = imread('welcomeImg.jpg');
imshow(welcomeFig);


%%%% add figure data to the share data handles
handles.currFig = welcomeFig;
handles.lastFig = welcomeFig;
handles.nextFig = welcomeFig;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IPSbyYuminLiu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%% --- Outputs from this function are returned to the command line.
function varargout = IPSbyYuminLiu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%% --- Executes on key release with focus on figure1 and none of its controls.
function figure1_KeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)

%% File menu
% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function openfile_Callback(hObject, eventdata, handles)
%%%% open an image file
[f, map] = imagein();
if(isempty(f))
    return;
end
imshow(f, map,'initialmagnification','fit');
%%%% update the current figure
handles.lastFig = handles.currFig;
handles.currFig = f;
handles.nextFig = f;
guidata(hObject,handles);

% --------------------------------------------------------------------
function savefile_Callback(hObject, eventdata, handles)
[filename, pathname, fmtindex] = uiputfile({'*.jpg','JPG Files(*.jpg)';...
                                            '*.png','PNG Files(*.png)';...
                                            '*.tif','TIF Files(*.tif)';...
                                            },'Save Image','currFig.jpg');                                      
saveImg = handles.currFig;                                      
if(pathname~=0)
    filename = strcat(pathname,filename);
    if(fmtindex==1)
        imwrite(saveImg,filename,'jpg');
    elseif(fmtindex==2)
        imwrite(saveImg,filename,'png');
    elseif(fmtindex==3)
        imwrite(saveImg,filename,'tif');
    end
end

% --------------------------------------------------------------------
function saveasfile_Callback(hObject, eventdata, handles)
[filename, pathname, fmtindex] = uiputfile({'*.jpg','JPG Files(*.jpg)';...
                                            '*.png','PNG Files(*.png)';...
                                            '*.tif','TIF Files(*.tif)';...
                                            },'Save Image','currFig.jpg');                                 
saveImg = handles.currFig;                                      
if(pathname~=0)
    filename = strcat(pathname,filename);
    if(fmtindex==1)
        imwrite(saveImg,filename,'jpg');
    elseif(fmtindex==2)
        imwrite(saveImg,filename,'png');
    elseif(fmtindex==3)
        imwrite(saveImg,filename,'tif');
    end
end

% --------------------------------------------------------------------
function closefile_Callback(hObject, eventdata, handles)
close(handles.figure1);


%% Edit menu
% --------------------------------------------------------------------
function edit_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function undo_Callback(hObject, eventdata, handles)
handles.nextFig = handles.currFig;
g = handles.lastFig;
imshow(g);
%%%% updatae shared data
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function redo_Callback(hObject, eventdata, handles)
g = handles.nextFig;
imshow(g);
handles.currFig = g;
guidata(hObject,handles);

% % --------------------------------------------------------------------
% function debug_Callback(hObject, eventdata, handles)
% % hObject    handle to debug (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% keyboard;


%% Filtering menu
% --------------------------------------------------------------------
function filtering_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function sharpening_Callback(hObject, eventdata, handles)
f = handles.currFig;
h = [1,1,1;1,-8,1;1,1,1];
g = imfilter(f,h,'same','symmetric');
g = imsubtract(f,g);
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function smoothing_Callback(hObject, eventdata, handles)
f = handles.currFig;
h = fspecial('average',3);
g = imfilter(f,h,'same','symmetric');
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);


%% Histogram Enhancement menu
% --------------------------------------------------------------------
function histoenhancement_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function stretching_Callback(hObject, eventdata, handles)
f = handles.currFig;
g = imadjust(f,stretchlim(f),[]);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function equalization_Callback(hObject, eventdata, handles)
f = handles.currFig;
g = f;
for depth = 1:length(f(1,1,:))
    g(:,:,depth) = histeq(f(:,:,depth));
end
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);


%% Tone menu
% --------------------------------------------------------------------
function tone_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function inverseColor_Callback(hObject, eventdata, handles)
f = handles.currFig;
if (ndims(f) ~= 3) || (size(f, 3) ~= 3)
   errordlg('Input image must be RGB.','Input Error');
   return;
end
g = 255 - f;
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function gloomy_Callback(hObject, eventdata, handles)
% hObject    handle to gloomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%% to get the current image displayed on the software
f = handles.currFig; 
if (ndims(f) ~= 3) || (size(f, 3) ~= 3)
   errordlg('Input image must be RGB.','Input Error');
   return;
end

%%%% start to process the image
f = double(f);
g = f;
for depth = 1:length(f(1,1,:))
    g(:,:,depth) = f(:,:,depth).^2 / 255;
end
g(g>255) = 255;
g = uint8(g);

%%%% show the image in the window of the software
imshow(g);

%%%% store the processed image to be shared between callback functions
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function vivid_Callback(hObject, eventdata, handles)
f = handles.currFig;
if(length(f(1,1,:))~=3)
    errordlg('Input image must be RGB.','Input Error');
    return;
end

%%%% popup hints
prompt = {'Enter Intensity Factor(>0):','Enter Saturation Factor(>0):'};
name = 'Input Vivid Factors';%% dialog title
numlines = 1;
defaultanswer = {'1.5','1.5'};%% show the default answer
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
%%%% display the dialog
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
%%%% get the user input
factorIn = str2double(answer{1});
factorSa = str2double(answer{2});
%%%% verify the input
if((~isnumeric(factorIn)) || (~isnumeric(factorSa))...
    || isnan(factorIn) || isnan(factorSa) ...
    || isinf(factorIn) || isinf(factorSa))
    errordlg('Input must be numeric.','Input Error');
    return;
elseif((factorIn<0) || (factorSa<0))
    errordlg('Input must be positive.','Input Error');
    return;
end

%factorIn = 1.5;
g = rgb2hsl(f);
A = g(:,:,2);
A = factorIn*A;
A(A>1) = 1;
g(:,:,2) = A;

B = g(:,:,2);
B = factorSa*B;
B(B>1) = 1;
g(:,:,2) = B;

g = hsl2rgb(g);
imshow(g);
%%%% update the shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function nostalgia_Callback(hObject, eventdata, handles)
f = handles.currFig;
if(length(f(1,1,:))==1)
    imshow(f);
    colormap('copper(256)');
    g = f;
elseif(length(f(1,1,:))==3)
    %A = [0.393,0.769,0.189;0.349,0.686,0.168;0.272,0.534,0.131]
    [M,N] = size(f(:,:,1));
    for ii = 1:M
        for jj = 1:N
            g(ii,jj,1) = 0.393*f(ii,jj,1)+0.769*f(ii,jj,2)+0.189*f(ii,jj,3);
            g(ii,jj,2) = 0.349*f(ii,jj,1)+0.686*f(ii,jj,2)+0.168*f(ii,jj,3);
            g(ii,jj,3) = 0.272*f(ii,jj,1)+0.534*f(ii,jj,2)+0.131*f(ii,jj,3);
        end
    end
end
       
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function virid_Callback(hObject, eventdata, handles)
f = handles.currFig;
f = double(f);
if (ndims(f) ~= 3) || (size(f, 3) ~= 3)
   errordlg('Input image must be RGB.','Input Error');
   return;
end

prompt = {'Bright offset(-255~255):'};
name = 'Input Virid Bright offset';
numlines = 1;
defaultanswer = {'50'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
offset = str2double(answer{1});
if((~isnumeric(offset)) || isnan(offset) || isinf(offset))
    errordlg('Input must be number.','Input Error');
    return;
elseif((offset<-255) ||(offset>255))
    errordlg('Input must be -255~255.','Input Error');
    return;
end

g = f;
r = f(:,:,1);
gg = f(:,:,2);
b = f(:,:,3);
%offset = 50;
g(:,:,1) = (gg-b).^2/128 + offset;
g(:,:,2) = (r-b).^2/128 + offset;
g(:,:,3) = (r-gg).^2/128 + offset;
g(g>255) = 255;
g(g<0) = 0;
g = uint8(g);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function frozen_Callback(hObject, eventdata, handles)
f = handles.currFig;
if (ndims(f) ~= 3) || (size(f, 3) ~= 3)
   errordlg('Input image must be RGB.','Input Error');
   return;
end
f = double(f);
g = f;
r = f(:,:,1);
gg = f(:,:,2);
b = f(:,:,3);
g(:,:,1) = abs(r-gg-b)*3/2;
g(:,:,1) = abs(gg-b-r)*3/2;
g(:,:,1) = abs(b-r-gg)*3/2;
g(g>255) = 255;
g = uint8(g);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function colorRender_Callback(hObject, eventdata, handles)
f = handles.currFig;
f = double(f);
if(ismatrix(f))
    gray = f;
elseif(length(f(1,1,:))==3)
    gray = 255*rgb2gray(f);
else
    errordlg('Input image must be RGB.','Input Error');
    return;
end
[M,N] = size(f(:,:,1));
g = zeros(M,N,3);

prompt = {'Enter R(0~255):','Enter G(0~255):','Enter B(0~255):'};
name = 'Input RGB Render Color';
numlines = 1;
defaultanswer = {'128','128','128'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
R = str2double(answer{1});
G = str2double(answer{2});
B = str2double(answer{3});
if((~isnumeric(R)) || (~isnumeric(G)) || (~isnumeric(B)) ...
   || isnan(R) || isnan(G) || isnan(B)...
   || isinf(R) || isinf(G) || isinf(B))
    errordlg('Input must be numeric.','Input Error');
    return;
end

if(ismatrix(f))
    g(:,:,1) = R*gray/255;
    g(:,:,2) = G*gray/255;
    g(:,:,3) = B*gray/255;
else
    g(:,:,1) = R*f(:,:,1)/255;
    g(:,:,2) = G*f(:,:,2)/255;
    g(:,:,3) = B*f(:,:,3)/255;
end
g(g<0) = 0;
g(g>255) = 255;
g = uint8(g);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

%% Distortion menu
% --------------------------------------------------------------------
function distortion_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function squeeze_Callback(hObject, eventdata, handles)
f = handles.currFig;
[M,N] = size(f(:,:,1));
midX = floor(N/2);
midY = floor(M/2);
maxBound = double(M/2)/nthroot(double(M^2+N^2)/4,4);
degree = max(10,maxBound);
g = f;
for ii = 1:M
    for jj = 1:N
        offsetX = double(jj-midX);
        offsetY = double(ii-midY);
        angle = atan2(offsetY,offsetX);%%% in randian [-pi,pi]
        radius = sqrt(offsetX^2 + offsetY^2);
        radius = sqrt(radius)*degree;
        X = floor(radius*cos(angle)) + midX;
        Y = floor(radius*sin(angle)) + midY;
        if(X<1)
            X = 1;
        elseif(X>N)
            X = N;
        end
        if(Y<1)
            Y = 1;
        elseif(Y>M)
            Y = M;
        end
        for depth = 1:length(f(1,1,:))
            g(ii,jj,depth) = f(Y,X,depth);
        end
    end
end
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);



% --------------------------------------------------------------------
function convexLen_Callback(hObject, eventdata, handles)
f = handles.currFig;
g = f;
[M,N] = size(f(:,:,1));
mc = round(M/2);
nc = round(N/2);
d1 = sqrt((1-mc)^2 + (1-nc)^2);
d2 = sqrt((1-mc)^2 + (N-nc)^2);
d3 = sqrt((M-mc)^2 + (1-nc)^2);
d4 = sqrt((M-mc)^2 + (N-nc)^2);
maxdis = max([d1,d2,d3,d4]);
for ii = 1:M
    for jj = 1:N
        dis = sqrt((ii-mc)^2 + (jj-nc)^2);
        %alpha = dis/maxdis;%%%% Convex Len factor
        alpha = max(dis/maxdis,0.4);%%%% Convex Len factor
        x = round(mc + alpha*(ii-mc));
        y = round(nc + alpha*(jj-nc));
        g(ii,jj,:) = f(x,y,:);
    end
end
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);


% --------------------------------------------------------------------
function whirlpool_Callback(hObject, eventdata, handles)
f = handles.currFig;
g = f;
[M,N] = size(f(:,:,1));
a = double(M/2); b = double(N/2); c = sqrt(a^2+b^2);
if(a<b)
    padNum = ceil(c - a) + 10;
else
    padNum = ceil(c - b) + 10;
end

f = padarray(f,[padNum,padNum],'replicate','both');
midY = floor(size(f(:,:,1),1)/2);
midX = floor(size(f(:,:,1),2)/2);
[M2,N2] = size(f(:,:,1));

prompt = {'Enter Rotation Factor(in radian):'};
name = 'Input Rotation Factor';
numlines = 1;
defaultanswer = {'3'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
maxdegree = str2double(answer{1});
if((~isnumeric(maxdegree)) || isnan(maxdegree) || isinf(maxdegree))
    errordlg('Input must be a number.','Input Error');
    return;
end

%maxdegree = pi;
for ii = padNum+1:padNum+M
    for jj = padNum+1:padNum+N
        offsetX = double(jj-midX);
        offsetY = double(ii-midY);
        angle = atan2(offsetY,offsetX);%%% in randian [-pi,pi]
        radius = sqrt(offsetX^2 + offsetY^2);
        %radius = sqrt(radius)*degree;
        angle = angle + maxdegree*nthroot(log(c/(radius+0.1)),3);
        %angle = angle + maxdegree*radius/c;
        X = floor(radius*cos(angle)) + midX;
        Y = floor(radius*sin(angle)) + midY;
        if(X<1)
            X = 1;
        elseif(X>N2)
            X = N2;
        end
        if(Y<1)
            Y = 1;
        elseif(Y>M2)
            Y = M2;
        end
        for depth = 1:length(f(1,1,:))
            g(ii-padNum,jj-padNum,depth) = f(Y,X,depth);
        end
    end
end

imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);
 
% --------------------------------------------------------------------
function waveForm_Callback(hObject, eventdata, handles)
prompt = {'Enter Manitude in Pixel(>0):','Enter Period in Pixel(>0):','Enter Rotation in Radian:'};
name = 'Input Wave Form Distortion Parameters';
numlines = 1;
defaultanswer = {'10','100','0.2'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
manitude = str2double(answer{1});
period = str2double(answer{2});
radian = str2double(answer{3});
if((~isnumeric(manitude)) || (~isnumeric(period)) || (~isnumeric(radian))...
   || isnan(manitude) || isnan(period) || isnan(radian)...
   || isinf(manitude) || isinf(period) || isinf(radian))
    errordlg('Input must be numeric.','Input Error');
    return;
end

f = handles.currFig;
f = double(f);
g = f;
[M,N]=size(f(:,:,1));
f = padarray(f,[manitude,manitude],'symmetric','both');
[M1,N1] = size(f(:,:,1));
for ii = 1:M
    for jj = 1:N
        y = round(ii - manitude - manitude*cos(2*pi/period*jj+radian)); 
        x = round(jj - manitude - manitude*cos(2*pi/period*ii+radian));
        for depth = 1:length(f(1,1,:))
            if((y>=1) && (y<=M1) && (x>=1) && (x<=N1))
                g(ii,jj,depth) = f(y,x,depth);
            end         
        end
    end
end
g = uint8(g);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function rightLeftMirror_Callback(hObject, eventdata, handles)
%%%% flip the image left and right
f = handles.currFig;
g = flip(f,2);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function upDownMirror_Callback(hObject, eventdata, handles)
%%%% flip the image upside down
f = handles.currFig;
g = flip(f,1);
imshow(g);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);



% --------------------------------------------------------------------
function stylization_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function lighting_Callback(hObject, eventdata, handles)
f = handles.currFig;
[M,N] = size(f(:,:,1));
Depth = length(f(1,1,:));
f = double(f);
g = f;

prompt = {'Enter Lighting Factor(0~255):'};
name = 'Input Lighting Factor';
numlines = 1;
defaultanswer = {'128'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
lightFactor = str2double(answer{1});
if((~isnumeric(lightFactor)) || isnan(lightFactor) || isinf(lightFactor)...
   || (lightFactor<0) || (lightFactor>255))
    errordlg('Input must be a number between 0 and 255.','Input Error');
    return;
end

%%%% select a region to light
uiwait(msgbox('Use leftkey to select vertices and rightkey to end selection',...
       'Select Lighting Areas','modal'));
[x,y] = getline(handles.axes1,'closed');
%lightFactor = 128;
lightCenter = [mean(x(1:end-1)),mean(y(1:end-1))];
maxDis = (x(1:end-1)-lightCenter(1)).^2 + (y(1:end-1)-lightCenter(2)).^2;
maxDis = sqrt(max(maxDis));
for ii = 1:M
    for jj = 1:N
        if(~inpolygon(ii,jj,y,x))
            continue;
        end
        dis = sqrt((ii-lightCenter(2))^2 + (jj-lightCenter(1))^2);
        addValue = lightFactor*(1-dis/maxDis);
        for depth = 1:Depth
            f(ii,jj,depth) = f(ii,jj,depth) + round(addValue);
            g(ii,jj,depth) = min(255,f(ii,jj,depth));
        end
    end
end
g(g>255) = 255;
g(g<0) = 0;
g = uint8(g);
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function paperCut_Callback(hObject, eventdata, handles)
f = handles.currFig;
[M,N] = size(f(:,:,1));
if(length(f(1,1,:))==3)
    gray = 0.3*f(:,:,1) + 0.59*f(:,:,2) + 0.1*f(:,:,3);
else
    gray = f;
end
gray = histeq(gray);
T = 128;%%%% threshold
rr = 255*ones(M,N);
gg = zeros(M,N);
bb = gg;
gg(gray>T) = 255;
bb(gray>T) = 255;
g = cat(3,rr,gg,bb);

imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function mosaic_Callback(hObject, eventdata, handles)
f = handles.currFig;
[M,N] = size(f(:,:,1));
Depth = length(f(1,1,:));
padNum = round(max([50,min([M/100,N/100])]));%10;

prompt = {'Enter Radius(>0):'};
name = sprintf('Input Mosaic Square Length(0<%d)',padNum);
numlines = 1;
defaultanswer = {'3'};
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);
if(isempty(answer))
    return;%% user cancle the input
end
rang = str2double(answer{1});
if((~isnumeric(rang)) || (rang<0) || isnan(rang) || isinf(rang))
    errordlg('Input must be positive number.','Input Error');
    return;
elseif(rang>padNum+1)
    errordlg(sprintf('Input must be smaller than %d',padNum),'Input Error');
    return;
end

%rang = round(padNum/2);
step = 2*rang+1;
g0 = padarray(f,[padNum,padNum],'symmetric');
g = g0;

for ii = padNum+1:step:padNum+M
    for jj = padNum+1:step:padNum+N
        neighbour = g0(ii-rang:ii+rang,jj-rang:jj+rang,:);
        for depth = 1:Depth
            meanValue = sum(sum(neighbour(:,:,depth)))/numel(neighbour(:,:,depth));
            g(ii-rang:ii+rang,jj-rang:jj+rang,depth) = meanValue;
        end
    end
end
g = g(padNum+1:padNum+M,padNum+1:padNum+N,:);
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function comic_Callback(hObject, eventdata, handles)
f = handles.currFig;
if(length(f(1,1,:))~=3)
    return;
end
f = double(f);
r = f(:,:,1);
gg = f(:,:,2);
b = f(:,:,3);
R = abs(2*gg-b+r).*r/256;
G = abs(2*b-gg+r).*r/256;
B = abs(2*b-gg+r).*gg/256;
gray = (R+G+B)/3;
R = gray+10;
R(R>255) = 255;
R = uint8(R);
G = R;
B = gray;
B(B>255) = 255;
B = uint8(B);
g = cat(3,R,G,B);
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function embossment_Callback(hObject, eventdata, handles)
f = handles.currFig;
[M,N] = size(f(:,:,1));

% prompt = {'Radian :'};
% name = 'Input Embossment radian';
% numlines = 1;
% defaultanswer = {'0.5'};
% options.Resize = 'on';
% options.WindowStyle = 'normal';
% options.Interpreter = 'tex';
% answer = inputdlg(prompt,name,numlines,defaultanswer,options);
% if(isempty(answer))
%     return;%% user cancle the input
% end
% radian = str2double(answer{1});
% if((~isnumeric(radian)) || isnan(radian) || isinf(radian))
%     errordlg('Input must be number.','Input Error');
%     return;
% end

radius = 1;
radian = pi/6;
offset = 127;
f = double(f);
g = f;
dKernel = [cos(radian+pi/4),-sin(radian),cos(radian+3*pi/4);...
           cos(radian),1,-1;...
           sin(radian),sin(radian),cos(radian-3*pi/4)];
for ii = radius+1:M-radius
    for jj = radius+1:N-radius
        for depth = 1:length(f(1,1,:))
            A = f(ii-radius:ii+radius,jj-radius:jj+radius,depth);
            sumValue = sum(sum(A .* dKernel));
            g(ii,jj,depth) = min(sumValue+offset,255);
        end
    end
end
g = uint8(g);
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function sketch_Callback(hObject, eventdata, handles)
f = handles.currFig;
if(ismatrix(f))
    ff = graygrad(f);
elseif((ndims(f)==3) && (size(f,3)==3))
    [~,~,ff] = colorgrad(f);
end
ff = im2uint8(ff);  
ff = 255 - ff; 
T=200; %%%% gray scale transformation threhold 
[M,N] = size(ff); 
g = zeros(M,N);  
for ii = 1:M  
    for jj = 1:N  
        if ff(ii,jj)<T  
            g(ii,jj)=0;  
        else  
            g(ii,jj)=235/(255-T)*(ff(ii,jj)-T);  
        end  
    end  
end  
g = uint8(g);
% % g = ff;
% % g = uint8(g);
% % g = histeq(g);
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function diffusion_Callback(hObject, eventdata, handles)
f = handles.currFig;
[M,N] = size(f(:,:,1));
Depth = length(f(1,1,:));
g = f;
padNum = round(max([15,min([M/100,N/100])]));%10;
g0 = padarray(f,[padNum,padNum],'symmetric');
for ii = padNum+1:padNum+M
    for jj = padNum+1:padNum+N
        xshift = randi(padNum,1)-round(padNum/2);
        yshift = randi(padNum,1)-round(padNum/2);
        for depth = 1:Depth
            g(ii-padNum,jj-padNum,depth) = g0(ii+xshift,jj+yshift,depth);
        end
    end
end
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function waterColor_Callback(hObject, eventdata, handles)
f = handles.currFig;
if (ndims(f) ~= 3) || (size(f, 3) ~= 3)
   errordlg('Input image must be RGB.','Input Error');
   return;
end
f = double(f);
[M,N] = size(f(:,:,1));
Depth = length(f(1,1,:));
g = f;
padNum = round(max([15,min([M/100,N/100])]));%10;
g0 = padarray(f,[padNum,padNum],'symmetric');
for ii = padNum+1:padNum+M
    for jj = padNum+1:padNum+N
        xshift = randi(padNum,1)-round(padNum/2);
        yshift = randi(padNum,1)-round(padNum/2);
        for depth = 1:Depth
            %g(ii-padNum,jj-padNum,depth) = g0(ii+xshift,jj+yshift,depth);
            g(ii-padNum,jj-padNum,depth) = (g0(ii+xshift,jj+yshift,depth)*g(ii-padNum,jj-padNum,depth))/255;
        end
    end
end
g = uint8(g);
colorLevel = 12;
g = rgb2hsl(g);
A = round(g(:,:,1)/colorLevel);
A = 360*(A/max(max(A)));
g(:,:,1) = A;
g = hsl2rgb(g);
h = fspecial('gaussian',[15,15],0.1);
g = imfilter(g,h);
imshow(g,[]);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function pencilDrawing_Callback(hObject, eventdata, handles)
f = handles.currFig;
if((ndims(f) ~= 3) || (size(f, 3) ~= 3))
   errordlg('Input image must be RGB.','Input Error');
   return;
end
ff = rgb2hsl(f);
gg = uint8(255*ff(:,:,3));
gg = graygrad(gg);

gg = im2uint8(gg);  
gg = 255 - gg; 
T=200; %%%% gray scale transformation threhold 
[M,N] = size(gg); 
ggg = zeros(M,N);  
for ii = 1:M  
    for jj = 1:N  
        if gg(ii,jj)<T  
            ggg(ii,jj)=0;  
        else  
            ggg(ii,jj)=235/(255-T)*(gg(ii,jj)-T);  
        end  
    end  
end 
ggg = double(ggg/255);
ff(:,:,3) = ggg;
g = hsl2rgb(ff);
g = uint8(g);
% % g = gg;
% % g = uint8(g);
% % g = histeq(g);
 
imshow(g);
%%%% updatae shared data
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);


%% Detection menu
% --------------------------------------------------------------------
function detection_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function edgeDetection_Callback(hObject, eventdata, handles)
%%%% find the edes
f = handles.currFig;
if(~ismatrix(f))
    f = rgb2gray(f);
end
g = edge(f,'canny');
imshow(g,[]);
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);

% --------------------------------------------------------------------
function faceDetection_Callback(hObject, eventdata, handles)
f = handles.currFig;
if(~isa(f,'uint8') && ~isa(f,'uint16') && ~isa(f,'double') ...
   && ~isa(f,'single') && ~isa(f,'int16'))
   errordlg('Input image must be uint8, uint16, double, single or int16.','Input Error');
   return;
end
% Create a detector object
faceDetector = vision.CascadeObjectDetector;   
% Detect faces
bbox = step(faceDetector, f); 
if(isempty(bbox))
    warndlg('No face detected!');
    %uiwait(msgbox('No face detected!','Warning'));
    return;
end
% Create a shape inserter object to draw bounding boxes around detections
shapeInserter = vision.ShapeInserter('BorderColor','Custom',...
                'CustomBorderColor',[255 0 0],'LineWidth',5); 
% Draw boxes around detected faces and display results              
g = step(shapeInserter, f, int32(bbox));    
imshow(g);  
handles.lastFig = handles.currFig;
handles.currFig = g;
guidata(hObject,handles);


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
CreateStruct.WindowStyle='modal';
CreateStruct.Interpreter='tex';
msgbox({'{Version 1.0}','Running Environment: MATLAB 2015b',...
        'Design & Created by Yumin Liu','All Rights Reserved'}...
        ,'About',CreateStruct);
