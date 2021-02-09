function varargout = DeteksiRetakan(varargin)
% DETEKSIRETAKAN MATLAB code for DeteksiRetakan.fig
%      DETEKSIRETAKAN, by itself, creates a new DETEKSIRETAKAN or raises the existing
%      singleton*.
%
%      H = DETEKSIRETAKAN returns the handle to a new DETEKSIRETAKAN or the handle to
%      the existing singleton*.
%
%      DETEKSIRETAKAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETEKSIRETAKAN.M with the given input arguments.
%
%      DETEKSIRETAKAN('Property','Value',...) creates a new DETEKSIRETAKAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DeteksiRetakan_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DeteksiRetakan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DeteksiRetakan

% Last Modified by GUIDE v2.5 09-Sep-2020 13:49:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DeteksiRetakan_OpeningFcn, ...
                   'gui_OutputFcn',  @DeteksiRetakan_OutputFcn, ...
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


% --- Executes just before DeteksiRetakan is made visible.
function DeteksiRetakan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DeteksiRetakan (see VARARGIN)

% Choose default command line output for DeteksiRetakan
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DeteksiRetakan wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DeteksiRetakan_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in input.
function input_Callback(hObject, eventdata, handles)
% hObject    handle to input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in inputbtn.
function inputbtn_Callback(hObject, eventdata, handles)
% hObject    handle to inputbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global im
axes(handles.axes1);
% upload file
[filename,pathname] = uigetfile('*.*','Choose the input image');
im = imread([pathname,filename]);

% set the image size
scale = 600/(max(size(im(:,:,1))));        
im = imresize(im,scale*size(im(:,:,1)));

% Image resize
[m,n,~] = size(im);

imshow(im);


% --- Executes on button press in detectbtn.
function detectbtn_Callback(hObject, eventdata, handles)
% hObject    handle to detectbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global im
axes(handles.axes2);

% Convert image RGB to gray scale
I = rgb2gray(im);

% Image enhancment
% low pass filter
[f1,f2]=freqspace(size(I),'meshgrid');
D=100/size(I,1);
LPF = ones(9); 
r=f1.^2+f2.^2;
for i=1:9
    for j=1:9
        t=r(i,j)/(D*D);
        LPF(i,j)=exp(-t);
    end
end
% applying filter
Y=fft2(double(I)); Y=fftshift(Y);
Y=convn(Y,LPF); Y=ifftshift(Y);
I_en=ifft2(Y);
% blurr image
I_en=imresize(I_en,size(I)); 
I_en=uint8(I_en);
I_en=imsubtract(I,I_en);
I_en=imadd(I_en,uint8(mean2(I)*ones(size(I))));

% Segmentation image
level = roundn(graythresh(I_en),-2); % Calculate threshold using  Otsu's method
BW = ~imbinarize(I_en,level);  % Convert image to binary image using threshold
BW = double(BW);

% Removing noise and conecting image
i = 25; BW1 = BW;
while 1
    BW2 = BW1; i = i + 1;
    BW_dilate = imdilate(BW1,strel('disk',i));  % dialate image
    BW_connect = bwmorph(BW_dilate,'bridge',inf);      % connecting close parts
    BW_fill = imfill(BW_connect,'holes');            % filling small spaces & holes
    BW_erode = imerode(BW_fill,strel('disk',i-1));   % erode image
    tmp = bwareafilt(BW_erode,1);              % get size of biggest connected shape
    tmp = fix(0.05*sum(sum(tmp)));        % size considered noise
    BW_open  = bwareaopen(BW_erode,tmp);           % remove isolated pixels
    CC = bwconncomp(BW_open);
    if CC.NumObjects<2,break;end          % break the loop at convergence
end
B = bwboundaries(BW_open); % Cracks boundaries

imshow(im);hold on
for i=1:length(B)
    tmp = B{i};
    plot(tmp(:,2),tmp(:,1),'b','LineWidth',2);
end
hold off;
warning on
