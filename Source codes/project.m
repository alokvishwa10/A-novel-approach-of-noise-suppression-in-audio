function varargout = project(varargin)
% PROJECT MATLAB code for project.fig
%      PROJECT, by itself, creates a new PROJECT or raises the existing
%      singleton*.
%
%      H = PROJECT returns the handle to a new PROJECT or the handle to
%      the existing singleton*.
%
%      PROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT.M with the given input arguments.
%
%      PROJECT('Property','Value',...) creates a new PROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project

% Last Modified by GUIDE v2.5 21-Jun-2021 22:31:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project_OpeningFcn, ...
                   'gui_OutputFcn',  @project_OutputFcn, ...
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


% --- Executes just before project is made visible.
function project_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project (see VARARGIN)

% Choose default command line output for project
handles.output = hObject;
fprintf('<------------------ Noise suppression project --------------> \n\n\n')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes project wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGui=guidata(handles.figure1);
[file,path] = uigetfile('.wav','Select an audio file');
global y Fe
[y,Fe] = audioread([path file]);
fprintf('-> Step 1: Loading Audio file ----------->');
% [y,Fe]=audioread('Noise.wav');
x=y(100000:end,1);  %remove the beginning of the sample
%Nx=length(x);
fprintf(' Done\n\n');
myGui.freqSam=Fe;
myGui.datasound=y;
myGui.x=x;
myGui.player=audioplayer(myGui.datasound,myGui.freqSam);
myGui.flag=2;
%handles.in=Fe;
%handles.in=y;
%guidata(hObject,handles);
guidata(handles.figure1,myGui);

time=(1:length(y))/Fe;  % Time vector on x-axis 
set(handles.text3,'String',[file]);
axes(handles.axes3);
plot(time,y(:,1));
guidata(handles.figure1,myGui);


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play
myGui=guidata(handles.figure1);
if(myGui.flag==2)
    myGui.flag=1;
    %disp('2');
    play(myGui.player);
else
    if(myGui.flag == 1)
        %disp('1');
        myGui.flag=0;
        pause(myGui.player);
        fprintf('-> Pause\n');
    else
        %disp('0');
        myGui.flag=1;
        resume(myGui.player)
        fprintf('-> Resume\n');
    end
end
guidata(handles.figure1,myGui);

function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop
myGui=guidata(handles.figure1);
myGui.flag=2;
stop(myGui.player);
guidata(handles.figure1,myGui);
fprintf('-> Stopped \n\n');


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
myGui=guidata(handles.figure1);
%data=get(handles.slider1,'value')
f=get(handles.slider1,'value'); %
set(myGui.player,'TimerFcn',{@ejecutar,f}) %f == handles.slider1 
play(myGui.player)

function []=ejecutar(~,~,h)
v=get(h,'Value');
set(h,'Value',v+0.005)



% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)
% hObject    handle to process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myGui=guidata(handles.figure1);
Fe=myGui.freqSam;
y=myGui.datasound;
x=myGui.x;
%algorithm parameters
apriori_SNR=1;  %select 0 for aposteriori SNR estimation and 1 for apriori (see [2])
alpha=0.05;      %only used if apriori_SNR=1
beta1=0.5;
beta2=1;
lambda=3;

%STFT parameters
NFFT=1024;
window_length=round(0.031 * Fe); 
window=hamming(window_length);
window = window(:);
overlap=floor(0.45*window_length); %number of windows samples without overlapping

%Signal parameters
t_min=0.4;    %interval for learning the noise
t_max=1.00;   %spectrum (in second)

%construct spectrogram 
[S,F,T] = spectrogram(x+i*eps,window,window_length-overlap,NFFT,Fe); %put a short imaginary part to obtain two-sided spectrogram
[Nf,Nw]=size(S);

%----------------------------%
%        noisy spectrum      %
%          extraction        %
%----------------------------%
fprintf('-> Step 2: Extracting noise spectrum ---->');
t_index=find(T>t_min & T<t_max);
absS_Noise=abs(S(:,t_index)).^2;
Noise_spectrum=mean(absS_Noise,2); %average spectrum of the Noise (assumed to be ergodic))
Noise_specgram=repmat(Noise_spectrum,1,Nw);
fprintf(' Done\n\n');

%---------------------------%
%       Estimate SNR        %
%---------------------------%
fprintf('-> Step 3: Estimating SNR --------------->');
absS=abs(S).^2;
SNR_est=max((absS./Noise_specgram)-1,0); % a posteriori SNR
if apriori_SNR==1
    SNR_est=filter((1-alpha),[1 -alpha],SNR_est);  %a priori SNR: see [2]
end    
fprintf(' Done\n\n');

%---------------------------%
%  Compute attenuation map  %
%---------------------------%
fprintf('-> Step 4: Compute TF attenuation map --->');
an_lk=max((1-lambda*((1./(SNR_est+1)).^beta1)).^beta2,0);  %an_l_k or anelka, sorry stupid french joke :)
STFT=an_lk.*S;
fprintf(' Done\n\n');

%--------------------------%
%   Compute Inverse STFT   %
%--------------------------%
fprintf('-> Step 5: Compute Inverse STFT --------->');
ind=mod((1:window_length)-1,Nf)+1;
output_signal=zeros((Nw-1)*overlap+window_length,1);

for indice=1:Nw %Overlapp add technique
    left_index=((indice-1)*overlap) ;
    index=left_index+[1:window_length];
    temp_ifft=real(ifft(STFT(:,indice),NFFT));
    output_signal(index)= output_signal(index)+temp_ifft(ind).*window;
end
fprintf(' Done\n\n\n');
audiowrite('output_audio.wav',output_signal,Fe);
[out,Freq]=audioread('output_audio.wav');
myGui.datasound1=out;
myGui.freqSamp1=Freq;
myGui.flag1=2;
myGui.player1=audioplayer(myGui.datasound1,myGui.freqSamp1);


%-----------------    Display Figure   ------------------------------------      

%show temporal signals
figure
subplot(2,1,1);
t_index=find(T>t_min & T<t_max);
plot([1:length(x)]/Fe,x);
xlabel('Time (s)');
ylabel('Amplitude');
hold on;
noise_interval=floor([T(t_index(1))*Fe:T(t_index(end))*Fe]);
plot(noise_interval/Fe,x(noise_interval),'r');
hold off;
legend('Original signal','Noise Only');
title('Original Sound');
%show denoised signal
subplot(2,1,2);
plot([1:length(output_signal)]/Fe,output_signal );
xlabel('Time (s)');
ylabel('Amplitude');
title('Sound without Noise');

%show spectrogram
t_epsilon=0.001;
figure
S_one_sided=max(S(1:length(F)/2,:),t_epsilon); %keep only the positive frequency
pcolor(T,F(1:end/2),10*log10(abs(S_one_sided))); 
shading interp;
colormap('hot');
title('Spectrogram: speech + Noise');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure
S_one_sided=max(STFT(1:length(F)/2,:),t_epsilon); %keep only the positive frequency
pcolor(T,F(1:end/2),10*log10(abs(S_one_sided))); 
shading interp;
colormap('hot');
title('Spectrogram: speech only');
xlabel('Time (s)');
ylabel('Frequency (Hz)');



time=(1:length(output_signal))/Fe;  % Time vector on x-axis 

axes(handles.axes4);
plot(time,output_signal(:,1));



guidata(handles.figure1,myGui);


% --- Executes on key press with focus on play and none of its controls.
function play_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in play2.
function play2_Callback(hObject, eventdata, handles)
% hObject    handle to play2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play2
myGui=guidata(handles.figure1);

if(myGui.flag1==2)
    myGui.flag1=1;
    %disp('2');
    play(myGui.player1);
else
    if(myGui.flag1 == 1)
        %disp('1');
        myGui.flag1=0;
        pause(myGui.player1);
        fprintf('-> Pause\n');
    else
        %disp('0');
        myGui.flag1=1;
        resume(myGui.player1)
        fprintf('-> Resume\n');
    end
end
guidata(handles.figure1,myGui);



% --- Executes on button press in stop2.
function stop2_Callback(hObject, eventdata, handles)
% hObject    handle to stop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop2
myGui=guidata(handles.figure1);
myGui.flag1=2;
stop(myGui.player1);
guidata(handles.figure1,myGui);
fprintf('-> Stopped \n\n');


% --- Executes on button press in stop.
