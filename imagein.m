function [I, map] = imagein(mypath)
%IMAGEIN Read image in from current-working or specified directory.
% I = IMAGEIN displays a window containing all the files in the
% current directory, and saves in I the image selected from the
% current directory.
% [I, MAP] = IMAGEIN variable MAP is required to be an output
% argument when the image being read is an indexed image.
% [ . . .] = IMAGEIN(’PATH’) is used when the image to be read
% resides in a specified directory. For example, the input
% argument ’C:\MY_WORK\MY_IMAGES’ opens a window showing
% the contents of directory MY_IMAGES. An image selected from
% that directory is read in as image I.
% Created by Yumin

%%%% if there is no input argument, then open the current-working directory
%%%% by default, otherwise specify the directory by the argument.
if nargin == 1
    currWD = cd; %% saving the current working directory to be restore
    cd(mypath);
end

%%%% open the file dialog box
% [fileName, pathName, filterIndex] = uigetfile({'*.jpg';'*.jpeg';'*.png';'*.jpeg';'*.gif';'*.bmp';'*.tif'}, 'select an image');
[fileName, pathName, filterIndex] = uigetfile({'*.jpg';'*.jpeg';'*.png';'*.bmp';'*.tif'}, 'select an image');
%%%% if the user presses Cancel, then return the function
if isnumeric(pathName) && isnumeric(fileName) && isnumeric(filterIndex)
    I = [];
    map = [];
    return;
end

%%%% combine the selected path and file name to be the reading path
readPath = [pathName,fileName];
%%%% read the selected image
[I map] = imread(readPath); 

%%%% restore the current working directory if it has been changed
if nargin == 1
    cd(currWD);
end


