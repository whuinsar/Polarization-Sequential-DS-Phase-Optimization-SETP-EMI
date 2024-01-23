function [data, count] = stamps_read(infile, lines, format)
% stamps_read  --  Read binary data file.
%
% modified from freadbk.m   by DJ 20210621



%% Handle input. (varargin?)

false=0; true=1; 

complextype=false; % defaults

if (nargin <  3) format='float32'; end %fall through
if (nargin <  2) lines=1; end %	fall through
if (nargin <  1)
    [infile, inpath] = uigetfile('*', 'Select binary input float file', 0,0);
    infile   = [inpath,infile];
    lines    = input('Enter number of lines in file: ');
    format = input('Format (enter between single quotes): ');% make a gui...
end

% Check format for complex type: 'cpx*'
if (~ischar(format)) 
    error('FREADBK: format must be string.'); 
end
if (~ischar(infile))   
    error('FREADBK: infile must be string.'); 
end
if (length(format)==3)
  if (format=='mph') 
    disp('changing mph format to cpxfloat32');
    format='cpxfloat32';
  end
  if (format=='hgt') 
    error('please use freadhgt for hgt format.');
  end
end

%% complex types defined as prepended 'cpx'
if (length(format)>6)
  if (format(1:3)=='cpx')
    complextype = true;
    format=format(4:length(format));
  end
end

%% Read from file in column vector.
fid = fopen(infile,'r');
if (fid<0)%	try one more time
  [infile, inpath] = uigetfile('*', 'Select binary input file', 0,0);
  infile = [inpath,infile];
  fid    = fopen(infile,'r');
  if (fid<0) 
      error(ferror(fid)); 
  end
end

[data,count]=fread(fid,format); % count is number of elements, not bytes...

fclose(fid);

%% real to complex
if (complextype==true)
  data = complex(data(1:2:count),data(2:2:count));
  count=count/2; % correction for complex types
end

%% convert to matrix
width = count/lines;
data  = reshape(data,width,lines).';



