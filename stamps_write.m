function count = stamps_write(data, outfile, format)
% stamps_write write matrix to a binary file 
% 
% modified from fwritebk.m   by DJ 20210621


%% Handle input.
false=0; true=1;
complextype=false;

if (nargin < 1) 
    error('writefloatfile: no data specified.');
end

if (nargin < 2) 
    outfile=[]; 
end

if (strcmp(outfile,'unknown')==1) 
    outfile=[]; 
end

if (isempty(outfile))
  [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
  outfile = [outpath,outfile];
end

if (nargin < 3)
  format = 'float32';% default
  disp('writing default float32 format.');
end

%% Check format for complex type: 'cpx*'
if (~ischar(format)) 
    error('FWRITEBK: format must be string.');
end

if (~ischar(outfile))  
    error('FWRITEBK: outfile must be string.'); 
end

if (length(format)>8)
  if (format(1:3)=='cpx')
    complextype = true;
    format=format(4:length(format));
  end
end


%% Write data to file in major row order.
fid = fopen(outfile,'w');
if (fid<0) %	try one more time.
  [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
  outfile = [outpath,outfile];
  fid = fopen(outfile,'w');
end

if (fid<0) 
    error('fwritebk: outfile could not be opened.'); 
end

data = data.';

if (complextype==true)
  data=[real(data), imag(data)];
  data=reshape(data,numel(data)/2,2).';
end

count=fwrite(fid,data,format); % write data in column order

fclose(fid);

if (complextype==true) 
    count=count/2; 
end








