function count = fwritebkj(data, outfile, bkformat,machinefmt)
% FWRITEBK write matrix to a binary file (major row order, pixel interleaved).
%   FWRITEBK (MAT) asks for a filename to write the specified matrix to in
%   float32 format.
%
%   COUNT = FWRITEBK (MAT) optionally returns the number of elements successfully
%   written.
%
%   FWRITEBK (MAT,FILENAME) writes matrix to specified file in float32 format.
%
%   FWRITEBK (MAT,FILENAME,BKFORMAT) uses bkformat to write matrix. bkformat
%   is either the same as the format flag of FREAD, or 'cpx' is prepended
%   for complex matrices:
%     'cpxfloat32'     complex floating point, 32 bits, store pixel interleaved.
%     'cpx...'         ...
%
%   FWRITEBKJ (MAT,FILENAME,BKFORMAT, MACHINEFMT) additionally 
%   specifies the order for reading or writing bytes or bits in the file using the machinefmt argument
%   MACHINEFMT:
%   'n' or 'native'         Your system byte ordering (default)
%   'b' or 'ieee-be'        Big-endian ordering (GAMMA software, for example)
%   'l' or 'ieee-le'        Little-endian ordering
%   's' or 'ieee-be.l64'    Big-endian ordering, 64-bit long data type
%   'a' or 'ieee-le.l64'    Little-endian ordering, 64-bit long data type
%
%   Use FWRITEBK (MAT,[],BKFORMAT) to ask for filename to write MAT to.
%   Use FWRITEBK (MAT,'unknown',BKFORMAT) to ask for filename.
%
%   See also FREADBK, FOPEN, FREAD, FWRITE, LOAD, SAVE, FPRINTF, FSEEK,
%   FWRITEHGT, FREADHGT, FSIZE
%

% $Revision: 1.7 $  $Date: 2001/05/04 16:37:30 $
% Bert Kampes, 4/3/00
%   ======================================================================
%   04/2009 MJ add machinefmt
%   ======================================================================

%%% Handle input.
false=0; true=1;
complextype=false;
%
if (nargin < 1) error('writefloatfile: no data specified.'); end;
if (nargin < 2) outfile=[]; end;
if (strcmp(outfile,'unknown')==1) outfile=[]; end;
if (isempty(outfile))%                        [] form
  [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
  outfile = [outpath,outfile];
end;
if (nargin < 3)
  bkformat = 'float32';%			default
  disp('writing default float32 format.');
%  if (isempty(outfile))%			[] form
%    [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
%    outfile = [outpath,outfile];
%  end;
end;
if (nargin < 4)
  machinefmt = 'b';%			default
  disp('writing default Big-endian ordering.');
%  if (isempty(outfile))%			[] form
%    [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
%    outfile = [outpath,outfile];
%  end;
end;
% Check bkformat for complex type: 'cpx*'
if (~ischar(bkformat)) error('FWRITEBK: bkformat must be string.'); end;
if (~ischar(outfile))  error('FWRITEBK: outfile must be string.'); end;
if (length(bkformat)>8)
  if (bkformat(1:3)=='cpx')
    complextype = true;
    bkformat=bkformat(4:length(bkformat));
  end;
end;


%%% Write data to file in major row order.
fid = fopen(outfile,'w',machinefmt);
if (fid<0)%					try one more time.
  [outfile, outpath] = uiputfile('*', 'Save outputfile as', 0,0);
  outfile = [outpath,outfile];
  fid = fopen(outfile,'w',machinefmt);
end;
if (fid<0) error('fwritebk: outfile could not be opened.'); end;


data = data.';
if (complextype==true)
  data=[real(data), imag(data)];
  data=reshape(data,prod(size(data))/2,2).';
end;

count=fwrite(fid,data,bkformat);%		write data in column order
fclose(fid);
if (complextype==true) count=count/2; end;

%%% EOF








