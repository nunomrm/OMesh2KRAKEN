function varargout = plot_TLgrid(filename)
%
% plot_TLgrid - a simple Matlab script to plot TL computed by RAM
%
% $Id: $

fid = fopen( filename, 'rb' );

if ( fid == -1 )
   error( 'No shade file with that name exists; you must run a model first' );
end

%
% read the header with the problem description
%   write( 3 ) freq,zs,zr,rmax,dr,ndr,zmax,dz,ndz,zmplt,c0,np,ns,rs
%
% Note that the FORTRAN standards do not specify precisely what gets
% written to a file when using unformatted write statements (so this
% is compiler-dependent in the most general sense).
%
% Often times, the record starts with an integer word that gives the
% data length (sometimes in bytes, sometimes in words). The code below
% is for gfortran, which starts AND ends the record with an integer
% word with the data length in units of bytes.


rec_len = fread(fid, 1, 'int32');

freq  = fread(fid, 1, 'float32');
zs    = fread(fid, 1, 'float32');
zr    = fread(fid, 1, 'float32');
rmax  = fread(fid, 1, 'float32');
dr    = fread(fid, 1, 'float32');
ndr   = fread(fid, 1,   'int32');
zmax  = fread(fid, 1, 'float32');
dz    = fread(fid, 1, 'float32');
ndz   = fread(fid, 1,   'int32');
zmplt = fread(fid, 1, 'float32');
c0    = fread(fid, 1, 'float32');
np    = fread(fid, 1,   'int32');
ns    = fread(fid, 1,   'int32');
rs    = fread(fid, 1, 'float32');
lz    = fread(fid, 1,   'int32');

rec_len = fread(fid, 1, 'int32');

% read the TL reciever matrix / grid

j = 1;
rr(1) = ndr*dr;
while (rr(j) < rmax),
  rr(j+1) = rr(j) + ndr*dr;
  j = j + 1;
end;

rd = ndz*dz*(1:lz);

Nrr = length(rr);
Nrd = length(rd);

TLgrid = zeros(Nrd, Nrr);

for j = 1:Nrr,
  rec_len = fread(fid,  1,   'int32');
  TLslice = fread(fid, lz, 'float32');
  rec_len = fread(fid,  1,   'int32');
  if ~(isempty(TLslice))        % Sometimes it gets [] on the last one.
    TLgrid(:,j) = TLslice;
  end
end;

% done reading data ...

fclose( fid );

% plot the real valued TL grid

imagesc( 1.0e-03*rr, rd, TLgrid );
colormap(flipud(jet(1024)));
if ispc; colormap(flipud(jet(256))); end;
caxis([min(min(TLgrid)) 130])
shading interp;
ha = gca();

% set(ha, 'YTick', []);

if ( nargout == 1 )
   varargout( 1 ) = { ha };   % return a handle to the figure
end
