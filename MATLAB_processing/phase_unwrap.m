function [ar]=phase_unwrap(arfidata,f,startframe,varargin)
%
% [ar]=phase_unwrap(arfidata,f,startframe)
%
% f is in Hz
% startframe is first frame to look for discontinuities
% scale
% iterations - zero to do while phase wraps exist

scale = 1;
iter = 1;

if (nargin>3);
  scale = varargin{1};
  if (nargin>4)
    iter = varargin{2};
  end;
end;


c=1540;
al=(c/(f*4*pi))*1e6;
ar = arfidata;

if (iter~=0)
  for (j=1:iter)
    s=(convn(ar,cat(3,1,-1)));
    s=-(s>al*pi*scale)+(s<-al*pi*scale);
    s(:,:,1:startframe)=0;
    s=cumsum(s,3)*al*pi;
    ar = ar+s(:,:,1:end-1)*2;
  end;
else
  s=cat(3,1,1);
  while(sum(abs(reshape(s(:,:,1:end-1),[],1)))~=0)
    s=(convn(ar,cat(3,1,-1)));
    s=-(s>al*pi*scale)+(s<-al*pi*scale);
    s(:,:,1:startframe)=0;
    s=cumsum(s,3)*al*pi;  
    ar = ar+s(:,:,1:end-1)*2;
  end;
end;
