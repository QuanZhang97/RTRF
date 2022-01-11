function [out] = radon_op(in,Param,operator);
%Forward and Adjoint operators for Radon RT in the time domain.
%  IN   	in:   		intput data 
%       	Param:  		parameters combination
%		Param.h:		offset
%		Param.v:		velocity
%		Param.nt:		number of samples
%		Param.dt:		time interval
%		Param.type:    1: linear 2: parablic 3: hyperbolic
%
%	   	operator: 
%% 			operator =  1 means impute is m(tau,v) and output is d(t,x) FORWARD  OP
%% 			operator = -1 means input is d(t,x) and output is m(tau,v)  ADJOINT  OP 
%      
%  OUT   out:  		output data
% 
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%
%  Example:
%  test/test_radon_recon_hyper.m
%  test/test_radon_recon_linear.m
%  test/test_radon_demul.m
%  
%  Dot test example:
%  nt=500;
%  nv=100;
%  nh=100;
%  h=1:nh;
%  dt=1;
%  type=1;
%  
%  Param.h=h;
%  Param.nt=nt;
%  Param.dt=dt;
%  Param.v=linspace(-5,10,nv);
%  Param.type=type;
% 
%  m1 = randn(nt,nv); 
% [d1 ] = radon_op(m1,Param,1);
% 
%  d2 = randn(nt,nh); 
% [m2 ] = radon_op(d2,Param,-1);
% 
% dot1 = sum(sum(d1.*d2))
% dot2 = sum(sum(m1.*m2))
% 
% % dot1 and dot2 should be equal within machine precision if the operators
% % were properly written. 
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
% 



h = Param.h;
v = Param.v;
nt = Param.nt;
dt = Param.dt;
type=Param.type;

 nh = length(h);
 nv = length(v);

 if operator == -1; m = zeros(nt,nv); end;
 if operator ==  1; d = zeros(nt,nh); end;

 if operator == -1; d = in; end; 
 if operator ==  1; m = in; end; 

 hmax=max(abs(h));
 
  for itau = 1:nt
 
    for ih = 1:nh
 
      for iv = 1:nv

%% This can also be replaced by Parabolic or linear integration 

    switch type
        case 1
             t = (itau-1)*dt + h(ih)/v(iv);	   	
            it = floor(t/dt)+1;            
        case 2
             t = (itau-1)*dt + h(ih)*h(ih)*v(iv)/hmax/hmax;   %curvature
            it = floor(t/dt)+1; 
            %if(it<=0) it=1;end
        case 3
            t = sqrt (((itau-1)*dt)^2 + (h(ih)/v(iv))^2 ) ;
            it = floor(t/dt)+1; 
        otherwise
             t = sqrt (((itau-1)*dt)^2 + (h(ih)/v(iv))^2 ) ;
            it = floor(t/dt)+1;   
    end

%% This code does not apply interpolation 

 if (it<=nt && it>0);
     
	if operator==-1;     m(itau,iv) = m(itau,iv) + d(it,  ih);  end
	if operator== 1;     d(it,  ih)   = d(it,ih) + m(itau,iv); end
 end

  end
  end
  end

 if operator == 1; out = d; end; 
 if operator ==-1; out = m; end; 

