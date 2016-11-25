%% Reading background problem info
function A=assemble_farfield(filename)
    fidff=fopen(filename,'r');
    nang=fscanf(fidff,'%f \n',1);
    nphi=fscanf(fidff,'%f \n',1);
    N=2*nphi*nphi*nang*nang;
    A=zeros(nphi*nang,6*nphi*nang);
    AA=fscanf(fidff,'%i %i %i %i %i %e %e %e %e %e %e\n',[11,N])';
    for i=1:N
        row=AA(i,:);
        %dhat
        ii=row(1)+1;
        jj=row(2)+1;
        %xhat
        kk=row(3)+1;
        mm=row(4)+1;
        %polarization
        ll=row(5);%0 or 1

        A(nang*(kk-1)+mm,nang*(ii-1)+6*(jj-1)+3*ll+1)= row(5)+1i*row(6);
        A(nang*(kk-1)+mm,nang*(ii-1)+6*(jj-1)+3*ll+2)= row(7)+1i*row(8);
        A(nang*(kk-1)+mm,nang*(ii-1)+6*(jj-1)+3*ll+3)= row(9)+1i*row(10);
    end

    return A
end
