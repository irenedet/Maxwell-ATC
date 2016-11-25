function rhs=assemble_rhs(filename)
    fidff=fopen(filename,'r');
    nang=fscanf(fidff,'%f \n',1);
    nphi=fscanf(fidff,'%f \n',1);
    N=2*nphi*nphi*nang*nang;
    rhs=zeros(nphi*nang,6*nphi*nang);
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
        ll=row(5);
        
        %Rotate the incidence angles
        if ii<= (nphi-1)/2
            iii = nphi-ii-1;
        else
            iii = nphi-ii+1;
        end

        if jj<= nang/2
            jjj=jj+nang/2;
        else
            jjj=jj-nang/2;
        end

        rhs(nang*(iii-1)+jjj,nang*(kk-1)+6*(mm-1)+3*ll+1)=row(6)+1i*row(7);
        rhs(nang*(iii-1)+jjj,nang*(kk-1)+6*(mm-1)+3*ll+2)=row(8)+1i*row(9);
        rhs(nang*(iii-1)+jjj,nang*(kk-1)+6*(mm-1)+3*ll+3)=row(10)+1i*row(11);
    end
end