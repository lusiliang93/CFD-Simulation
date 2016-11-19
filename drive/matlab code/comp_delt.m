if k==1
    dt=0.02;
else
    delta=1/dx^2+1/dy^2;
    first=Re*Pr/2/delta;
    second=Re/2/delta;
    third=dx/abs(max(max(u)));
    fourth=dy/abs(max(max(v)));
    dt=tau*min([first,second,third,fourth]);
end